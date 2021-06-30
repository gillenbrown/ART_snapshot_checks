from pathlib import Path

import numpy as np
import yt
from yt.extensions.astro_analysis.halo_analysis.api import HaloCatalog

yt.funcs.mylog.setLevel(50)  # ignore yt's output

from . import plot_utils

# ======================================================================================
#
# Setup for precalculation of key quantities
#
# ======================================================================================
class Galaxy(object):
    def __init__(self, sphere, rank):
        self.sphere = sphere
        self.rank = rank
        self.ds = self.sphere.ds

        # then get the mask showing which clusters are done forming
        self.mask_done_forming = self.sphere[("STAR", "age")] > 15 * yt.units.Myr

        # have a place for some things (like the CIMF) to be stored as a user
        # calculates them.
        self.precalculated = dict()

    def __getitem__(self, property):
        """
        Automatically apply the cut to only pick formed clusters
        """
        quantity = self.sphere[property]
        if property[0] == "STAR":
            quantity = quantity[self.mask_done_forming]
        return quantity

    def prop_all_clusters(self, property):
        """
        Like __getitem__, but does not apply the restriction that clusters must be
        done forming
        """
        return self.sphere[property]


class Simulation(object):
    def __init__(self, ds_path, sphere_radius_kpc):
        # get the dataset and corresponding halo file
        run_dir = ds_path.parent.parent
        halo_path = run_dir / "halos"
        halo_name = ds_path.name.replace("continuous_", "halos_")
        halo_name = halo_name.replace(".art", ".0.bin")

        self.ds = yt.load(str(ds_path))
        self.halos_ds = yt.load(str(halo_path / halo_name))

        # get the axis names and other stuff
        self.name = plot_utils.names[run_dir]
        self.axes = plot_utils.axes[self.name]
        self.color = plot_utils.colors[self.name]

        # have the place to store some precalculated things, particularly those that
        # include all galaxies in the simulation
        self.precalculated = dict()

        # if we are the old IC set, we have one galaxy, otherwise two
        # check what kind of particles are present
        if ("N-BODY_0", "MASS") in self.ds.derived_field_list:
            self.n_galaxies = 2
        else:
            self.n_galaxies = 1

        # create halo catalog object
        self.hc = HaloCatalog(halos_ds=self.halos_ds, data_ds=self.ds)
        self.hc.create(save_halos=True, save_catalog=False)

        # Do some parsing of the halo catalogs
        # We get the indices that sort it. The reversing there makes the biggest
        # halos first, like we want.
        halo_masses = yt.YTArray(
            [halo.quantities["particle_mass"] for halo in self.hc.halo_list]
        )
        rank_idxs = np.argsort(halo_masses)[::-1]
        # Add this info to each halo object, and put the halos into a new sorted list,
        # with the highest mass (lowest rank) halos first). We only keep the number
        # the user requested
        self.galaxies = []
        for rank, idx in enumerate(rank_idxs[: self.n_galaxies], start=1):
            halo = self.hc.halo_list[idx]
            radius = min(
                halo.quantities["virial_radius"], sphere_radius_kpc * yt.units.kpc
            )
            center = [
                halo.quantities["particle_position_x"],
                halo.quantities["particle_position_y"],
                halo.quantities["particle_position_z"],
            ]
            sphere = self.ds.sphere(center=center, radius=radius)
            self.galaxies.append(Galaxy(sphere, rank))


# ======================================================================================
#
# Loading datasets
#
# ======================================================================================
# Get the datasets and halo catalogs. When doing these we need to be a bit
# careful about the datasets. We will make one set of comparisons at the last
# common output of all simulations, then one with the last output of each
# simulation. Those all need to be stored separately.
def filename_to_scale_factor(filename):
    return float(filename[-10:-4])


def get_outputs_in_dir(sim_dir):
    directory = Path(sim_dir)
    if directory not in plot_utils.names:
        print(f"Skipping {directory}")
        return []

    out_dir = directory / "out"
    return [
        file
        for file in out_dir.iterdir()
        if file.is_file()
        and str(file.name).endswith(".art")
        and str(file.name).startswith("continuous_a")
    ]


def get_simulations_last(sim_dirs, sphere_radius_kpc=30):
    sims_last = []
    for directory in sim_dirs:
        all_outputs = get_outputs_in_dir(directory)
        if len(all_outputs) > 0:
            last_output = sorted(all_outputs)[-1]
            sims_last.append(Simulation(last_output, sphere_radius_kpc))

    return sims_last


def get_common_scale_factor(sim_dirs, z_max):
    # Start by getting the last common output among the production runs
    last_outputs = []
    for directory in sim_dirs:
        directory = Path(directory)
        if directory not in plot_utils.names or "stampede2/production" not in str(
            directory
        ):
            continue

        all_outputs = get_outputs_in_dir(directory)
        # restrict to be a reasonable redshift
        a_min = 1 / (1 + z_max)
        this_last_output = sorted(all_outputs)[-1]
        if filename_to_scale_factor(this_last_output.name) > a_min:
            last_outputs.append(this_last_output)

    earliest_last_output = sorted(last_outputs)[0]
    # include fudge factor for scale comparisons (so 0.1801 and 0.1802 match)
    return filename_to_scale_factor(earliest_last_output.name) + 0.001


def get_simulations_common(sim_dirs, z_max=5, sphere_radius_kpc=30):
    common_scale = get_common_scale_factor(sim_dirs, z_max)

    sims_common = []
    for directory in sim_dirs:
        all_outputs = get_outputs_in_dir(directory)

        # get the last one that's in common with the other simulations
        all_common_outputs = [
            file
            for file in all_outputs
            if filename_to_scale_factor(file.name) <= common_scale
            and abs(filename_to_scale_factor(file.name) - common_scale) < 0.02
        ]
        # if there are no snapshots early enough for this, don't add them
        if len(all_common_outputs) > 0:
            last_common_snapshot = sorted(all_common_outputs)[-1]
            sims_common.append(Simulation(last_common_snapshot, sphere_radius_kpc))

    return sims_common, common_scale
