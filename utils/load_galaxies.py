from pathlib import Path

import numpy as np
import yt
from astropy import table

from . import run_attributes

yt.funcs.mylog.setLevel(50)  # ignore yt's output


def time_units(ds, array):
    return ds._handle.tphys_from_tcode_array(array) * yt.units.year


def duration(region):
    end_time = time_units(region.ds, region[("STAR", "TERMINATION_TIME")])
    return end_time - region[("STAR", "creation_time")]


# ======================================================================================
#
# Setup for precalculation of key quantities
#
# ======================================================================================
class Galaxy(object):
    def __init__(
        self,
        ds,
        center=None,
        sphere_radius=None,
        m_vir=None,
        r_vir=None,
        rank=None,
        name=None,
    ):
        self.ds = ds
        self.center = center
        self.m_vir = m_vir
        self.r_vir = r_vir
        self.rank = rank
        self.name = name

        self.mask_done_forming = None

        # only create the sphere if the user wants to
        if sphere_radius is not None:
            self.sphere = self.ds.sphere(center=self.center, radius=sphere_radius)
        else:
            self.sphere = None

        # have a place for some things (like the CIMF) to be stored as a user
        # calculates them.
        self.precalculated = dict()

    def __getitem__(self, property):
        """
        Automatically apply the cut to only pick formed clusters
        """
        if self.sphere is None:
            raise ValueError("No sphere set on galaxy initialization")

        # then get the mask showing which clusters are done forming
        if self.mask_done_forming is None and property[0] == "STAR":
            self.mask_done_forming = self.sphere[("STAR", "age")] > 15 * yt.units.Myr

        quantity = self.sphere[property]
        if property[0] == "STAR":
            quantity = quantity[self.mask_done_forming]
        return quantity

    def prop_all_clusters(self, property):
        """
        Like __getitem__, but does not apply the restriction that clusters must be
        done forming
        """
        if self.sphere is None:
            raise ValueError("No sphere set on galaxy initialization")

        return self.sphere[property]


class Simulation(object):
    def __init__(
        self, ds_path, sphere_radius_kpc=None, min_virial=False, n_galaxies=None
    ):
        """
        ds_path must be a Path object
        """
        # get the dataset and corresponding halo file
        self.run_dir = ds_path.parent.parent
        halo_path = self.run_dir / "halos"
        halo_name = ds_path.name.replace("continuous_", "out_")
        halo_name = halo_name.replace(".art", ".list")

        self.ds = yt.load(str(ds_path))

        self.z = self.ds.current_redshift
        self.scale_factor = 1 / (1 + self.z)

        # get the axis names and other stuff. The dictionaries are defaultdicts, so
        # there is no need to worry about key errors
        self.names = run_attributes.names[self.run_dir]
        self.axes = list(self.names.keys())
        self.color = run_attributes.colors[self.run_dir]
        self.marker = run_attributes.markers[self.run_dir]
        self.ls = run_attributes.lss[self.run_dir]

        # have the place to store some precalculated things, particularly those that
        # include all galaxies in the simulation
        self.precalculated = dict()

        # if we are the old IC set, we have one galaxy, otherwise two
        # check what kind of particles are present
        if n_galaxies is None:
            if ("N-BODY_0", "MASS") in self.ds.derived_field_list:
                self.n_galaxies = 2
            else:
                self.n_galaxies = 1
        else:
            self.n_galaxies = n_galaxies

        # load halo catalogs
        cat = table.Table.read(halo_path / halo_name, format="ascii.commented_header")
        # delete some unneeded quantities. If I want to inlcude these later, I'll need
        # to verify I'm doing the units correctly.
        for col in cat.colnames:
            if col not in ["Mvir", "Vmax", "Vrms", "Rvir", "X", "Y", "Z"]:
                del cat[col]
        # To modify the units of things, we need to know little h. It's in the
        # file header. We also want a, to turn things into physical units.
        with open(halo_path / halo_name, "r") as in_file:
            line_num = 1
            for line in in_file:
                if line_num == 2:
                    a = float(line.split()[-1])
                elif line_num == 3:
                    h = float(line.split()[-1])
                    break

                line_num += 1
        assert 0 < a < 1.01  # slight buffer for last output
        assert 0.6 < h < 0.8

        # Masses are in Msun / h
        cat["Mvir"] = cat["Mvir"] / h
        # Positions in Mpc / h (comoving)
        for col in ["X", "Y", "Z"]:
            cat[col] = cat[col] / h
        # Halo Distances, Lengths, and Radii in kpc / h (comoving)
        cat["Rvir"] = cat["Rvir"] / h
        # Velocities in km / s (physical, peculiar) -- no change needed

        # add units to names
        cat.rename_column("Mvir", "Mvir_msun")
        cat.rename_column("X", "X_mpccm")
        cat.rename_column("Y", "Y_mpccm")
        cat.rename_column("Z", "Z_mpccm")
        cat.rename_column("Rvir", "Rvir_kpccm")
        cat.rename_column("Vmax", "Vmax_kms")
        cat.rename_column("Vrms", "Vrms_kms")

        # Do some parsing of the halo catalogs
        # We get the indices that sort it. The reversing there makes the biggest
        # halos first, like we want.
        rank_idxs = np.argsort(cat["Mvir_msun"])[::-1]
        # Add this info to each halo object, and put the halos into a new sorted list,
        # with the highest mass (lowest rank) halos first). We only keep the number
        # the user requested
        self.galaxies = []
        for rank, idx in enumerate(rank_idxs[: self.n_galaxies], start=1):
            row = cat[idx]

            r_vir = self.ds.arr(row["Rvir_kpccm"], "kpccm")
            M_vir = self.ds.arr(row["Mvir_msun"], "Msun")
            center = self.ds.arr(
                [row["X_mpccm"], row["Y_mpccm"], row["Z_mpccm"]], "Mpccm"
            )

            if sphere_radius_kpc is None:
                radius = None
            else:
                if min_virial:
                    radius = self.ds.arr(
                        min(r_vir.to("kpc").value, sphere_radius_kpc), "kpc"
                    )
                else:
                    radius = self.ds.arr(sphere_radius_kpc, "kpc")

            self.galaxies.append(Galaxy(self.ds, center, radius, M_vir, r_vir, rank))

        # # we need to be careful with a couple of runs. Runs with epsff=1% or runs
        # # with 10% and fboost=1 are unreliable above 10^5 Msun
        # run_dir_str = str(self.run_dir)
        # if "sfe001" in run_dir_str or (
        #     "sfe010" in run_dir_str and "fboost1" in run_dir_str
        # ):
        #     self.unreliable_mass = 1e5
        #     self.reliable = False
        # else:
        #     self.unreliable_mass = np.inf
        #     self.reliable = True
        self._determine_failed()

    def _determine_failed(self):
        # probably change this at some point
        if (
            "sfe010" not in str(self.run_dir) and "sfe001" not in str(self.run_dir)
        ) or "rj" in str(self.run_dir):
            self.reliable = True
            self.unreliable_mass = np.inf
            return
        # Figure out at what mass the durations reach a median of 14 Myr
        durations = self.func_all_galaxies(lambda g: duration(g).to("Myr").value)
        masses = self.func_all_galaxies(
            lambda g: g[("STAR", "INITIAL_MASS")].to("Msun").value
        )
        # assume I'm reliable, then break otherwise
        self.reliable = True
        self.unreliable_mass = np.inf
        dm = 0.25
        m_min = 3
        while m_min < 7:
            m_max = m_min + dm
            good_idx = np.logical_and(masses > 10 ** m_min, masses < 10 ** m_max)
            # check that there aren't too few clusters
            if np.sum(good_idx) < 10:
                m_min += dm
                continue

            this_durations = durations[good_idx]
            median = np.median(this_durations)
            if median > 14:
                self.unreliable_mass = 10 ** m_min
                self.reliable = False
                break

            m_min += dm

    def __repr__(self):
        return str(self.run_dir)

    def func_all_galaxies(self, func):
        """
        Apply one function to all galaxies, and append the results to each other.

        This is likely used for something like getting the stellar masses of all stars
        in all galaxies, or their bound fractions.
        """
        return np.concatenate([func(galaxy) for galaxy in self.galaxies])


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
    # This is useful to avoid loading things without names that won't appear on plots,
    # but it makes this code less useful elsewhere when I may want to use it to load
    # all kinds of simulations
    # if directory not in run_attributes.names:
    #     print(f"Skipping {directory}")
    #     return []

    out_dir = directory / "out"
    return [
        file
        for file in out_dir.iterdir()
        if file.is_file()
        and str(file.name).endswith(".art")
        and str(file.name).startswith("continuous_a")
    ]


def get_simulations_last(sim_dirs, sphere_radius_kpc=30, min_virial=True):
    sims_last = []
    for directory in sorted(sim_dirs):
        all_outputs = get_outputs_in_dir(directory)
        if len(all_outputs) > 0:
            last_output = sorted(all_outputs)[-1]
            sims_last.append(Simulation(last_output, sphere_radius_kpc, min_virial))

    return sims_last


def get_common_scale_factor(sim_dirs, z_max):
    # Start by getting the last common output among the provided runs
    last_outputs = []
    for directory in sim_dirs:
        directory = Path(directory)
        if directory not in run_attributes.names:
            continue

        all_outputs = get_outputs_in_dir(directory)
        if len(all_outputs) == 0:
            print(f"This has no outputs: {directory}")
            continue

        # restrict to be a reasonable redshift
        a_min = 1 / (1 + z_max)
        this_last_output = sorted(all_outputs)[-1]
        if filename_to_scale_factor(this_last_output.name) > a_min:
            last_outputs.append(filename_to_scale_factor(this_last_output.name))

    # include fudge factor for scale comparisons (so 0.1801 and 0.1802 match)
    return min(last_outputs) + 0.001


def get_simulations_common(sim_dirs, z_max=5, sphere_radius_kpc=30, min_virial=True):
    common_scale = get_common_scale_factor(sim_dirs, z_max)

    sims_common = []
    for directory in sorted(sim_dirs):
        all_outputs = get_outputs_in_dir(directory)
        if len(all_outputs) == 0:
            # happens when skipped by get_outputs_in_dir
            continue

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
            sims_common.append(
                Simulation(last_common_snapshot, sphere_radius_kpc, min_virial)
            )

    return sims_common, common_scale


def get_simulations_same_scale(
    sim_dirs, desired_z, z_tolerance=0.05, sphere_radius_kpc=30, min_virial=True
):
    sims_common = []
    for directory in sorted(sim_dirs):
        all_outputs = get_outputs_in_dir(directory)
        if len(all_outputs) == 0:
            # happens when skipped by get_outputs_in_dir
            continue

        closest_z = np.inf
        closest_snapshot = None
        for out in all_outputs:
            a = filename_to_scale_factor(out.name)
            z = (1 / a) - 1

            if abs(z - desired_z) < abs(closest_z - desired_z):
                closest_z = z
                closest_snapshot = out

        # then validate that we got close enough.
        if abs(closest_z - desired_z) / desired_z < z_tolerance:
            sims_common.append(
                Simulation(closest_snapshot, sphere_radius_kpc, min_virial)
            )
        else:
            print(f"No outputs close to z={desired_z} in {directory}")
    return sims_common


def get_plot_names(sims):
    all_plots = []
    for sim in sims:
        for ax in sim.axes:
            all_plots.append(ax)
    return list(set(all_plots))


def get_plot_names_dirs(dirs):
    all_plots = []
    for d in dirs:
        for ax in run_attributes.names[d]:
            all_plots.append(ax)
    return list(set(all_plots))
