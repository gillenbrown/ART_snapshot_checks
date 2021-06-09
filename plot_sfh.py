"""
plots_sfh.py

Creates a plot showing the comparative SFH of galaxies in an output

The parameters passed to this script must be directories with the 
"""
import sys
from pathlib import Path
from collections import defaultdict

import numpy as np
from yt.extensions.astro_analysis.halo_analysis.api import HaloCatalog
from astropy import cosmology
from astropy import units as u

import abundance_matching

um = abundance_matching.UniverseMachine()

import betterplotlib as bpl

bpl.set_style()

import yt

yt.funcs.mylog.setLevel(50)  # ignore yt's output

import plot_utils

sfh_sentinel = Path(sys.argv[1])

# ======================================================================================
#
# Setup for precalculation of key quantities
#
# ======================================================================================
class Galaxy(object):
    def __init__(self, sphere, rank):
        self.sphere = sphere
        self.rank = rank

        # storing calculated quantities
        self._precalc_sfh = None
        self._precalc_growth = None

    def sfh(self):
        if self._precalc_sfh is None:
            self._precalc_sfh = self._create_sfh()
        return self._precalc_sfh

    def mass_growth(self):
        if self._precalc_growth is None:
            self._precalc_growth = self._create_cumulative_mass()
        return self._precalc_growth

    def _create_sfh(self):
        masses = self.sphere[("STAR", "INITIAL_MASS")].in_units("msun")
        creation_times = self.sphere[("STAR", "creation_time")]

        # have the bins start with the last time, so that the final bin is not
        # strange due to only being incomplete
        dt = 300 * yt.units.Myr
        max_bin = self.sphere.ds.current_time.to("Myr")
        bin_edges = [max_bin]
        while bin_edges[-1] > 0:
            bin_edges.append(bin_edges[-1] - dt)
        # then reverse it, so the biggest bin is at the end
        bin_edges = bin_edges[::-1]

        bin_centers = []
        sfr_values = []
        # go through each bin, seeing which stars were born there
        for left_idx in range(len(bin_edges) - 1):
            # get the age range
            min_age = bin_edges[left_idx]
            max_age = bin_edges[left_idx + 1]

            # then get the stars in that range
            age_more_idx = np.where(creation_times >= min_age)
            age_less_idx = np.where(creation_times < max_age)
            age_good_idx = np.intersect1d(age_more_idx, age_less_idx)

            # and sum their masses
            this_formed_mass = np.sum(masses[age_good_idx])

            sfr_values.append(this_formed_mass / dt)
            bin_centers.append((min_age + max_age) / 2.0)

        return yt.YTArray(bin_centers), yt.YTArray(sfr_values)

    def _create_cumulative_mass(self):
        """
        Create the cumulative stellar mass of the stars within some region.

        This will create many timesteps, then for each timestep record the mass
        that formed earlier than this time.

        :param data_obj: region of the simulation that stars will be selected from.
                         Can be something like a sphere, or even all_data()
        """
        masses = self.sphere[("STAR", "INITIAL_MASS")].in_units("msun")
        creation_times = self.sphere[("STAR", "creation_time")]

        sort_idxs = np.argsort(creation_times.to("Myr").value)

        times = creation_times[sort_idxs]
        mass_in_order = masses[sort_idxs]
        mass_cumulative = np.cumsum(mass_in_order)

        return times.in_units("Gyr"), yt.YTArray(mass_cumulative)


class Simulation(object):
    def __init__(self, ds_path):
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

        # if we are the old IC set, we have one galaxy, otherwise two
        # check what kind of particles are present
        if ("N-BODY_0", "MASS") in self.ds.derived_field_list:
            self.n_galaxies_each = 2
        else:
            self.n_galaxies_each = 1

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
        for rank, idx in enumerate(rank_idxs[: self.n_galaxies_each], start=1):
            halo = self.hc.halo_list[idx]
            radius = min(halo.quantities["virial_radius"], 30 * yt.units.kpc)
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
# make dictionary to store the resulting datasets
sims = []
for directory in sys.argv[1:]:
    directory = Path(directory)
    if directory not in plot_utils.names:
        print(f"Skipping {directory}")
        continue

    out_dir = directory / "out"
    halos_dir = directory / "halos"

    last_out = sorted(
        [
            file
            for file in out_dir.iterdir()
            if file.is_file()
            and str(file.name).endswith(".art")
            and str(file.name).startswith("continuous_a")
        ]
    )[-1]

    sims.append(Simulation(last_out))

# ======================================================================================
#
# Set up cosmology
#
# ======================================================================================
temp_ds = sims[0].ds
# Initialize the cosmology object, used to put a redshift scale on the plots
h = temp_ds.artio_parameters["hubble"][0]
H_0 = h * 100 * u.km / (u.Mpc * u.second)
omega_matter = temp_ds.artio_parameters["OmegaM"][0]
cosmo = cosmology.FlatLambdaCDM(H0=H_0, Om0=omega_matter, Tcmb0=2.725)
# and functions to transfer between redshift and age
def age_to_z(age):
    try:
        return cosmology.z_at_value(cosmo.age, age)
    except cosmology.funcs.CosmologyError:  # happens at very high z
        return 1000


def z_to_age_Gyr(z):
    # this is needed since the twin axis function always works with a scalar
    return z_to_age(z).to("Gyr").value


def z_to_age(z):
    return cosmo.age(z).to("Gyr")


# ======================================================================================
#
# Make the SFH plot
#
# ======================================================================================
label_redshifts = [10, 5, 3, 2, 1, 0.5, 0.3, 0.2, 0.1]


def plot_sfh(axis_name):
    fig, ax = bpl.subplots()

    # store data about times
    max_time = 0

    # then go through the sims
    for sim in sims:
        if axis_name not in sim.axes:
            continue
        for galaxy in sim.galaxies:
            times, sfh_values = galaxy.sfh()
            # don't plot halos without few points
            if len(times) < 2:
                continue

            plot_times = times.to("Gyr").value
            plot_sfh = sfh_values.to("msun/yr").value
            dt = plot_times[1] - plot_times[0]
            if galaxy.rank == 1:
                label = sim.name
            else:
                label = None
            ax.errorbar(
                plot_times,
                plot_sfh,
                xerr=0.5 * dt,
                markersize=8,
                c=sim.color,
                label=label,
            )
            ax.plot(plot_times, plot_sfh, lw=1.0, c=sim.color)

            # figure out the max time to use for the plot limit
            if max(plot_times) > max_time:
                max_time = max(plot_times)

    # compare to Milky Way prediction
    zs, sfhs, hi_lim, lo_lim = um.get_sfh("halo", 0, 1e12)
    ages = [z_to_age(z).to("Gyr").value for z in zs]
    ax.fill_between(
        x=ages, y1=lo_lim, y2=hi_lim, alpha=0.4, lw=0, color="0.3", label="MW-like"
    )

    ax.legend(loc=2, fontsize=10)
    ax.set_yscale("log")
    ax.set_limits(0, 1.05 * max_time, 0.1, 20)
    ax.add_labels("Time [Gyr]", "SFR  [$M_\odot$/yr]")

    # then add the redshift axis. The process of selecting the labels raises
    # warnings, so we can ignore that
    # with warnings.catch_warnings():
    #     warnings.simplefilter('ignore', UserWarning)
    ax.twin_axis("x", label_redshifts, "Redshift", new_to_old_func=z_to_age_Gyr)

    fig.savefig(sfh_sentinel.parent / f"sfh_{axis_name}.pdf")


# ======================================================================================
#
# Make the cumulative growth plot
#
# ======================================================================================
def plot_cumulative_growth(axis_name):
    fig, ax = bpl.subplots()
    # store data about times
    max_time = 0

    # iterate through simulations
    for sim in sims:
        if axis_name not in sim.axes:
            continue
        for galaxy in sim.galaxies:
            times, cumulative_mass = galaxy.mass_growth()
            # handle halos with few points
            if len(times) == 1:
                times = np.concatenate([times, times])
                cumulative_mass = np.concatenate([cumulative_mass, cumulative_mass])
            elif len(times) == 0:
                times = [0, 0, 0] * yt.units.Gyr
                cumulative_mass = [0, 0, 0] * yt.units.Msun

            plot_times = times.to("Gyr").value
            cumulative_mass = cumulative_mass.to("msun").value
            if galaxy.rank == 1:
                label = f"{sim.name}: z = {1 / sim.ds.scale_factor - 1:.1f}"
            else:
                label = None
            ax.plot(plot_times, cumulative_mass, c=sim.color, label=label)

            # figure out the max time to use for the plot limit
            if max(plot_times) > max_time:
                max_time = max(plot_times)

    ax.legend(fontsize=10)
    ax.set_yscale("log")
    ax.set_limits(0, 1.05 * max_time, 2e7, 3e10)
    ax.add_labels("Time [Gyr]", "Stellar Mass  [$M_\odot$]")

    # then add the redshift axis. The process of selecting the labels raises
    # warnings, so we can ignore that
    # with warnings.catch_warnings():
    #     warnings.simplefilter('ignore', UserWarning)
    ax.twin_axis("x", label_redshifts, "Redshift", new_to_old_func=z_to_age_Gyr)

    fig.savefig(sfh_sentinel.parent / f"mass_growth_{axis_name}.pdf")


# ======================================================================================
#
# then actually do the plotting
#
# ======================================================================================
for plot_name in plot_utils.get_plot_names([sim.name for sim in sims]):
    plot_sfh(plot_name)
    plot_cumulative_growth(plot_name)

sfh_sentinel.touch()
