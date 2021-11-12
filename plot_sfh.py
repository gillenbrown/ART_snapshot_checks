"""
plots_sfh.py

Creates a plot showing the comparative SFH of galaxies in an output

The parameters passed to this script must be directories with the 
"""
import sys
from pathlib import Path

import numpy as np
import yt
from astropy import cosmology
from astropy import units as u
from scipy import integrate
import betterplotlib as bpl
import abundance_matching

from utils import load_galaxies, plot_utils

bpl.set_style()
yt.funcs.mylog.setLevel(50)  # ignore yt's output
um = abundance_matching.UniverseMachine()

sfh_sentinel = Path(sys.argv[1])

# ======================================================================================
#
# load the simulations
#
# ======================================================================================
sims = load_galaxies.get_simulations_last(sys.argv[2:])


# ======================================================================================
#
# functions to calculate key quantities
#
# ======================================================================================
def sfh(galaxy):
    if "sfh" not in galaxy.precalculated:
        masses = galaxy.prop_all_clusters(("STAR", "INITIAL_MASS")).in_units("msun")
        creation_times = galaxy.prop_all_clusters(("STAR", "creation_time"))

        # have the bins start with the last time, so that the final bin is not
        # strange due to only being incomplete
        dt = 300 * yt.units.Myr
        max_bin = galaxy.ds.current_time.to("Myr")
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

        galaxy.precalculated["sfh"] = yt.YTArray(bin_centers), yt.YTArray(sfr_values)

    return galaxy.precalculated["sfh"]


def cumulative_growth(galaxy):
    """
    Create the cumulative stellar mass of the stars within some region.

    This will create many timesteps, then for each timestep record the mass
    that formed earlier than this time.

    :param galaxy: galaxy object
    """
    if "mass_growth" not in galaxy.precalculated:
        masses = galaxy.prop_all_clusters(("STAR", "INITIAL_MASS")).in_units("msun")
        creation_times = galaxy.prop_all_clusters(("STAR", "creation_time"))

        sort_idxs = np.argsort(creation_times.to("Myr").value)

        times = creation_times[sort_idxs]
        mass_cumulative = np.cumsum(masses[sort_idxs])

        # add the current time to this, duplicating the last mass entry. This is to make
        # the plot go all the way to the current simulation epoch, even if no stars
        # have formed lately. I have to do this mess with units since I'm using an
        # intermediate numpy array
        times = (
            np.concatenate(
                [
                    times.in_units("Gyr").value,
                    [galaxy.ds.current_time.in_units("Gyr").value],
                ]
            )
            * yt.units.Gyr
        )
        mass_cumulative = (
            np.concatenate(
                [
                    mass_cumulative.in_units("Msun").value,
                    [mass_cumulative[-1].in_units("Msun").value],
                ]
            )
            * yt.units.Msun
        )

        galaxy.precalculated["mass_growth"] = times.in_units("Gyr"), yt.YTArray(
            mass_cumulative
        )

    return galaxy.precalculated["mass_growth"]


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
            times, sfh_values = sfh(galaxy)
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
            ax.scatter(
                plot_times,
                plot_sfh,
                s=100,
                c=sim.color,
                alpha=1,
                label=label,
                zorder=10,
            )
            ax.plot(plot_times, plot_sfh, lw=1.0, c=sim.color, zorder=9)

            # figure out the max time to use for the plot limit
            if max(plot_times) > max_time:
                max_time = max(plot_times)

    # compare to Milky Way prediction
    zs, sfhs, hi_lim, lo_lim = um.get_sfh("halo", 0, 1e12)
    ages = [z_to_age(z).to("Gyr").value for z in zs]
    ax.fill_between(
        x=ages,
        y1=lo_lim,
        y2=hi_lim,
        alpha=0.4,
        lw=0,
        color="0.3",
        zorder=0,
        label="MW-like (Universe Machine)",
    )

    plot_utils.add_legend(ax, loc=2, frameon=False, fontsize=10)
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
    for sim in sorted(sims, key=lambda x: x.name):
        if axis_name not in sim.axes:
            continue

        if axis_name == "old_ic_code" and sim.name == "Discrete $\\alpha<10$":
            name = "ART 2.0 Adiabatic"
        elif axis_name == "old_ic_code" and sim.name == "ART 2.1 Entropy $f_{boost}=5$":
            name = "ART 2.1 Entropy"
        else:
            name = sim.name
        for galaxy in sim.galaxies:
            times, cumulative_mass = cumulative_growth(galaxy)
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
                label = f"{name}: z = {1 / sim.ds.scale_factor - 1:.1f}"
            else:
                label = None
            ax.plot(plot_times, cumulative_mass, c=sim.color, label=label)

            # figure out the max time to use for the plot limit
            if max(plot_times) > max_time:
                max_time = max(plot_times)

    # compare to Milky Way prediction
    zs, sfhs, hi_lim, lo_lim = um.get_sfh("halo", 0, 1e12)
    um_ages_gyr = [z_to_age_Gyr(z) for z in zs]
    um_ages_yr = [a * 1e9 for a in um_ages_gyr]
    # integrate the cumulative mass as a function of time
    um_mass = integrate.cumulative_trapezoid(y=sfhs, x=um_ages_yr, initial=0)
    ax.plot(
        um_ages_gyr,
        um_mass,
        c=bpl.almost_black,
        ls="--",
        zorder=0,
        label="MW-like (Universe Machine)",
    )

    plot_utils.add_legend(ax, loc=4, fontsize=10, frameon=False)
    ax.set_yscale("log")
    ax.set_limits(0, 1.05 * max_time, 2e7, 1e10)
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
for plot_name in load_galaxies.get_plot_names([sim.name for sim in sims]):
    plot_sfh(plot_name)
    plot_cumulative_growth(plot_name)

sfh_sentinel.touch()
