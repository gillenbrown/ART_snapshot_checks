"""
plot_bound_fraction.py

Plot the cluster bound fraction and some related quantities
"""
from pathlib import Path

import sys
import numpy as np
from scipy import special
from matplotlib import pyplot as plt

from utils import load_galaxies, plot_utils

import betterplotlib as bpl

bpl.set_style()

sentinel = Path(sys.argv[1])
plot_dir = sentinel.parent

# ======================================================================================
#
# load the simulations
#
# ======================================================================================
sims_last = load_galaxies.get_simulations_last(sys.argv[2:])
common_z = 4.0
sims_common = load_galaxies.get_simulations_same_scale(sys.argv[2:], common_z)

# ======================================================================================
#
# analysis functions
#
# ======================================================================================
def f_bound(eps_int):
    # Li et al 2019: https://ui.adsabs.harvard.edu/abs/2019MNRAS.487..364L/abstract
    # equation 17
    alpha_star = 0.48
    f_sat = 0.94
    term_a = special.erf(np.sqrt(3 * eps_int / alpha_star))
    term_b = np.sqrt(12 * eps_int / (np.pi * alpha_star))
    term_c = np.exp(-3 * eps_int / alpha_star)
    return (term_a - (term_b * term_c)) * f_sat


def get_eps_int(galaxy):
    star_initial_mass = get_initial_mass_msun(galaxy)
    # the variable named INITIAL_BOUND_FRACTION is not the initial_bound fraction,
    # it's actually the accumulated mass nearby through the course of accretion, in
    # code masses. This is used to calculate the formation efficiency, which is then
    # used to get the bound fraction.
    star_accumulated_mass = galaxy[("STAR", "INITIAL_BOUND_FRACTION")].to("").value
    star_accumulated_mass *= galaxy.ds.mass_unit
    star_accumulated_mass = star_accumulated_mass.to("Msun").value

    return star_initial_mass / star_accumulated_mass


def get_initial_mass_msun(galaxy):
    return galaxy[("STAR", "INITIAL_MASS")].to("Msun").value


def get_initial_bound_fraction(galaxy):
    return f_bound(get_eps_int(galaxy))


def get_all_galaxies(sim, func):
    return np.concatenate([func(galaxy) for galaxy in sim.galaxies])


# ======================================================================================
#
# plotting functions
#
# ======================================================================================
def plot_bound_fraction(axis_name, sim_share_type):
    """
    Plot the cluster initial bound fraction as a function of mass.

    :param axis_name: The name of the axis, so we can select which simulations go
                      on this plot.
    :param sim_share_type: Either "last" or "common". This will determine if
                           the redshift is labeled as common to all, or
                           individually per simulation

    """
    # validate options given
    if sim_share_type not in ["last", "common"]:
        raise ValueError("bad sim_share_type")

    # choose which simulations we're using
    if sim_share_type == "common":
        sims = sims_common
    else:
        sims = sims_last

    fig, ax = bpl.subplots(figsize=[9, 7])

    for sim in sims:
        if axis_name not in sim.axes:
            continue

        mass = get_all_galaxies(sim, get_initial_mass_msun)
        f_b = get_all_galaxies(sim, get_initial_bound_fraction)

        plot_utils.shaded_region(
            ax,
            mass,
            f_b,
            sim.color,
            p_lo=25,
            p_hi=75,
            log_x=True,
            label=plot_utils.plot_label(sim, sim_share_type, axis_name),
        )

    # format axes
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_limits(1e2, 1e7, 0.001, 1)
    ax.add_labels("Initial Particle Mass [$M_\odot$]", "Initial Bound Fraction")
    ax.legend(loc=4)
    plot_utils.nice_log_axis(ax, "y")

    # if there is a common redshift, annotate it
    if sim_share_type == "common":
        ax.easy_add_text(f"z = {common_z:.1f}", "upper left")

    fig.savefig(sentinel.parent / f"bound_fraction_{axis_name}_{sim_share_type}.png")
    plt.close(fig)  # to save memory


def plot_eps_int(axis_name, sim_share_type):
    """
    Plot kde histograms of the cluster integrated star formation efficiency

    :param axis_name: The name of the axis, so we can select which simulations go
                      on this plot.
    :param sim_share_type: Either "last" or "common". This will determine if
                           the redshift is labeled as common to all, or
                           individually per simulation

    """
    # validate options given
    if sim_share_type not in ["last", "common"]:
        raise ValueError("bad sim_share_type")

    # choose which simulations we're using
    if sim_share_type == "common":
        sims = sims_common
    else:
        sims = sims_last

    fig, ax = bpl.subplots(figsize=[9, 7])
    x_values = np.logspace(-3, 0, 301)
    # I need to adjust the limits of the plot depending on the SFE of the sims present
    # set a default value which will be adjusted
    x_min = 1e-2
    for sim in sims:
        if axis_name not in sim.axes:
            continue

        mass = get_all_galaxies(sim, get_initial_mass_msun)
        eps_int = get_all_galaxies(sim, get_eps_int)

        plot_y = plot_utils.kde(x_values, eps_int, width=0.05, weights=mass, log=True)
        ax.plot(
            x_values,
            plot_y,
            c=sim.color,
            label=plot_utils.plot_label(sim, sim_share_type, axis_name),
        )

        # if there is something other than sfe100, update the limit
        if "sfe010" in str(sim.run_dir) or "sfe001" in str(sim.run_dir):
            x_min = 1e-3

    # format axes. limits depend on which runs are shown
    ax.set_xscale("log")
    ax.set_limits(x_min, 1, 0)
    ax.legend(loc=2, fontsize=16, frameon=False)
    plot_utils.nice_log_axis(ax, "x")
    ax.add_labels(
        "Integrated Star Formation Efficiency $\epsilon_{int}$",
        "Mass Weighted KDE Density",
    )

    # if there is a common redshift, annotate it
    if sim_share_type == "common":
        ax.easy_add_text(f"z = {common_z:.1f}", "upper right")

    fig.savefig(sentinel.parent / f"eps_int_{axis_name}_{sim_share_type}.pdf")
    plt.close(fig)  # to save memory


# ======================================================================================
#
# make the plots
#
# ======================================================================================
for plot_name in load_galaxies.get_plot_names(sims_last):
    for share_type in ["last", "common"]:
        plot_bound_fraction(plot_name, share_type)
        plot_eps_int(plot_name, share_type)

# touch the sentinel once done
sentinel.touch()
