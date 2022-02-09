"""
plot_age_spread.py

Creates a plot showing the comparative cluster age spreads of clusters in 
different outputs. This code is very similar to plot_cimf.py, some of it was
copied. 
"""
import sys
from pathlib import Path

import numpy as np
import yt
from tqdm import tqdm
import betterplotlib as bpl
from matplotlib import pyplot as plt

from utils import load_galaxies, plot_utils
from analysis_functions import age_spreads

bpl.set_style()
yt.funcs.mylog.setLevel(50)  # ignore yt's output

sentinel = Path(sys.argv[1]).resolve()
plot_dir = sentinel.parent

# ======================================================================================
#
# load the simulations
#
# ======================================================================================
sims_last = load_galaxies.get_simulations_last(sys.argv[2:])
sims_common, common_scale = load_galaxies.get_simulations_common(sys.argv[2:])


def time_cumulative_hist(sim, time_func, mask_name):
    """
    Make the cluster initial mass function.

    :param time_func: the function to be used to calculate the specific
                      implementation of the age spread
    :returns: Two lists. The first is the list of age spreads, in sorted order.
              The second is the fraction of clusters with age spreads less or
              equal to this value.
    """
    # check mask_name
    if mask_name not in ["lo", "hi"]:
        raise RuntimeError("Bad mask name")

    # check if precalculated
    quantity_name = f"{time_func.__name__}_{mask_name}"
    if quantity_name not in sim.precalculated:
        spreads = []
        for galaxy in sim.galaxies:
            # make the mask. Note that the cut that the cluster must be done forming
            # is automatically included in the galaxy object
            cluster_masses = galaxy[("STAR", "initial_mass")]
            cluster_cut = 1e5 * yt.units.Msun
            if mask_name == "lo":
                mask = cluster_masses <= cluster_cut
            else:
                mask = cluster_masses > cluster_cut

            # then calculate the age spreads, if we have objects in this category
            if np.sum(mask) > 0:
                this_spreads = time_func(galaxy).to("Myr").value[mask]
            else:
                this_spreads = []

            spreads = np.concatenate([this_spreads, spreads])

        # sort the age spreads, so we can use them for the cumulative distribution
        spreads = np.sort(spreads)

        if len(spreads) > 0:
            ranks = np.linspace(1 / len(spreads), 1, len(spreads))
        else:
            ranks = []

        sim.precalculated[quantity_name] = spreads, ranks

    return sim.precalculated[quantity_name]


# ======================================================================================
#
# plotting
#
# ======================================================================================
def plot_age_growth_base(
    ax,
    axis_name,
    age_quantity,
    sim_share_type,
    mass_side,
    legend,
    both,
    ls=None,
    label_mass=True,
):
    if sim_share_type == "last":
        sims = sims_last
    else:
        sims = sims_common

    if mass_side not in ["lo", "hi"]:
        raise ValueError("bad mass_side")
    high = mass_side == "hi"

    funcs = {
        "Duration": age_spreads.duration,
        "Average Age": age_spreads.ave_time,
        "Age Spread": age_spreads.age_spread,
    }
    func = funcs[age_quantity]

    # add the labels here
    ax.add_labels(age_quantity + " [Myr]", "Cumulative Fraction")
    if high and label_mass:
        ax.easy_add_text("M > $10^5 M_\odot$", "upper left")
    elif label_mass:
        ax.easy_add_text("M < $10^5 M_\odot$", "upper left")

    for sim in sims:
        if axis_name not in sim.axes:
            continue

        # make the plot legend
        if legend:
            label = sim.names[axis_name]
            z = 1 / sim.ds.scale_factor - 1
            if sim_share_type == "last" and not 1.49 < z < 1.51:
                label += f": z = {z:.1f}"
        else:
            label = None

        if ls is None:
            ls = sim.ls

        ages, fractions = time_cumulative_hist(sim, func, mass_side)
        ax.plot(ages, fractions, c=sim.color, ls=ls, label=label)

    # add limits appropriately. This dictionary holds the maximum x value for
    # low and high mass, respectively. If we're plotting both, use the high mass for
    # both panels.
    limits = {"Duration": (5, 15), "Average Age": (3, 6), "Age Spread": (2, 4)}
    if both or high:
        idx = 1
    else:
        idx = 0
    ax.set_limits(0, limits[age_quantity][idx], 0, 1)

    if legend:
        plot_utils.add_legend(ax, loc=4, fontsize=16)
        # if there is a common redshift, annotate it
        if sim_share_type == "common":
            ax.easy_add_text(f"z = {1/common_scale - 1:.1f}", "upper right")


def plot_age_growth(axis_name, age_quantity, sim_share_type, which_mass):
    """
    Plot the cluster age spread quantities

    :param axis_name: The name of the sim group to put on this plot
    :param sim_share_type: Either "last" or "common", depending on which
                           dictionaries are passed in. This will determine if
                           the redshift is labeled as common to all, or
                           individually per simulation
    """
    if sim_share_type not in ["last", "common"]:
        raise ValueError("bad plot_name_suffix")

    if which_mass not in ["lo", "hi", "both_split", "both_share"]:
        raise ValueError("bad which_mass")

    if which_mass == "both_split":
        fig, axs = bpl.subplots(figsize=[20, 7], ncols=2)
        plot_age_growth_base(
            axs[0], axis_name, age_quantity, sim_share_type, "lo", True, True
        )
        plot_age_growth_base(
            axs[1], axis_name, age_quantity, sim_share_type, "hi", False, True
        )
    elif which_mass == "both_share":
        fig, ax = bpl.subplots()
        # add dummy legend for line styles. I do this here since the
        # plot_age_growth_base function adds the legend.
        ax.plot([0, 0], [2, 2], c="0.5", ls="--", label="M < $10^5 M_\odot$")
        ax.plot([0, 0], [2, 2], c="0.5", ls="-", label="M > $10^5 M_\odot$")

        plot_age_growth_base(
            ax,
            axis_name,
            age_quantity,
            sim_share_type,
            "lo",
            legend=False,
            both=True,
            ls="--",
            label_mass=False,
        )
        plot_age_growth_base(
            ax,
            axis_name,
            age_quantity,
            sim_share_type,
            "hi",
            legend=True,
            both=True,
            ls="-",
            label_mass=False,
        )

    else:
        fig, ax = bpl.subplots()
        plot_age_growth_base(
            ax, axis_name, age_quantity, sim_share_type, which_mass, True, False
        )

    plot_name = f"age_{age_quantity.lower().replace(' ', '_')}"
    plot_name += f"_{axis_name}_{sim_share_type}_{which_mass}.pdf"
    fig.savefig(plot_dir / plot_name)
    # then remove figure for memory purposes
    plt.close(fig)


def plot_age_mass(axis_name, age_quantity, sim_share_type):
    if sim_share_type == "last":
        sims = sims_last
    elif sim_share_type == "common":
        sims = sims_common
    else:
        raise ValueError("bad sim_share_type")

    funcs = {
        "Duration": age_spreads.duration,
        "Average Age": age_spreads.ave_time,
        "Age Spread": age_spreads.age_spread,
    }
    func = funcs[age_quantity]

    fig, ax = bpl.subplots()

    for sim in sims:
        if axis_name not in sim.axes:
            continue

        spreads = sim.func_all_galaxies(lambda g: func(g).to("Myr").value)
        masses = sim.func_all_galaxies(
            lambda g: g[("STAR", "INITIAL_MASS")].to("Msun").value
        )

        plot_utils.shaded_region(
            ax,
            masses,
            spreads,
            sim.color,
            p_lo=25,
            p_hi=75,
            log_x=True,
            label=plot_utils.plot_label(sim, sim_share_type, axis_name),
        )

    # format axis
    ax.add_labels("Initial Particle Mass [$M_\odot$]", age_quantity + " [Myr]")
    ax.set_xscale("log")
    ax.legend(loc=2)
    plot_utils.nice_log_axis(ax, "x")

    # add limits appropriately. This dictionary holds the maximum x value for
    # low and high mass, respectively. If we're plotting both, use the high mass for
    # both panels.
    limits = {"Duration": 15, "Average Age": 6, "Age Spread": 4}
    ax.set_limits(1e2, 1e7, 0, limits[age_quantity])

    plot_name = f"age_mass_{age_quantity.lower().replace(' ', '_')}"
    plot_name += f"_{axis_name}_{sim_share_type}.pdf"
    fig.savefig(plot_dir / plot_name)
    # then remove figure for memory purposes
    plt.close(fig)


# ======================================================================================
#
# actually make the plots
#
# ======================================================================================
for plot_name in tqdm(load_galaxies.get_plot_names(sims_last)):
    for share_type in ["common", "last"]:
        for age_type in ["Duration", "Average Age", "Age Spread"]:
            plot_age_mass(plot_name, age_type, share_type)
            for which_mass in ["both_split", "both_share"]:  # ["lo", "hi", ]
                plot_age_growth(plot_name, age_type, share_type, which_mass)

sentinel.touch()
