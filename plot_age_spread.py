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
common_redshift = 4
sims_common = load_galaxies.get_simulations_same_scale(sys.argv[2:], common_redshift)


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


def find_median_of_cumulative(xs, cumulative_fraction):
    for idx in range(len(xs)):
        if cumulative_fraction[idx] > 0.5:
            # find the exact place it crosses, assuming a straight line between the
            # point before median cross and after cross
            x1 = xs[idx - 1]
            x2 = xs[idx]
            y1 = cumulative_fraction[idx - 1]
            y2 = cumulative_fraction[idx]
            # get the slope
            slope = (y2 - y1) / (x2 - x1)
            # then find the value at which y = 0.5
            return x1 + (0.5 - y1) / slope


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
    plot_median=True,
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

    # I'll calculate the median of the latest run, the add that to the figure
    last_median = -1
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

        # store the median if it's the last one
        if plot_median:
            this_median = find_median_of_cumulative(ages, fractions)
            if this_median > last_median:
                last_median = this_median

    # then plot this median
    if plot_median:
        ax.plot(
            [0, last_median, last_median],
            [0.5, 0.5, 0],
            ls=":",
            lw=2,
            c=bpl.almost_black,
            zorder=0,
        )

    # add limits appropriately. This dictionary holds the maximum x value for
    # low and high mass, respectively. If we're plotting both, use the high mass for
    # both panels.
    limits = {"Duration": (5, 15), "Average Age": (3, 6), "Age Spread": (2, 15)}
    if both or high:
        idx = 1
    else:
        idx = 0
    ax.set_limits(0, limits[age_quantity][idx], 0, 1)

    if legend:
        # if there is a common redshift, annotate it
        if sim_share_type == "common":
            title = f"z = {common_redshift:.1f}"
        else:
            title = None
        plot_utils.add_legend(ax, loc=4, fontsize=16, title=title)


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
        fig, axs = bpl.subplots(figsize=[18, 7], ncols=2)
        plot_age_growth_base(
            axs[0],
            axis_name,
            age_quantity,
            sim_share_type,
            "lo",
            True,
            True,
        )
        plot_age_growth_base(
            axs[1],
            axis_name,
            age_quantity,
            sim_share_type,
            "hi",
            False,
            True,
        )
    elif which_mass == "both_share":
        fig, ax = bpl.subplots(figsize=[9, 7])
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
            plot_median=False,
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
            plot_median=True,
        )

    else:
        fig, ax = bpl.subplots(figsize=[9, 7])
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
    limits = {"Duration": 15, "Average Age": 10, "Age Spread": 10}
    ax.set_limits(1e2, 1e7, 0, limits[age_quantity])

    plot_name = f"age_mass_{age_quantity.lower().replace(' ', '_')}"
    plot_name += f"_{axis_name}_{sim_share_type}.pdf"
    fig.savefig(plot_dir / plot_name)
    # then remove figure for memory purposes
    plt.close(fig)


def plot_spread_vs_duration(sim):
    fig, ax = bpl.subplots(figsize=[8, 8])
    ax.equal_scale()

    # plot a scatterplot of duration vs age spread
    duration = sim.func_all_galaxies(lambda g: age_spreads.duration(g).to("Myr").value)
    spread = sim.func_all_galaxies(lambda g: age_spreads.age_spread(g).to("Myr").value)
    ax.scatter(duration, spread, s=0.1, marker=",", alpha=1, edgecolor="")

    # plot guiding lines for different accretion histories/
    # This uses equation 12 from Li et al paper 2
    def t_spread_power_law(t_dur, alpha):
        return t_dur * (2 * alpha + 1) / (alpha + 1) ** 2

    toy_duration = np.arange(0, 15, 0.01)
    ax.plot(
        toy_duration,
        t_spread_power_law(toy_duration, 0),
        c=bpl.almost_black,
        ls="--",
        label="$\dot{M}=const$",
    )
    ax.plot(
        toy_duration,
        t_spread_power_law(toy_duration, 1),
        c=bpl.almost_black,
        ls=":",
        label="$\dot{M}\propto t$",
    )
    ax.set_limits(0, 15, 0, 15)
    ax.add_labels("Duration [Myr]", "Age Spread [Myr]")

    ax.legend(loc=2)

    plot_name = f"age_comparison_{plot_utils.get_sim_dirname(sim.run_dir)}.png"
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
for sim in tqdm(sims_last):
    plot_spread_vs_duration(sim)

sentinel.touch()
