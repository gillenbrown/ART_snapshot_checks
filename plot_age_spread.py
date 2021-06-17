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
import betterplotlib as bpl

from utils import load_galaxies, plot_utils

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


# ======================================================================================
#
# functions to calculate key quantities
#
# ======================================================================================
def time_units(ds, array):
    return ds._handle.tphys_from_tcode_array(array) * yt.units.year


def duration(region, mask):
    end_time = time_units(region.ds, region[("STAR", "TERMINATION_TIME")][mask])
    return end_time - region[("STAR", "creation_time")][mask]


def ave_time(region, mask):
    art_units_ave_age = region[("STAR", "AVERAGE_AGE")][mask]
    art_units_birth = region[("STAR", "BIRTH_TIME")][mask]

    ave_time = time_units(region.ds, art_units_birth + art_units_ave_age)
    return ave_time - region[("STAR", "creation_time")][mask]


def age_spread(region, mask):
    initial_mass = region[("STAR", "initial_mass")][mask]
    age_spread = region[("STAR", "AGE_SPREAD")][mask] * region.ds.arr(1, "code_mass**2")
    birth_time = region[("STAR", "BIRTH_TIME")][mask]
    creation_time = region[("STAR", "creation_time")][mask]

    time = time_units(region.ds, (initial_mass ** 2 / age_spread) + birth_time)

    return time - creation_time


def time_cumulative_hist(galaxy, time_func, mask_name):
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
    if quantity_name not in galaxy.precalculated:
        sphere = galaxy.sphere

        # make the mask. First is stars that are done forming
        mask_done_forming = sphere[("STAR", "age")] > 15 * yt.units.Myr
        cluster_masses = sphere[("STAR", "initial_mass")]
        cluster_cut = 1e5 * yt.units.Msun
        if mask_name == "lo":
            mask = mask_done_forming & (cluster_masses <= cluster_cut)
        else:
            mask = mask_done_forming & (cluster_masses > cluster_cut)

        # then calculate the age spreads, if we have objects in this category
        if np.sum(mask) > 0:
            spreads = np.sort(time_func(sphere, mask).to("Myr").value)
            ranks = np.linspace(1 / len(spreads), 1, len(spreads))
        else:
            spreads, ranks = [], []
        galaxy.precalculated[quantity_name] = spreads, ranks

    return galaxy.precalculated[quantity_name]


# ======================================================================================
#
# plotting
#
# ======================================================================================
def plot_age_growth(axis_name, sim_share_type):
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

    if sim_share_type == "last":
        sims = sims_last
    else:
        sims = sims_common

    fig, axs = bpl.subplots(figsize=[13, 17], ncols=2, nrows=3)

    funcs = [duration, ave_time, age_spread]
    names = ["Duration", "Average Age", "Age Spread"]

    for ax_row, func, name in zip(axs, funcs, names):
        # add the labels here
        ax_row[0].add_labels(
            name + " [Myr]", "Cumulative Fraction", "M < $10^5 M_\odot$"
        )
        ax_row[1].add_labels(
            name + " [Myr]", "Cumulative Fraction", "M > $10^5 M_\odot$"
        )
        ax_row[0].set_limits(0, 5, 0, 1)
        ax_row[1].set_limits(0, 8, 0, 1)
        for sim in sims:
            if axis_name not in sim.axes:
                continue
            for galaxy in sim.galaxies:
                # make the label only for the biggest halo
                if galaxy.rank == 1:
                    # and include the redshift if it's different for each sim
                    if sim_share_type == "last":
                        label = f"{sim.name}: z = {1/sim.ds.scale_factor - 1:.1f}"
                    else:
                        label = sim.name
                else:
                    label = None

                ages_lo, fractions_lo = time_cumulative_hist(galaxy, func, "lo")
                ages_hi, fractions_hi = time_cumulative_hist(galaxy, func, "hi")
                ax_row[1].plot(ages_hi, fractions_hi, c=sim.color, lw=2, label=label)
                ax_row[0].plot(ages_lo, fractions_lo, c=sim.color, lW=2, label=label)

    # put a legend in the top left panel
    axs[0][0].legend(loc=4, fontsize=14)

    # if there is a common redshift, annotate it
    if sim_share_type == "common":
        axs[0][0].easy_add_text(f"z = {1/common_scale - 1:.1f}", "upper left")

    plot_name = f"age_spread_{axis_name.replace('_', '')}_{sim_share_type}.pdf"
    fig.savefig(plot_dir / plot_name)


# then actually call this function to build the plots
for plot_name in plot_utils.get_plot_names([sim.name for sim in sims_last]):
    for share_type in ["common", "last"]:
        plot_age_growth(plot_name, share_type)

sentinel.touch()
