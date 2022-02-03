"""
plot_age_metallicity.py

Plot the age-metallicity relation, as well as the hypernova fraction for runs
with nonzero hypernovae
"""
from pathlib import Path

import sys
import numpy as np
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
sims = load_galaxies.get_simulations_last(sys.argv[2:])

# functions to get the simulation name and hypernova fraction
def get_sim_name(sim_loc):
    sim_path = Path(sim_loc)
    if sim_path.name == "run":
        sim_path = sim_path.parent

    # now I should have the directory
    return sim_path.name


def get_f_hn0(sim_name):
    # find it in the name.
    for chunk in sim_name.split("_"):
        if chunk.startswith("hn"):
            return float(chunk[2:]) / 100
    # if it's not specified, use 0
    return 0.0


# ======================================================================================
#
# plotting functions, including the hypernova fraction at a given metallicity
#
# ======================================================================================
f_hn_scale = -0.001


def f_hn_at_z(metallicity, f_hn0):
    return f_hn0 * np.exp(metallicity / f_hn_scale)
    # return max(f_hn, 0.001)
    # Ignore the minimum value, I'll just only put ticks in the good range.


def z_at_f_hn(hn_fraction, f_hn0):
    return f_hn_scale * np.log(hn_fraction / f_hn0)


def age_metallicity_plot_base(ax, sim, add_hn_axis=False, label=None):
    sim_name = get_sim_name(sim.run_dir)
    f_hn0 = get_f_hn0(sim_name)
    for gal in sim.galaxies:
        metallicity = (
            gal[("STAR", "METALLICITY_SNII")] + gal[("STAR", "METALLICITY_SNIa")]
        )
        creation_time = gal[("STAR", "creation_time")].to("Gyr").value
        plot_utils.shaded_region(ax, creation_time, metallicity, sim.color, label=label)

    # format axes
    ax.set_yscale("log")
    ax.set_limits(0, 5, 1e-4, 0.03)
    ax.add_labels("Time of Cluster Formation [Gyr]", "Metallicity")

    def z_at_f_hn_wrapper(z):
        return z_at_f_hn(z, f_hn0)

    # if we have a nonzero hn fraction, add the second axis
    if f_hn0 > 0 and add_hn_axis:
        ax.twin_axis(
            "y",
            [0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5],
            "$f_{HN}$",
            new_to_old_func=z_at_f_hn_wrapper,
        )


def age_metallicity_plot(sim):
    # get basic info
    sim_name = get_sim_name(sim.run_dir)

    fig, ax = bpl.subplots()
    age_metallicity_plot_base(ax, sim)

    # then save the plot
    fig.savefig(plot_dir / f"cluster_age_metallicity_{sim_name}.pdf")
    # then remove figure for memory purposes
    plt.close(fig)


def mass_metallicity_plot(sim):
    # get basic info
    sim_name = get_sim_name(sim.run_dir)

    n_galaxies = len(sim.galaxies)
    fig, axs = bpl.subplots(ncols=n_galaxies, figsize=[10 * n_galaxies, 7])
    if n_galaxies == 1:
        axs = [axs]  # so we can iterate
    for gal, ax in zip(sim.galaxies, axs):
        metallicity = (
            gal[("STAR", "METALLICITY_SNII")] + gal[("STAR", "METALLICITY_SNIa")]
        )
        mass = gal[("STAR", "INITIAL_MASS")].to("Msun").value
        ax.shaded_density(
            metallicity,
            mass,
            log_hist=True,
            log_xy=True,
            smoothing=0,
            bin_size=(0.2, 0.5),
        )
        plot_utils.shaded_region(
            ax, metallicity, mass, sim.color, p_lo=50, p_hi=50, log_x=True, dx=0.2
        )

        # format axes
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_limits(1e-4, 0.03, 1e2, 1e7)
        ax.add_labels("Metallicity", "Cluster Initial Mass [$M_\odot$]")

    # then save the plot
    fig.savefig(plot_dir / f"cluster_metallicity_mass_{sim_name}.pdf")
    # then remove figure for memory purposes
    plt.close(fig)


def mass_age_plot(sim):
    # get basic info
    sim_name = get_sim_name(sim.run_dir)

    n_galaxies = len(sim.galaxies)
    fig, axs = bpl.subplots(ncols=n_galaxies, figsize=[10 * n_galaxies, 7])
    if n_galaxies == 1:
        axs = [axs]  # so we can iterate
    for gal, ax in zip(sim.galaxies, axs):
        creation_time = gal[("STAR", "creation_time")].to("Gyr").value
        mass = gal[("STAR", "INITIAL_MASS")].to("Msun").value
        ax.shaded_density(
            creation_time,
            mass,
            log_hist=True,
            log_xy=(False, True),
            smoothing=0,
            bin_size=(0.2, 0.5),
        )
        plot_utils.shaded_region(
            ax, creation_time, mass, sim.color, p_lo=50, p_hi=50, log_x=False, dx=0.2
        )

        # format axes
        ax.set_yscale("log")
        ax.set_limits(0, 5, 1e2, 1e7)
        ax.add_labels(
            "Time of Cluster Formation [Gyr]", "Cluster Initial Mass [$M_\odot$]"
        )

    # then save the plot
    fig.savefig(plot_dir / f"cluster_age_mass_{sim_name}.pdf")
    # then remove figure for memory purposes
    plt.close(fig)


# ======================================================================================
#
# make the plots
#
# ======================================================================================
for sim in sims[:3]:
    age_metallicity_plot(sim)
    mass_metallicity_plot(sim)
    mass_age_plot(sim)

# ======================================================================================
#
# then make one plot with multiple sims
#
# ======================================================================================
# fig, ax = bpl.subplots()
# plot_sims = [
#     "discrete_hn50_virial10_entropy_fboost1",
#     "discrete_hn00_virial10_entropy_fboost1",
# ]
# for sim in sims:
#     if get_sim_name(sim.run_dir) in plot_sims:
#         age_metallicity_plot_base(
#             ax, sim, add_hn_axis=False, label=sim.names["old_ic_sn_feedback"]
#         )
#
# # then save the plot
# fig.savefig(plot_dir / f"cluster_age_metallicity_comparison.pdf")
# # then remove figure for memory purposes
# plt.close(fig)

# touch the sentinel once done
sentinel.touch()
