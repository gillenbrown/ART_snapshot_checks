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


def f_hn_at_z(metallicity):
    return f_hn0 * np.exp(metallicity / f_hn_scale)
    # return max(f_hn, 0.001)
    # Ignore the minimum value, I'll just only put ticks in the good range.


def z_at_f_hn(hn_fraction):
    return f_hn_scale * np.log(hn_fraction / f_hn0)


# ======================================================================================
#
# make the plots
#
# ======================================================================================
for sim in sims:
    # get basic info
    sim_name = get_sim_name(sim.run_dir)
    f_hn0 = get_f_hn0(sim_name)

    gal = sim.galaxies[0]

    metallicity = gal[("STAR", "METALLICITY_SNII")] + gal[("STAR", "METALLICITY_SNIa")]
    creation_time = gal[("STAR", "creation_time")].to("Gyr").value

    fig, ax = bpl.subplots()
    plot_utils.shaded_region(ax, creation_time, metallicity, sim.color)
    ax.set_yscale("log")
    ax.set_limits(0, 5, 1e-4, 0.03)
    ax.add_labels("Time of Cluster Formation [Gyr]", "Metallicity")

    # if we have a nonzero hn fraction, add the second axis
    if f_hn0 > 0:
        ax.twin_axis(
            "y",
            [0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5],
            "$f_{HN}$",
            new_to_old_func=z_at_f_hn,
        )

    # then save the plot
    fig.savefig(plot_dir / f"age_metallicity_{sim_name}.pdf")
    # then remove figure for memory purposes
    plt.close(fig)

# touch the sentinel once done
sentinel.touch()
