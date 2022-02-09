"""
plot_cimf_evolution.py

Plot the evolution of the CIMF for one simulation.

This takes two parameters:
- the location to save the plot to. The name of this plot should tell the directory
  where the sim lives
"""
import sys
from pathlib import Path

import numpy as np
from scipy import special
import yt
import betterplotlib as bpl

from utils import load_galaxies, plot_utils, run_attributes

bpl.set_style()
yt.funcs.mylog.setLevel(50)  # ignore yt's output

plot = Path(sys.argv[1])
# then get the simulation directory. This is based on reverse engineering the function
# in the makefile to get the plot name.
sim_dir = Path(
    plot.name.replace("cimf_zevolution_", "").replace("__", "/").replace(".pdf", "")
)

# ======================================================================================
#
# CIMF calculations borrowed from plot_cimf.py and simplified
#
# ======================================================================================
# Then the functions to calculate the CIMF. Here we need to do some analysis
# of the bound fraction.
def f_bound(eps_int):
    # Li et al 2019: https://ui.adsabs.harvard.edu/abs/2019MNRAS.487..364L/abstract
    # equation 17
    alpha_star = 0.48
    f_sat = 0.94
    term_a = special.erf(np.sqrt(3 * eps_int / alpha_star))
    term_b = np.sqrt(12 * eps_int / (np.pi * alpha_star))
    term_c = np.exp(-3 * eps_int / alpha_star)
    return (term_a - (term_b * term_c)) * f_sat


def get_initial_bound_fraction(galaxy):
    star_initial_mass = galaxy[("STAR", "INITIAL_MASS")].to("Msun").value
    # the variable named INITIAL_BOUND_FRACTION is not the initial_bound fraction,
    # it's actually the accumulated mass nearby through the course of accretion, in
    # code masses. This is used to calculate the formation efficiency, which is then
    # used to get the bound fraction.
    star_accumulated_mass = galaxy[("STAR", "INITIAL_BOUND_FRACTION")].to("").value
    star_accumulated_mass *= galaxy.ds.mass_unit
    star_accumulated_mass = star_accumulated_mass.to("Msun").value

    eps_int = star_initial_mass / star_accumulated_mass

    return f_bound(eps_int)


def cimf(sim):

    mass = []
    for galaxy in sim.galaxies:
        # use mass that accounts for stellar evolution
        raw_mass = galaxy[("STAR", "MASS")].to("Msun").value
        star_initial_bound = get_initial_bound_fraction(galaxy)
        tidal_bound_fraction = galaxy[("STAR", "BOUND_FRACTION")].value
        this_mass = raw_mass * star_initial_bound * tidal_bound_fraction

        mass = np.concatenate([this_mass, mass])

    # create bins with spacing of 0.16 dex
    bin_width = 0.16  # dex
    m_boundaries_log = np.arange(3 - 0.5 * bin_width, 7, 0.16)
    m_centers_log = [
        np.mean([m_boundaries_log[idx], m_boundaries_log[idx + 1]])
        for idx in range(len(m_boundaries_log) - 1)
    ]

    m_boundaries = 10 ** m_boundaries_log
    m_centers = 10 ** np.array(m_centers_log)

    # then make the histogram showing how many there are per bin
    hist, edges = np.histogram(mass, bins=m_boundaries)
    assert np.array_equiv(m_boundaries, edges)

    # We have dN, make it per dLogM
    hist = np.array(hist) / (bin_width * np.log(10))

    # Also normalize for the number of galaxies
    hist = hist / sim.n_galaxies

    return m_centers, hist


# ======================================================================================
#
# redshift evolution of a single run
#
# ======================================================================================
# This is pretty simple. We load the simulation outputs at a range of redshift and
# plot the observable mass function at that redshift.
fig, ax = bpl.subplots()
# colors = ["#E996CA", "#9E729F", "#534D74", "#082949"]
colors = [
    run_attributes.h(*hsv)
    for hsv in [
        (0.90, 0.20, 0.90),
        (0.80, 0.35, 0.75),
        (0.70, 0.50, 0.60),
        (0.60, 0.80, 0.10),
    ]
]
for z, c in zip([7, 6, 5, 4], colors):
    sim = load_galaxies.get_simulations_same_scale([sim_dir], z)
    if len(sim) == 0:  # no sim at this redshift found
        continue
    assert len(sim) == 1
    plot_masses, dn_dlogM = cimf(sim[0])
    ax.plot(plot_masses, dn_dlogM, c=c, label=f"z={z}")

# format axis
plot_utils.add_legend(ax, loc=1, fontsize=18)
ax.set_yscale("log")
ax.set_xscale("log")
ax.set_limits(1e3, 1e7, 10, 1e4)
ax.add_labels("$f_b$M [$M_\odot$]", "dN/dlogM")
fig.savefig(plot)
