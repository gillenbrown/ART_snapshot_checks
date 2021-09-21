"""
plot_cimf.py

Creates a plot showing the comparative cluster initial mass function of 
galaxies in an output
"""
import sys
from pathlib import Path

import numpy as np
from scipy import special
import yt
import betterplotlib as bpl

from utils import load_galaxies

bpl.set_style()
yt.funcs.mylog.setLevel(50)  # ignore yt's output

cimf_sentinel = Path(sys.argv[1])

# ======================================================================================
#
# load the simulations
#
# ======================================================================================
sims_last = load_galaxies.get_simulations_last(sys.argv[2:])
sims_common, common_scale = load_galaxies.get_simulations_common(sys.argv[2:])

# ======================================================================================
#
# CIMF calculations
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


def cimf(sim, mass_type, max_age_myr):
    """
    Make the cluster initial mass function.

    Notes on the different mass variables:
    ('STAR', 'INITIAL_MASS') - initial stellar mass of the star particle
    ('STAR', 'MASS') - current stellar mass of the star particle, accounting for
        stellar evolution
    ('STAR', 'INITIAL_BOUND_FRACTION') - NOT the actual initial_bound fraction. See
        the `get_initial_bound_fraction` function above for more on how to use this,
        but this variable is the accumulated mass near the cluster over the course
        of accretion. This is used to calculate formation efficiency, which is then
        used to get the actual initial_bound fraction
    ('STAR', 'BOUND_FRACTION') - This is the actual bound fraction at the current
        time, but NOT accounting for the proper initial_bound fraction

    :param sim: simulation object object
    :param mass_type: String encoding which mass to get here. The options are:
                      "initial" - just the initial stellar masses
                      "initial_bound" - initial masses including initial_bound
                                        fraction
                      "current" - current bound mass, accounting for the initial
                                  bound fraction, tidal disruption, and stellar
                                  death
    :param include_initial_bound: whether to incorporate the initial bound
                                  fraction of clusters, or just get the
                                  distribution of initial particle masses
    :param max_age_myr: The maximum age to restrict the plot to. Is infinity as the
                        default, which plots all stars.
    :returns: Two lists. The first is f_i * M, representing the initial
              bound mass, of M if include_initial_bound=False. This will be
              binned values suitable to plot. The second is dN/dLogM for each
              of the bins in the first list.
    """
    id_string = f"{mass_type}_{max_age_myr}"

    if id_string not in sim.precalculated:
        mass = []
        for galaxy in sim.galaxies:
            if mass_type == "initial":
                this_mass = galaxy[("STAR", "INITIAL_MASS")].to("Msun").value
            elif mass_type == "initial_bound":
                initial_mass = galaxy[("STAR", "INITIAL_MASS")].to("Msun").value
                star_initial_bound = get_initial_bound_fraction(galaxy)
                this_mass = initial_mass * star_initial_bound
            elif mass_type == "current":
                raw_mass = galaxy[("STAR", "INITIAL_MASS")].to("Msun").value
                star_initial_bound = get_initial_bound_fraction(galaxy)
                tidal_bound_fraction = galaxy[("STAR", "BOUND_FRACTION")].value
                this_mass = raw_mass * star_initial_bound * tidal_bound_fraction
            else:
                raise ValueError("Mass not recognized")

            # then restrict to recently formed clusters. This can be set to infinity,
            # which plots everything.
            mask = galaxy[("STAR", "age")] < max_age_myr * yt.units.Myr
            this_mass = this_mass[mask]

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

        sim.precalculated[id_string] = m_centers, hist
    return sim.precalculated[id_string]


# ======================================================================================
#
# plotting
#
# ======================================================================================
def plot_power_law(ax, slope, x1, x2, y1):
    # Here the slope is dN/dM \propto M^slope
    # Since we plot this in dN/dlogM space we
    # have to add 1 to the slope
    # dlogM \propto dM/M
    y2 = y1 * (x2 / x1) ** (slope + 1)

    ax.plot([x1, x2], [y1, y2], c=bpl.almost_black, lw=1, ls="--")
    ax.add_text(1.1 * x2, y2, text=slope, va="center", ha="left", fontsize=18)


def plot_cimf(axis_name, sim_share_type, masses_to_plot, max_age_myr=np.inf):
    """
    Plot various versions of the cluster mass function.

    :param axis_name: The name of the axis, so we can select which simulations go
                      on this plot.
    :param sim_share_type: Either "last" or "common". This will determine if
                           the redshift is labeled as common to all, or
                           individually per simulation
    :param masses_to_plot: A lit of the masses to plot, see the `cimf` function for the
                           allowed options
    :param max_age_myr: The maximum age to restrict the plot to. Is infinity as the
                        default, which plots all stars.
    """
    # validate options given
    if sim_share_type not in ["last", "common"]:
        raise ValueError("bad sim_share_type")
    # if we plot current, it should be alone, with nothing else
    if "current" in masses_to_plot and len(masses_to_plot) > 1:
        raise RuntimeError("Current masses must be plotted alone.")

    # make the name of the plot
    plot_name = f"cimf_{axis_name.replace('_', '')}_{sim_share_type}"
    if "current" in masses_to_plot:
        plot_name += "_current"
    elif (
        len(masses_to_plot) == 2
        and "initial" in masses_to_plot
        and "initial_bound" in masses_to_plot
    ):
        plot_name += "_initial"
    else:
        raise ValueError("This plot not yet named:", masses_to_plot)
    # add the age if it's not infinity
    if not np.isinf(max_age_myr):
        plot_name += f"_{max_age_myr}myr"
    plot_name += ".pdf"

    # choose which simulations we're using
    if sim_share_type == "common":
        sims = sims_common
    else:
        sims = sims_last

    fig, ax = bpl.subplots(figsize=[9, 7])

    for sim in sims:
        if axis_name not in sim.axes:
            continue
        for mass_type in masses_to_plot:
            mass_plot, dn_dlogM = cimf(sim, mass_type, max_age_myr)

            # make the label only for the biggest halo, and not for initial only
            if mass_type != "initial":
                # and include the redshift if it's different for each sim
                if sim_share_type == "last":
                    label = f"{sim.name}: z = {1/sim.ds.scale_factor - 1:.1f}"
                else:
                    label = sim.name
            else:
                label = None

            # have different line styles
            lss = {"initial": ":", "initial_bound": "-", "current": "-"}

            ax.plot(mass_plot, dn_dlogM, c=sim.color, ls=lss[mass_type], label=label)

    # formax axes
    ax.legend(loc=1)
    ax.set_yscale("log")
    ax.set_xscale("log")
    # have different y limits for different versions of the plot
    # the minimum value in the plot is 1 / (0.16 * ln(10) = 2.5
    # put the plot limit just above that, to make it cleaner and stop
    # weird vertical lines.
    if not np.isinf(max_age_myr):
        # small timeframe, don't need to show much.
        y_min = 3
        y_max = 1e4
    elif "current" in masses_to_plot:
        y_min = 10
        y_max = 1e4
    else:
        y_min = 10
        y_max = 1e5
    ax.set_limits(1e3, 1e7, y_min, y_max)

    # plot the guiding lines
    log_space = 0.4 * (np.log10(y_max) - np.log10(y_min))
    y_guide = 10 ** (np.log10(y_min) + log_space)
    plot_power_law(ax, -2, 1e6, 3e6, y_guide)
    plot_power_law(ax, -3, 1e6, 3e6, y_guide)

    # if there is a common redshift, annotate it
    if sim_share_type == "common":
        ax.easy_add_text(f"z = {1/common_scale - 1:.1f}", "upper left")

    if "current" in masses_to_plot:
        ax.add_labels("$f_b$M [$M_\odot$]", "dN/dlogM")
    else:
        ax.add_labels("$f_i M_i$ [$M_\odot$]", "dN/dlogM")

    fig.savefig(cimf_sentinel.parent / plot_name)


# ======================================================================================
#
# loop to actually make the plots
#
# ======================================================================================
# then actually call this function to build the plots
for plot_name in load_galaxies.get_plot_names([sim.name for sim in sims_last]):
    for share_type in ["common", "last"]:
        # plot main CIMF and unbound CIMF
        plot_cimf(plot_name, share_type, ["initial_bound", "initial"])
        # plot recently formed clusters
        plot_cimf(plot_name, share_type, ["initial_bound", "initial"], 300)
        # plot surviving clusters
        plot_cimf(plot_name, share_type, ["current"])

cimf_sentinel.touch()
