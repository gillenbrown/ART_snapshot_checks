"""
plot_cimf.py

Creates a plot showing the comparative cluster initial mass function of 
galaxies in an output
"""
import sys
from pathlib import Path

import numpy as np
from scipy import special
from astropy import cosmology
from astropy import units as u
import yt
from matplotlib import pyplot as plt
import betterplotlib as bpl

from tqdm import tqdm

from utils import load_galaxies, plot_utils

bpl.set_style()
yt.funcs.mylog.setLevel(50)  # ignore yt's output

cimf_sentinel = Path(sys.argv[1])

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


# ======================================================================================
#
# tidal evolution to z=0
#
# ======================================================================================
dt = 100 * yt.units.Myr


def t_tidal(M):
    omega_tid = 50 / yt.units.Gyr

    term_1 = 10 * yt.units.Gyr
    term_2 = (M / (2e5 * yt.units.Msun)) ** (2 / 3)
    term_3 = 100 / (omega_tid * yt.units.Gyr)
    total = term_1 * term_2 * term_3
    return total.to("Myr")


def get_n_steps(t_now, t_z_0):
    return int(np.ceil((t_z_0 - t_now).to("Myr") / dt))


precalculated_disruption = dict()


def precalculate_disruption(n_steps):
    if n_steps in precalculated_disruption:
        return
    else:
        precalculated_disruption[n_steps] = dict()

    # otherwise, we need to build this
    print(f"\n\n\nprecalculating\n{n_steps}\n\n\n")
    for log_m in tqdm(np.arange(2, 8, 0.01)):
        key_m = str(round(log_m, 2))

        m = 10 ** log_m * yt.units.Msun
        for _ in range(n_steps):
            if m < 100 * yt.units.Msun:
                m = 0
                break
            m *= np.exp(-dt / t_tidal(m))
        precalculated_disruption[n_steps][key_m] = m


def get_evolution_single_cluster(n_steps, log_m):
    try:
        key_m = str(round(log_m, 2))
        return precalculated_disruption[n_steps][key_m]
    except KeyError:  # only happens for initial M below 100, which will be fully
        # disrupted by z=0
        return 0


def evolve_cluster_population(galaxy):
    # get cosmology to get times
    # Initialize the cosmology object, used to put a redshift scale on the plots
    H_0 = galaxy.ds.artio_parameters["hubble"][0] * 100 * u.km / (u.Mpc * u.second)
    omega_matter = galaxy.ds.artio_parameters["OmegaM"][0]
    cosmo = cosmology.FlatLambdaCDM(H0=H_0, Om0=omega_matter, Tcmb0=2.725)

    t_now = galaxy.ds.current_time
    # need to convert the astropy units into yt units.
    t_z_0 = cosmo.age(0).to("Gyr").value * yt.units.Gyr
    n_steps = get_n_steps(t_now, t_z_0)

    # precalculate the disruption. This won't do anything if it already exists.
    precalculate_disruption(n_steps)

    # then get star masses (note that this ignores stellar evolution).
    raw_mass = galaxy[("STAR", "INITIAL_MASS")].to("Msun")
    star_initial_bound = get_initial_bound_fraction(galaxy)
    tidal_bound_fraction = galaxy[("STAR", "BOUND_FRACTION")].value
    cluster_masses = raw_mass * star_initial_bound * tidal_bound_fraction
    # get the log here. This makes using it in precalculation easier, as this is
    # vectorized. Set a minimum of 0.1 to avoid log of zero errors.
    log_cluster_masses_msun = np.log10(np.maximum(0.1, cluster_masses.to("Msun").value))

    evolved_masses = [
        get_evolution_single_cluster(n_steps, log_m)
        for log_m in tqdm(log_cluster_masses_msun)
    ]
    return np.array(evolved_masses)


# ======================================================================================
#
# CIMF itself
#
# ======================================================================================
def cimf(sim, mass_type, max_age_myr, max_z):
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
    :param max_z: The maximum metallicity to restrict the plot to. Is infinity as the
                  default, which plots all stars.
    :returns: Two lists. The first is f_i * M, representing the initial
              bound mass, of M if include_initial_bound=False. This will be
              binned values suitable to plot. The second is dN/dLogM for each
              of the bins in the first list.
    """
    id_string = f"{mass_type}_{max_age_myr}_{max_z}"

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
            elif mass_type == "evolved":
                this_mass = evolve_cluster_population(galaxy)
            else:
                raise ValueError("Mass not recognized")

            # then restrict the metallicity and age
            metallicity = galaxy[("STAR", "METALLICITY_SNII")].value
            metallicity += galaxy[("STAR", "METALLICITY_SNIa")].value
            if ("STAR", "METALLICITY_AGB") in sim.ds.field_list:
                metallicity += galaxy[("STAR", "METALLICITY_AGB")].value
            mask_z = metallicity < max_z
            mask_age = galaxy[("STAR", "age")] < max_age_myr * yt.units.Myr
            this_mass = this_mass[np.logical_and(mask_z, mask_age)]

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


def plot_cimf(
    axis_name, sim_share_type, masses_to_plot, max_age_myr=np.inf, max_z=np.inf
):
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
    if ("current" in masses_to_plot or "evolved" in masses_to_plot) and len(
        masses_to_plot
    ) > 1:
        raise RuntimeError("Current or evolved masses must be plotted alone.")

    # make the name of the plot
    plot_name = f"cimf_{axis_name}_{sim_share_type}"
    if "current" in masses_to_plot:
        plot_name += "_current"
    elif "evolved" in masses_to_plot:
        plot_name += "_evolved"
    elif len(masses_to_plot) == 1 and "initial" in masses_to_plot:
        plot_name += "_initial"
    elif (
        len(masses_to_plot) == 2
        and "initial" in masses_to_plot
        and "initial_bound" in masses_to_plot
    ):
        plot_name += "_initial_with_bound"
    else:
        raise ValueError("This plot not yet named:", masses_to_plot)
    # add the age if it's not infinity
    if not np.isinf(max_age_myr):
        plot_name += f"_{max_age_myr}myr"
    if not np.isinf(max_z):
        plot_name += f"_Z{max_z:.1e}"
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
            mass_plot, dn_dlogM = cimf(sim, mass_type, max_age_myr, max_z)

            # make the label only for the biggest halo, and not for initial only
            if len(masses_to_plot) == 1 or mass_type != "initial":
                label = plot_utils.plot_label(sim, sim_share_type, axis_name)
            else:
                label = None

            # have different line styles
            if len(masses_to_plot) == 1:
                ls = "-"
            else:
                ls = {
                    "initial": ":",
                    "initial_bound": "-",
                    "current": "-",
                    "evolved": "-",
                }[mass_type]

            ax.plot(mass_plot, dn_dlogM, c=sim.color, ls=ls, label=label)

    # format axes
    plot_utils.add_legend(ax, loc=1, fontsize=18)
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
    elif sim_share_type == "common" or axis_name == "lg_sfe":
        y_min = 10
        y_max = 1e5
    elif "current" in masses_to_plot:
        y_min = 10
        y_max = 1e4
    else:
        y_min = 10
        y_max = 1e6
    ax.set_limits(1e3, 1e7, y_min, y_max)

    # plot the guiding lines
    # log_space = 0.4 * (np.log10(y_max) - np.log10(y_min))
    # y_guide = 10 ** (np.log10(y_min) + log_space)
    # plot_power_law(ax, -2, 1e6, 3e6, y_guide)
    # plot_power_law(ax, -3, 1e6, 3e6, y_guide)

    # if there is a common redshift, annotate it
    if sim_share_type == "common":
        ax.easy_add_text(f"z = {common_z:.1f}", "upper left")

    if "current" in masses_to_plot or "evolved" in masses_to_plot:
        ax.add_labels("$f_b$M [$M_\odot$]", "dN/dlogM")
    elif len(masses_to_plot) == 1 and "initial" in masses_to_plot:
        ax.add_labels("$M_i$ [$M_\odot$]", "dN/dlogM")
    else:
        ax.add_labels("$f_i M_i$ [$M_\odot$]", "dN/dlogM")

    fig.savefig(cimf_sentinel.parent / plot_name)
    # then remove figure for memory purposes
    plt.close(fig)


# ======================================================================================
#
# loop to actually make the plots
#
# ======================================================================================
# then actually call this function to build the plots
for plot_name in load_galaxies.get_plot_names(sims_last):
    for share_type in ["common", "last"]:
        # plot particle masses
        plot_cimf(plot_name, share_type, ["initial"])
        # plot main CIMF and unbound CIMF
        plot_cimf(plot_name, share_type, ["initial_bound", "initial"])
        # plot recently formed clusters
        plot_cimf(plot_name, share_type, ["initial_bound", "initial"], max_age_myr=300)
        # plot surviving clusters
        plot_cimf(plot_name, share_type, ["current"])
    # plot low metallicity clusters. Only do this on the last output
    plot_cimf(plot_name, "last", ["initial"], max_z=0.001)
    # plot cluster population evolved to z=0. Only use the last output for this
    # plot_cimf(plot_name, "last", ["evolved"])

cimf_sentinel.touch()
