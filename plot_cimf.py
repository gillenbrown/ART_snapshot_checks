"""
plot_cimf.py

Creates a plot showing the comparative cluster initial mass function of 
galaxies in an output
"""
import sys
from pathlib import Path

import numpy as np
from scipy import optimize

import yt
from matplotlib import pyplot as plt
import betterplotlib as bpl

from utils import load_galaxies, plot_utils
from analysis_functions import cimf

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
# fitting power law to CIMF
#
# ======================================================================================
# define the minium mass to use
min_fit_m = 1e5
max_fit_m = 1e7
m_pivot = 1e5


def power_law(xs, slope, log_norm):
    return 10 ** log_norm * ((xs / m_pivot) ** slope)


def fit_power_law_base(masses, dn_dlogm):
    """
    Fit a power law slope to the CIMF.

    Note that this fit's the slope of dN/dlogM and reports that slope. Typically one
    reports the slope dN/dM, so to correct for the one must subtract 1 from the value
    returned here.

    Proof for myself since I forget this every time (treat equals as proportional):
    dN/dM = M^alpha

    but dlogM = dM / M
    so dM = M dlogM

    dN/dM = dN/(M dlogM) = M^alpha
    dN/dlogM = M^(alpha + 1) = M^k
    where k is what I'm fitting here, the slope in dN/dlogM.
    alpha + 1 = k
    alpha = k - 1
    """
    # get errors. Assume Poisson errors. Turn dn_dlogm to N first
    dlogm = np.log(10) * (np.log10(masses[1]) - np.log10(masses[0]))

    dn = dn_dlogm * dlogm
    dn_err = np.sqrt(dn)
    # the Poisson error for a bin with zero clusters is technically zero, but make it
    # 0.5 so as to not break the fitting, and be less than the error for one cluster.
    dn_err[np.isclose(dn_err, 0)] = 0.5

    # restrict the fitting to the mass range. Note that I don't restrict that there are
    # nonzero clusters in a bin, but I do change the errors to account for that.
    idxs = np.logical_and(masses < max_fit_m, masses > min_fit_m)
    fit_masses = masses[idxs]
    fit_dn = dn[idxs]
    fit_dn_err = dn_err[idxs]

    if len(fit_masses) == 1:
        return -99, -99

    # do the fitting
    def to_minimize(params):
        slope, log_norm = params
        predicted_n = power_law(fit_masses, slope, log_norm)
        return np.sum((predicted_n - fit_dn) ** 2 / fit_dn_err)

    fit_result = optimize.minimize(to_minimize, x0=(-2, 4), bounds=((-5, 0), (0, 6)))
    if fit_result.success:
        # correct log_norm for the offset for the log10 in the normalization of the y
        # value. We want to fit in raw numbers, but then plot in dN/dlogM, so we need to
        # add this factor back. We take another log since we return the log of the
        # normalization
        return fit_result.x[0], fit_result.x[1] + np.log10(np.log(10)), dn_err / dlogm
    else:
        return -99, -99, dn_err / dlogm


def fit_power_law(sim_share_type, max_age_myr=np.inf, max_z=np.inf):
    # choose which simulations we're using
    if sim_share_type == "common":
        sims = sims_common
    else:
        sims = sims_last

    with open(
        cimf_sentinel.parent / f"cimf_slopes_{sim_share_type}.txt", "w"
    ) as out_file:
        for sim in sims:
            sim_name = plot_utils.get_sim_dirname(sim.run_dir)
            mass_plot, dn_dlogM = cimf.cimf(sim, "initial", max_age_myr, max_z)
            slope, log_norm, yerr = fit_power_law_base(mass_plot, dn_dlogM)

            # make a plot
            fig, ax = bpl.subplots()
            ax.plot(mass_plot, dn_dlogM, c=sim.color)
            # include shaded region. Artificially adjust the shaded region for visual
            # purposes
            visual_factor = 2
            ax.fill_between(
                mass_plot,
                dn_dlogM + visual_factor * yerr,
                dn_dlogM - visual_factor * yerr,
                color=sim.color,
                alpha=0.5,
            )

            # power law slope lines
            fit_m = np.logspace(np.log10(min_fit_m), np.log10(max_fit_m), 100)
            fit_n = power_law(fit_m, slope, log_norm)
            ax.plot(fit_m, fit_n, ls="--", lw=2, c=bpl.almost_black)
            ax.easy_add_text("$\\alpha=" + f"{slope - 1:.2f}" + "$", "upper right")

            # format axes
            ax.set_yscale("log")
            ax.set_xscale("log")
            # annotate limits used in fit
            ax.axvline(min_fit_m, ls=":", lw=1, c=bpl.almost_black)
            ax.axvline(max_fit_m, ls=":", lw=1, c=bpl.almost_black)
            ax.set_limits(1e3, 1e7, 1, 1e6)
            ax.add_labels("$f_i M_i$ [$M_\odot$]", "dN/dlogM")
            ax.easy_add_text(f"z={sim.z:.1f}", "upper left")
            fig.savefig(
                cimf_sentinel.parent / f"cimf_slopefit_{sim_name}_{sim_share_type}.pdf"
            )
            # then remove figure for memory purposes
            plt.close(fig)

            # write to output file
            out_file.write(f"{str(sim_name)}: {slope - 1:.2f}\n")


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
    axis_name,
    sim_share_type,
    masses_to_plot,
    max_age_myr=np.inf,
    max_z=np.inf,
    guiding_lines=False,
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

    # note that guiding lines are added
    if guiding_lines:
        plot_name += "_slopeguide"

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
            mass_plot, dn_dlogM = cimf.cimf(sim, mass_type, max_age_myr, max_z)

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

            # for runs that are unreliable, plot them with solid lines up to 10^5, then
            # dashed lines above that
            plot_utils.plot_line_with_cut(
                ax,
                mass_plot,
                dn_dlogM,
                cut_x=sim.unreliable_mass,
                c=sim.color,
                ls1=ls,
                ls2="--",
                label=label,
            )

    # if we're plotting evolved masses, include the MW GCs
    if "evolved" in masses_to_plot:
        ax.plot(
            *cimf.harric_gc_mass_function(), ls="--", c=bpl.almost_black, label="MW GCs"
        )
    # format axes
    plot_utils.add_legend(ax, loc=1, fontsize=18)
    ax.set_yscale("log")
    ax.set_xscale("log")
    # have different y limits for different versions of the plot
    # the minimum value in the plot is 1 / (0.16 * ln(10) = 2.5
    # put the plot limit just below that to show the maximum cluster mass
    y_min = 2
    if not np.isinf(max_age_myr):
        # small timeframe, don't need to show much.
        y_max = 1e4
    elif "current" in masses_to_plot:
        y_max = 1e4
    elif (
        sim_share_type == "common"
        and axis_name == "lg_sfe"
        and "current" in masses_to_plot
    ):
        y_max = 1e6
    elif sim_share_type == "common" or axis_name == "lg_sfe":
        y_max = 1e5
    else:
        y_max = 1e6
    ax.set_limits(1e3, 1e7, y_min, y_max)

    # plot the guiding lines
    if guiding_lines:
        log_space = 0.5 * (np.log10(y_max) - np.log10(y_min))
        y_guide = 10 ** (np.log10(y_min) + log_space)
        plot_power_law(ax, -2, 2e6, 5e6, y_guide)
        plot_power_law(ax, -3, 2e6, 5e6, y_guide)

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
        # plot particle masses with and without guiding lines.
        plot_cimf(plot_name, share_type, ["initial"], guiding_lines=False)
        plot_cimf(plot_name, share_type, ["initial"], guiding_lines=True)
        # plot main CIMF and unbound CIMF. Guiding lines off by default
        plot_cimf(plot_name, share_type, ["initial_bound", "initial"])
        # plot recently formed clusters
        plot_cimf(plot_name, share_type, ["initial_bound", "initial"], max_age_myr=300)
        # plot surviving clusters, with and without guiding lines
        plot_cimf(plot_name, share_type, ["current"], guiding_lines=False)
        plot_cimf(plot_name, share_type, ["current"], guiding_lines=True)
    # plot low metallicity clusters. Only do this on the last output
    plot_cimf(plot_name, "last", ["initial"], max_z=0.001)
    # plot cluster population evolved to z=0. Only use the last output for this
    # plot_cimf(plot_name, "last", ["evolved"])
# fit the power law slope
fit_power_law("last")
fit_power_law("common")

cimf_sentinel.touch()
