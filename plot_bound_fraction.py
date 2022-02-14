"""
plot_bound_fraction.py

Plot the cluster bound fraction and some related quantities
"""
from pathlib import Path

import sys
import numpy as np
from scipy import special, optimize
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

# and set up files to write mean and std to
files = dict()
for share_type in ["common", "last"]:
    for scaled in [True, False]:
        key = (share_type, scaled)
        if scaled:
            filename = f"fit_eps_int_ff_hist_{share_type}.txt"
        else:
            filename = f"fit_eps_int_hist_{share_type}.txt"
        files[key] = open(sentinel.parent / filename, "w")

        files[key].write("Format of the file is as follows:\n")
        files[key].write("Simulations are organized by the plots they are on.\n")
        files[key].write("A given simulation can appear multiple times.\n")
        files[key].write("Each row is log(mean), mean, std[dex]\n")
        files[key].write("\n\n")

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


def get_initial_bound_mass_msun(galaxy):
    return get_initial_mass_msun(galaxy) * get_initial_bound_fraction(galaxy)


def get_ages_gyr(galaxy):
    return galaxy[("STAR", "age")].to("Gyr").value


def get_dynamical_bound_fraction(galaxy):
    return galaxy[("STAR", "BOUND_FRACTION")].value


# ======================================================================================
#
# histograms of eps_int
#
# ======================================================================================
def get_eps_ff(sim_dir):
    sim_path = Path(sim_dir)
    if sim_path.name == "run":
        sim_path = sim_path.parent

    # now I should have the directory
    if "sfe001" in sim_path.name:
        return 0.01
    elif "sfe010" in sim_path.name:
        return 0.1
    else:
        return 1


def fit_gaussian(xs, ys, log=True):
    if log:
        xs = np.log10(xs)

    def to_minimize(params):
        norm = 10 ** params[2]
        test_ys = norm * plot_utils.gaussian(xs, params[0], params[1] ** 2)
        return np.sum((test_ys - ys) ** 2)

    result = optimize.minimize(
        to_minimize,
        x0=np.array([-1, 0.5, 0]),
        bounds=((-2, 1), (0.001, 0.5), (-3, 1)),
        tol=1e-10,
    )
    assert result.success
    return result.x


def plot_eps_int_histogram(axis_name, sim_share_type, scale_by_epsff=False):
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
    x_values = np.logspace(-3, 2, 501)
    # I need to adjust the limits of the plot depending on the SFE of the sims present
    # set a default value which will be adjusted
    x_min = 1e-2
    x_max = 1

    # set up the file to write out the gaussian fits
    out_file = files[(sim_share_type, scale_by_epsff)]
    out_file.write(f"{axis_name}\n")
    for sim in sims:
        if axis_name not in sim.axes:
            continue

        mass = sim.func_all_galaxies(get_initial_mass_msun)
        eps_int = sim.func_all_galaxies(get_eps_int)

        if scale_by_epsff:
            eps_ff = get_eps_ff(sim.run_dir)
            eps_int = eps_int / eps_ff

        plot_y = plot_utils.kde(x_values, eps_int, width=0.05, weights=mass, log=True)
        if sim.reliable:
            ls = "-"
        else:
            ls = "--"
        ax.plot(
            x_values,
            plot_y,
            c=sim.color,
            ls=ls,
            label=plot_utils.plot_label(sim, sim_share_type, axis_name),
        )

        # also fit a Gaussian to this distribution
        log_mean, std_dex, log_norm = fit_gaussian(x_values, plot_y, log=True)
        # gauss_ys = 10 ** log_norm * plot_utils.gaussian(
        #     np.log10(x_values), log_mean, std_dex ** 2
        # )
        # ax.plot(x_values, gauss_ys, ls=":", c=sim.color)

        # write them to the correct file
        out_file.write(
            f"{sim.names[axis_name]} {log_mean:.2f} {10**log_mean:.4f} {std_dex:.2f}\n"
        )

        # if there is something other than sfe100, update the limit
        if "sfe010" in str(sim.run_dir) or "sfe001" in str(sim.run_dir):
            if scale_by_epsff:
                x_max = 10
            else:
                x_min = 1e-3

    # format axes. limits depend on which runs are shown
    ax.set_xscale("log")
    ax.set_limits(x_min, x_max, 0)
    plot_utils.add_legend(ax, loc=2, fontsize=15, frameon=False)
    plot_utils.nice_log_axis(ax, "x")
    if scale_by_epsff:
        ax.add_labels(
            "$\epsilon_{int} / \epsilon_{ff}$",
            "Mass Weighted KDE Density",
        )
    else:
        ax.add_labels(
            "Integrated Star Formation Efficiency $\epsilon_{int}$",
            "Mass Weighted KDE Density",
        )

    # if there is a common redshift, annotate it
    if sim_share_type == "common":
        ax.easy_add_text(f"z = {common_z:.1f}", "upper right")

    if scale_by_epsff:
        fig.savefig(
            sentinel.parent / f"eps_int_ff_hist_{axis_name}_{sim_share_type}.pdf"
        )

    else:
        fig.savefig(sentinel.parent / f"eps_int_hist_{axis_name}_{sim_share_type}.pdf")
    plt.close(fig)  # to save memory
    out_file.write("\n\n")


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

        mass = sim.func_all_galaxies(get_initial_mass_msun)
        f_b = sim.func_all_galaxies(get_initial_bound_fraction)

        plot_utils.shaded_region(
            ax,
            mass,
            f_b,
            sim.color,
            p_lo=25,
            p_hi=75,
            log_x=True,
            cut_x=sim.unreliable_mass,
            label=plot_utils.plot_label(sim, sim_share_type, axis_name),
        )

    # format axes
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_limits(1e2, 1e7, 0.001, 1)
    ax.add_labels("Initial Particle Mass [$M_\odot$]", "Initial Bound Fraction")
    ax.legend(loc=4, frameon=False)
    plot_utils.nice_log_axis(ax, "y")

    # if there is a common redshift, annotate it
    if sim_share_type == "common":
        ax.easy_add_text(f"z = {common_z:.1f}", "upper left")

    fig.savefig(sentinel.parent / f"bound_fraction_{axis_name}_{sim_share_type}.pdf")
    plt.close(fig)  # to save memory


def plot_eps_int(axis_name, sim_share_type):
    """
    Plot the cluster integrated star formation efficiency as a function of mass.

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

        mass = sim.func_all_galaxies(get_initial_mass_msun)
        eps_int = sim.func_all_galaxies(get_eps_int)

        plot_utils.shaded_region(
            ax,
            mass,
            eps_int,
            sim.color,
            p_lo=25,
            p_hi=75,
            log_x=True,
            cut_x=sim.unreliable_mass,
            label=plot_utils.plot_label(sim, sim_share_type, axis_name),
        )

    # format axes
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_limits(1e2, 1e7, 1e-3, 1)
    ax.add_labels(
        "Initial Particle Mass [$M_\odot$]",
        "Integrated Star Formation Efficiency $\epsilon_{int}$",
    )
    ax.legend(loc=4, frameon=False)
    plot_utils.nice_log_axis(ax, "y")

    # if there is a common redshift, annotate it
    if sim_share_type == "common":
        ax.easy_add_text(f"z = {common_z:.1f}", "upper left")

    fig.savefig(sentinel.parent / f"eps_int_{axis_name}_{sim_share_type}.pdf")
    plt.close(fig)  # to save memory


def plot_dynamical_bound(sim):
    """
    Plot the dynamical bound fraction for a given simulation
    """
    m_initial = sim.func_all_galaxies(get_initial_bound_mass_msun)
    age = sim.func_all_galaxies(get_ages_gyr)
    f_bound_dyn = sim.func_all_galaxies(get_dynamical_bound_fraction)

    # then make a plot showing the band with age for different mass ranges
    fig, ax = bpl.subplots()

    dlogm = 1
    logm_max = 6
    c_idx = 0
    colors = ["#03A09B", "#76BD81", "#E9DA67"]
    while logm_max > 3.5:
        logm_min = logm_max - dlogm
        assert np.isclose(logm_min, int(logm_min))  # only integers for the legend
        label = (
            "$10^{" + f"{logm_min:.0f}" + "} - 10^{" + f"{logm_max:.0f}" + "} M_\odot$"
        )

        # select the clusters in this mass range.
        idx = np.logical_and(m_initial > 10 ** logm_min, m_initial < 10 ** logm_max)
        if np.sum(idx) < 2:  # not enough clusters in this range
            continue

        # figure out if this mass range is reliable. To do this we assign the whole
        # age range to be unreliable at this mass range
        if logm_min >= np.log10(sim.unreliable_mass):
            cut_age = 0
        else:
            cut_age = np.inf

        plot_utils.shaded_region(
            ax,
            age[idx],
            f_bound_dyn[idx],
            colors[c_idx],
            p_lo=25,
            p_hi=75,
            cut_x=cut_age,
            dx={3: 0.1, 4: 0.3, 5: 0.3}[logm_min],
            log_x=False,
            label=label,
        )

        c_idx += 1
        logm_max -= dlogm

    ax.add_labels("Cluster Age [Gyr]", "Dynamical Bound Fraction")
    # ax.set_xscale("log")
    # ax.set_limits(0.01, 5, 0, 1)
    ax.set_limits(0, 4.5, 0, 1)
    ax.legend(loc=1, frameon=False, fontsize=14)
    fig.savefig(
        sentinel.parent / f"dynamical_{plot_utils.get_sim_dirname(sim.run_dir)}.pdf"
    )
    plt.close(fig)  # to save memory


# ======================================================================================
#
# make the plots
#
# ======================================================================================
for plot_name in load_galaxies.get_plot_names(sims_last):
    for share_type in ["last", "common"]:
        plot_bound_fraction(plot_name, share_type)
        plot_eps_int_histogram(plot_name, share_type, scale_by_epsff=False)
        plot_eps_int_histogram(plot_name, share_type, scale_by_epsff=True)
        plot_eps_int(plot_name, share_type)

for sim in sims_last:
    plot_dynamical_bound(sim)

# close the files
for f in files.values():
    f.close()

# touch the sentinel once done
sentinel.touch()
