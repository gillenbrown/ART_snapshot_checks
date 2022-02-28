import sys
from pathlib import Path
import numpy as np
import betterplotlib as bpl
from utils import run_attributes, plot_utils

bpl.set_style()

plot_name = Path(sys.argv[1]).resolve()

# ======================================================================================
#
# Set up data to show on plot
#
# ======================================================================================
# determine what sims to show and at what scale factor
direc_a = {
    run_attributes.old_ic("discrete_hn00_virial10_entropy_fboost1_sfe001"): "0.3901",
    run_attributes.old_ic("discrete_hn00_virial10_entropy_fboost1_sfe010"): "0.3906",
    run_attributes.old_ic("discrete_hn00_virial10_entropy_fboost1"): "0.3905",
    run_attributes.old_ic("discrete_hn00_virial10_entropy_fboost2_sfe001"): "0.2900",
    run_attributes.old_ic("discrete_hn00_virial10_entropy_fboost2_sfe010"): "0.2802",
    # run_attributes.old_ic("discrete_hn00_virial10_entropy_fboost2"): "0.4004",
    # run_attributes.production("tl_sfe001_hn20"): "0.2313",
    # run_attributes.production("tl_sfe010_hn20"): "0.2980",
    # run_attributes.production("tl_sfe100_hn20"): "0.2667",
    # run_attributes.production("tl_sfe100_hn05"): "0.3491",
    # run_attributes.production("tl_sfe100_hn00"): "0.3528",
    # run_attributes.production("tl_sfe100_hn00_fboost1"): "0.2418",
    # run_attributes.production("tl_sfe100_hn00_fboost3"): "0.2726",
    # run_attributes.production("rj_sfe100_hn20"): "0.3480",
    # run_attributes.production("rj_sfe010_hn20"): "0.2645",
}

# see here for more on linestyles:
# https://matplotlib.org/stable/gallery/lines_bars_and_markers/linestyles.html
success_ls = "-"
fail_ls = (0, (5, 1))
lss = {
    run_attributes.old_ic("discrete_hn00_virial10_entropy_fboost1_sfe001"): fail_ls,
    run_attributes.old_ic("discrete_hn00_virial10_entropy_fboost1_sfe010"): fail_ls,
    run_attributes.old_ic("discrete_hn00_virial10_entropy_fboost1"): success_ls,
    run_attributes.old_ic("discrete_hn00_virial10_entropy_fboost2_sfe001"): fail_ls,
    run_attributes.old_ic("discrete_hn00_virial10_entropy_fboost2_sfe010"): success_ls,
    # run_attributes.old_ic("discrete_hn00_virial10_entropy_fboost2"): success_ls,
    # run_attributes.production("tl_sfe001_hn20"): "--",
    # run_attributes.production("tl_sfe010_hn20"): success_ls,
    # run_attributes.production("tl_sfe100_hn20"): success_ls,
    # run_attributes.production("tl_sfe100_hn05"): success_ls,
    # run_attributes.production("tl_sfe100_hn00"): success_ls,
    # run_attributes.production("tl_sfe100_hn00_fboost1"): success_ls,
    # run_attributes.production("tl_sfe100_hn00_fboost3"): success_ls,
    # run_attributes.production("rj_sfe100_hn20"): success_ls,
    # run_attributes.production("rj_sfe010_hn20"): success_ls,
}

# define my own names for this
labels = {
    run_attributes.old_ic(
        "discrete_hn00_virial10_entropy_fboost1_sfe001"
    ): "$f_{boost}=1$, $\epsilon_{ff}=1$%",
    run_attributes.old_ic(
        "discrete_hn00_virial10_entropy_fboost1_sfe010"
    ): "$f_{boost}=1$, $\epsilon_{ff}=10$%",
    run_attributes.old_ic(
        "discrete_hn00_virial10_entropy_fboost1"
    ): "$f_{boost}=1$, $\epsilon_{ff}=100$%",
    run_attributes.old_ic(
        "discrete_hn00_virial10_entropy_fboost2_sfe001"
    ): "$f_{boost}=2$, $\epsilon_{ff}=1$%",
    run_attributes.old_ic(
        "discrete_hn00_virial10_entropy_fboost2_sfe010"
    ): "$f_{boost}=2$, $\epsilon_{ff}=10$%",
    run_attributes.old_ic(
        "discrete_hn00_virial10_entropy_fboost2"
    ): "$f_{boost}=2$, $\epsilon_{ff}=100$%",
}
# labels = {
#     run_attributes.old_ic(
#         "discrete_hn00_virial10_entropy_fboost1_sfe001"
#     ): "Isolated MW: $f_{boost}=1$, $\epsilon_{ff}=1$%, $f_{hn0}=0$",
#     run_attributes.old_ic(
#         "discrete_hn00_virial10_entropy_fboost1_sfe010"
#     ): "Isolated MW: $f_{boost}=1$, $\epsilon_{ff}=10$%, $f_{hn0}=0$",
#     run_attributes.old_ic(
#         "discrete_hn00_virial10_entropy_fboost1"
#     ): "Isolated MW: $f_{boost}=1$, $\epsilon_{ff}=100$%, $f_{hn0}=0$",
#     run_attributes.old_ic(
#         "discrete_hn00_virial10_entropy_fboost2_sfe001"
#     ): "Isolated MW: $f_{boost}=2$, $\epsilon_{ff}=1$%, $f_{hn0}=0$",
#     run_attributes.old_ic(
#         "discrete_hn00_virial10_entropy_fboost2_sfe010"
#     ): "Isolated MW: $f_{boost}=2$, $\epsilon_{ff}=10$%, $f_{hn0}=0$",
#     run_attributes.old_ic(
#         "discrete_hn00_virial10_entropy_fboost2"
#     ): "Isolated MW: $f_{boost}=2$, $\epsilon_{ff}=100$%",
#     run_attributes.production(
#         "tl_sfe001_hn20"
#     ): "Thelma & Louise: $f_{boost}=5$, $\epsilon_{ff}=1$%, $f_{hn0}=20$%",
#     run_attributes.production(
#         "tl_sfe010_hn20"
#     ): "Thelma & Louise: $f_{boost}=5$, $\epsilon_{ff}=10$%, $f_{hn0}=20$%",
#     run_attributes.production(
#         "tl_sfe100_hn20"
#     ): "Thelma & Louise: $f_{boost}=5$, $\epsilon_{ff}=100$%, $f_{hn0}=20$%",
#     run_attributes.production(
#         "tl_sfe100_hn05"
#     ): "Thelma & Louise: $f_{boost}=5$, $\epsilon_{ff}=100$%, $f_{hn0}=5$%",
#     run_attributes.production(
#         "tl_sfe100_hn00"
#     ): "Thelma & Louise: $f_{boost}=5$, $\epsilon_{ff}=100$%, $f_{hn0}=0$",
#     run_attributes.production(
#         "tl_sfe100_hn00_fboost1"
#     ): "Thelma & Louise: $f_{boost}=1$, $\epsilon_{ff}=100$%, $f_{hn0}=0$",
#     run_attributes.production(
#         "tl_sfe100_hn00_fboost3"
#     ): "Thelma & Louise: $f_{boost}=3$, $\epsilon_{ff}=100$%, $f_{hn0}=0$",
#     run_attributes.production(
#         "rj_sfe100_hn20"
#     ): "Romeo & Juliet: $f_{boost}=5$, $\epsilon_{ff}=100$%, $f_{hn0}=20$%",
#     run_attributes.production(
#         "rj_sfe010_hn20"
#     ): "Romeo & Juliet: $f_{boost}=5$, $\epsilon_{ff}=10$%, $f_{hn0}=20$%",
# }
# ======================================================================================
#
# analysis functions
#
# ======================================================================================
def read_gas_file(file_loc):
    n_h2 = []
    m_h2 = []
    for line in open(file_loc, "r"):
        if line.startswith("#"):
            continue
        this_n, this_m = line.split()
        n_h2.append(float(this_n))
        m_h2.append(float(this_m))
    return np.array(n_h2), np.array(m_h2)


def plot_pdf(ax, values, weights, x_min, x_max, bin_size, **kwargs):
    boundaries_log = np.arange(x_min - 0.5 * bin_size, x_max + 0.5 * bin_size, bin_size)
    centers_log = [
        np.mean([boundaries_log[idx], boundaries_log[idx + 1]])
        for idx in range(len(boundaries_log) - 1)
    ]

    boundaries = 10 ** boundaries_log
    centers = 10 ** np.array(centers_log)

    dx = np.histogram(a=values, bins=boundaries, weights=weights)[0]
    # make dx per log rho
    dx = dx / (bin_size * np.log(10))
    ax.plot(centers, dx, **kwargs)
    ax.axvline(500, ls=":", c=bpl.almost_black)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.add_labels("H$_2$ Number Density [$cm^{-3}$]", "dM/dlogn [$M_\odot$]")
    ax.set_limits(1, 2e5, 1e5, 1e9)


def plot_cumulative(ax, values, weights, normalize, **kwargs):
    sort_idxs = np.argsort(values)
    values = values[sort_idxs]
    weights = weights[sort_idxs]
    cumulative_weight = np.cumsum(weights)
    if normalize:
        cumulative_weight = cumulative_weight / cumulative_weight[-1]

    ax.plot(values, cumulative_weight, **kwargs)
    ax.set_xscale("log")
    if normalize:
        ax.add_labels("H$_2$ Number Density [$cm^{-3}$]", "Cumulative Mass Fraction")
        ax.set_limits(1e-3, 1e6, 0, 1)
    else:
        ax.add_labels("H$_2$ Number Density [$cm^{-3}$]", "Cumulative Mass [$M_\odot$]")
        ax.set_limits(1e-3, 1e6, 0, 5e8)


def format_exp(value):
    exp = int(np.floor(np.log10(value)))
    factor = value / 10 ** exp
    return f"{factor:.2f}" + "$\\times 10^{" + str(exp) + "}$"


# ======================================================================================
#
# Then make the plot
#
# ======================================================================================
fig, ax = bpl.subplots()
for direc, a in direc_a.items():
    n_h2, m_h2 = read_gas_file(direc / "checks" / f"gas_pdf_a{a}.txt")

    color = run_attributes.colors[direc]
    label = labels[direc] + "$M_{H_2}$ = " + format_exp(np.sum(m_h2)) + "$M_\odot$"
    ls = lss[direc]

    plot_pdf(ax, n_h2, m_h2, -4, 7, 0.2, color=color, ls=ls, label=label)
    # plot_cumulative(axs[1], n_h2, m_h2, False, color, label=None)
    # plot_cumulative(axs[2], n_h2, m_h2, True, color, label=label)

legend = plot_utils.add_legend(ax, loc=1, fontsize=14)
# legend = plot_utils.add_legend(
#     ax, bbox_to_anchor=(1.04, 1), loc="upper left", fontsize=10
# )
# legend._legend_box.align = "left"

fig.savefig(plot_name)
