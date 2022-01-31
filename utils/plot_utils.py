from pathlib import Path
from collections import defaultdict

import numpy as np
from matplotlib import cm, ticker
from matplotlib import colors as mpl_col

import betterplotlib as bpl

# ======================================================================================
#
# Set up plot info
#
# ======================================================================================
# I have to hardcode some labels to make this easier, parsing them won't work
# nearly as well
base_dir = Path.home() / "art_runs" / "runs"


def hui(suffix):
    return base_dir / "shangrila" / "hui" / suffix


def old_ic(suffix):
    return (
        base_dir / "stampede2" / "old_ic_comparison_production_analog" / suffix / "run"
    )


def production(suffix):
    return base_dir / "stampede2" / "production" / suffix / "run"


def rj_nbody(suffix):
    return base_dir / "stampede2" / "rj_nbody" / suffix / "run"


def stampede2_analysis(suffix):
    analysis_dir = Path("/scratch/06912/tg862118/art_runs/analysis/production")
    return analysis_dir / suffix / "run"


def prod_fmt(ic, eps_ff, f_hn):
    name = f"{ic} "
    name += "$\epsilon_{ff} = $" + f"{eps_ff:}%, "
    # add spaces to pad lower percents:
    if eps_ff < 10:
        name += " " * 4
    elif eps_ff < 100:
        name += " " * 2
    name += "$f_{HN} = $" + f"{f_hn}%"
    return name


# if the simulation is not included here, it is not used in any plots
names = defaultdict(
    dict,
    {
        hui("sfe_100"): {
            "adi_adv": "L18 Hydro - L18 Feedback",
        },
        old_ic("continuous_hn00_virial10_entropy_fboost1"): {
            "old_ic_discreteness": "Continuous"
        },
        old_ic("continuoushui_hn00_novirial"): {
            "adi_adv": "New Hydro - L18 Feedback",
            "old_ic_sn_timing": "L18",
        },
        old_ic("discrete_hn00_novirial_entropy_fboost1"): {
            "old_ic_virial": "No Virial Criterion"
        },
        old_ic("discrete_hn00_virial10"): {
            "adi_adv": "New Hydro - New Feedback",
        },
        old_ic("discrete_hn00_virial10_advect"): {
            "adi_adv": "L18 Hydro - New Feedback",
        },
        old_ic("discrete_hn00_virial10_entropy"): {
            "old_ic_sn_feedback": "$f_{boost}=5$, $f_{HN,0}=0$",
        },
        old_ic("discrete_hn00_virial10_entropy_fboost3"): {
            "old_ic_sn_feedback": "$f_{boost}=3$, $f_{HN,0}=0$",
        },
        old_ic("discrete_hn00_virial10_entropy_fboost2"): {
            "old_ic_sn_feedback": "$f_{boost}=2$, $f_{HN,0}=0$",
        },
        old_ic("discrete_hn00_virial10_entropy_fboost1"): {
            "old_ic_sn_feedback": "$f_{boost}=1$, $f_{HN,0}=0$",
            "old_ic_sfe": "$\epsilon_{ff}=100$%",
            "old_ic_virial": "Virial Criterion",
            "old_ic_discreteness": "Discrete SN",
            "old_ic_molecular": "$c_\\rho$=10",
        },
        old_ic("discrete_hn00_virial10_entropy_fboost1_crho03"): {
            "old_ic_molecular": "$c_\\rho$=3",
        },
        old_ic("discrete_hn00_virial10_entropy_fboost1_crho30"): {
            "old_ic_molecular": "$c_\\rho$=30",
        },
        old_ic("discrete_hn00_virial10_entropy_fboost1_sfe001"): {
            "old_ic_sfe": "$\epsilon_{ff}=1$%",
        },
        old_ic("discrete_hn00_virial10_entropy_fboost1_sfe010"): {
            "old_ic_sfe": "$\epsilon_{ff}=10$%",
        },
        old_ic("discrete_hn00_virial10_entropy_molvadim_fboost1"): {
            "old_ic_molecular": "$c_\\rho$=3, GK11 Shielding",
        },
        old_ic("discrete_hn00_virial10_entropy_newagediff"): {
            "old_ic_sn_timing": "Hybrid Approach",
        },
        old_ic("discrete_hn00_virial10_entropy_newagediffallave"): {
            "old_ic_sn_timing": "Average Approach",
        },
        old_ic("discrete_hn00_virial10_entropy_newagediffallbirth"): {
            "old_ic_sn_timing": "Birth Approach",
        },
        old_ic("discrete_hn50_virial10_entropy_fboost1"): {
            "old_ic_sn_feedback": "$f_{boost}=1$, $f_{HN,0}=50$%",
        },
        production("tl_sfe001_hn20"): {
            "lg_sfe": "T&L, $\epsilon_{ff}=1$%",
        },
        production("tl_sfe010_hn20"): {
            "lg_sfe": "T&L, $\epsilon_{ff}=10$%",
        },
        production("tl_sfe100_hn20"): {
            "lg_sn_feedback": "$f_{boost}=5$, $f_{HN,0}=20$%",
            "lg_hn_fraction": "$f_{HN,0}=20$%",
            "lg_sfe": "T&L, $\epsilon_{ff}=100$%",
        },
        production("tl_sfe100_hn05"): {
            "lg_sn_feedback": "$f_{boost}=5$, $f_{HN,0}=5$%",
            "lg_hn_fraction": "$f_{HN,0}=5$%",
        },
        production("tl_sfe100_hn00"): {
            "lg_sn_feedback": "$f_{boost}=5$, $f_{HN,0}=0$",
            "lg_fboost": "$f_{boost}=5$",
            "lg_hn_fraction": "$f_{HN,0}=0$",
        },
        production("tl_sfe100_hn00_fboost1"): {
            "lg_sn_feedback": "$f_{boost}=1$, $f_{HN,0}=0$",
            "lg_fboost": "$f_{boost}=1$",
        },
        production("tl_sfe100_hn00_fboost3"): {
            "lg_sn_feedback": "$f_{boost}=3$, $f_{HN,0}=0$",
            "lg_fboost": "$f_{boost}=3$",
        },
        production("rj_sfe010_hn20"): {
            "lg_sfe": "R&J, $\epsilon_{ff}=10$%",
        },
        production("rj_sfe100_hn20"): {
            "lg_sfe": "R&J, $\epsilon_{ff}=100$%",
        },
    },
)


def h(h, s, v):  # stands for hsv to hex
    return mpl_col.to_hex(mpl_col.hsv_to_rgb([h, s, v]))


cmap_rj_collisionless = cm.Reds


def default_color():
    return "#bafc03"  # an ugly color that's noticeable, to let me know to replace it


# I have a bit of a scheme to follow here. Each IC gets it's own distinct color:
# T&L is blue, R&J is green, and the old IC is red.
# Within each IC, a decrease in the SFE fades the color by increasing value and
# decreasing saturation (for old IC I'll just change value)
# f_boost is just a change in saturation for the old IC

colors = defaultdict(
    default_color,
    {
        hui("sfe_100"): "#AAAAAA",
        old_ic("discrete_hn00_virial10_entropy_fboost1"): h(0.00, 0.70, 0.5),
        old_ic("discrete_hn00_virial10_entropy_fboost2"): h(0.05, 0.70, 0.75),
        old_ic("discrete_hn00_virial10_entropy_fboost3"): h(0.08, 0.70, 0.90),
        old_ic("discrete_hn00_virial10_entropy"): h(0.13, 0.75, 0.90),
        old_ic("continuous_hn00_virial10_entropy_fboost1"): bpl.color_cycle[0],
        old_ic("discrete_hn50_virial10_entropy_fboost1"): h(0.80, 0.4, 0.6),
        old_ic("discrete_hn00_virial10_entropy_fboost1_sfe010"): h(0.00, 0.50, 0.60),
        old_ic("discrete_hn00_virial10_entropy_fboost1_sfe001"): h(0.00, 0.20, 0.80),
        old_ic("discrete_hn00_virial10_advect"): bpl.almost_black,
        old_ic("discrete_hn00_virial10"): h(0.06, 0.70, 0.80),
        old_ic("continuoushui_hn00_novirial"): h(0.06, 0.30, 0.85),
        old_ic("discrete_hn00_novirial_entropy_fboost1"): bpl.color_cycle[1],
        old_ic("discrete_hn00_virial10_entropy_fboost1_crho03"): h(0.05, 0.50, 0.60),
        old_ic("discrete_hn00_virial10_entropy_fboost1_crho30"): h(0.85, 0.30, 0.70),
        old_ic("discrete_hn00_virial10_entropy_molvadim_fboost1"): h(0.95, 0.4, 0.3),
        old_ic("discrete_hn00_virial10_entropy_newagediff"): h(0.55, 0.35, 0.45),
        old_ic("discrete_hn00_virial10_entropy_newagediffallave"): h(0.35, 0.45, 0.75),
        old_ic("discrete_hn00_virial10_entropy_newagediffallbirth"): h(0.2, 0.6, 0.8),
        # These colors are very carefully chosen to avoid colorblindness issues. The hue
        # changes between the SFE variations (blue) to the HN variations (purple), with
        # the shared SFE 100 HN 20 run in the middle. The blues are more saturated,
        # while the purples are less saturated. I found this essential to making the
        # colors distinguishable to those with colorblindess.
        # If you ever change these, use davidmathlogic.com/colorblind to doublecheck.
        production("tl_sfe001_hn20"): h(0.60, 0.45, 0.85),
        production("tl_sfe010_hn20"): h(0.60, 0.70, 0.75),
        production("tl_sfe100_hn20"): h(0.65, 0.80, 0.50),
        production("tl_sfe100_hn05"): h(0.70, 0.30, 0.65),
        production("tl_sfe100_hn00"): h(0.70, 0.15, 0.85),
        production("tl_sfe100_hn00_fboost1"): h(0.70, 0.75, 0.90),
        production("tl_sfe100_hn00_fboost3"): h(0.70, 0.55, 0.85),
        production("rj_sfe010_hn20"): h(0.35, 0.20, 0.70),
        production("rj_sfe100_hn20"): h(0.35, 0.30, 0.50),
        rj_nbody("original_92.48mpc_level07"): cmap_rj_collisionless(1.0),
        rj_nbody("hybrid_46.24mpc_level08"): cmap_rj_collisionless(0.7),
        rj_nbody("hybrid_23.12mpc_level08"): cmap_rj_collisionless(0.4),
    },
)

if len(list(colors.values())) != len(set(colors.values())):
    raise ValueError("Not all colors are unique! Fix please!")

# add duplicates for Stampede2 analysis
production_list = [
    "tl_sfe001_hn20",
    "tl_sfe010_hn20",
    "tl_sfe100_hn20",
    "tl_sfe100_hn05",
    "tl_sfe100_hn00",
    "tl_sfe100_hn00_fboost1",
    "tl_sfe100_hn00_fboost3",
    "rj_sfe010_hn20",
    "rj_sfe100_hn20",
]
for suffix in production_list:
    names[stampede2_analysis(suffix)] = names[production(suffix)]
    colors[stampede2_analysis(suffix)] = colors[production(suffix)]


# ======================================================================================
#
# functions for commonly used plotting ideas
#
# ======================================================================================
def add_legend(ax, loc=0, fontsize=12, frameon=False, **kwargs):
    # make a default legend
    ax.legend()
    # then get everything from it
    handles, labels = ax.get_legend_handles_labels()

    if len(handles) == 0:
        return

    # then figure out how to sort them
    sorts = dict()
    for h, l in zip(handles, labels):
        # get rid of the redshift part of the label
        if ": z = " in l:
            sort_label = l[: l.find(":")]
        else:
            sort_label = l
        # add extras to manually change the sort for labels I want at the front
        # or back
        if "NBm" in l or "$f_{boost}=1$, $f_{HN,0}=50$%" in l:
            # ! is the first real ASCII character
            sorts[h, l] = "!!!" + sort_label
        elif "Universe Machine" in l:
            sorts[h, l] = "zzz" + sort_label
        else:
            sorts[h, l] = sort_label

    # then actually do the sorting
    handles, labels = zip(
        *sorted([hl for hl in zip(handles, labels)], key=lambda hl: sorts[hl])
    )

    return ax.legend(
        handles=handles,
        labels=labels,
        loc=loc,
        fontsize=fontsize,
        frameon=frameon,
        **kwargs,
    )


def moving_percentiles(xs, ys, percentile, dx):
    bins = np.arange(min(xs), max(xs) + dx, dx)

    bin_centers = []
    radii_percentiles = []
    for idx in range(len(bins) - 1):
        lower = bins[idx]
        upper = bins[idx + 1]

        # then find all clusters in this mass range
        mask_above = xs > lower
        mask_below = xs < upper
        mask_good = np.logical_and(mask_above, mask_below)

        good_ys = ys[mask_good]
        if len(good_ys) > 1:
            radii_percentiles.append(np.percentile(good_ys, percentile))
            # the bin centers will be the mean in log space
            bin_centers.append(np.mean([lower, upper]))

    return np.array(bin_centers), np.array(radii_percentiles)


def shaded_region(ax, xs, ys, color, p_lo=10, p_hi=90, dx=0.2, log_x=False, label=None):
    if log_x:
        xs = np.log10(xs)

    c_lo, p_lo = moving_percentiles(xs, ys, p_lo, dx)
    c_50, p_50 = moving_percentiles(xs, ys, 50, dx)
    c_hi, p_hi = moving_percentiles(xs, ys, p_hi, dx)

    if log_x:
        c_lo = 10 ** c_lo
        c_50 = 10 ** c_50
        c_hi = 10 ** c_hi

    # The X values should be the same for all percentiles
    assert np.allclose(c_lo, c_50)
    assert np.allclose(c_hi, c_50)

    ax.plot(c_50, p_50, c=color, zorder=1, label=label)
    ax.fill_between(x=c_lo, y1=p_lo, y2=p_hi, color=color, alpha=0.5, zorder=0)


# Function to use to set the ticks
@ticker.FuncFormatter
def nice_log_formatter(x, pos):
    exp = np.log10(x)
    # this only works for labels that are factors of 10. Other values will produce
    # misleading results, so check this assumption.
    assert np.isclose(exp, int(exp))

    # for values between 0.01 and 100, just use that value.
    # Otherwise use the log.
    if abs(exp) < 3:
        return f"{x:g}"
    else:
        return "$10^{" + f"{exp:.0f}" + "}$"


def nice_log_axis(ax, which):
    if which not in ["x", "y", "both"]:
        raise ValueError("error specifying axis in nice_log_axis")
    if which in ["x", "both"]:
        ax.xaxis.set_major_formatter(nice_log_formatter)
    if which in ["y", "both"]:
        ax.yaxis.set_major_formatter(nice_log_formatter)


# ======================================================================================
#
# functions for KDE histograms
#
# ======================================================================================
def gaussian(x, mean, variance):
    """
    Normalized Gaussian Function at a given value.

    Is normalized to integrate to 1.

    :param x: value to calculate the Gaussian at
    :param mean: mean value of the Gaussian
    :param variance: Variance of the Gaussian distribution
    :return: log of the likelihood at x
    """
    exp_term = np.exp(-((x - mean) ** 2) / (2 * variance))
    normalization = 1.0 / np.sqrt(2 * np.pi * variance)
    return exp_term * normalization


def kde(x_grid, x_values, width, weights, log=False):
    ys = np.zeros(x_grid.size)

    if log:
        x_grid = np.log10(x_grid)
        x_values = np.log10(x_values)

    for x, weight in zip(x_values, weights):
        ys += weight * gaussian(x_grid, x, width ** 2)

    # normalize the y value
    ys = np.array(ys)
    ys = 10 * ys / np.sum(ys)
    return ys
