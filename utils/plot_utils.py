from pathlib import Path
from collections import defaultdict

import numpy as np
from matplotlib import cm
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
            "old_ic_sn_timing": "L18: Continuous",
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
            "old_ic_molecular": "$c_\\rho$=30, Other Changes",
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
            "lg_hn_fraction": "$f_{HN,0}=20$%",
            "lg_sfe": "T&L, $\epsilon_{ff}=100$%",
        },
        production("tl_sfe100_hn05"): {
            "lg_hn_fraction": "$f_{HN,0}=5$%",
        },
        production("tl_sfe100_hn00"): {
            "lg_fboost": "$f_{boost}=5$",
            "lg_hn_fraction": "$f_{HN,0}=0$",
        },
        production("tl_sfe100_hn00_fboost1"): {
            "lg_fboost": "$f_{boost}=1$",
        },
        production("tl_sfe100_hn00_fboost3"): {
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
        old_ic("discrete_hn00_virial10_entropy_fboost1"): h(0.00, 0.60, 0.65),
        old_ic("discrete_hn00_virial10_entropy_fboost2"): h(0.05, 0.70, 0.75),
        old_ic("discrete_hn00_virial10_entropy_fboost3"): h(0.08, 0.70, 0.90),
        old_ic("discrete_hn00_virial10_entropy"): h(0.13, 0.75, 0.90),
        old_ic("continuous_hn00_virial10_entropy_fboost1"): bpl.color_cycle[1],
        old_ic("discrete_hn50_virial10_entropy_fboost1"): h(0.80, 0.4, 0.6),
        old_ic("discrete_hn00_virial10_entropy_fboost1_sfe010"): h(0.00, 0.45, 0.70),
        old_ic("discrete_hn00_virial10_entropy_fboost1_sfe001"): h(0.00, 0.20, 0.80),
        old_ic("discrete_hn00_virial10_advect"): bpl.almost_black,
        old_ic("discrete_hn00_virial10"): h(0.06, 0.70, 0.80),
        old_ic("continuoushui_hn00_novirial"): h(0.06, 0.30, 0.85),
        old_ic("discrete_hn00_novirial_entropy_fboost1"): bpl.color_cycle[0],
        old_ic("discrete_hn00_virial10_entropy_fboost1_crho03"): h(0.05, 0.50, 0.60),
        old_ic("discrete_hn00_virial10_entropy_fboost1_crho30"): h(0.85, 0.30, 0.60),
        old_ic("discrete_hn00_virial10_entropy_molvadim_fboost1"): h(0.95, 0.4, 0.4),
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
        production("tl_sfe100_hn00_fboost1"): h(0.70, 0.55, 0.85),
        production("tl_sfe100_hn00_fboost3"): h(0.70, 0.35, 0.85),
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
