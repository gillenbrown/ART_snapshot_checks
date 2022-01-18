from pathlib import Path

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


names = {
    base_dir
    / "pleiades/nbody/new_ic_trim_25mpc/root_08/run/outputs/vel_offset/": "T&L Collisionless",
    hui("sfe_10"): "NBm SFE10",
    hui("sfe_100"): "NBm SFE100",
    old_ic("continuous_hn00_novirial"): "Continuous",
    old_ic(
        "continuous_hn00_virial10_entropy_fboost1"
    ): "ART 2.1 Entropy $f_{boost}=1$ Continuous",
    old_ic("continuoushui_hn00_novirial"): "Continuous Hui",
    # old_ic("continuouspopmcluster_hn00_novirial"): "Continuous PopM",
    old_ic("continuoussnr_hn00_novirial"): "Continuous SNR",
    old_ic("discrete_hn00_novirial"): "Discrete",
    old_ic(
        "discrete_hn00_novirial_entropy_fboost1"
    ): "ART 2.1 Entropy $f_{boost}=1$ No Virial",
    old_ic("discrete_hn00_virial10"): "Discrete $\\alpha<10$",
    old_ic("discrete_hn00_virial10_fboost3"): "Discrete $\\alpha<10, f_{boost}=3$",
    old_ic("discrete_hn20_virial10"): "Discrete $\\alpha<10$ HN20",
    old_ic("discrete_hn00_virial10_19"): "ART 1.9 Adiabatic",
    old_ic("discrete_hn00_virial10_19_advect"): "ART 1.9 Advect",
    old_ic("discrete_hn00_virial10_advect"): "ART 2.0 Advect",
    old_ic("discrete_hn00_virial10_elements"): "ART 2.0 Adiabatic All Elements",
    old_ic("discrete_hn00_virial10_nostars"): "ART 2.0 Adiabatic No Stars",
    old_ic("discrete_hn00_virial10_19_old"): "ART 1.9 Hui's version Adiabatic",
    old_ic("discrete_hn00_virial10_advect_nostars"): "ART 2.0 Advect No Stars",
    old_ic("discrete_hn00_virial10_noadvoradia_nostars"): "ART 2.0 No Flags",
    old_ic("discrete_hn00_virial10_entropy"): "ART 2.1 Entropy $f_{boost}=5$",
    old_ic("discrete_hn00_virial10_entropy_fboost3"): "ART 2.1 Entropy $f_{boost}=3$",
    old_ic("discrete_hn00_virial10_entropy_fboost2"): "ART 2.1 Entropy $f_{boost}=2$",
    old_ic("discrete_hn00_virial10_entropy_fboost1"): "ART 2.1 Entropy $f_{boost}=1$",
    old_ic(
        "discrete_hn00_virial10_entropy_fboost1_crho03"
    ): "ART 2.1 Entropy $f_{boost}=1$ $C_\\rho=3",
    old_ic(
        "discrete_hn00_virial10_entropy_fboost1_crho30"
    ): "ART 2.1 Entropy $f_{boost}=1$ $C_\\rho=30",
    old_ic(
        "discrete_hn00_virial10_entropy_fboost1_sfe001"
    ): "ART 2.1 Entropy $f_{boost}=1$ $\epsilon_{ff}=1$%",
    old_ic(
        "discrete_hn00_virial10_entropy_fboost1_sfe010"
    ): "ART 2.1 Entropy $f_{boost}=1$ $\epsilon_{ff}=10$%",
    old_ic(
        "discrete_hn00_virial10_entropy_fboost3_nosnia"
    ): "ART 2.1 Entropy $f_{boost}=3$ No SNIa",
    old_ic(
        "discrete_hn00_virial10_entropy_molvadim"
    ): "ART 2.1 Entropy Molecular Changes $f_{boost}=5$",
    old_ic(
        "discrete_hn00_virial10_entropy_molvadim_fboost1"
    ): "ART 2.1 Entropy Molecular Changes $f_{boost}=1$",
    old_ic(
        "discrete_hn00_virial10_entropy_molvadim_fboost2"
    ): "ART 2.1 Entropy Molecular Changes $f_{boost}=2$",
    old_ic(
        "discrete_hn00_virial10_entropy_molvadim_fboost3"
    ): "ART 2.1 Entropy Molecular Changes $f_{boost}=3$",
    old_ic(
        "discrete_hn00_virial10_entropy_newagediff"
    ): "ART 2.1 Entropy SN Timing Hybrid",
    old_ic(
        "discrete_hn00_virial10_entropy_newagediffallave"
    ): "ART 2.1 Entropy SN Timing Average",
    old_ic(
        "discrete_hn00_virial10_entropy_newagediffallbirth"
    ): "ART 2.1 Entropy SN Timing Birth",
    old_ic(
        "discrete_hn00_virial10_entropy_noagediff"
    ): "ART 2.1 Entropy $f_{boost}=5$ No Age Diff",
    old_ic(
        "discrete_hn00_virial10_entropy_hybridagediff"
    ): "ART 2.1 Entropy $f_{boost}=5$ Hybrid Age Diff",
    old_ic(
        "discrete_hn50_virial10_entropy_fboost1"
    ): "ART 2.1 Entropy $f_{boost}=1$ HN50",
    old_ic("discrete_hn00_virial10_entropy_nosync"): "ART 2.1 No Sync",
    old_ic("discrete_hn00_virial10_noturb_adi"): "ART 2.0 No Turbulence Adiabatic",
    old_ic("discrete_hn00_virial10_noturb_adv"): "ART 2.0 No Turbulence Advect",
    production("tl_sfe001_hn20"): prod_fmt("T&L", 1, 20),
    production("tl_sfe010_hn20"): prod_fmt("T&L", 10, 20),
    production("tl_sfe100_hn20"): prod_fmt("T&L", 100, 20),
    production("tl_sfe100_hn05"): prod_fmt("T&L", 100, 5),
    production("tl_sfe100_hn00"): prod_fmt("T&L", 100, 0),
    production("tl_sfe100_hn00_fboost1"): prod_fmt("T&L", 100, 0) + " $f_{boost}=1",
    production("tl_sfe100_hn00_fboost3"): prod_fmt("T&L", 100, 0) + " $f_{boost}=3",
    production("rj_sfe010_hn20"): prod_fmt("R&J", 10, 20),
    production("rj_sfe100_hn20"): prod_fmt("R&J", 100, 20),
    rj_nbody("original_92.48mpc_level07"): "R&J Collisionless Original",
    rj_nbody("hybrid_46.24mpc_level08"): "R&J Collisionless 2x Trim",
    rj_nbody("hybrid_23.12mpc_level08"): "R&J Collisionless 4x Trim",
}

cmap_nbm = cm.Greys
cmap_rj_collisionless = cm.Reds


def hsv_to_hex(h, s, v):
    return mpl_col.to_hex(mpl_col.hsv_to_rgb([h, s, v]))


colors = {
    "NBm SFE10": cmap_nbm(0.35),
    "NBm SFE100": cmap_nbm(0.65),
    "Continuous": bpl.color_cycle[0],
    "ART 2.1 Entropy $f_{boost}=1$ Continuous": bpl.color_cycle[0],
    "Continuous Hui": bpl.color_cycle[3],
    "Continuous PopM": "lightblue",
    "Continuous SNR": bpl.color_cycle[3],
    "Discrete": bpl.color_cycle[1],
    "ART 2.1 Entropy $f_{boost}=1$ No Virial": bpl.color_cycle[0],
    "Discrete $\\alpha<10$": bpl.color_cycle[5],
    "Discrete $\\alpha<10, f_{boost}=3$": bpl.color_cycle[6],
    "Discrete $\\alpha<10$ HN20": bpl.color_cycle[7],
    "ART 1.9": "red",
    "ART 1.9 Advect": "green",
    "ART 1.9 Adiabatic": "blue",
    "ART 2.0 Adiabatic All Elements": "orange",
    "ART 2.0 Adiabatic No Stars": "brown",
    "ART 1.9 Hui's version Adiabatic": "yellow",
    "ART 2.0 Advect No Stars": "purple",
    "ART 2.0 No Flags": "cyan",
    "ART 2.0 Advect": bpl.color_cycle[0],
    "ART 2.0 No Turbulence Adiabatic": bpl.color_cycle[4],
    "ART 2.0 No Turbulence Advect": bpl.color_cycle[6],
    "ART 2.1 Entropy $f_{boost}=5$": bpl.color_cycle[1],
    "ART 2.1 Entropy $f_{boost}=3$": bpl.color_cycle[0],
    "ART 2.1 Entropy $f_{boost}=2$": bpl.color_cycle[3],
    "ART 2.1 Entropy $f_{boost}=1$": bpl.color_cycle[4],
    "ART 2.1 Entropy $f_{boost}=1$ $C_\\rho=3": bpl.color_cycle[1],
    "ART 2.1 Entropy $f_{boost}=1$ $C_\\rho=30": bpl.color_cycle[2],
    "ART 2.1 Entropy $f_{boost}=1$ $\epsilon_{ff}=1$%": bpl.color_cycle[0],
    "ART 2.1 Entropy $f_{boost}=1$ $\epsilon_{ff}=10$%": bpl.color_cycle[1],
    "ART 2.1 Entropy $f_{boost}=3$ No SNIa": bpl.color_cycle[5],
    "ART 2.1 Entropy $f_{boost}=5$ No Age Diff": bpl.color_cycle[0],
    "ART 2.1 Entropy $f_{boost}=5$ Hybrid Age Diff": bpl.color_cycle[3],
    "ART 2.1 Entropy Molecular Changes $f_{boost}=5$": bpl.color_cycle[0],
    "ART 2.1 Entropy Molecular Changes $f_{boost}=3$": bpl.color_cycle[3],
    "ART 2.1 Entropy Molecular Changes $f_{boost}=2$": bpl.color_cycle[5],
    "ART 2.1 Entropy Molecular Changes $f_{boost}=1$": bpl.color_cycle[6],
    "ART 2.1 Entropy SN Timing Hybrid": bpl.color_cycle[0],
    "ART 2.1 Entropy SN Timing Average": bpl.color_cycle[3],
    "ART 2.1 Entropy SN Timing Birth": bpl.color_cycle[4],
    "ART 2.1 Entropy $f_{boost}=1$ HN50": bpl.color_cycle[6],
    "ART 2.1 No Sync": bpl.color_cycle[3],
    # These colors are very carefully chosen to avoid colorblindness issues. The hue
    # changes between the SFE variations (blue) to the HN variations (purple), with
    # the shared SFE 100 HN 20 run in the middle. The blues are more saturated, while
    # the purples are less saturated. I found this essential to making the colors
    # distinguishable to those with colorblindess.
    # If you ever change these, use davidmathlogic.com/colorblind to doublecheck.
    prod_fmt("T&L", 1, 20): hsv_to_hex(0.60, 0.45, 0.85),
    prod_fmt("T&L", 10, 20): hsv_to_hex(0.60, 0.70, 0.75),
    prod_fmt("T&L", 100, 20): hsv_to_hex(0.65, 0.80, 0.50),
    prod_fmt("T&L", 100, 5): hsv_to_hex(0.70, 0.30, 0.65),
    prod_fmt("T&L", 100, 0): hsv_to_hex(0.70, 0.15, 0.85),
    prod_fmt("T&L", 100, 0) + " $f_{boost}=3": hsv_to_hex(0.90, 0.15, 0.85),
    prod_fmt("T&L", 100, 0) + " $f_{boost}=1": hsv_to_hex(0.95, 0.15, 0.85),
    prod_fmt("R&J", 10, 20): hsv_to_hex(0.35, 0.20, 0.70),
    prod_fmt("R&J", 100, 20): hsv_to_hex(0.35, 0.30, 0.50),
    "R&J Collisionless Original": cmap_rj_collisionless(1.0),
    "R&J Collisionless 2x Trim": cmap_rj_collisionless(0.7),
    "R&J Collisionless 4x Trim": cmap_rj_collisionless(0.4),
    "T&L Collisionless": bpl.almost_black,
}

axes = {
    "NBm SFE10": [],
    "NBm SFE100": [
        "adi_adv",
    ],
    "Continuous": [],
    "ART 2.1 Entropy $f_{boost}=1$ Continuous": ["old_ic_discreteness"],
    "Continuous Hui": ["adi_adv", "old_ic_sn_timing"],
    "Continuous PopM": [],
    "Continuous SNR": [],
    "Discrete": [],
    "ART 2.1 Entropy $f_{boost}=1$ No Virial": ["old_ic_virial"],
    "Discrete $\\alpha<10$": ["adi_adv"],
    "Discrete $\\alpha<10, f_{boost}=3$": [],
    "Discrete $\\alpha<10$ HN20": [],
    "ART 1.9": [],
    "ART 1.9 Advect": [],
    "ART 1.9 Adiabatic": [],
    "ART 2.0 Adiabatic All Elements": [],
    "ART 2.0 Adiabatic No Stars": [],
    "ART 1.9 Hui's version Adiabatic": [],
    "ART 2.0 Advect No Stars": [],
    "ART 2.0 No Flags": [],
    "ART 2.0 Advect": ["adi_adv"],
    "ART 2.0 No Turbulence Adiabatic": [],
    "ART 2.0 No Turbulence Advect": [],
    "ART 2.1 Entropy $f_{boost}=5$": ["old_ic_fboost"],
    "ART 2.1 Entropy $f_{boost}=3$": ["old_ic_fboost"],
    "ART 2.1 Entropy $f_{boost}=2$": ["old_ic_fboost"],
    "ART 2.1 Entropy $f_{boost}=1$": [
        "old_ic_fboost",
        "old_ic_hn_fraction",
        "old_ic_sfe",
        "old_ic_virial",
        "old_ic_discreteness",
        "old_ic_molecular",
    ],
    "ART 2.1 Entropy $f_{boost}=1$ $C_\\rho=3": ["old_ic_molecular"],
    "ART 2.1 Entropy $f_{boost}=1$ $C_\\rho=30": ["old_ic_molecular"],
    "ART 2.1 Entropy $f_{boost}=1$ $\epsilon_{ff}=1$%": ["old_ic_sfe"],
    "ART 2.1 Entropy $f_{boost}=1$ $\epsilon_{ff}=10$%": ["old_ic_sfe"],
    "ART 2.1 Entropy $f_{boost}=3$ No SNIa": [],
    "ART 2.1 Entropy $f_{boost}=5$ No Age Diff": [],
    "ART 2.1 Entropy $f_{boost}=5$ Hybrid Age Diff": [],
    "ART 2.1 Entropy Molecular Changes $f_{boost}=1$": ["old_ic_molecular"],
    "ART 2.1 Entropy Molecular Changes $f_{boost}=2$": [],
    "ART 2.1 Entropy Molecular Changes $f_{boost}=3$": [],
    "ART 2.1 Entropy Molecular Changes $f_{boost}=5$": [],
    "ART 2.1 Entropy SN Timing Hybrid": ["old_ic_sn_timing"],
    "ART 2.1 Entropy SN Timing Average": ["old_ic_sn_timing"],
    "ART 2.1 Entropy SN Timing Birth": ["old_ic_sn_timing"],
    "ART 2.1 Entropy $f_{boost}=1$ HN50": ["old_ic_hn_fraction"],
    "ART 2.1 No Sync": [],
    prod_fmt("T&L", 1, 20): ["lg_sfe"],
    prod_fmt("T&L", 10, 20): ["lg_sfe"],
    prod_fmt("T&L", 100, 20): ["lg_hn_fraction", "lg_sfe"],
    prod_fmt("T&L", 100, 5): ["lg_hn_fraction"],
    prod_fmt("T&L", 100, 0): ["lg_fboost", "lg_hn_fraction"],
    prod_fmt("T&L", 100, 0) + " $f_{boost}=1": ["lg_fboost"],
    prod_fmt("T&L", 100, 0) + " $f_{boost}=3": ["lg_fboost"],
    prod_fmt("R&J", 10, 20): ["lg_sfe"],
    prod_fmt("R&J", 100, 20): ["lg_sfe"],
    "R&J Collisionless Original": [],
    "R&J Collisionless 2x Trim": [],
    "R&J Collisionless 4x Trim": [],
    "T&L Collisionless": [],
}


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
        if "NBm" in l:
            sorts[h, l] = "000" + sort_label
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
