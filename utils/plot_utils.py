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
names = {
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
        "old_ic_fboost": "$f_{boost}=5$",
    },
    old_ic("discrete_hn00_virial10_entropy_fboost3"): {
        "old_ic_fboost": "$f_{boost}=3$",
    },
    old_ic("discrete_hn00_virial10_entropy_fboost2"): {
        "old_ic_fboost": "$f_{boost}=2$",
    },
    old_ic("discrete_hn00_virial10_entropy_fboost1"): {
        "old_ic_fboost": "$f_{boost}=1$",
        "old_ic_hn_fraction": "$f_{HN,0}=0$",
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
        "old_ic_hn_fraction": "$f_{HN,0}=50$%",
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
}


def create_color_cmap(hex_color, min_saturation=0.1, max_value=0.8):
    """
    Create a colormap that fades from one color to nearly white.

    This is done by converting the color to HSV, then decreasing the saturation while
    increasing the value (which makes it closer to white)

    :param hex_color: Original starting color, must be in hex format
    :param min_saturation: The saturation of the point farthest from the original color
    :param max_value: The value of the point farthest from the original color
    :return: A matplotilb colormap. Calling it with 0 returns the color specififed
             by `min_saturation` and `max_value` while keeping the same hue, while
             1 will return the original color.
    """
    # convert to HSV (rgb required as an intermediate)
    base_color_rgb = mpl_col.hex2color(hex_color)
    h, s, v = mpl_col.rgb_to_hsv(base_color_rgb)
    N = 256  # number of points in final colormap
    # check that this color is within the range specified
    assert s > min_saturation
    assert v < max_value
    # reduce the saturation and up the brightness. Start from the outer values, as these
    # will correspond to 0, while the original color will be 1
    saturations = np.linspace(min_saturation, s, N)
    values = np.linspace(max_value, v, N)
    out_xs = np.linspace(0, 1, N)

    # set up the weird format required by LinearSegmentedColormap
    cmap_dict = {"red": [], "blue": [], "green": []}
    for idx in range(N):
        r, g, b = mpl_col.hsv_to_rgb((h, saturations[idx], values[idx]))
        out_x = out_xs[idx]
        # LinearSegmentedColormap requires a weird format. I don't think the difference
        # in the last two values matters, it seems to work fine without it.
        cmap_dict["red"].append((out_x, r, r))
        cmap_dict["green"].append((out_x, g, g))
        cmap_dict["blue"].append((out_x, b, b))
    return mpl_col.LinearSegmentedColormap(hex_color, cmap_dict, N=256)


def hsv_to_hex(h, s, v):
    return mpl_col.to_hex(mpl_col.hsv_to_rgb([h, s, v]))


cmap_old_ic_sfe = create_color_cmap(bpl.color_cycle[0])


def default_color():
    return "#bafc03"  # an ugly color that's noticeable, to let me know to replace it


colors = defaultdict(
    default_color,
    {
        hui("sfe_100"): "#AAAAAA",
        old_ic("discrete_hn00_virial10_entropy_fboost1"): bpl.color_cycle[0],
        old_ic("continuous_hn00_virial10_entropy_fboost1"): bpl.color_cycle[1],
        old_ic("discrete_hn50_virial10_entropy_fboost1"): bpl.color_cycle[4],
        old_ic("discrete_hn00_virial10_entropy_fboost1_sfe010"): cmap_old_ic_sfe(0.65),
        old_ic("discrete_hn00_virial10_entropy_fboost1_sfe001"): cmap_old_ic_sfe(0.3),
        old_ic("discrete_hn00_virial10_advect"): bpl.almost_black,
        old_ic("discrete_hn00_virial10"): bpl.color_cycle[3],
        # These colors are very carefully chosen to avoid colorblindness issues. The hue
        # changes between the SFE variations (blue) to the HN variations (purple), with
        # the shared SFE 100 HN 20 run in the middle. The blues are more saturated,
        # while the purples are less saturated. I found this essential to making the
        # colors distinguishable to those with colorblindess.
        # If you ever change these, use davidmathlogic.com/colorblind to doublecheck.
        production("tl_sfe001_hn20"): hsv_to_hex(0.60, 0.45, 0.85),
        production("tl_sfe010_hn20"): hsv_to_hex(0.60, 0.70, 0.75),
        production("tl_sfe100_hn20"): hsv_to_hex(0.65, 0.80, 0.50),
        production("tl_sfe100_hn05"): hsv_to_hex(0.70, 0.30, 0.65),
        production("tl_sfe100_hn00"): hsv_to_hex(0.70, 0.15, 0.85),
        production("tl_sfe100_hn00_fboost1"): hsv_to_hex(0.90, 0.15, 0.85),
        production("tl_sfe100_hn00_fboost3"): hsv_to_hex(0.80, 0.15, 0.85),
        production("rj_sfe010_hn20"): hsv_to_hex(0.35, 0.20, 0.70),
        production("rj_sfe100_hn20"): hsv_to_hex(0.35, 0.30, 0.50),
        rj_nbody("original_92.48mpc_level07"): cmap_rj_collisionless(1.0),
        rj_nbody("hybrid_46.24mpc_level08"): cmap_rj_collisionless(0.7),
        rj_nbody("hybrid_23.12mpc_level08"): cmap_rj_collisionless(0.4),
    },
)

if len(list(colors.values())) != len(set(colors.values())):
    raise ValueError("Not all colors are unique! Fix please!")


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
