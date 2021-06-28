from pathlib import Path

from matplotlib import cm
from matplotlib import colors as mpl_col

# I have to hardcode some labels to make this easier, parsing them won't work
# nearly as well
def full_dir(partial_path):
    base_dir = Path.home() / "art_runs" / "runs"
    return base_dir / partial_path


def prod_fmt(ic, eps_ff, f_hn):
    name = f"{ic} "
    name += "$\epsilon_{ff} = $" + f"{eps_ff:<3}%, "
    name += "$f_{HN} = $" + f"{f_hn}%"
    return name


names = {
    full_dir(
        "pleiades/nbody/new_ic_trim_25mpc/root_08/run/outputs/vel_offset/"
    ): "T&L Collisionless",
    full_dir("shangrila/hui/sfe_10"): "NBm SFE10",
    full_dir("shangrila/hui/sfe_100"): "NBm SFE100",
    full_dir("stampede2/old_ic_comparison/cap5000kms_hn00/run"): "Old IC SFE100 0% HN",
    full_dir("stampede2/old_ic_comparison/cap5000kms_hn01/run"): "Old IC SFE100 1% HN",
    full_dir("stampede2/old_ic_comparison/cap5000kms_hn10/run"): "Old IC SFE100 10% HN",
    full_dir("stampede2/old_ic_comparison/cap5000kms_hn20/run"): "Old IC SFE100 20% HN",
    full_dir(
        "stampede2/old_ic_comparison/cap5000kms_hn50_v1/run"
    ): "Old IC SFE100 50% HN",
    full_dir(
        "stampede2/old_ic_comparison/cap5000kms_hn50/run"
    ): "Old IC SFE100 50% HN new",
    full_dir("stampede2/production/tl_sfe001_hn20/run"): prod_fmt("T&L", 1, 20),
    full_dir("stampede2/production/tl_sfe010_hn20/run"): prod_fmt("T&L", 10, 20),
    full_dir("stampede2/production/tl_sfe100_hn20/run"): prod_fmt("T&L", 100, 20),
    full_dir("stampede2/production/tl_sfe100_hn05/run"): prod_fmt("T&L", 100, 5),
    full_dir("stampede2/production/tl_sfe100_hn00/run"): prod_fmt("T&L", 100, 0),
    full_dir("stampede2/production/rj_sfe010_hn20/run"): prod_fmt("R&J", 10, 20),
    full_dir("stampede2/production/rj_sfe100_hn20/run"): prod_fmt("R&J", 100, 20),
    full_dir(
        "stampede2/rj_nbody/original_92.48mpc_level07/run"
    ): "R&J Collisionless Original",
    full_dir(
        "stampede2/rj_nbody/hybrid_46.24mpc_level08/run"
    ): "R&J Collisionless 2x Trim",
    full_dir(
        "stampede2/rj_nbody/hybrid_23.12mpc_level08/run"
    ): "R&J Collisionless 4x Trim",
}

cmap_nbm = cm.Greys
cmap_rerun = cm.Greens
cmap_rj_collisionless = cm.Reds


def hsv_to_hex(h, s, v):
    return mpl_col.to_hex(mpl_col.hsv_to_rgb([h, s, v]))


colors = {
    "NBm SFE10": cmap_nbm(0.35),
    "NBm SFE100": cmap_nbm(0.65),
    "Old IC SFE100 0% HN": cmap_rerun(0.3),
    "Old IC SFE100 1% HN": cmap_rerun(0.4),
    "Old IC SFE100 10% HN": cmap_rerun(0.5),
    "Old IC SFE100 20% HN": cmap_rerun(0.6),
    "Old IC SFE100 50% HN": cmap_rerun(0.7),
    "Old IC SFE100 50% HN new": cmap_rerun(0.8),
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
    prod_fmt("R&J", 10, 20): hsv_to_hex(0.35, 0.20, 0.70),
    prod_fmt("R&J", 100, 20): hsv_to_hex(0.35, 0.30, 0.50),
    "R&J Collisionless Original": cmap_rj_collisionless(1.0),
    "R&J Collisionless 2x Trim": cmap_rj_collisionless(0.7),
    "R&J Collisionless 4x Trim": cmap_rj_collisionless(0.4),
}

axes = {
    "NBm SFE10": ["old_ic", "all"],
    "NBm SFE100": ["old_ic", "all"],
    "Old IC SFE100 0% HN": ["old_ic", "all"],
    "Old IC SFE100 1% HN": ["old_ic", "all"],
    "Old IC SFE100 10% HN": ["old_ic", "all"],
    "Old IC SFE100 20% HN": ["old_ic", "all"],
    "Old IC SFE100 50% HN": ["old_ic", "all"],
    "Old IC SFE100 50% HN new": ["old_ic", "all"],
    prod_fmt("T&L", 1, 20): ["tl", "lg", "all"],
    prod_fmt("T&L", 10, 20): ["tl", "lg", "all"],
    prod_fmt("T&L", 100, 20): ["tl", "lg", "all"],
    prod_fmt("T&L", 100, 5): ["tl", "lg", "all"],
    prod_fmt("T&L", 100, 0): ["tl", "lg", "all"],
    prod_fmt("R&J", 10, 20): ["rj", "lg", "all"],
    prod_fmt("R&J", 100, 20): ["rj", "lg", "all"],
    "R&J Collisionless Original": ["rj", "lg", "nbody", "all"],
    "R&J Collisionless 2x Trim": ["rj", "lg", "nbody", "all"],
    "R&J Collisionless 4x Trim": ["rj", "lg", "nbody", "all"],
}


def get_plot_names(sim_names):
    all_plots = []
    for name in sim_names:
        for plot in axes[name]:
            all_plots.append(plot)
    return list(set(all_plots))
