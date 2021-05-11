from pathlib import Path

from matplotlib import cm

# I have to hardcode some labels to make this easier, parsing them won't work
# nearly as well
def full_dir(partial_path):
    base_dir = Path.home() / "art_runs" / "runs"
    return base_dir / partial_path


names = {
    full_dir("pleiades/nbody/new_ic_trim_25mpc/root_08/run/outputs/vel_offset/"): "T&L Collisionless",
    full_dir("shangrila/hui/sfe_10"): "NBm SFE10",
    full_dir("shangrila/hui/sfe_100"): "NBm SFE100",
    full_dir("stampede2/old_ic_comparison/cap5000kms_hn00/run"): "Old IC SFE100 0% HN",
    full_dir("stampede2/old_ic_comparison/cap5000kms_hn01/run"): "Old IC SFE100 1% HN",
    full_dir("stampede2/old_ic_comparison/cap5000kms_hn10/run"): "Old IC SFE100 10% HN",
    full_dir("stampede2/old_ic_comparison/cap5000kms_hn20/run"): "Old IC SFE100 20% HN",
    full_dir("stampede2/old_ic_comparison/cap5000kms_hn50_v1/run"): "Old IC SFE100 50% HN",
    full_dir("stampede2/old_ic_comparison/cap5000kms_hn50/run"): "Old IC SFE100 50% HN new",
    full_dir("stampede2/production/tl_sfe001_hn20/run"): "T&L SFE1 20% HN",
    full_dir("stampede2/production/tl_sfe010_hn20/run"): "T&L SFE10 20% HN",
    full_dir("stampede2/production/tl_sfe100_hn20/run"): "T&L SFE100 20% HN",
    full_dir("stampede2/production/tl_sfe100_hn05/run"): "T&L SFE100 5% HN",
    full_dir("stampede2/production/tl_sfe100_hn00/run"): "T&L SFE100 0% HN",
    full_dir("stampede2/production/rj_sfe010_hn20/run"): "R&J SFE10 20% HN",
    full_dir("stampede2/production/rj_sfe100_hn20/run"): "R&J SFE100 20% HN",
    full_dir("stampede2/rj_nbody/original_92.48mpc_level07/run"): "R&J Collisionless Original",
    full_dir("stampede2/rj_nbody/hybrid_46.24mpc_level08/run"): "R&J Collisionless 2x Trim",
    full_dir("stampede2/rj_nbody/hybrid_23.12mpc_level08/run"): "R&J Collisionless 4x Trim",
}

cmap_tl_sfe = cm.Blues
cmap_tl_hnf = cm.Purples
cmap_rj_sfe = cm.Oranges
# note that cmap_lg_sfe(1.0) should approximately equal cmap_lg_hnf(1.0)
cmap_nbm = cm.Greys
cmap_rerun = cm.Greens

colors = {
    "NBm SFE10": cmap_nbm(0.35),
    "NBm SFE100": cmap_nbm(0.65),
    "Old IC SFE100 0% HN": cmap_rerun(0.3),
    "Old IC SFE100 1% HN": cmap_rerun(0.4),
    "Old IC SFE100 10% HN": cmap_rerun(0.5),
    "Old IC SFE100 20% HN": cmap_rerun(0.6),
    "Old IC SFE100 50% HN": cmap_rerun(0.7),
    "Old IC SFE100 50% HN new": cmap_rerun(0.8),
    "T&L SFE1 20% HN": cmap_tl_sfe(0.4),
    "T&L SFE10 20% HN": cmap_tl_sfe(0.7),
    "T&L SFE100 20% HN": cmap_tl_hnf(1.0),
    "T&L SFE100 5% HN": cmap_tl_hnf(0.7),
    "T&L SFE100 0% HN": cmap_tl_hnf(0.4),
    "R&J SFE10 20% HN": cmap_rj_sfe(0.5),
    "R&J SFE100 20% HN": cmap_rj_sfe(1.0),
}

axis_number = {
    "NBm SFE10": 0,
    "NBm SFE100": 0,
    "Old IC SFE100 0% HN": 0,
    "Old IC SFE100 1% HN": 0,
    "Old IC SFE100 10% HN": 0,
    "Old IC SFE100 20% HN": 0,
    "Old IC SFE100 50% HN": 0,
    "Old IC SFE100 50% HN new": 0,
    "T&L SFE1 20% HN": 1,
    "T&L SFE10 20% HN": 1,
    "T&L SFE100 20% HN": 1,
    "T&L SFE100 5% HN": 1,
    "T&L SFE100 0% HN": 1,
    "R&J SFE10 20% HN": 2,
    "R&J SFE100 20% HN": 2,
}