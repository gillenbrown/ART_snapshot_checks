from pathlib import Path

from matplotlib import cm

# I have to hardcode some labels to make this easier, parsing them won't work
# nearly as well
def full_dir(partial_path):
    base_dir = Path.home() / "art_runs" / "runs"
    return base_dir / partial_path


names = {
    full_dir("shangrila/hui/sfe_10"): "NBm SFE10",
    full_dir("shangrila/hui/sfe_50"): "NBm SFE50",
    full_dir("shangrila/hui/sfe_100"): "NBm SFE100",
    full_dir("shangrila/hui/sfe_200"): "NBm SFE200",
    full_dir("stampede2/old_ic_comparison/cap5000kms_hn00/run"): "Old IC SFE100 0% HN",
    full_dir("stampede2/old_ic_comparison/cap5000kms_hn05/run"): "Old IC SFE100 5% HN",
    full_dir("stampede2/old_ic_comparison/cap5000kms_hn20/run"): "Old IC SFE100 20% HN",
    full_dir("stampede2/old_ic_comparison/cap5000kms_hn50_v1/run"): "Old IC SFE100 50% HN",
    full_dir("stampede2/old_ic_comparison/cap5000kms_hn50/run"): "Old IC SFE100 50% HN new",
    full_dir("stampede2/production/sfe001_hn20/run"): "LG SFE1 20% HN",
    full_dir("stampede2/production/sfe010_hn20/run"): "LG SFE10 20% HN",
    full_dir("stampede2/production/sfe100_hn20/run"): "LG SFE100 20% HN",
    full_dir("stampede2/production/sfe100_hn05/run"): "LG SFE100 5% HN",
    full_dir("stampede2/production/sfe100_hn00/run"): "LG SFE100 0% HN",
}

cmap_lg_sfe = cm.Blues
cmap_lg_hnf = cm.Purples
# note that cmap_lg_sfe(1.0) should approximately equal cmap_lg_hnf(1.0)
cmap_nbm = cm.Greys
cmap_rerun = cm.Greens

colors = {
    "NBm SFE10": cmap_nbm(0.2),
    "NBm SFE50": cmap_nbm(0.35),
    "NBm SFE100": cmap_nbm(0.5),
    "NBm SFE200": cmap_nbm(0.65),
    "Old IC SFE100 0% HN": cmap_rerun(0.2),
    "Old IC SFE100 5% HN": cmap_rerun(0.35),
    "Old IC SFE100 20% HN": cmap_rerun(0.5),
    "Old IC SFE100 50% HN": cmap_rerun(0.65),
    "Old IC SFE100 50% HN new": cmap_rerun(0.8),
    "LG SFE1 20% HN": cmap_lg_sfe(0.4),
    "LG SFE10 20% HN": cmap_lg_sfe(0.7),
    "LG SFE100 20% HN": cmap_lg_hnf(1.0),
    "LG SFE100 5% HN": cmap_lg_hnf(0.7),
    "LG SFE100 0% HN": cmap_lg_hnf(0.4),
}
