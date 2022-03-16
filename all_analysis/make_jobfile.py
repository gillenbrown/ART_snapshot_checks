from pathlib import Path

out_file = open("./job_list.txt", "w")

scratch_dir = Path("/scratch/06912/tg862118/art_runs/analysis/production")

directories = [
    scratch_dir / "rj_sfe010_hn20",
    scratch_dir / "rj_sfe100_hn20",
    scratch_dir / "tl_sfe001_hn20",
    scratch_dir / "tl_sfe010_hn20",
    scratch_dir / "tl_sfe100_hn00",
    scratch_dir / "tl_sfe100_hn00_fboost1",
    scratch_dir / "tl_sfe100_hn00_fboost3",
    scratch_dir / "tl_sfe100_hn05",
    scratch_dir / "tl_sfe100_hn20",
]

# ======================================================================================
#
# Then actually doing things
#
# ======================================================================================
for d in directories:
    out_dir = d / "run" / "out"
    art_files = sorted([f for f in out_dir.iterdir() if f.suffix == ".art"])
    for output in art_files:
        out_file.write(f"python -u all_analysis.py {str(output)}\n")

out_file.close()
