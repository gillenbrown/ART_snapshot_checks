"""
merger_prep.py

This handles taking the out.list files in form of out_a0.xxxx.list and turns
them into the out_N.list where N is an integer, such that Consistent Trees
can use them for merger trees. This assumes that all such out.list files are
in the halos directory.

This only takes one argument:
the path to the sentinel
"""
import sys
from pathlib import Path
import shutil

# get the location of the sentinel
sentinel = Path(sys.argv[1])
# then turn that into other directories we'll need later
rockstar_dir = sentinel.parent
halos_dir = rockstar_dir.parent / "halos"

# first get rid of all the out.list files in rockstar dir. We'll replace those
# with new ones later
for f in rockstar_dir.iterdir():
    if f.name.startswith("out_") and f.name.endswith(".list"):
        f.unlink()

# then go through them all in halos and copy them to rockstar halos with their
# name changed
all_halos = [
    f
    for f in halos_dir.iterdir()
    if f.name.startswith("out_") and f.name.endswith(".list")
]
for n, f in enumerate(sorted(all_halos)):
    new_path = rockstar_dir / f"out_{n}.list"
    # copy the file
    shutil.copy2(f, new_path)

# I also need to modify the rockstar file. Specifically, for the production runs I
# need to change the Stampede2 directory into the shangrila directory.
rockstar_config_file = rockstar_dir / "rockstar.cfg"
rockstar_config_file_temp = rockstar_dir / "rockstar.cfg.temp"
with open(rockstar_config_file, "r") as cfg_old:
    with open(rockstar_config_file_temp, "w") as cfg_new:
        for line in cfg_old:
            if "tg862118" in line:
                line = line.replace(
                    "/scratch/06912/tg862118/art_runs/runs/",
                    "/u/home/gillenb/art_runs/runs/stampede2/",
                )
                print(line)
            # also update the number of files:
            if line.startswith("NUM_SNAPS"):
                value = int(line.split()[-1])
                if value == n + 1:
                    pass
                elif value == 2:
                    line = line.replace(str(value), str(n + 1))
                else:
                    raise RuntimeError("NUM_SNAPS doesn't make sense")
            cfg_new.write(line)
# then replace the old with the new
rockstar_config_file_temp.replace(rockstar_config_file)

# then touch the sentinel to finish the work
sentinel.touch()
