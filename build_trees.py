"""
build_trees.py

This handles taking the out.list files in form of out_a0.xxxx.list and turns
them into the out_N.list where N is an integer, such that Consistent Trees
can use them for merger trees. This assumes that all such out.list files are
in the halos directory. It then builds the trees.

This takes three arguments:
1. The path to the tree that will be created.
2. The path to the script used to build the tree config file
3. The directory where consistent_trees is located. 
"""
import sys
import subprocess
import os
from pathlib import Path
import shutil

# get the location of the sentinel
tree = Path(sys.argv[1]).resolve()
tree_config_script = Path(sys.argv[2]).resolve()
consistent_trees_dir = Path(sys.argv[3]).resolve()
# then turn that into other directories we'll need later
rockstar_dir = tree.parent.parent
halos_dir = rockstar_dir.parent / "halos"


def run_command(command):
    subprocess.call(command, shell=True)


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
# need to change the Stampede2 directory into the shangrila directory, and replace
# the number of snapshots with the actual value.
rockstar_config_file = rockstar_dir / "rockstar.cfg"
rockstar_config_file_temp = rockstar_dir / "rockstar.cfg.temp"
num_halos = len(all_halos)
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
                line = line.replace(str(value), str(num_halos))
            cfg_new.write(line)
# then replace the old with the new
rockstar_config_file_temp.replace(rockstar_config_file)

# Then build the config files
# perl $(tree_config_script) $(call tree_cfg_to_rockstar_cfg,$@)
rockstar_config_file = rockstar_dir / "rockstar.cfg"
run_command(f"perl {str(tree_config_script)} {rockstar_config_file}")

# then build the actual halo files. I need to change to the consistent trees
# directory to do this. I'll store where I was, then move back once done.
original_dir = os.getcwd()
# I need to change to the home dire
os.chdir(consistent_trees_dir)
config_file = rockstar_dir / "outputs" / "merger_tree.cfg"
run_command(f"perl do_merger_tree.pl {str(config_file)}")
# move back.
os.chdir(original_dir)

# then cleanup all the halo out files we copied at the beginning.
for f in rockstar_dir.iterdir():
    if f.name.startswith("out_") and f.name.endswith(".list"):
        f.unlink()
