"""
merger_prep.py - Parse the files needed to get the merger trees ready to be made

This involves moving the out_0.list files from all halo subdirectories into 
the main directory, but renamed to put them in the right order. I also copy
the rockstar config file from the last one into the main directory to be used
by the merger config script.

This takes one command line argument: the name of the final rockstar config
file to be created.
"""
import sys
import os
from pathlib import Path 
import shutil

final_rockstar_file = Path(sys.argv[1]).absolute()
rockstar_dir = final_rockstar_file.parent

# remove all the existing out*.list files, since we'll recreate them soon
for item in rockstar_dir.iterdir():
    if not item.is_dir() and item.name.startswith("out_") and item.name.endswith(".list"):
        os.remove(item)

# get all the subirectories with halo information
subdirs = []
for item in rockstar_dir.iterdir():
    if item.is_dir() and item.name.startswith("halos_a"):
        subdirs.append(item)

# We want the files in order, so sort this
subdirs = sorted(subdirs)

# then go through and move the out_0.list files as appropriate
for idx, directory in enumerate(subdirs):
    original_filename = "out_0.list"
    new_filename = original_filename.replace("0", str(idx))

    shutil.copy2(directory / original_filename, rockstar_dir / new_filename)

# then copy the final rockstar config file to this directory too. I can't quite
# copy it, since there are a few things I need to modify
cfg_name = "rockstar.cfg"
with open(subdirs[-1] / cfg_name, "r") as subdir_cfg:
    with open(rockstar_dir / cfg_name, "w") as out_cfg:
        for line in subdir_cfg:
            # modify the directory name
            if line.split()[0] == "OUTBASE":
                # We want to remove the subdirectory name. To make this easier
                # we first throw out the quotes, then we can parse the name
                # directly
                clean_line = line.strip().replace('"', "")
                old_dir = Path(clean_line.split()[-1])
                line = line.replace("/" + old_dir.name, "")
            elif line.split()[0] == "NUM_SNAPS":
                old_snaps = line.split()[-1]
                new_snaps = str(len(subdirs))
                line = line.replace(old_snaps, new_snaps)
            # write either the unmodified line or the modified outbase line
            out_cfg.write(line)
    