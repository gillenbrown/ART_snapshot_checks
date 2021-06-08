"""
manage_halos.py - Get the halo directory ready on Stampede2

We do the following things in this script:
- Remove all the old files from the rockstar halos directory. We want to end up
  with only the last output from the previous run 
- Copy this last output over to the directory, renaming it appropriately
- Change the paramters in the rockstar restart.cfg file. We modify restart_snap
  to be 1, since we'll be starting from the second output next (the one after 
  the one we copied into this directory). We need to do this so the descendant
  information is still included when running rockstar. 

This script takes two parameters: 
- The name of the sentinel file to be created once this is done
- The hostname. We only do something if this is stampede2
"""
import sys
from pathlib import Path
import shutil
import re

if len(sys.argv) != 3:
    raise ValueError("Bad parameters to manage_halos.py")

sentinel = Path(sys.argv[1])
hostname = sys.argv[2]

if hostname != "stampede2":
    sentinel.touch()
    exit()

# ==============================================================================
# 
# Remove everything from the rockstar directory
# 
# ==============================================================================
rockstar_dir = sentinel.parent
for file in rockstar_dir.iterdir():
    if file.name not in ["restart.cfg", "sentinel.txt"] and not file.is_dir():
        file.unlink()
    elif file.name == "profiling" and file.is_dir():
        shutil.rmtree(file)

# ==============================================================================
# 
# Find the last output so we can copy it over
# 
# ==============================================================================
max_scale = "a0.0000"  # keep as string to avoid parsing. comparison still work
human_dir = rockstar_dir.parent / "halos"
# use a regular expression to get what we need to compare to
regex = re.compile("[a][0-1]\.[\d]{4}")
for file in human_dir.iterdir():
    match = regex.search(file.name)
    if match:
        scale = match.group(0)
        if scale > max_scale:
            max_scale = scale

# we can then find all the files that have this and copy them over
for file in human_dir.iterdir():
    if max_scale in file.name:
        new_file = rockstar_dir / file.name.replace(max_scale, "0")
        shutil.copy2(file, new_file)

# ==============================================================================
# 
# Modify the parameter file
# 
# ==============================================================================
restart_old = rockstar_dir / "restart.cfg"
restart_new = rockstar_dir / "restart.cfg.temp"

with open(restart_old, "r") as old:
    with open(restart_new, "w") as new:
        for line in old:
            split = line.split()
            key = split[0]
            value = split[-1]
            # we need to change the restart number
            if key == "RESTART_SNAP":
                line = line.replace(value, "1")
            # and the scale factor
            if key == "SCALE_NOW":
                # get rid of the "a" in the scale factor string we have above
                line = line.replace(value, max_scale.replace("a", ""))
            new.write(line)

# then copy the file over
shutil.move(restart_new, restart_old)

# then touch the sentinel file
sentinel.touch()