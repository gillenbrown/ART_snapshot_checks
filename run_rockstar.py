import sys
from pathlib import Path
import os
import shutil
import yt
yt.funcs.mylog.setLevel(0)  # ignore yt's output

# format of sys.argv:
# idx 0: name of script
# ids 1: directory of simulation outputs to run halos on
# idx 2: directory to store halo finding outputs in
# idx 3: directory to store the nicely formatted outputs in

# print function that always flushes
import functools
print = functools.partial(print, flush=True)

# do error checking on root only.
if yt.is_root():
    if len(sys.argv) != 5:
        raise ValueError("Please provide the proper command line argument.")

# turn the directories the user passes into the absolute path
sim_dir = Path(sys.argv[1]).resolve()
rockstar_dir = Path(sys.argv[2]).resolve()
halos_dir = Path(sys.argv[3]).resolve()
machine = sys.argv[4]
# Create the subdirectory where we'll do the actual halo finding
temp_dir = sim_dir / "temp_output_dir_for_halo_finding"
if not temp_dir.is_dir() and yt.is_root():
    temp_dir.mkdir()

print(f"Reading simulations from: {sim_dir}")
print(f"Writing intermediate catalogs to: {rockstar_dir}")
print(f"Writing final catalogs to: {halos_dir}")

# ==============================================================================
#
# Key functions to be used later
#
# ==============================================================================
def move_all_simulation_files(art_file, old_dir, new_dir):
    for file in old_dir.iterdir():
        if file.stem == art_file.stem:
            file.replace(new_dir / file.name)

def scale_factor_from_art_file(art_file):
    return art_file.stem[-6:]

def make_human_readable_halo_files(rockstar_idx, scale_factor, method="copy"):
    """
    This moves one set of halo files to the halos directory with nicer names

    :param rockstar_idx: The index Rockstar has assigned to this set of outputs
    :param scale_factor: the scale factor to use in the nicely formatted name
    :param method: Whether to copy or move the files. It must be either "copy" or "move"
    """
    if method not in ["copy", "move"]:
        raise ValueError("Incorrect method in rename_halo_files()")

    old_halo_prefix = f"halos_{rockstar_idx}."
    new_halo_prefix = f"halos_a{scale_factor}."
    old_list = f"out_{rockstar_idx}.list"
    new_list = f"out_a{scale_factor}.list"
    for file in rockstar_dir.iterdir():
        if file.name.startswith(old_halo_prefix):
            new_name = file.name.replace(old_halo_prefix, new_halo_prefix)
            new_path = halos_dir / new_name
            if method == "move":
                file.replace(new_path)
            else:
                shutil.copy2(file, new_path)
        if file.name == old_list:
            new_path = halos_dir / new_list
            if method == "move":
                file.replace(new_path)
            else:
                shutil.copy2(file, new_path)

def shift_halo_output_index():
    """
    This moves the set of halo files with index 1 to index 0

    This lets them be the start of the next output
    """
    old_halo_prefix = "halos_1."
    new_halo_prefix = "halos_0."
    old_list = "out_1.list"
    new_list = "out_0.list"
    for file in rockstar_dir.iterdir():
        if file.name.startswith(old_halo_prefix):
            new_name = file.name.replace(old_halo_prefix, new_halo_prefix)
            new_path = rockstar_dir / new_name
            # check that this file isn't already there - it should never be!
            assert not new_path.is_file()
            file.replace(new_path)
        if file.name == old_list:
            new_path = rockstar_dir / new_list
            # check that this file isn't already there - it should never be!
            assert not new_path.is_file()
            file.replace(new_path)

def modify_parameter_file(scale_of_last_completed_output):
    """
    Change some values in the parameter file.

    Specifically, we change the restart number and the scale factor of the output. I'm
    not sure whether the scale factor change matters, but the restart number definitely
    does. In this case it should already be 1, since when we do it 2 at a time that's
    the index of the second one, but we double check.
    """
    restart_old = rockstar_dir / "restart.cfg"
    restart_new = rockstar_dir / "restart.cfg.temp"

    with open(restart_old, "r") as old:
        with open(restart_new, "w") as new:
            for line in old:
                split = line.split()
                key = split[0]
                value = split[-1]
                # we need to double check the restart number
                if key == "RESTART_SNAP":
                    assert value == "1"
                # and the scale factor
                if key == "SCALE_NOW":
                    # get rid of the "a" in the scale factor string we have above
                    line = line.replace(value, scale_of_last_completed_output)
                new.write(line)

    # then copy the file over
    restart_new.replace(restart_old)

def clean_up_rockstar_files():
    # Remove some of the files not needed to restart. They'll be regenerated
    # in the next loop of the rockstar finding
    for file in rockstar_dir.iterdir():
        if (
            file.name not in ["restart.cfg", "sentinel.txt"] and
            not file.name.startswith("halos_0.") and
            not (file.name.startswith("out_") and file.suffix == ".list") and
            not file.is_dir()
        ):
            file.unlink()
        elif file.name == "profiling" and file.is_dir():
            shutil.rmtree(file)

# ==============================================================================
#
# Do the main loop of moving files to the temporary directory for analysis
#
# ==============================================================================
# first get all the output files
art_files = sorted([file for file in sim_dir.iterdir() if file.suffix == ".art"])
# start by moving the first file there
move_all_simulation_files(art_files[0], sim_dir, temp_dir)
# Then loop through all the rest of the files
for art_file_idx_second in range(1, len(art_files)):
    # move this file to the temporary directory
    move_all_simulation_files(art_files[art_file_idx_second], sim_dir, temp_dir)
    # Do the halo analysis. I know using os in this way is "wrong", but yt and rockstar
    # have issues with MPI ranks when doing multiple runs of rockstar, as well as
    # memory issues when running rockstar on many outputs, so completely separating it
    # out will make sure that gets taken care of correctly.
    if machine == "stampede2":
        os.system(f"ibrun -n 47 -o 1 python ./halo_finding_rockstar.py {temp_dir} {rockstar_dir} 47")
    elif machine == "shangrila":
        os.system(f"mpiexec -np 4 python ./halo_finding_rockstar.py {temp_dir} {rockstar_dir} 4")

    # Then need to handle the existing files to prepare for the next iteration
    # First we get the oldest outputs out to the other directory. These have index 0
    # be definition
    first_scale = scale_factor_from_art_file(art_files[art_file_idx_second - 1])
    make_human_readable_halo_files(0, first_scale, "move")
    # Then move the ones with index 1 to be index 0, since those will be the starting
    # point for the next iteration of halo finding
    shift_halo_output_index()
    # then modify the parameter file
    modify_parameter_file(scale_factor_from_art_file(art_files[art_file_idx_second]))
    # and clean up the directory
    clean_up_rockstar_files()

    # then move the first file out. We'll move the next file in at the start
    # of the next loop
    move_all_simulation_files(art_files[art_file_idx_second - 1], temp_dir, sim_dir)

# ==============================================================================
#
# Final cleanup
#
# ==============================================================================
# copy the outputs of the last halo run to the human directory. Copy, not move, them
# so they can be the start of the next set of halo finding
make_human_readable_halo_files(0, scale_factor_from_art_file(art_files[-1]), "copy")

# move the files out of the temporary directory, then delete it
for out_file in temp_dir.iterdir():
    new_file = sim_dir / out_file.name
    out_file.replace(new_file)
# make sure temp_dir is empty before deleting it
assert len([item for item in temp_dir.iterdir()]) == 0
temp_dir.rmdir()

# update the sentinel file
(rockstar_dir / "sentinel.txt").touch()