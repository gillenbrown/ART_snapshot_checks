import sys
from pathlib import Path
import os
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

    # Then need to handle the existing files to prepare for the next iteration
    # handle_halo_files()

    # then move the first file out. We'll move the next file in at the start
    # of the next loop
    move_all_simulation_files(art_files[art_file_idx_second - 1], temp_dir, sim_dir)

# ==============================================================================
#
# Final cleanup
#
# ==============================================================================
# move the files out of the temporary directory, then delete it
for out_file in temp_dir.iterdir():
    new_file = sim_dir / out_file.name
    out_file.replace(new_file)
temp_dir.rmdir()

#     # update the sentinel file
#     (rockstar_dir / "sentinel.txt").touch()