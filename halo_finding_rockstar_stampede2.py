import sys
import os
import gc
from pathlib import Path
import yt
from yt.extensions.astro_analysis.halo_finding.rockstar.api import RockstarHaloFinder
yt.funcs.mylog.setLevel(0)  # ignore yt's output
yt.enable_parallelism()

# format of sys.argv:
# idx 0: name of script
# ids 1: directory of simulation outputs to run halos on
# idx 2: directory to store halo finding outputs in

# print function that always flushes
import functools
print = functools.partial(print, flush=True)

# do error checking on root only.
if yt.is_root():
    if len(sys.argv) != 4:
        raise ValueError("Please provide the proper command line argument.")

# turn the directories the user passes into the absolute path
sim_dir = Path(sys.argv[1]).resolve()
rockstar_dir = Path(sys.argv[2]).resolve()
cores_to_use = int(sys.argv[3])
# Create the subdirectory where we'll do the actual halo finding
temp_dir = sim_dir / "temp_output_dir_for_halo_finding"
if not temp_dir.is_dir() and yt.is_root():
    temp_dir.mkdir()

if yt.is_root():
    print("Reading simulations from: {}".format(sim_dir))
    print("Writing halo catalogs to: {}".format(rockstar_dir))

# determine how many readers and writers to use, depending on what machine we're on. We
# have one master process, plus the number of readers and writers. Reading is quicker
# then analysis, so we typically have fewer readers
processing_cores = cores_to_use - 1
if processing_cores < 5:
    readers = 1
elif processing_cores > 30:
    readers = 16
else:
    raise NotImplementedError("Choose writers for 5 < cores < 30")
writers = processing_cores - readers
if yt.is_root():
    print("Rockstar will be running with:")
    print(f"\t- {cores_to_use} total cores")
    print(f"\t- 1 master process")
    print(f"\t- {readers} readers")
    print(f"\t- {writers} writers")

def move_all_simulation_files(art_file, old_dir, new_dir):
    if yt.is_root():
        for file in old_dir.iterdir():
            if file.stem == art_file.stem:
                file.replace(new_dir / file.name)

# ==============================================================================
#
# Key functions to be used later
#
# ==============================================================================
def rockstar_iteration():
    """
    Run rockstar on two outputs in the directory. This is the function

    We do this in a function to make sure that there are not many datasets loaded at
    once, which will overload the memory
    """
    # check to see if there is a currently existing halo catalog already here
    # to restart from.
    if (rockstar_dir / "restart.cfg").is_file():
        restart = True
    else:
        restart = False
    print("=== able to restart:", restart)

    ts = yt.load(str(temp_dir) + '/continuous_a?.????.art')

    # check what kind of particles are present
    if ('N-BODY_0', 'MASS') in ts[0].derived_field_list:
        particle_type = "N-BODY_0"
    else:
        particle_type = "N-BODY"

    rh = RockstarHaloFinder(ts, num_readers=readers, num_writers=writers, outbase=out_dir,
                            particle_type=particle_type)
    rh.run(restart=restart)

    # then make sure to clean up the memory. We do this explicitly just to be sure
    del ts
    gc.collect()

# ==============================================================================
#
# Do the main loop of moving files to the temporary directory for analysis
#
# ==============================================================================
def print_temp_dir():
    if yt.is_root():
        print("\n===\nState of temp dir:")
        for file in temp_dir.iterdir():
            print(file)
        print("===")

# first get all the output files
art_files = sorted([file for file in sim_dir.iterdir() if file.suffix == ".art"])
if yt.is_root():
    print(art_files)
# start by moving the first file there
move_all_simulation_files(art_files[0], sim_dir, temp_dir)
# Then loop through all the rest of the files
for art_file_idx_second in range(1, len(art_files)):
    if yt.is_root():
        print("\n", art_file_idx_second, art_files[art_file_idx_second], "\n")
    # move this file to the temporary directory
    move_all_simulation_files(art_files[art_file_idx_second], sim_dir, temp_dir)
    print_temp_dir()
    # Do the halo analysis
    rockstar_iteration()
    # Then need to handle the existing files to prepare for the next iteration
    # handle_halo_files()

    # then move the first file out. We'll move the next file in at the start
    # of the next loop
    move_all_simulation_files(art_files[art_file_idx_second - 1], temp_dir, sim_dir)
    print_temp_dir()

# ==============================================================================
#
# Final cleanup
#
# ==============================================================================
# move the files out of the temporary directory, then delete it
if yt.is_root():
    for out_file in temp_dir.iterdir():
        new_file = sim_dir / out_file.name
        out_file.replace(new_file)
    temp_dir.rmdir()

#     # update the sentinel file
#     (rockstar_dir / "sentinel.txt").touch()
