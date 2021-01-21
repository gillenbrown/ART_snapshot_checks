import sys
from pathlib import Path
import os
import shutil
import yt
from tqdm import tqdm
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
temp_dir = rockstar_dir.parent / "temp_output_dir_for_halo_finding"
if not temp_dir.is_dir() and yt.is_root():
    temp_dir.mkdir()

# ==============================================================================
#
# Key functions to be used later
#
# ==============================================================================
def already_done(art_file):
    scale_factor = scale_factor_from_art_file(art_file)
    halo_file = f"halos_a{scale_factor}.0.bin"
    list_file = f"out_a{scale_factor}.list"
    all_output_files = [item.name for item in halos_dir.glob('*')]
    return halo_file in all_output_files and list_file in all_output_files

def move_all_simulation_files(art_file, old_dir, new_dir):
    for file in old_dir.iterdir():
        if file.stem == art_file.stem:
            file.replace(new_dir / file.name)

def scale_factor_from_art_file(art_file):
    return art_file.stem[-6:]

def art_file_from_scale_factor(directory, scale_factor):
    return directory / f"continuous_a{scale_factor}.art"

def move_restart_halo_catalogs(scale_factor):
    """
    This moves one set of halo catalogs from the directory with nicer names to the
    rockstar directory, with index 0, as we will restart from this.
    """
    new_halo_prefix = "halos_0."
    old_halo_prefix = f"halos_a{scale_factor}."
    new_list = "out_0.list"
    old_list = f"out_a{scale_factor}.list"
    for file in halos_dir.iterdir():
        if file.name.startswith(old_halo_prefix):
            new_name = file.name.replace(old_halo_prefix, new_halo_prefix)
            new_path = rockstar_dir / new_name
            file.replace(new_path)
        if file.name == old_list:
            new_path = rockstar_dir / new_list
            file.replace(new_path)

def make_human_readable_halo_files(rockstar_idx, scale_factor):
    """
    This moves one set of halo files to the halos directory with nicer names

    :param rockstar_idx: The index Rockstar has assigned to this set of outputs
    :param scale_factor: the scale factor to use in the nicely formatted name
    """
    old_halo_prefix = f"halos_{rockstar_idx}."
    new_halo_prefix = f"halos_a{scale_factor}."
    old_list = f"out_{rockstar_idx}.list"
    new_list = f"out_a{scale_factor}.list"
    for file in rockstar_dir.iterdir():
        if file.name.startswith(old_halo_prefix):
            new_name = file.name.replace(old_halo_prefix, new_halo_prefix)
            new_path = halos_dir / new_name
            file.replace(new_path)
        if file.name == old_list:
            new_path = halos_dir / new_list
            file.replace(new_path)

# I have this set up to use the restart file to tell whether or not we need to restart.
# I do this because without the restart file, it will NOT restart properly.
# This means I do need to be careful not to delete that restart file. But I do have a
# check that should guard against this if things go wrong.
def has_restart():
    restart_file_loc = rockstar_dir / "restart.cfg"
    return restart_file_loc.is_file()
    # n_halo_files = len(list(halos_dir.iterdir()))
    # return n_halo_files > 0

# def get_last_scale():
#     """
#     Get the scale factor of the last completed output.
#     """
#     max_scale_str = "-1"
#     for item in halos_dir.iterdir():
#         if item.name.startswith("out_") and item.name.endswith(".list"):
#             this_scale_str = item.name[5:11]
#             if float(this_scale_str) > float(max_scale_str):
#                 max_scale_str = this_scale_str
#
#     assert float(max_scale_str) > 0
#     return max_scale_str

def get_last_scale_from_parameter_file():
    """
    Get the scale factor of the last completed output.
    """
    restart_file_loc = rockstar_dir / "restart.cfg"
    assert restart_file_loc.is_file()

    with open(restart_file_loc, "r") as restart_file:
        for line in restart_file:
            split = line.split()
            key = split[0]
            value = split[-1]

            if key == "SCALE_NOW":
                return value

def modify_parameter_file(scale_of_last_completed_output, restarted):
    """
    Change some values in the parameter file.

    Specifically, we change the restart number and the scale factor of the output. I'm
    not sure whether the scale factor change matters, but the restart number definitely
    does. In this case it should already be 1, since when we do it 2 at a time that's
    the index of the second one, but we double check.

    Note that the only reason I need to modify the scale factor is to tell which output
    to use next. Normally the scale factor is a much longer string, but setting it
    manually to the same number of sig figs as the ART output file makes sure I can
    read it.
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
                    # if we started fresh, the value should be zero. We want to make
                    # sure it is 1, so we can restart from the next output next time
                    if not restarted:
                        assert line.strip() == "RESTART_SNAP = 0"
                        line = line.replace("RESTART_SNAP = 0", "RESTART_SNAP = 1")
                    else:
                        assert value == "1"
                # and the scale factor
                if key == "SCALE_NOW":
                    line = line.replace(value, scale_of_last_completed_output)
                new.write(line)

    # then copy the file over
    restart_new.replace(restart_old)

def clean_up_rockstar_files():
    # Remove some of the files not needed to restart. They'll be regenerated
    # in the next loop of the rockstar finding. Note that this removes any residual
    # halos files too, so be careful to move those out before using this. Note that I
    # do note remove the rockstar.cfg file because it's needed by the merger trees. It
    # gets overwritten each time the halo finder runs, so I don't need to worry about
    # updating it manually.
    for file in rockstar_dir.iterdir():
        if (
            file.name not in ["restart.cfg", "sentinel.txt", "rockstar.cfg"] and
            not file.is_dir()
        ):
            file.unlink()
        elif file.name == "profiling" and file.is_dir():
            shutil.rmtree(file)

def halo_finding_successfull(this_scale_factor, restarted):
    """
    Check if the halo finding worked. This is determined by checking the number of lines
    in the outputs file. If it has more than the lines in the header, we say it worked
    """
    # if it's early, there may not be halos, so let that pass.
    if this_scale_factor < 0.05:
        return True
    # see which rockstar index to check
    if restarted:
        rockstar_idx = 1
    else:
        rockstar_idx = 0
    # then count the number of lines in that file.
    halo_file_loc = rockstar_dir / f"halos_{rockstar_idx}.0.ascii"
    with open(halo_file_loc, "r") as halo_file:
        n_lines = len([1 for line in halo_file if len(line.strip()) > 0])
    return n_lines > 20  # there are 20 header lines

# ==============================================================================
#
# Do the main loop of moving files to the temporary directory for analysis
#
# ==============================================================================
# first clean up any leftover files in the rockstar directory. This includes the
# out*.list files that are made at the end of this script.
clean_up_rockstar_files()
# get all the output files, but exclude the restart files (which are just
# 'continuous_0.art'
all_art_files = [file for file in sim_dir.iterdir()
                 if file.name.startswith('continuous_a') and file.suffix == ".art"]
# see which ones need halos to be made for
art_files = sorted([art_file for art_file in all_art_files
                    if not already_done(art_file)])

# I want to do a check that the files were made correctly (since sometimes I find bugs
# on shangrila that result in zero halos being found). I'll redo a set of halos if it
# doesn't work. To do this, I'll use an index, then iterate it if the halo finding is
# successfull.
idx = 0
while idx < len(art_files):
    art_file = art_files[idx]
    # first check whether we have a file to restart from. If we do, we'll restart from
    # it. Otherwise, we'll do a run with one simulation output to create a restart file
    # which we can then build off of.
    restart = has_restart()

    # first compare the scale factor to the scale factor given in the config file. The
    # scale factor of this item should be larger than the previously done one.
    this_scale_str = scale_factor_from_art_file(art_file)
    if restart:
        restart_scale_str = get_last_scale_from_parameter_file()
        if float(this_scale_str) <= float(restart_scale_str):
            out_str = "Error in halo finding: " \
                      f"Last done output is a={restart_scale_str}, " \
                      f"yet we're trying to do a={this_scale_str} now"
            raise RuntimeError(out_str)

    # If there's no restart, make sure the scale factor for this is reasonable. This
    # is just a double check to make sure it's not accidentally missing
    if not restart and float(this_scale_str) > 0.15:
        raise RuntimeError(f"No restart file found for a={this_scale_str}, are you sure?")

    # print what we're about to do
    if restart:
        print("\n" + "=" * 80 + "\n" +
              f"Restarting from a={restart_scale_str} to do a={this_scale_str}\n" +
              f"{str(sim_dir)}\n" +
              "=" * 80 + "\n")
    else:
        print("\n" + "=" * 80 + "\n" +
              f"Starting fresh from a={this_scale_str}\n" +
              f"{str(sim_dir)}\n" +
              "=" * 80 + "\n")

    # Then move the output files corresponding to these scales to the temp directory.
    move_all_simulation_files(art_file, sim_dir, temp_dir)
    # I need to get the appropriate ART restart files. If there's not a restart file, we
    # don't do this
    if restart:
        restart_file = art_file_from_scale_factor(sim_dir, restart_scale_str)
        if not restart_file.is_file():
            raise RuntimeError(f"ART output for a={restart_scale_str} not found in {sim_dir}")
        move_all_simulation_files(restart_file, sim_dir, temp_dir)


    # Do the same thing with the rockstar output files. I only need to move the set
    # we're restarting from, as of course the second one hasn't been made yet. We only
    # do this if we're restarting
    if restart:
        move_restart_halo_catalogs(restart_scale_str)

    # Do the halo analysis. I know using os in this way is "wrong", but yt and rockstar
    # have issues with MPI ranks when doing multiple runs of rockstar, as well as
    # memory issues when running rockstar on many outputs, so completely separating it
    # out will make sure that gets taken care of correctly.
    if machine == "stampede2":
        os.system(f"ibrun -n 47 -o 1 python ./halo_finding_rockstar.py {temp_dir} {rockstar_dir} 47")
    elif machine == "shangrila":
        os.system(f"mpiexec -np 4 python ./halo_finding_rockstar.py {temp_dir} {rockstar_dir} 4")

    # if successfull, continue on. Otherwise, delete the halo catalogs and try again.
    if halo_finding_successfull(float(this_scale_str), restart):
        # move the resulting halo catalogs to the directory for human readable files,
        # formatting them nicely as we do so. The indices are different if we restart, so
        # we need to be careful.
        if restart:
            make_human_readable_halo_files(0, restart_scale_str)
            make_human_readable_halo_files(1, this_scale_str)
        else:
            make_human_readable_halo_files(0, this_scale_str)

        # reset the rockstar restart file
        modify_parameter_file(this_scale_str, restart)
        # move to next set of halos
        idx += 1
    else:
        print("\n" + "=" * 80 + "\n" +
              f"Halo catalogs for a={this_scale_str} failed! Trying again...\n" +
              f"{str(sim_dir)}\n" +
              "=" * 80 + "\n")
        # if we are restarting, move back the restart halo catalog, which was
        # successful. When we clean up the rockstar files below, it will remove the
        # unsuccessfull halo catalogs that are left in this directory
        if restart:
            make_human_readable_halo_files(0, restart_scale_str)
            # we modify the parameter file to go back to the last output
            modify_parameter_file(restart_scale_str, restart)
        else:
            # remove the restart file
            restart_file_loc = rockstar_dir / "restart.cfg"
            restart_file_loc.unlink()

    # and clean out other rockstar files, including any failed halo catalogs
    clean_up_rockstar_files()

    # move simulation output files back to where they belong
    move_all_simulation_files(art_file, temp_dir, sim_dir)
    if restart:
        move_all_simulation_files(restart_file, temp_dir, sim_dir)

# ==============================================================================
#
# Formatting out_*.list files for consistent trees
#
# ==============================================================================
# I need to move all the out_*.list files to the rockstar directory, but without their
# nice scale factor names, so that consistent trees knows what to do with them. Note
# that I delete these files at the beginning of the script so they don't interfere with
# the actual creation of the files
out_files = sorted([halo_file for halo_file in halos_dir.iterdir()
                   if halo_file.name.startswith("out_")
                   and halo_file.name.endswith(".list")])

for idx, out_file in tqdm(enumerate(out_files)):
    new_out = f"out_{idx}.list"
    new_path = rockstar_dir / new_out
    shutil.copy2(out_file, new_path)

# I also need to modify the rockstar.cfg file to have the correct number of output
# files, as this is what tells it how many to include in the tree.
config_old = rockstar_dir / "rockstar.cfg"
config_new = rockstar_dir / "rockstar.cfg.temp"

with open(config_old, "r") as old:
    with open(config_new, "w") as new:
        for line in old:
            split = line.split()
            key = split[0]
            value = split[-1]

            if key == "NUM_SNAPS":
                line = line.replace(value, str(len(out_files)))
            new.write(line)

# then copy the file over
config_new.replace(config_old)

# ==============================================================================
#
# Final cleanup
#
# ==============================================================================
# make sure temp_dir is empty before deleting it
assert len([item for item in temp_dir.iterdir()]) == 0
temp_dir.rmdir()

# update the sentinel file
(rockstar_dir / "sentinel.txt").touch()
