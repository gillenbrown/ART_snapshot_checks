# this script takes one parameter: The name of the desired (renamed) halo
# catalog. For example, this could be ~/dirname/halos_a0.0078.0.bin. All
# that matters is the scale factor contained therein, which must be in the
# 12th to last to 7th to last characters.

# This script finds the original rockstar halo catalog (which has unhelpful
# names like halos_0.0.bin) and copies it to the desired location with the
# helpful name with the scale factor in it. It does this for all of the original
# rockstar outputs at the same scale factor (so both ascii and binary, plus
# all the individual files for one output written by different Rockstar threads)

import sys
from pathlib import Path
import shutil

# get the representative halo catalog, and parse it to get the directory to
# store our catalogs and the directory of the original rockstar outputs.
halo_output_file = Path(sys.argv[1]).resolve()
halo_output_dir = halo_output_file.parent
rockstar_dir = Path(str(halo_output_dir).replace("/halos", "/rockstar_halos"))

# get the scale factor of the desired halo catalog
target_a_str = str(halo_output_file)[-12:-6]

# create function to check a given halo catalog against the expected value
def file_matches_scale_factor(input_file, target_a_str):
    with open(rockstar_dir / file, "r") as halo:
        # first line is not needed
        halo.readline()
        # second line has scale factor
        this_a_str = halo.readline().split()[-1]

        # python format automatically rounds so we can compare it to the
        # scale factor in the target str
        this_a_str = "{:.4f}".format(float(this_a_str))

        return this_a_str == target_a_str


# go through all the rockstar outputs to find the desired one. In the end we'll
# store the prefix for those files, so for now we'll set a sentinel value
# telling us we didn't find it.
old_prefix = "not found"
for file in rockstar_dir.iterdir():
    # Each output has multiple halo files. Get only one of them to read in
    # and test. The ASCII files have the scale factor right in them
    if "halos_" in file.name and ".0.ascii" in file.name:
        if file_matches_scale_factor(file.name, target_a_str):
            # get everything up to the period. This will contain the file
            # index that matches what we want. We do want the period, though, so
            # that we don't include halos_80.0.ascii when we found
            # halos_8.0.ascii as the match
            old_prefix = file.name.split(".")[0] + "."
            break

# there is one pathological case where the rounding breaks.
# I think this is because the halo output file already reports a
# rounded value, and rounding that can case inaccuracies
# i.e. a=0.20004999 should be rounded to a=0.2000, but the halo file
# rounds to more digits, reports a=0.200050, which then gets rounded
# to a=0.2001. There's also the 0.5 rounding to even thing that
# I'm not sure if that happens in ART, so I'll check one digit on
# either side to see if those match too
if old_prefix == "not found":
    print(" - {} is a weird case".format(halo_output_file))
    temp_target_str_up = "{:.4f}".format(float(target_a_str) + 0.0001)
    temp_target_str_down = "{:.4f}".format(float(target_a_str) - 0.0001)
    # then do the same looking we did before
    for file in rockstar_dir.iterdir():
        if "halos_" in file.name and ".0.ascii" in file.name:
            if file_matches_scale_factor(
                file.name, temp_target_str_up
            ) or file_matches_scale_factor(file.name, temp_target_str_down):

                old_prefix = file.name.split(".")[0] + "."
                break

# if we still can't find it, exit
if old_prefix == "not found":
    raise ValueError(
        "Rockstar halo file for\n{}\nwas not found in "
        "directory\n{}".format(halo_output_file, rockstar_dir)
    )

# now we have the right index and can copy all the files we need
old_halo_prefix = old_prefix
new_halo_prefix = "halos_a{}.".format(target_a_str)
old_list = old_prefix.replace("halos_", "out_") + "list"
new_list = "out_a{}.list".format(target_a_str)
for file in rockstar_dir.iterdir():
    old_path = rockstar_dir / file
    if file.name.startswith(old_halo_prefix):
        new_name = file.name.replace(old_halo_prefix, new_halo_prefix)
        new_path = halo_output_dir / new_name
        shutil.copy2(old_path, new_path)
        # copy2 keeps modification time, touch it
        new_path.touch()
    if file.name == old_list:
        new_path = halo_output_dir / new_list
        shutil.copy2(old_path, new_path)
        # copy2 keeps modification time, touch it
        new_path.touch()
