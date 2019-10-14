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
import os
import shutil

# get the representative halo catalog, and parse it to get the directory to 
# store our catalogs and the directory of the original rockstar outputs.
halo_output_file = os.path.abspath(sys.argv[1])
halo_output_dir = os.path.dirname(halo_output_file)
rockstar_dir = halo_output_dir.replace("/halos", "/rockstar_halos")
if not rockstar_dir.endswith(os.sep):
    rockstar_dir += os.sep
if not halo_output_dir.endswith(os.sep):
    halo_output_dir += os.sep

# get the scale factor of the desired halo catalog
target_a_str = halo_output_file[-12:-6]

# go through all the rockstar outputs to find the desired one. In the end we'll 
# store the prefix for those files, so for now we'll set a sentinel value 
# telling us we didn't find it.
old_prefix = "not found"
for file in os.listdir(rockstar_dir):
    # Each output has multiple halo files. Get only one of them to read in
    # and test. The ASCII files have the scale factor right in them
    if "halos_" in file and ".0.ascii" in file:
        with open(rockstar_dir + file, "r") as halo:
            # first line is not needed
            halo.readline()
            # second line has scale factor
            this_a_str = halo.readline().split()[-1]

            # python format automatically rounds so we can compare it to the
            # scale factor in the target str
            # there is one pathological case where the rounding breaks 
            # because of floating point inaccuracies
            if this_a_str == "0.068650":
                this_a_str = "0.0686"
            else:
                this_a_str = "{:.4f}".format(float(this_a_str))

            if this_a_str == target_a_str:
                # get everything up to the period. This will contain the file
                # index that matches what we want
                old_prefix = file.split(".")[0]
                break

# if we didn't find it, exit
if old_prefix == "not found":
    raise ValueError("Rockstar halo file for\n{}\nwas not found in "
                     "directory\n{}".format(halo_output_file, rockstar_dir))

# now we have the right index and can copy all the files we need
new_prefix = "halos_a{}".format(target_a_str)
for file in os.listdir(rockstar_dir):
    if file.startswith(old_prefix):
        old_path = rockstar_dir + os.sep + file
        new_name = file.replace(old_prefix, new_prefix)
        new_path = halo_output_dir + os.sep + new_name
        shutil.copy(old_path, new_path)
