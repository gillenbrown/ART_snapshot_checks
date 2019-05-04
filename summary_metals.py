"""
global_properties.py

Reports some properties related to the metals in the simulation

Takes 1 required and 2 optional parameters.
1 - Location of the simulation output. Can be relative to the working directory
    where this code was called from, or an absolute path.
2 - Optional argument. Must be "clobber" if included. This will bypass the 
    check whether we overwrite a previously existing output file.
3 - Optional argument. Must be "silent" if included. Will print info and 
    write it to the file if this is not included. Will only write to file
    if this is included.
"""

import sys
import os

import yt
import numpy as np

import betterplotlib as bpl

yt.funcs.mylog.setLevel(50)  # ignore yt's output
bpl.presentation_style()
bpl.presentation_style()  # for some reason this needs to be there twice

# Check that the third argument is correct
if len(sys.argv) == 3 and sys.argv[2] not in ["clobber", "silent"]:
    raise ValueError("Argument 2 (if used), must be 'clobber' or 'silent'")
if "silent" in sys.argv:
    silent = True
else:
    silent = False

def print_and_write(info, file_obj):
    if file_obj is not None:
        file_obj.write(info + "\n")
    if not silent:
        print(info)

ds_loc = os.path.abspath(sys.argv[1])
scale_factor = ds_loc[-10:-4]
ds = yt.load(ds_loc)
ad = ds.all_data()

# get the location of where to write the file.
sim_dir = os.path.dirname(ds_loc) + os.sep
file_dir = sim_dir.replace("/out/", "/checks/")
file_path = file_dir + "summary_metals_a" + scale_factor + ".txt"
plots_dir = sim_dir.replace("/out/", "/plots/")

print_and_write("Output being written to:", None)
print_and_write(file_path, None)

# see if there is an existing file here that we don't want to replace.
if "clobber" not in sys.argv:
    if os.path.isfile(file_path):
        good = False
        while not good:
            choice = input("This output file already exists. "
                           "Overwrite? (y/n): ")
            if choice.lower() in ["y", "n"]:
                good = True
            else:
                print("Choose 'y' or 'n'")

        if choice.lower() == "n":
            # don't overwrite, so the code ends
            print("File won't be overwritten, so nothing will be done.")
            exit()
        # Don't need to do anything here inside this block if the user does 
        # want to overwrite, since the file created later will automatically 
        # overwrite any existing file.
# open the file
out_file = open(file_path, "w")

# =============================================================================
#         
# Basic information
# 
# =============================================================================
print_and_write("", None)  # no write since this is just for console formatting
print_and_write("Simulation location:", out_file)
print_and_write(ds_loc, out_file)

a = ds.scale_factor
z = 1.0 / a - 1.0
print_and_write("\na = {:.4f}".format(a), out_file)
print_and_write("z = {:.4f}".format(z), out_file)

# =============================================================================
#         
# Gas cells
# 
# =============================================================================
# See some statistics about the metal properties of the gas at various levels


# =========================================================================
#         
# Plots
# 
# =========================================================================


print_and_write("\nPlots will be saved to:", out_file)
print_and_write("Replace with name of plot", out_file)

out_file.close()
