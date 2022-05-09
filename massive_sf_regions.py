"""
massive_sf_regions.py

Create a text file with the most massive SF region in each output
"""
from pathlib import Path

import sys

import numpy as np
import betterplotlib as bpl

from utils import plot_utils

bpl.set_style()

out_file_loc = Path(sys.argv[1])
sim_dirs = [Path(d) for d in sys.argv[2:]]

# =========================================================================
#
# Go through the files
#
# =========================================================================
def get_scale(filename):
    """ Filename is assumed to be for a forming_clusters output """
    if not filename.startswith("forming_clusters"):
        raise RuntimeError("Should not happen")
    return filename[-10:-4]


def find_max_in_column(file_location, col_idx):
    """ Returns the index of the row with the highest value in the desired column """
    max_value = -np.inf
    max_idx = None
    with open(file_location, "r") as in_file:
        for idx, line in enumerate(in_file):
            if line.startswith("#") or len(line.strip()) == 0:
                continue
            # try to get this index, but check to make sure the line is long enough
            try:
                if float(line[col_idx]) > float(max_value):
                    max_value = line[col_idx]
                    max_idx = idx
            except IndexError:
                return None
    return max_idx


def get_info_from_given_line(file_location, row_idx):
    """ Get the position, mass, SFR, and average age of the cluster complex desired """
    with open(file_location, "r") as in_file:
        for idx, line in enumerate(in_file):
            if idx == row_idx:
                split_line = line.split()
                return [
                    split_line[1],  # x
                    split_line[2],  # y
                    split_line[3],  # z
                    split_line[9],  # total mass within 30pc
                    split_line[10],  # star formation rate
                    split_line[11],  # average age
                ]


out_file = open(out_file_loc, "w")
# define a shorthand function that handles the parameters to out,
# so we don't have to duplicate that each time
def out(info):
    out_file.write(info + "\n")


# print header for this output file
out("# This catalog contains the most massive star forming regions in each snapshot.")
out("# The columns have the following info: ")
out("# - Simulation Name")
out("# - Scale Factor")
out("# - Position x [kpc]")
out("# - Position y [kpc]")
out("# - Position z [kpc]")
out("# - Mass of clusters younger than 15 Myr within 30 pc of this location [Msun]")
out("# - Star formation rate within 30 pc of this cluster [Msun/year]")
out("# - Ave age of clusters younger than 15 Myr within 30 pc of this cluster [Myr]")

for dir in sim_dirs:
    check_files = sorted(
        [
            item
            for item in (dir / "checks").iterdir()
            if item.name.startswith("forming_clusters")
        ]
    )

    for file in check_files:
        # for each file, print the cluster complex with the highest mass and highest
        # SFR. Those are columns 9 and 10
        for col_idx in [9, 10]:
            idx = find_max_in_column(file, col_idx)
            if idx is not None:  # if the value was found
                x, y, z, mass, sfr, age = get_info_from_given_line(file, idx)
                # note that these quantities are strings
                out(
                    f"{plot_utils.get_sim_dirname(dir)} "
                    f"{get_scale(file.name)} "
                    f"{x} "
                    f"{y} "
                    f"{z} "
                    f"{mass} "
                    f"{sfr} "
                    f"{age} "
                )

out_file.close()
