"""
massive_sf_regions.py

Create a text file with the most massive SF region in each output
"""
from pathlib import Path

import sys

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


def parse_forming_clusters(file_location):
    max_mass = 0
    x, y, z = 0, 0, 0
    with open(file_location, "r") as in_file:
        for line in in_file:
            if line.startswith("#") or len(line.strip()) == 0:
                continue

            # we have a data line. There is an old version of this file without the
            # appropriate info that we need to check against
            split = line.split()
            if len(split) == 6:
                return 0, 0, 0, 0
            elif len(split) != 10:
                raise RuntimeError("Shouldn't happen.")

            # then get the info and check
            _, this_x, this_y, this_z, _, _, _, _, _, this_mass = split
            if float(this_mass) > float(max_mass):
                max_mass = this_mass
                x = this_x
                y = this_y
                z = this_z

    return x, y, z, max_mass


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

for dir in sim_dirs:
    check_files = sorted(
        [
            item
            for item in (dir / "checks").iterdir()
            if item.name.startswith("forming_clusters")
        ]
    )

    for file in check_files:
        x, y, z, mass = parse_forming_clusters(file)
        # note that these quantities are strings
        out(
            f"{plot_utils.get_sim_dirname(dir)} "
            f"{get_scale(file.name)} "
            f"{x} "
            f"{y} "
            f"{z} "
            f"{mass} "
        )

out_file.close()
