"""
print_galaxy_masses.py

Makes a simple text file with the log stellar mass of the main galaxy in all
simulations. This is for easy copying into my spreadsheet

This script takes the output file name, then the list of all summary files.
"""

import sys
from pathlib import Path
import numpy as np
import betterplotlib as bpl

bpl.set_style()

output_file = Path(sys.argv[1])

# ======================================================================================
#
# Find the last output for each simulation
#
# ======================================================================================
last_summary = dict()
for summary_path in sys.argv[2:]:
    s_path = Path(summary_path)
    home_dir = s_path.parent.parent  # summary dir, then that parent
    # may have one further directory to go up if I used the run dir
    if home_dir.name == "run":
        home_dir = home_dir.parent

    sim_name = home_dir.name
    try:
        if s_path.name > last_summary[sim_name].name:
            last_summary[sim_name] = s_path
    except KeyError:
        # this simulation is not present.
        last_summary[sim_name] = s_path

# ======================================================================================
#
# Then read these summaries and output the data
#
# ======================================================================================
def read_summary_file(file_loc, n_to_read):
    m_star = []
    with open(file_loc, "r") as summary:
        for line in summary:
            # the galaxies are in order of rank in the file
            if line.startswith("Stellar Mass within 30 kpc"):
                m_star.append(float(line.split()[-2]))

            if len(m_star) == n_to_read:
                return m_star


with open(output_file, "w") as out:
    for sim_name, summary_path in last_summary.items():
        if "/production/" in str(summary_path):
            n_to_read = 2
        else:
            n_to_read = 1

        this_stellar_masses = read_summary_file(summary_path, n_to_read)
        # then format this to be written
        out_str = f"{sim_name:<50}"
        for m in this_stellar_masses:
            out_str += f" {np.log10(m):>5.2f}"

        out.write(out_str + "\n")
