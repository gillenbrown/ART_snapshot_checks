""""
print_sf_levels.py

Print the distribution of levels for cells that meet the criteria for star formation
"""
import sys
from pathlib import Path

import numpy as np

from utils import load_galaxies, plot_utils
from analysis_functions import gas_pdfs

out_file_path = Path(sys.argv[1])
out_file = open(out_file_path, "w")

# ======================================================================================
#
# load the simulations
#
# ======================================================================================
sims_last = load_galaxies.get_simulations_last(sys.argv[2:])
common_z = 4.0
sims_common = load_galaxies.get_simulations_same_scale(sys.argv[2:], common_z)

# ======================================================================================
#
# analysis functions
#
# ======================================================================================
def find_sf_cells(sim):
    density = sim.func_all_galaxies(
        lambda g: gas_pdfs.get_gas_number_density(g).to("cm**(-3)").value
    )
    temp = sim.func_all_galaxies(lambda g: gas_pdfs.get_gas_temp(g).to("K").value)
    alpha = sim.func_all_galaxies(gas_pdfs.get_gas_virial_criterion)
    f_h2 = sim.func_all_galaxies(gas_pdfs.get_h2_frac)

    idx_1 = density > 1000
    idx_2 = temp < 1e4
    idx_3 = alpha < 10
    idx_4 = f_h2 > 0.5

    return np.logical_and.reduce((idx_1, idx_2, idx_3, idx_4))


def print_level(sim):
    out_file.write(f"{plot_utils.get_sim_dirname(sim.run_dir)}\n")
    out_file.write(f"z={sim.z:.4f}\n")

    idx_sf = find_sf_cells(sim)

    levels = sim.func_all_galaxies(gas_pdfs.get_gas_level)
    sizes = sim.func_all_galaxies(lambda g: gas_pdfs.get_cell_size_pc(g).to("pc").value)

    # match levels to cell size
    level_size_dict = {l: s for l, s in zip(np.unique(levels), np.unique(sizes)[::-1])}

    # then count the number of cells at the levels that have sf
    good_levels = levels[idx_sf]

    uniques, counts = np.unique(good_levels, return_counts=True)
    for idx in range(len(uniques)):
        out_file.write(
            f"level {uniques[idx]:>2.0f} = "
            f"{level_size_dict[uniques[idx]]:>5.2f} pc: "
            f"{counts[idx]:>5.0f} cells \n"
        )

    out_file.write("\n")


# ======================================================================================
#
# actual running
#
# ======================================================================================
for idx in range(len(sims_last)):
    # print the same run at different redshifts together
    print_level(sims_last[idx])
    print_level(sims_common[idx])


out_file.close()
