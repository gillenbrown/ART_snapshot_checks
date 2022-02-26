"""
write_gas_pdf.py

Writes out the gas pdf.

Takes 2 required parameters.
1 - Location of the output file
2 - Location of the simulation output. Can be relative to the working directory
    where this code was called from, or an absolute path.

"""
from utils import load_galaxies as lg

import sys
from pathlib import Path

import yt
from analysis_functions import gas_pdfs

yt.funcs.mylog.setLevel(50)  # ignore yt's output

# =========================================================================
#
# load sims and prep output file
#
# =========================================================================
out_file = open(sys.argv[1], "w")
ds_loc = Path(sys.argv[2]).resolve()
sim = lg.Simulation(ds_loc, sphere_radius_kpc=10, min_virial=False)

# define a shorthand function that handles the parameters to out,
# so we don't have to duplicate that each time
def out(info):
    out_file.write(info + "\n")


# print header for this output file
out("# This output contains the cell properties for all cells with nonzero H2: ")
out("# - Cell H2 number density [cm^-3]")
out("# - Cell H2 mass [Msun]")
out(f"# This simulation has {sim.n_galaxies} galaxies.")

# =========================================================================
#
# Analysis
#
# =========================================================================
for galaxy in sim.galaxies:
    n_h2 = gas_pdfs.get_h2_number_density(galaxy.sphere)
    m_h2 = gas_pdfs.get_h2_mass(galaxy.sphere)

    # clean up units
    n_h2 = n_h2.to("cm**(-3)").value
    m_h2 = m_h2.to("Msun").value

    # get only the cells that have significant amounts of H2. Note that this threshold
    # was empirically determined by looking at pdfs and picking a value at which the
    # cumulative pdf is essentially zero
    good_idxs = n_h2 > 1e-3
    n_h2 = n_h2[good_idxs]
    m_h2 = m_h2[good_idxs]

    # then print info
    for i in range(len(n_h2)):
        out(f"{n_h2[i]:.2e} {m_h2[i]:.2e} ")


out_file.close()
