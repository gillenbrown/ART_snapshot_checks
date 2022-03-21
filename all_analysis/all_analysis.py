import sys
from pathlib import Path
import numpy as np
import yt

# add path to sys to import utils
sys.path.append(str(Path("..").resolve() / "utils"))
sys.path.append(str(Path("..").resolve() / "analysis_functions"))
import load_galaxies
import gas_pdfs

yt.funcs.mylog.setLevel(50)  # ignore yt's output

output_loc = sys.argv[1]


def write(output):
    print(output, flush=True)


write(output_loc)
# ======================================================================================
#
# Then actually doing things
#
# ======================================================================================
try:
    sim = load_galaxies.Simulation(Path(output_loc), sphere_radius_kpc=30, n_galaxies=1)
except:
    write(f"Error in creation! {output_loc}")
    exit()

# then do analysis
try:
    ad = sim.ds.all_data()
except:
    write(f"Error in ad! {output_loc}")
    exit()

try:
    densities = gas_pdfs.get_gas_number_density(ad).to("cm**(-3)").value
    write(f"max number density = {np.max(densities):.2e} cm^-3 {output_loc}")
except:
    write(f"Error in densities! {output_loc}")
