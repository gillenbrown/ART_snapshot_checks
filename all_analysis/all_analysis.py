import sys
from pathlib import Path
import gc
import numpy as np
import yt

# add path to sys to import utils
sys.path.append(str(Path("..").resolve() / "utils"))
sys.path.append(str(Path("..").resolve() / "analysis_functions"))
import load_galaxies
import gas_pdfs

yt.funcs.mylog.setLevel(50)  # ignore yt's output

# out_file = open("./molecular_gas.txt", "w")

scratch_dir = Path("/scratch/06912/tg862118/art_runs/analysis/production")

# directories = [
#     scratch_dir / "rj_sfe010_hn20",
#     scratch_dir / "rj_sfe100_hn20",
#     scratch_dir / "tl_sfe001_hn20",
#     scratch_dir / "tl_sfe010_hn20",
#     scratch_dir / "tl_sfe100_hn00",
#     scratch_dir / "tl_sfe100_hn00_fboost3",
#     scratch_dir / "tl_sfe100_hn05",
#     scratch_dir / "tl_sfe100_hn20",
# ]
directories = [scratch_dir / d for d in sys.argv[1:]]

# ======================================================================================
#
# analysis functions
#
# ======================================================================================
def write(out_str):
    # I got rid of writing to the file for now, since when parallelizing with launcher
    # I don't want the output of different runs interleaved. Just printing to stdout is
    # fine with me. Note that if you do want to write to a file, add a newline to the
    # end of the string to print.
    # out_file.write(out_str)
    # out_file.flush()  # this and next line ensure writing happens now
    # os.fsync(out_file.fileno())
    print(out_str, flush=True)


# ======================================================================================
#
# Then actually doing things
#
# ======================================================================================
for d in directories:
    out_dir = d / "run" / "out"
    art_files = sorted([f for f in out_dir.iterdir() if f.suffix == ".art"])
    run_name = d.name
    write(f"{run_name}\n")

    # then parse this to get the halo files
    for output in art_files:
        gc.collect()
        try:
            sim = load_galaxies.Simulation(output, sphere_radius_kpc=30, n_galaxies=1)
        except:
            write(f"Error in creation! {output.name}")
            gc.collect()
            continue
        # then do analysis
        try:
            ad = sim.ds.all_data()
        except:
            write(f"{sim.z:.4f} Error in ad!")
            gc.collect()
            continue

        try:
            densities = gas_pdfs.get_gas_number_density(ad).to("cm**(-3)").value
            write(f"z={sim.z:.4f} max number density = {np.max(densities):.2e} cm^-3")
            del densities
        except:
            write(f"{sim.z:.4f} Error in densities!")

        del ad
        del sim
        gc.collect()


# out_file.close()
