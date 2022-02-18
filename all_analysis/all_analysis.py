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


def find_sf_cells(ad):
    density = gas_pdfs.get_gas_number_density(ad).to("cm**(-3)").value
    temp = gas_pdfs.get_gas_temp(ad).to("K").value
    alpha = gas_pdfs.get_gas_virial_criterion(ad)
    f_h2 = gas_pdfs.get_h2_frac(ad)

    idx_1 = density > 1000
    idx_2 = temp < 1e4
    idx_3 = alpha < 10
    idx_4 = f_h2 > 0.5

    return np.logical_and.reduce((idx_1, idx_2, idx_3, idx_4))


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
        # try:
        sim = load_galaxies.Simulation(output, sphere_radius_kpc=30, n_galaxies=1)
        # except:
        #     print("", flush=True)
        #     print("FAILED - sim creation failed", flush=True)
        #     print(output, flush=True)
        #     print("", flush=True)
        #     continue

        try:
            ad = sim.ds.all_data()
            idx_sf = find_sf_cells(ad)
            levels = gas_pdfs.get_gas_level(ad)
            sizes = gas_pdfs.get_cell_size_pc(ad).to("pc").value

            # match levels to cell size
            level_size_dict = {
                l: s for l, s in zip(np.unique(levels), np.unique(sizes)[::-1])
            }

            # then count the number of cells at the levels that have sf
            good_levels = levels[idx_sf]
            uniques, counts = np.unique(good_levels, return_counts=True)

            # then write this all. I'll do each output on one line
            out_str = f" z={sim.z:6.4f}"
            for idx in range(len(uniques)):
                out_str += (
                    f" - level {uniques[idx]:>2.0f} = "
                    f"{level_size_dict[uniques[idx]]:>5.2f} pc = "
                    f"{counts[idx]} cells"
                )
            write(out_str)

            del ad
            del levels
            del sizes
            del level_size_dict
            del good_levels
            del uniques
            del counts

        except:
            write(f"{run_name} {sim.z:.4f} Error!")

        del sim
        gc.collect()


# out_file.close()
