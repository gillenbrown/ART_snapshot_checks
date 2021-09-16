from pathlib import Path
import os
import numpy as np
import yt

from utils import load_galaxies

yt.funcs.mylog.setLevel(50)  # ignore yt's output

out_file = open("./molecular_gas.txt", "w")

scratch_dir = Path("/scratch/06912/tg862118/art_runs/analysis/production")

directories = [
    scratch_dir / "rj_sfe010_hn20",
    scratch_dir / "rj_sfe100_hn20",
    scratch_dir / "tl_sfe001_hn20",
    scratch_dir / "tl_sfe010_hn20",
    scratch_dir / "tl_sfe100_hn00",
    scratch_dir / "tl_sfe100_hn05",
    scratch_dir / "tl_sfe100_hn20",
]


def write(run_name, redshift, value1, value2):
    this_str = f"{run_name} {redshift:.4f} {value1:.3e} {value2:.3e}\n"
    out_file.write(this_str)
    out_file.flush()  # this and next line ensure writing happens now
    os.fsync(out_file.fileno())
    print(this_str, flush=True)


for d in directories:
    out_dir = d / "run" / "out"
    art_files = sorted([f for f in out_dir.iterdir() if f.suffix == ".art"])

    print("", flush=True)
    print(d, flush=True)
    print("", flush=True)

    run_name = d.name

    # then parse this to get the halo files
    for output in art_files:
        print("", flush=True)
        print(output, flush=True)
        print("", flush=True)

        # try:
        sim = load_galaxies.Simulation(output, sphere_radius_kpc=30, n_galaxies=1)
        # except:
        #     print("", flush=True)
        #     print("FAILED - sim creation failed", flush=True)
        #     print(output, flush=True)
        #     print("", flush=True)
        #     continue

        ad = sim.ds.all_data()

        # Then do analysis. This is copied from another analysis script
        # check that we have any halos at all. If not, we can exit. This can happen
        # for early outputs where nothing has collapsed yet.
        if len(sim.galaxies) == 0:
            write(run_name, sim.z, np.nan, np.nan)
            continue

        # only work on the biggest halo
        for gal in sim.galaxies:
            if gal.rank == 1:
                break

        assert gal.rank == 1

        # Print the molecular gas mass
        gas_mass_H2_totals = []
        for container in [gal.sphere, ad]:
            try:
                cell_volumes = container[("index", "cell_volume")]
                gas_mass_H2 = (
                    2
                    * container[("artio", "RT_HVAR_H2")]
                    * sim.ds.arr(1, "code_mass/code_length**3")
                    * cell_volumes
                )
                gas_mass_H2_totals.append(np.sum(gas_mass_H2.to("Msun").value))
            except:
                gas_mass_H2_totals.append(np.nan)

        write(run_name, sim.z, gas_mass_H2_totals[0], gas_mass_H2_totals[1])

out_file.close()
