from pathlib import Path

import numpy as np
import yt
from astropy import table
from matplotlib import cm
from matplotlib import colors as mpl_col

import betterplotlib as bpl

yt.funcs.mylog.setLevel(50)  # ignore yt's output

# ======================================================================================
#
# Set up plot info
#
# ======================================================================================
# I have to hardcode some labels to make this easier, parsing them won't work
# nearly as well
base_dir = Path.home() / "art_runs" / "runs"


def hui(suffix):
    return base_dir / "shangrila" / "hui" / suffix


def old_ic(suffix):
    return (
        base_dir / "stampede2" / "old_ic_comparison_production_analog" / suffix / "run"
    )


def production(suffix):
    return base_dir / "stampede2" / "production" / suffix / "run"


def rj_nbody(suffix):
    return base_dir / "stampede2" / "rj_nbody" / suffix / "run"


def prod_fmt(ic, eps_ff, f_hn):
    name = f"{ic} "
    name += "$\epsilon_{ff} = $" + f"{eps_ff:}%, "
    # add spaces to pad lower percents:
    if eps_ff < 10:
        name += " " * 4
    elif eps_ff < 100:
        name += " " * 2
    name += "$f_{HN} = $" + f"{f_hn}%"
    return name


names = {
    base_dir
    / "pleiades/nbody/new_ic_trim_25mpc/root_08/run/outputs/vel_offset/": "T&L Collisionless",
    hui("sfe_10"): "NBm SFE10",
    hui("sfe_100"): "NBm SFE100",
    old_ic("continuous_hn00_novirial"): "Continuous",
    old_ic(
        "continuous_hn00_virial10_entropy_fboost1"
    ): "ART 2.1 Entropy $f_{boost}=1$ Continuous",
    old_ic("continuoushui_hn00_novirial"): "Continuous Hui",
    # old_ic("continuouspopmcluster_hn00_novirial"): "Continuous PopM",
    old_ic("continuoussnr_hn00_novirial"): "Continuous SNR",
    old_ic("discrete_hn00_novirial"): "Discrete",
    old_ic(
        "discrete_hn00_novirial_entropy_fboost1"
    ): "ART 2.1 Entropy $f_{boost}=1$ No Virial",
    old_ic("discrete_hn00_virial10"): "Discrete $\\alpha<10$",
    old_ic("discrete_hn00_virial10_fboost3"): "Discrete $\\alpha<10, f_{boost}=3$",
    old_ic("discrete_hn20_virial10"): "Discrete $\\alpha<10$ HN20",
    old_ic("discrete_hn00_virial10_19"): "ART 1.9 Adiabatic",
    old_ic("discrete_hn00_virial10_19_advect"): "ART 1.9 Advect",
    old_ic("discrete_hn00_virial10_advect"): "ART 2.0 Advect",
    old_ic("discrete_hn00_virial10_elements"): "ART 2.0 Adiabatic All Elements",
    old_ic("discrete_hn00_virial10_nostars"): "ART 2.0 Adiabatic No Stars",
    old_ic("discrete_hn00_virial10_19_old"): "ART 1.9 Hui's version Adiabatic",
    old_ic("discrete_hn00_virial10_advect_nostars"): "ART 2.0 Advect No Stars",
    old_ic("discrete_hn00_virial10_noadvoradia_nostars"): "ART 2.0 No Flags",
    old_ic("discrete_hn00_virial10_entropy"): "ART 2.1 Entropy $f_{boost}=5$",
    old_ic("discrete_hn00_virial10_entropy_fboost3"): "ART 2.1 Entropy $f_{boost}=3$",
    old_ic("discrete_hn00_virial10_entropy_fboost2"): "ART 2.1 Entropy $f_{boost}=2$",
    old_ic("discrete_hn00_virial10_entropy_fboost1"): "ART 2.1 Entropy $f_{boost}=1$",
    old_ic(
        "discrete_hn00_virial10_entropy_fboost1_crho03"
    ): "ART 2.1 Entropy $f_{boost}=1$ $C_\\rho=3",
    old_ic(
        "discrete_hn00_virial10_entropy_fboost1_crho30"
    ): "ART 2.1 Entropy $f_{boost}=1$ $C_\\rho=30",
    old_ic(
        "discrete_hn00_virial10_entropy_fboost1_sfe001"
    ): "ART 2.1 Entropy $f_{boost}=1$ $\eps_{ff}=1$%",
    old_ic(
        "discrete_hn00_virial10_entropy_fboost1_sfe010"
    ): "ART 2.1 Entropy $f_{boost}=1$ $\eps_{ff}=10$%",
    old_ic(
        "discrete_hn00_virial10_entropy_fboost3_nosnia"
    ): "ART 2.1 Entropy $f_{boost}=3$ No SNIa",
    old_ic(
        "discrete_hn00_virial10_entropy_molvadim"
    ): "ART 2.1 Entropy Molecular Changes $f_{boost}=5$",
    old_ic(
        "discrete_hn00_virial10_entropy_molvadim_fboost1"
    ): "ART 2.1 Entropy Molecular Changes $f_{boost}=1$",
    old_ic(
        "discrete_hn00_virial10_entropy_molvadim_fboost2"
    ): "ART 2.1 Entropy Molecular Changes $f_{boost}=2$",
    old_ic(
        "discrete_hn00_virial10_entropy_molvadim_fboost3"
    ): "ART 2.1 Entropy Molecular Changes $f_{boost}=3$",
    old_ic(
        "discrete_hn00_virial10_entropy_newagediff"
    ): "ART 2.1 Entropy SN Timing Hybrid",
    old_ic(
        "discrete_hn00_virial10_entropy_newagediffallave"
    ): "ART 2.1 Entropy SN Timing Average",
    old_ic(
        "discrete_hn00_virial10_entropy_newagediffallbirth"
    ): "ART 2.1 Entropy SN Timing Birth",
    old_ic(
        "discrete_hn00_virial10_entropy_noagediff"
    ): "ART 2.1 Entropy $f_{boost}=5$ No Age Diff",
    old_ic(
        "discrete_hn00_virial10_entropy_hybridagediff"
    ): "ART 2.1 Entropy $f_{boost}=5$ Hybrid Age Diff",
    old_ic(
        "discrete_hn50_virial10_entropy_fboost1"
    ): "ART 2.1 Entropy $f_{boost}=1$ HN50",
    old_ic("discrete_hn00_virial10_entropy_nosync"): "ART 2.1 No Sync",
    old_ic("discrete_hn00_virial10_noturb_adi"): "ART 2.0 No Turbulence Adiabatic",
    old_ic("discrete_hn00_virial10_noturb_adv"): "ART 2.0 No Turbulence Advect",
    production("tl_sfe001_hn20"): prod_fmt("T&L", 1, 20),
    production("tl_sfe010_hn20"): prod_fmt("T&L", 10, 20),
    production("tl_sfe100_hn20"): prod_fmt("T&L", 100, 20),
    production("tl_sfe100_hn05"): prod_fmt("T&L", 100, 5),
    production("tl_sfe100_hn00"): prod_fmt("T&L", 100, 0),
    production("tl_sfe100_hn00_fboost1"): prod_fmt("T&L", 100, 0) + " $f_{boost}=1",
    production("tl_sfe100_hn00_fboost3"): prod_fmt("T&L", 100, 0) + " $f_{boost}=3",
    production("rj_sfe010_hn20"): prod_fmt("R&J", 10, 20),
    production("rj_sfe100_hn20"): prod_fmt("R&J", 100, 20),
    rj_nbody("original_92.48mpc_level07"): "R&J Collisionless Original",
    rj_nbody("hybrid_46.24mpc_level08"): "R&J Collisionless 2x Trim",
    rj_nbody("hybrid_23.12mpc_level08"): "R&J Collisionless 4x Trim",
}

cmap_nbm = cm.Greys
cmap_rj_collisionless = cm.Reds


def hsv_to_hex(h, s, v):
    return mpl_col.to_hex(mpl_col.hsv_to_rgb([h, s, v]))


colors = {
    "NBm SFE10": cmap_nbm(0.35),
    "NBm SFE100": cmap_nbm(0.65),
    "Continuous": bpl.color_cycle[0],
    "ART 2.1 Entropy $f_{boost}=1$ Continuous": bpl.color_cycle[0],
    "Continuous Hui": bpl.color_cycle[3],
    "Continuous PopM": "lightblue",
    "Continuous SNR": bpl.color_cycle[3],
    "Discrete": bpl.color_cycle[1],
    "ART 2.1 Entropy $f_{boost}=1$ No Virial": bpl.color_cycle[0],
    "Discrete $\\alpha<10$": bpl.color_cycle[5],
    "Discrete $\\alpha<10, f_{boost}=3$": bpl.color_cycle[6],
    "Discrete $\\alpha<10$ HN20": bpl.color_cycle[7],
    "ART 1.9": "red",
    "ART 1.9 Advect": "green",
    "ART 1.9 Adiabatic": "blue",
    "ART 2.0 Adiabatic All Elements": "orange",
    "ART 2.0 Adiabatic No Stars": "brown",
    "ART 1.9 Hui's version Adiabatic": "yellow",
    "ART 2.0 Advect No Stars": "purple",
    "ART 2.0 No Flags": "cyan",
    "ART 2.0 Advect": bpl.color_cycle[0],
    "ART 2.0 No Turbulence Adiabatic": bpl.color_cycle[4],
    "ART 2.0 No Turbulence Advect": bpl.color_cycle[6],
    "ART 2.1 Entropy $f_{boost}=5$": bpl.color_cycle[1],
    "ART 2.1 Entropy $f_{boost}=3$": bpl.color_cycle[0],
    "ART 2.1 Entropy $f_{boost}=2$": bpl.color_cycle[3],
    "ART 2.1 Entropy $f_{boost}=1$": bpl.color_cycle[4],
    "ART 2.1 Entropy $f_{boost}=1$ $C_\\rho=3": bpl.color_cycle[1],
    "ART 2.1 Entropy $f_{boost}=1$ $C_\\rho=30": bpl.color_cycle[2],
    "ART 2.1 Entropy $f_{boost}=1$ $\eps_{ff}=1$%": bpl.color_cycle[0],
    "ART 2.1 Entropy $f_{boost}=1$ $\eps_{ff}=10$%": bpl.color_cycle[1],
    "ART 2.1 Entropy $f_{boost}=3$ No SNIa": bpl.color_cycle[5],
    "ART 2.1 Entropy $f_{boost}=5$ No Age Diff": bpl.color_cycle[0],
    "ART 2.1 Entropy $f_{boost}=5$ Hybrid Age Diff": bpl.color_cycle[3],
    "ART 2.1 Entropy Molecular Changes $f_{boost}=5$": bpl.color_cycle[0],
    "ART 2.1 Entropy Molecular Changes $f_{boost}=3$": bpl.color_cycle[3],
    "ART 2.1 Entropy Molecular Changes $f_{boost}=2$": bpl.color_cycle[5],
    "ART 2.1 Entropy Molecular Changes $f_{boost}=1$": bpl.color_cycle[6],
    "ART 2.1 Entropy SN Timing Hybrid": bpl.color_cycle[0],
    "ART 2.1 Entropy SN Timing Average": bpl.color_cycle[3],
    "ART 2.1 Entropy SN Timing Birth": bpl.color_cycle[4],
    "ART 2.1 Entropy $f_{boost}=1$ HN50": bpl.color_cycle[6],
    "ART 2.1 No Sync": bpl.color_cycle[3],
    # These colors are very carefully chosen to avoid colorblindness issues. The hue
    # changes between the SFE variations (blue) to the HN variations (purple), with
    # the shared SFE 100 HN 20 run in the middle. The blues are more saturated, while
    # the purples are less saturated. I found this essential to making the colors
    # distinguishable to those with colorblindess.
    # If you ever change these, use davidmathlogic.com/colorblind to doublecheck.
    prod_fmt("T&L", 1, 20): hsv_to_hex(0.60, 0.45, 0.85),
    prod_fmt("T&L", 10, 20): hsv_to_hex(0.60, 0.70, 0.75),
    prod_fmt("T&L", 100, 20): hsv_to_hex(0.65, 0.80, 0.50),
    prod_fmt("T&L", 100, 5): hsv_to_hex(0.70, 0.30, 0.65),
    prod_fmt("T&L", 100, 0): hsv_to_hex(0.70, 0.15, 0.85),
    prod_fmt("T&L", 100, 0) + " $f_{boost}=3": hsv_to_hex(0.90, 0.15, 0.85),
    prod_fmt("T&L", 100, 0) + " $f_{boost}=1": hsv_to_hex(0.95, 0.15, 0.85),
    prod_fmt("R&J", 10, 20): hsv_to_hex(0.35, 0.20, 0.70),
    prod_fmt("R&J", 100, 20): hsv_to_hex(0.35, 0.30, 0.50),
    "R&J Collisionless Original": cmap_rj_collisionless(1.0),
    "R&J Collisionless 2x Trim": cmap_rj_collisionless(0.7),
    "R&J Collisionless 4x Trim": cmap_rj_collisionless(0.4),
    "T&L Collisionless": bpl.almost_black,
}

axes = {
    "NBm SFE10": [],
    "NBm SFE100": [
        "adi_adv",
    ],
    "Continuous": [],
    "ART 2.1 Entropy $f_{boost}=1$ Continuous": ["old_ic_discreteness"],
    "Continuous Hui": ["adi_adv", "old_ic_sn_timing"],
    "Continuous PopM": [],
    "Continuous SNR": [],
    "Discrete": [],
    "ART 2.1 Entropy $f_{boost}=1$ No Virial": ["old_ic_virial"],
    "Discrete $\\alpha<10$": ["adi_adv"],
    "Discrete $\\alpha<10, f_{boost}=3$": [],
    "Discrete $\\alpha<10$ HN20": [],
    "ART 1.9": [],
    "ART 1.9 Advect": [],
    "ART 1.9 Adiabatic": [],
    "ART 2.0 Adiabatic All Elements": [],
    "ART 2.0 Adiabatic No Stars": [],
    "ART 1.9 Hui's version Adiabatic": [],
    "ART 2.0 Advect No Stars": [],
    "ART 2.0 No Flags": [],
    "ART 2.0 Advect": ["adi_adv"],
    "ART 2.0 No Turbulence Adiabatic": [],
    "ART 2.0 No Turbulence Advect": [],
    "ART 2.1 Entropy $f_{boost}=5$": ["old_ic_fboost"],
    "ART 2.1 Entropy $f_{boost}=3$": ["old_ic_fboost"],
    "ART 2.1 Entropy $f_{boost}=2$": ["old_ic_fboost"],
    "ART 2.1 Entropy $f_{boost}=1$": [
        "old_ic_fboost",
        "old_ic_hn_fraction",
        "old_ic_sfe",
        "old_ic_virial",
        "old_ic_discreteness",
        "old_ic_molecular",
    ],
    "ART 2.1 Entropy $f_{boost}=1$ $C_\\rho=3": ["old_ic_molecular"],
    "ART 2.1 Entropy $f_{boost}=1$ $C_\\rho=30": ["old_ic_molecular"],
    "ART 2.1 Entropy $f_{boost}=1$ $\eps_{ff}=1$%": ["old_ic_sfe"],
    "ART 2.1 Entropy $f_{boost}=1$ $\eps_{ff}=10$%": ["old_ic_sfe"],
    "ART 2.1 Entropy $f_{boost}=3$ No SNIa": [],
    "ART 2.1 Entropy $f_{boost}=5$ No Age Diff": [],
    "ART 2.1 Entropy $f_{boost}=5$ Hybrid Age Diff": [],
    "ART 2.1 Entropy Molecular Changes $f_{boost}=1$": ["old_ic_molecular"],
    "ART 2.1 Entropy Molecular Changes $f_{boost}=2$": [],
    "ART 2.1 Entropy Molecular Changes $f_{boost}=3$": [],
    "ART 2.1 Entropy Molecular Changes $f_{boost}=5$": [],
    "ART 2.1 Entropy SN Timing Hybrid": ["old_ic_sn_timing"],
    "ART 2.1 Entropy SN Timing Average": ["old_ic_sn_timing"],
    "ART 2.1 Entropy SN Timing Birth": ["old_ic_sn_timing"],
    "ART 2.1 Entropy $f_{boost}=1$ HN50": ["old_ic_hn_fraction"],
    "ART 2.1 No Sync": [],
    prod_fmt("T&L", 1, 20): ["lg_sfe"],
    prod_fmt("T&L", 10, 20): ["lg_sfe"],
    prod_fmt("T&L", 100, 20): ["lg_hn_fraction", "lg_sfe"],
    prod_fmt("T&L", 100, 5): ["lg_hn_fraction"],
    prod_fmt("T&L", 100, 0): ["lg_fboost", "lg_hn_fraction"],
    prod_fmt("T&L", 100, 0) + " $f_{boost}=1": ["lg_fboost"],
    prod_fmt("T&L", 100, 0) + " $f_{boost}=3": ["lg_fboost"],
    prod_fmt("R&J", 10, 20): ["lg_sfe"],
    prod_fmt("R&J", 100, 20): ["lg_sfe"],
    "R&J Collisionless Original": [],
    "R&J Collisionless 2x Trim": [],
    "R&J Collisionless 4x Trim": [],
    "T&L Collisionless": [],
}


def get_plot_names(sim_names):
    all_plots = []
    for name in sim_names:
        for plot in axes[name]:
            all_plots.append(plot)
    return list(set(all_plots))


# ======================================================================================
#
# Setup for precalculation of key quantities
#
# ======================================================================================
class Galaxy(object):
    def __init__(
        self,
        ds,
        center=None,
        sphere_radius=None,
        m_vir=None,
        r_vir=None,
        rank=None,
        name=None,
    ):
        self.ds = ds
        self.center = center
        self.m_vir = m_vir
        self.r_vir = r_vir
        self.rank = rank
        self.name = name

        self.mask_done_forming = None

        # only create the sphere if the user wants to
        if sphere_radius is not None:
            self.sphere = self.ds.sphere(center=self.center, radius=sphere_radius)
        else:
            self.sphere = None

        # have a place for some things (like the CIMF) to be stored as a user
        # calculates them.
        self.precalculated = dict()

    def __getitem__(self, property):
        """
        Automatically apply the cut to only pick formed clusters
        """
        if self.sphere is None:
            raise ValueError("No sphere set on galaxy initialization")

        # then get the mask showing which clusters are done forming
        if self.mask_done_forming is None and property[0] == "STAR":
            self.mask_done_forming = self.sphere[("STAR", "age")] > 15 * yt.units.Myr

        quantity = self.sphere[property]
        if property[0] == "STAR":
            quantity = quantity[self.mask_done_forming]
        return quantity

    def prop_all_clusters(self, property):
        """
        Like __getitem__, but does not apply the restriction that clusters must be
        done forming
        """
        if self.sphere is None:
            raise ValueError("No sphere set on galaxy initialization")

        return self.sphere[property]


class Simulation(object):
    def __init__(
        self, ds_path, sphere_radius_kpc=None, min_virial=False, n_galaxies=None
    ):
        """
        ds_path must be a Path object
        """
        # get the dataset and corresponding halo file
        run_dir = ds_path.parent.parent
        halo_path = run_dir / "halos"
        halo_name = ds_path.name.replace("continuous_", "out_")
        halo_name = halo_name.replace(".art", ".list")

        self.ds = yt.load(str(ds_path))

        self.z = self.ds.current_redshift
        self.scale_factor = 1 / (1 + self.z)

        # get the axis names and other stuff. There are times where these aren't
        # defined, but when it's not defined it's not needed. Using .get returns none
        # as the default value.
        self.name = names.get(run_dir)
        self.axes = axes.get(self.name)
        self.color = colors.get(self.name)

        # have the place to store some precalculated things, particularly those that
        # include all galaxies in the simulation
        self.precalculated = dict()

        # if we are the old IC set, we have one galaxy, otherwise two
        # check what kind of particles are present
        if n_galaxies is None:
            if ("N-BODY_0", "MASS") in self.ds.derived_field_list:
                self.n_galaxies = 2
            else:
                self.n_galaxies = 1
        else:
            self.n_galaxies = n_galaxies

        # load halo catalogs
        cat = table.Table.read(halo_path / halo_name, format="ascii.commented_header")
        # delete some unneeded quantities. If I want to inlcude these later, I'll need
        # to verify I'm doing the units correctly.
        for col in cat.colnames:
            if col not in ["Mvir", "Vmax", "Vrms", "Rvir", "X", "Y", "Z"]:
                del cat[col]
        # To modify the units of things, we need to know little h. It's in the
        # file header. We also want a, to turn things into physical units.
        with open(halo_path / halo_name, "r") as in_file:
            line_num = 1
            for line in in_file:
                if line_num == 2:
                    a = float(line.split()[-1])
                elif line_num == 3:
                    h = float(line.split()[-1])
                    break

                line_num += 1
        assert 0 < a < 1.01  # slight buffer for last output
        assert 0.6 < h < 0.8

        # Masses are in Msun / h
        cat["Mvir"] = cat["Mvir"] / h
        # Positions in Mpc / h (comoving)
        for col in ["X", "Y", "Z"]:
            cat[col] = cat[col] / h
        # Halo Distances, Lengths, and Radii in kpc / h (comoving)
        cat["Rvir"] = cat["Rvir"] / h
        # Velocities in km / s (physical, peculiar) -- no change needed

        # add units to names
        cat.rename_column("Mvir", "Mvir_msun")
        cat.rename_column("X", "X_mpccm")
        cat.rename_column("Y", "Y_mpccm")
        cat.rename_column("Z", "Z_mpccm")
        cat.rename_column("Rvir", "Rvir_kpccm")
        cat.rename_column("Vmax", "Vmax_kms")
        cat.rename_column("Vrms", "Vrms_kms")

        # Do some parsing of the halo catalogs
        # We get the indices that sort it. The reversing there makes the biggest
        # halos first, like we want.
        rank_idxs = np.argsort(cat["Mvir_msun"])[::-1]
        # Add this info to each halo object, and put the halos into a new sorted list,
        # with the highest mass (lowest rank) halos first). We only keep the number
        # the user requested
        self.galaxies = []
        for rank, idx in enumerate(rank_idxs[: self.n_galaxies], start=1):
            row = cat[idx]

            r_vir = self.ds.arr(row["Rvir_kpccm"], "kpccm")
            M_vir = self.ds.arr(row["Mvir_msun"], "Msun")
            center = self.ds.arr(
                [row["X_mpccm"], row["Y_mpccm"], row["Z_mpccm"]], "Mpccm"
            )

            if sphere_radius_kpc is None:
                radius = None
            else:
                if min_virial:
                    radius = self.ds.arr(
                        min(r_vir.to("kpc").value, sphere_radius_kpc), "kpc"
                    )
                else:
                    radius = self.ds.arr(sphere_radius_kpc, "kpc")

            self.galaxies.append(Galaxy(self.ds, center, radius, M_vir, r_vir, rank))

    def __repr__(self):
        return self.name


# ======================================================================================
#
# Loading datasets
#
# ======================================================================================
# Get the datasets and halo catalogs. When doing these we need to be a bit
# careful about the datasets. We will make one set of comparisons at the last
# common output of all simulations, then one with the last output of each
# simulation. Those all need to be stored separately.
def filename_to_scale_factor(filename):
    return float(filename[-10:-4])


def get_outputs_in_dir(sim_dir):
    directory = Path(sim_dir)
    if directory not in names:
        print(f"Skipping {directory}")
        return []

    out_dir = directory / "out"
    return [
        file
        for file in out_dir.iterdir()
        if file.is_file()
        and str(file.name).endswith(".art")
        and str(file.name).startswith("continuous_a")
    ]


def get_simulations_last(sim_dirs, sphere_radius_kpc=30, min_virial=True):
    sims_last = []
    for directory in sorted(sim_dirs):
        all_outputs = get_outputs_in_dir(directory)
        if len(all_outputs) > 0:
            last_output = sorted(all_outputs)[-1]
            sims_last.append(Simulation(last_output, sphere_radius_kpc, min_virial))

    return sims_last


def get_common_scale_factor(sim_dirs, z_max):
    # Start by getting the last common output among the provided runs
    last_outputs = []
    for directory in sim_dirs:
        directory = Path(directory)
        if directory not in names:
            continue

        all_outputs = get_outputs_in_dir(directory)
        if len(all_outputs) == 0:
            print(f"This has no outputs: {directory}")
            continue

        # restrict to be a reasonable redshift
        a_min = 1 / (1 + z_max)
        this_last_output = sorted(all_outputs)[-1]
        if filename_to_scale_factor(this_last_output.name) > a_min:
            last_outputs.append(filename_to_scale_factor(this_last_output.name))

    # include fudge factor for scale comparisons (so 0.1801 and 0.1802 match)
    return min(last_outputs) + 0.001


def get_simulations_common(sim_dirs, z_max=5, sphere_radius_kpc=30, min_virial=True):
    common_scale = get_common_scale_factor(sim_dirs, z_max)

    sims_common = []
    for directory in sorted(sim_dirs):
        all_outputs = get_outputs_in_dir(directory)
        if len(all_outputs) == 0:
            # happens when skipped by get_outputs_in_dir
            continue

        # get the last one that's in common with the other simulations
        all_common_outputs = [
            file
            for file in all_outputs
            if filename_to_scale_factor(file.name) <= common_scale
            and abs(filename_to_scale_factor(file.name) - common_scale) < 0.02
        ]
        # if there are no snapshots early enough for this, don't add them
        if len(all_common_outputs) > 0:
            last_common_snapshot = sorted(all_common_outputs)[-1]
            sims_common.append(
                Simulation(last_common_snapshot, sphere_radius_kpc, min_virial)
            )

    return sims_common, common_scale


def get_simulations_same_scale(
    sim_dirs, desired_z, z_tolerance=0.05, sphere_radius_kpc=30, min_virial=True
):
    sims_common = []
    for directory in sorted(sim_dirs):
        all_outputs = get_outputs_in_dir(directory)
        if len(all_outputs) == 0:
            # happens when skipped by get_outputs_in_dir
            continue

        closest_z = np.inf
        closest_snapshot = None
        for out in all_outputs:
            a = filename_to_scale_factor(out.name)
            z = (1 / a) - 1

            if abs(z - desired_z) < abs(closest_z - desired_z):
                closest_z = z
                closest_snapshot = out

        # then validate that we got close enough.
        if abs(closest_z - desired_z) / desired_z < z_tolerance:
            sims_common.append(
                Simulation(closest_snapshot, sphere_radius_kpc, min_virial)
            )
        else:
            print(f"No outputs close to z={desired_z} in {directory}")
    return sims_common
