"""
plot_z0_age_metallicity.py

Plot the age-metallicity relation for surviving clusters at z=0
"""
from pathlib import Path

import sys
import numpy as np
from matplotlib import pyplot as plt
from astropy import cosmology, table
from astropy import units as u

from utils import load_galaxies, plot_utils
from analysis_functions import cimf

import betterplotlib as bpl

bpl.set_style()

sentinel = Path(sys.argv[1])
plot_dir = sentinel.parent

# ======================================================================================
#
# load the simulations and observational data
#
# ======================================================================================
sims = load_galaxies.get_simulations_last(sys.argv[2:])

amr_loc = Path(__file__).parent / "data" / "vandenberg_13_amr.txt"
amr_table = table.Table.read(str(amr_loc), format="ascii")
gc_age = amr_table["Age"]
gc_age_err = amr_table["Age_err"]
gc_feh = amr_table["[Fe/H]"]

# ======================================================================================
#
# analysis functions
#
# ======================================================================================
def get_metallicity(galaxy):
    z = galaxy[("STAR", "METALLICITY_SNII")] + galaxy[("STAR", "METALLICITY_SNIa")]
    try:
        z += galaxy[("STAR", "METALLICITY_AGB")]
    except:
        pass

    return z


def get_feh(galaxy):
    try:
        star_fe = galaxy[("STAR", "METALLICITY_Fe")]
        sun_fe = 1.23e-03  # Grevesse & Sauval 1998, SSR, 85, 161
        return np.log10(star_fe / sun_fe)
    except:
        # no Fe present. Simply use raw metallicity.
        # calibrated by examining the distribution of Z vs [Fe/H] for runs that
        # actually have [Fe/H]
        return np.log10(get_metallicity(galaxy) / 0.03)


def get_birth_time(galaxy):
    return galaxy[("STAR", "creation_time")].to("Gyr").value


def plot_age_metallicity(sim):
    fig, ax = bpl.subplots(figsize=[8, 7])

    # get cosmology to get age of the universe at z=0
    H_0 = sim.ds.artio_parameters["hubble"][0] * 100 * u.km / (u.Mpc * u.second)
    omega_matter = sim.ds.artio_parameters["OmegaM"][0]
    cosmo = cosmology.FlatLambdaCDM(H0=H_0, Om0=omega_matter, Tcmb0=2.725)
    t_z_0 = cosmo.age(0).to("Gyr").value

    feh = sim.func_all_galaxies(get_feh)
    birth_time = sim.func_all_galaxies(get_birth_time)
    ages = t_z_0 - birth_time
    evolved_masses = sim.func_all_galaxies(cimf.evolve_cluster_population)

    good_idx = evolved_masses > 1e4

    if np.sum(good_idx) > 0:  # only plot if we have clusters!
        # ax.scatter(feh[good_idx], ages[good_idx], alpha=1, label="Simulations")
        ax.shaded_density(
            feh[good_idx],
            ages[good_idx],
            bin_size=0.05,
            smoothing=0.2,
            cmap="Greys",
            log_hist=False,
        )
        ax.density_contour(
            feh[good_idx],
            ages[good_idx],
            bin_size=0.05,
            percent_levels=[0.5, 0.9],
            smoothing=0.2,
            colors=bpl.color_cycle[2],
            labels=True,
        )

    # plot observational data
    ax.errorbar(
        x=gc_feh,
        y=gc_age,
        yerr=gc_age_err,
        xerr=0.15,
        c=bpl.color_cycle[0],
        label="MW GCs",
        zorder=100,
    )

    # format plot
    ax.add_labels("[Fe/H]", "Age [Gyr]")
    ax.set_limits(-3, 0.5, 10, 14)
    ax.legend()

    # then save the plot
    sim_name = plot_utils.get_sim_dirname(sim.run_dir)
    fig.savefig(plot_dir / f"z0_age_metallicity_{sim_name}.pdf")
    # then remove figure for memory purposes
    plt.close(fig)


# ======================================================================================
#
# do the work
#
# ======================================================================================
for sim in sims:
    plot_age_metallicity(sim)

sentinel.touch()
