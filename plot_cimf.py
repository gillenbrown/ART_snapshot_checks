"""
plot_cimf.py

Creates a plot showing the comparative cluster initial mass function of 
galaxies in an output
"""
import sys
from pathlib import Path

import numpy as np
from scipy import special
import yt
from yt.extensions.astro_analysis.halo_analysis.api import HaloCatalog

import betterplotlib as bpl

bpl.set_style()
yt.funcs.mylog.setLevel(50)  # ignore yt's output

import plot_utils

cimf_sentinel = Path(sys.argv[1])

# ======================================================================================
#
# Setup for precalculation of key quantities
#
# ======================================================================================
class Galaxy(object):
    def __init__(self, sphere, rank):
        self.sphere = sphere
        self.rank = rank

        # storing calculated cimfs
        self._precalc_cimf = dict()

    # Then the functions to calculate the CIMF. Here we need to do some analysis
    # of the bound fraction.
    @staticmethod
    def _f_bound(eps_int):
        # Li et al 2019: https://ui.adsabs.harvard.edu/abs/2019MNRAS.487..364L/abstract
        # equation 17
        alpha_star = 0.48
        f_sat = 0.94
        term_a = special.erf(np.sqrt(3 * eps_int / alpha_star))
        term_b = np.sqrt(12 * eps_int / (np.pi * alpha_star))
        term_c = np.exp(-3 * eps_int / alpha_star)
        return (term_a - (term_b * term_c)) * f_sat

    def _get_initial_bound_fraction(self):
        star_initial_mass = self.sphere[("STAR", "INITIAL_MASS")].to("Msun").value
        # the variable named INITIAL_BOUND_FRACTION is not the initial_bound fraction,
        # it's actually the accumulated mass nearby through the course of accretion, in
        # code masses. This is used to calculate the formation efficiency, which is then
        # used to get the bound fraction.
        star_accumulated_mass = (
            self.sphere[("STAR", "INITIAL_BOUND_FRACTION")].to("").value
        )
        star_accumulated_mass *= self.sphere.ds.mass_unit
        star_accumulated_mass = star_accumulated_mass.to("Msun").value

        eps_int = star_initial_mass / star_accumulated_mass

        return self._f_bound(eps_int)

    def cimf(self, mass_type, max_age_myr):
        id_string = f"{mass_type}_{max_age_myr}"

        if id_string not in self._precalc_cimf:
            self._precalc_cimf[id_string] = self._create_cimf(mass_type, max_age_myr)

        return self._precalc_cimf[id_string]

    def _create_cimf(self, mass_type, max_age_myr):
        """
        Make the cluster initial mass function.

        Notes on the different mass variables:
        ('STAR', 'INITIAL_MASS') - initial stellar mass of the star particle
        ('STAR', 'MASS') - current stellar mass of the star particle, accounting for
            stellar evolution
        ('STAR', 'INITIAL_BOUND_FRACTION') - NOT the actual initial_bound fraction. See
            the `get_initial_bound_fraction` function above for more on how to use this,
            but this variable is the accumulated mass near the cluster over the course
            of accretion. This is used to calculate formation efficiency, which is then
            used to get the actual initial_bound fraction
        ('STAR', 'BOUND_FRACTION') - This is the actual bound fraction at the current
            time, but NOT accounting for the proper initial_bound fraction

        :param data_obj: Sphere object representing a galaxy
        :param mass_type: String encoding which mass to get here. The options are:
                          "initial" - just the initial stellar masses
                          "initial_bound" - initial masses including initial_bound
                                            fraction
                          "current" - current bound mass, accounting for the initial
                                      bound fraction, tidal disruption, and stellar
                                      death
        :param include_initial_bound: whether to incorporate the initial bound
                                      fraction of clusters, or just get the
                                      distribution of initial particle masses
        :param max_age_myr: The maximum age to restrict the plot to. Is infinity as the
                            default, which plots all stars.
        :returns: Two lists. The first is f_i * M, representing the initial
                  bound mass, of M if include_initial_bound=False. This will be
                  binned values suitable to plot. The second is dN/dLogM for each
                  of the bins in the first list.
        """
        if mass_type == "initial":
            mass = self.sphere[("STAR", "INITIAL_MASS")].to("Msun").value
        elif mass_type == "initial_bound":
            initial_mass = self.sphere[("STAR", "INITIAL_MASS")].to("Msun").value
            star_initial_bound = self._get_initial_bound_fraction()
            mass = initial_mass * star_initial_bound
        elif mass_type == "current":
            raw_mass = self.sphere[("STAR", "INITIAL_MASS")].to("Msun").value
            star_initial_bound = self._get_initial_bound_fraction()
            tidal_bound_fraction = self.sphere[("STAR", "BOUND_FRACTION")].value
            mass = raw_mass * star_initial_bound * tidal_bound_fraction

        # then restrict to recently formed clusters. This can be set to infinity, which
        # plots everything
        max_age = max_age_myr * yt.units.Myr
        mask = self.sphere[("STAR", "age")] < max_age
        mass = mass[mask]

        # create bins with spacing of 0.16 dex
        bin_width = 0.16  # dex
        m_boundaries_log = np.arange(3 - 0.5 * bin_width, 7, 0.16)
        m_centers_log = [
            np.mean([m_boundaries_log[idx], m_boundaries_log[idx + 1]])
            for idx in range(len(m_boundaries_log) - 1)
        ]

        m_boundaries = 10 ** m_boundaries_log
        m_centers = 10 ** np.array(m_centers_log)

        # then make the histogram showing how many there are per bin
        hist, edges = np.histogram(mass, bins=m_boundaries)
        assert np.array_equiv(m_boundaries, edges)

        # We have dN, make it per dLogM
        hist = np.array(hist) / (bin_width * np.log(10))

        return m_centers, hist


class Simulation(object):
    def __init__(self, ds_path):
        # get the dataset and corresponding halo file
        run_dir = ds_path.parent.parent
        halo_path = run_dir / "halos"
        halo_name = ds_path.name.replace("continuous_", "halos_")
        halo_name = halo_name.replace(".art", ".0.bin")

        self.ds = yt.load(str(ds_path))
        self.halos_ds = yt.load(str(halo_path / halo_name))

        # get the axis names and other stuff
        self.name = plot_utils.names[run_dir]
        self.axes = plot_utils.axes[self.name]
        self.color = plot_utils.colors[self.name]

        # if we are the old IC set, we have one galaxy, otherwise two
        # check what kind of particles are present
        if ("N-BODY_0", "MASS") in self.ds.derived_field_list:
            self.n_galaxies_each = 2
        else:
            self.n_galaxies_each = 1

        # create halo catalog object
        self.hc = HaloCatalog(halos_ds=self.halos_ds, data_ds=self.ds)
        self.hc.create(save_halos=True, save_catalog=False)

        # Do some parsing of the halo catalogs
        # We get the indices that sort it. The reversing there makes the biggest
        # halos first, like we want.
        halo_masses = yt.YTArray(
            [halo.quantities["particle_mass"] for halo in self.hc.halo_list]
        )
        rank_idxs = np.argsort(halo_masses)[::-1]
        # Add this info to each halo object, and put the halos into a new sorted list,
        # with the highest mass (lowest rank) halos first). We only keep the number
        # the user requested
        self.galaxies = []
        for rank, idx in enumerate(rank_idxs[: self.n_galaxies_each], start=1):
            halo = self.hc.halo_list[idx]
            radius = min(halo.quantities["virial_radius"], 30 * yt.units.kpc)
            center = [
                halo.quantities["particle_position_x"],
                halo.quantities["particle_position_y"],
                halo.quantities["particle_position_z"],
            ]
            sphere = self.ds.sphere(center=center, radius=radius)
            self.galaxies.append(Galaxy(sphere, rank))


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


# Start by getting the last common output among the production runs
last_snapshots = []
for directory in sys.argv[2:]:
    directory = Path(directory)
    if directory not in plot_utils.names or "stampede2/production" not in str(
        directory
    ):
        continue

    out_dir = directory / "out"

    all_snapshots = [
        file.name
        for file in out_dir.iterdir()
        if file.is_file()
        and str(file.name).endswith(".art")
        and str(file.name).startswith("continuous_a")
    ]
    # restrict to be a reasonable redshift
    this_last_snapshot = sorted(all_snapshots)[-1]
    if filename_to_scale_factor(this_last_snapshot) > 0.15:
        last_snapshots.append(this_last_snapshot)

earliest_last_snapshot = sorted(last_snapshots)[0]
common_scale = filename_to_scale_factor(earliest_last_snapshot) + 0.001
# include fudge factor for scale comparisons (so 0.1801 and 0.1802 match)

# set up the dictionaries where we will store the datasets and halos
sims_common = []
sims_last = []

for directory in sys.argv[2:]:
    directory = Path(directory)
    if directory not in plot_utils.names:
        print(f"Skipping {directory}")
        continue

    out_dir = directory / "out"

    all_snapshots = [
        file
        for file in out_dir.iterdir()
        if file.is_file()
        and str(file.name).endswith(".art")
        and str(file.name).startswith("continuous_a")
    ]
    # get the actual last snapshot
    last_snapshot = sorted(all_snapshots)[-1]
    sims_last.append(Simulation(last_snapshot))

    # then the last one that's in common with the other simulations
    all_common_snapshots = [
        file
        for file in all_snapshots
        if filename_to_scale_factor(file.name) <= common_scale
        and abs(filename_to_scale_factor(file.name) - common_scale) < 0.02
    ]
    # if there are no snapshots early enough for this, don't add them
    if len(all_common_snapshots) > 0:
        last_common_snapshot = sorted(all_common_snapshots)[-1]
        sims_common.append(Simulation(last_common_snapshot))

# ======================================================================================
#
# plotting
#
# ======================================================================================
def plot_power_law(ax, slope, x1, x2, y1):
    # Here the slope is dN/dM \propto M^slope
    # Since we plot this in dN/dlogM space we
    # have to add 1 to the slope
    # dlogM \propto dM/M
    y2 = y1 * (x2 / x1) ** (slope + 1)

    ax.plot([x1, x2], [y1, y2], c=bpl.almost_black, lw=1, ls="--")
    ax.add_text(1.1 * x2, y2, text=slope, va="center", ha="left", fontsize=18)


def plot_cimf(axis_name, sim_share_type, masses_to_plot, max_age_myr=np.inf):
    """
    Plot various versions of the cluster mass function.

    :param axis_name: The name of the axis, so we can select which simulations go
                      on this plot.
    :param sim_share_type: Either "last" or "common". This will determine if
                           the redshift is labeled as common to all, or
                           individually per simulation
    :param masses_to_plot: A lit of the masses to plot, see the `cimf` function for the
                           allowed options
    :param max_age_myr: The maximum age to restrict the plot to. Is infinity as the
                        default, which plots all stars.
    """
    # validate options given
    if sim_share_type not in ["last", "common"]:
        raise ValueError("bad sim_share_type")
    # if we plot current, it should be alone, with nothing else
    if "current" in masses_to_plot and len(masses_to_plot) > 1:
        raise RuntimeError("Current masses must be plotted alone.")

    # make the name of the plot
    plot_name = f"cimf_{axis_name.replace('_', '')}_{sim_share_type}"
    if "current" in masses_to_plot:
        plot_name += "_current"
    elif (
        len(masses_to_plot) == 2
        and "initial" in masses_to_plot
        and "initial_bound" in masses_to_plot
    ):
        plot_name += "_initial"
    else:
        raise ValueError("This plot not yet named:", masses_to_plot)
    # add the age if it's not infinity
    if not np.isinf(max_age_myr):
        plot_name += f"_{max_age_myr}myr"
    plot_name += ".pdf"

    # choose which simulations we're using
    if sim_share_type == "common":
        sims = sims_common
    else:
        sims = sims_last

    fig, ax = bpl.subplots(figsize=[9, 7])

    for sim in sims:
        if axis_name not in sim.axes:
            continue
        for galaxy in sim.galaxies:
            for mass_type in masses_to_plot:
                mass_plot, dn_dlogM = galaxy.cimf(mass_type, max_age_myr)

                # make the label only for the biggest halo, and not for initial only
                if galaxy.rank == 1 and mass_type != "initial":
                    # and include the redshift if it's different for each sim
                    if sim_share_type == "last":
                        label = f"{sim.name}: z = {1/sim.ds.scale_factor - 1:.1f}"
                    else:
                        label = sim.name
                else:
                    label = None

                # have different line styles
                lss = {"initial": ":", "initial_bound": "-", "current": "-"}

                ax.plot(
                    mass_plot, dn_dlogM, c=sim.color, ls=lss[mass_type], label=label
                )

    # formax axes
    ax.legend(loc=1, fontsize=10)
    ax.set_yscale("log")
    ax.set_xscale("log")
    # have different y limits for different versions of the plot
    # the minimum value in the plot is 1 / (0.16 * ln(10) = 2.5
    # put the plot limit just above that, to make it cleaner and stop
    # weird vertical lines.
    if not np.isinf(max_age_myr):
        # small timeframe, don't need to show much.
        y_min = 3
        y_max = 1e5
    elif masses_to_plot == "current":
        y_min = 10
        y_max = 1e5
    else:
        y_min = 10
        y_max = 1e6
    ax.set_limits(1e3, 1e7, y_min, y_max)

    # plot the guiding lines
    log_space = 0.4 * (np.log10(y_max) - np.log10(y_min))
    y_guide = 10 ** (np.log10(y_min) + log_space)
    plot_power_law(ax, -2, 1e6, 3e6, y_guide)
    plot_power_law(ax, -3, 1e6, 3e6, y_guide)

    # if there is a common redshift, annotate it
    if sim_share_type == "common":
        ax.easy_add_text(f"z = {1/common_scale - 1:.1f}", "upper left")

    if "current" in masses_to_plot:
        ax.add_labels("$f_b$M [$M_\odot$]", "dN/dlogM")
    else:
        ax.add_labels("$f_i M_i$ [$M_\odot$]", "dN/dlogM")

    fig.savefig(cimf_sentinel.parent / plot_name)


# ======================================================================================
#
# loop to actually make the plots
#
# ======================================================================================
# then actually call this function to build the plots
for plot_name in plot_utils.get_plot_names([sim.name for sim in sims_last]):
    for share_type in ["common", "last"]:
        # plot main CIMF and unbound CIMF
        plot_cimf(plot_name, share_type, ["initial_bound", "initial"])
        # plot recently formed clusters
        plot_cimf(plot_name, share_type, ["initial_bound", "initial"], 300)
        # plot surviving clusters
        plot_cimf(plot_name, share_type, ["current"])

cimf_sentinel.touch()
