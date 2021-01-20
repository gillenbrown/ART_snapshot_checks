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

from plot_utils import names, colors, axis_number
    
def filename_to_scale_factor(filename):
    return float(filename[-10:-4]) 

def get_ds_and_halos(ds_path):
    """ get the dataset and corresponding halo file """
    halo_path = ds_path.replace("out/continuous_", "halos/halos_").replace(".art", ".0.bin")

    ds = yt.load(ds_path)
    halos_ds = yt.load(halo_path)

    # if we are the old IC set, we have one galaxy, otherwise two
    # check what kind of particles are present
    if ('N-BODY_0', 'MASS') in ds.derived_field_list:
        n_galaxies_each = 2
    else:
        n_galaxies_each = 1

    # create halo catalog object
    hc = HaloCatalog(halos_ds=halos_ds, data_ds=ds)
    hc.create(save_halos=True, save_catalog=False)

    # Do some parsing of the halo catalogs
    # We get the indices that sort it. The reversing there makes the biggest 
    # halos first, like we want.
    halo_masses = yt.YTArray([halo.quantities["particle_mass"] 
                              for halo in hc.halo_list])
    rank_idxs = np.argsort(halo_masses)[::-1]
    # Add this info to each halo object, and put the halos into a new sorted list,
    # with the highest mass (lowest rank) halos first). We only keep the number
    # the user requested
    halos = []
    for rank, idx in enumerate(rank_idxs[:n_galaxies_each], start=1):
        halo = hc.halo_list[idx]
        halo.quantities["rank"] = rank
        halos.append(halo)

    return ds, halos

# Get the datasets and halo catalogs. When doing these we need to be a bit 
# careful about the datasets. We will make one set of comparisons at the last
# common output of all simulations, then one with the last output of each 
# simulation. Those all need to be stored separately.

# Start by getting the last common output
last_snapshots = []
for directory in sys.argv[1:]:
    directory = Path(directory)
    if directory not in names:
        print(f"Skipping {directory}")
        continue

    out_dir = directory / "out"
    
    all_snapshots = [file.name for file in out_dir.iterdir()
                     if file.is_file()
                     and str(file.name).endswith(".art")
                     and str(file.name).startswith("continuous_a")]
    # restrict to be a reasonable redshift
    this_last_snapshot = sorted(all_snapshots)[-1]
    if filename_to_scale_factor(this_last_snapshot) > 0.15:
        last_snapshots.append(this_last_snapshot)
    
earliest_last_snapshot = sorted(last_snapshots)[0]
common_scale = filename_to_scale_factor(earliest_last_snapshot) + 0.001  
# include fudge factor for scale comparisons (so 0.1801 and 0.1802 match)

# set up the dictionaries where we will store the datasets and halos
common_ds = dict()
last_ds = dict()
common_halos = dict()
last_halos = dict()

for directory in sys.argv[1:]:
    directory = Path(directory)
    if directory not in names:
        continue

    out_dir = directory / "out"
    
    all_snapshots = [file for file in out_dir.iterdir()
                     if file.is_file()
                     and str(file.name).endswith(".art")
                     and str(file.name).startswith("continuous_a")]
    # get the actual last snapshot
    last_snapshot = sorted(all_snapshots)[-1]

    ds_last, halos_last = get_ds_and_halos(str(last_snapshot))
    last_ds[names[directory]] = ds_last
    last_halos[names[directory]] = halos_last
    
    # then the last one that's in common with the other simulations
    all_common_snapshots = [file for file in all_snapshots
                            if filename_to_scale_factor(file.name) <= common_scale
                            and abs(filename_to_scale_factor(file.name) - common_scale) < 0.02]
    # if there are no snapshots early enough for this, don't add them
    if len(all_common_snapshots) > 0:
        last_common_snapshot = sorted(all_common_snapshots)[-1]

        # then get the dataset and halo objects
        ds_common, halos_common = get_ds_and_halos(str(last_common_snapshot))
        common_ds[names[directory]] = ds_common
        common_halos[names[directory]] = halos_common
    

# Then the functions to calculate the CIMF. Here we need to do some analysis
# of the bound fraction. 
def f_bound(eps_int):
    # Li et al 2019: https://ui.adsabs.harvard.edu/abs/2019MNRAS.487..364L/abstract
    # equation 17
    alpha_star = 0.48
    f_sat = 0.94
    term_a = special.erf(np.sqrt(3 * eps_int / alpha_star))
    term_b = np.sqrt(12 * eps_int / (np.pi * alpha_star))
    term_c = np.exp(-3 * eps_int / alpha_star)
    return (term_a - (term_b * term_c)) * f_sat

def get_initial_bound_fraction(data_obj):
    star_initial_mass = data_obj[('STAR', 'INITIAL_MASS')].to("Msun").value
    # the variable named INITIAL_BOUND_FRACTION is not the initial bound fraction, 
    # it's actually the accumulated mass nearby through the course of accretion, in
    # code masses. This is used to calculate the formation efficiency, which is then
    # used to get the bound fraction.
    star_accumulated_mass = data_obj[('STAR', 'INITIAL_BOUND_FRACTION')].to("").value
    star_accumulated_mass *= data_obj.ds.mass_unit
    star_accumulated_mass = star_accumulated_mass.to("Msun").value
    
    eps_int = star_initial_mass / star_accumulated_mass
    
    return f_bound(eps_int)

# Then the actual plotting code
def plot_power_law(ax, slope, x1, x2, y1):
    # Here the slope is dN/dM \propto M^slope
    # Since we plot this in dN/dlogM space we 
    # have to add 1 to the slope
    # dlogM \propto dM/M
    y2 = y1 * (x2/x1)**(slope+1)
    
    ax.plot([x1, x2], [y1, y2], c=bpl.almost_black, lw=1, ls="--")
    ax.add_text(1.1*x2, y2, text=slope, va="center", ha="left", fontsize=18)

def cimf(data_obj, mass_type, max_age_myr):
    """
    Make the cluster initial mass function. \

    Notes on the different mass variables:
    ('STAR', 'INITIAL_MASS') - initial stellar mass of the star particle
    ('STAR', 'MASS') - current stellar mass of the star particle, accounting for
        stellar evolution
    ('STAR', 'INITIAL_BOUND_FRACTION') - NOT the actual initial bound fraction. See
        the `get_initial_bound_fraction` function above for more on how to use this,
        but this variable is the accumulated mass near the cluster over the course of
        accretion. This is used to calculate formation efficiency, which is then used
        to get the actual initial bound fraction
    ('STAR', 'BOUND_FRACTION') - This is the actual bound fraction at the current time,
        but NOT accounting for the proper initial bound fraction
    
    :param data_obj: Sphere object representing a galaxy
    :param mass_type: String encoding which mass to get here. The options are:
                      "initial" - just the initial stellar masses
                      "initial bound" - initial masses including initial bound fraction
                      "current" - current bound mass, accounting for the initial bound
                                  fraction, tidal disruption, and stellar death
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
        mass = data_obj[('STAR', 'INITIAL_MASS')].to("Msun").value
    elif mass_type == "initial bound":
        initial_mass = data_obj[('STAR', 'INITIAL_MASS')].to("Msun").value
        star_initial_bound = get_initial_bound_fraction(data_obj)
        mass = initial_mass * star_initial_bound
    elif mass_type == "current":
        raw_mass = data_obj[('STAR', 'INITIAL_MASS')].to("Msun").value
        star_initial_bound = get_initial_bound_fraction(data_obj)
        tidal_bound_fraction = data_obj[('STAR', 'BOUND_FRACTION')].value
        mass = raw_mass * star_initial_bound * tidal_bound_fraction

    # then restrict to recently formed clusters. This can be set to infinity, which
    # plots everything
    max_age = max_age_myr * yt.units.Myr
    mask = data_obj[('STAR', 'age')] < max_age
    mass = mass[mask]
    
    # create bins with spacing of 0.16 dex
    bin_width = 0.16  # dex
    m_boundaries_log = np.arange(3-0.5*bin_width, 7, 0.16)
    m_centers_log = [np.mean([m_boundaries_log[idx], m_boundaries_log[idx+1]])
                     for idx in range(len(m_boundaries_log)-1)]
    
    m_boundaries = 10**m_boundaries_log
    m_centers = 10**np.array(m_centers_log)
    
    # then make the histogram showing how many there are per bin
    hist, edges = np.histogram(mass, bins=m_boundaries)
    assert np.array_equiv(m_boundaries, edges)
    
    # We have dN, make it per dLogM
    hist = np.array(hist) / (bin_width * np.log(10))
    
    return m_centers, hist

def plot_cimf(ds_dict, halos_dict, plot_name_suffix, masses_to_plot, max_age_myr=np.inf):
    """
    Plot the initial cluster mass function.

    :param ds_dict: dictionary with keys of simulation names and values of the
                    snapshots for this simulations
    :param halos_dict: similar to ds_dict, except it contains the list of halo
                       objects to be plotted for this simulation
    :param plot_name_suffix: Either "last" or "common", depending on which
                             dictionaries are passed in. This will determine if
                             the redshift is labeled as common to all, or 
                             individually per simulation
    :param masses_to_plot: A lit of the masses to plot, see the `cimf` function for the
                           allowed options
    :param max_age_myr: The maximum age to restrict the plot to. Is infinity as the
                        default, which plots all stars.
    """
    if plot_name_suffix not in ["last", "common"]:
        raise ValueError("bad plot_name_suffix")
    # add the age if it's not infinity
    if not np.isinf(max_age_myr):
        plot_name_suffix += f"{max_age_myr}myr"

    for split in [True, False]:
        if split:
            fig, axs = bpl.subplots(figsize=[18, 7], ncols=2)
        else:
            fig, ax = bpl.subplots(figsize=[9, 7])
            axs = [ax]

        for idx, name in enumerate(halos_dict):
            c = colors[name]
            if split:
                ax = axs[axis_number[name]]
            for halo in halos_dict[name]:
                center = [halo.quantities["particle_position_x"],
                          halo.quantities["particle_position_y"],
                          halo.quantities["particle_position_z"]]

                sphere = ds_dict[name].sphere(center=center, radius=(100, "kpc"))

                for mass_type in masses_to_plot:
                    mass_plot, dn_dlogM = cimf(sphere, mass_type, max_age_myr)

                    # make the label only for the biggest halo, and not for initial only
                    if halo.quantities["rank"] == 1 and mass_type != "initial":
                        # and include the redshift if it's different for each sim
                        if "last" in plot_name_suffix:
                            label = f"{name}: z = {1/ds_dict[name].scale_factor - 1:.1f}"
                        else:
                            label = name
                    else:
                        label = None

                    # have different line styles
                    lss = {"initial": ":", "initial bound": "-", "current": "-"}

                    ax.plot(mass_plot, dn_dlogM, c=c, ls=lss[mass_type], label=label)

        for ax in axs:
            # plot the guiding lines
            plot_power_law(ax, -2, 1E6, 3E6, 1E4)
            plot_power_law(ax, -3, 1E6, 3E6, 1E4)

            ax.legend(loc=1, fontsize=10)
            ax.set_yscale("log")
            ax.set_xscale("log")
            ax.set_limits(1E3, 1E7, 10, 1E7)

            # if there is a common redshift, annotate it
            if "common" in plot_name_suffix:
                ax.easy_add_text(f"z = {1/common_scale - 1:.1f}", "upper left")

            ax.add_labels("$f_i$M [$M_\odot$]", "dN/dlogM")

        name = f"./comparison_plots/cimf_{plot_name_suffix}"
        if 'current' in masses_to_plot:
            name += "_current"
        if split:
            name += "_split"
        name += ".png"
        fig.savefig(name)

# then actually call this function to build the plots
plot_cimf(common_ds, common_halos, "common", ["initial bound", "initial"])
plot_cimf(last_ds, last_halos, "last", ["initial bound", "initial"])
plot_cimf(common_ds, common_halos, "common", ["initial bound", "initial"], 100)
plot_cimf(last_ds, last_halos, "last", ["initial bound", "initial"], 100)
plot_cimf(common_ds, common_halos, "common", ["current"])
plot_cimf(last_ds, last_halos, "last", ["current"])
