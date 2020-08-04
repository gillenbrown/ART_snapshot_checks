"""
plot_cimf.py

Creates a plot showing the comparative cluster initial mass function of 
galaxies in an output
"""
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy import special
import yt
from yt.extensions.astro_analysis.halo_analysis.api import HaloCatalog

import betterplotlib as bpl

bpl.set_style()
yt.funcs.mylog.setLevel(50)  # ignore yt's output

# I have to hardcode some labels to make this easier, parsing them won't work
# nearly as well
def full_dir(partial_path):
    base_dir = Path("/u/home/gillenb/art_runs/runs/").absolute()
    return base_dir / partial_path

names = {full_dir("shangrila/old_ic_comparison/default/run"): "ART 2.0 SFE100",
         full_dir("shangrila/old_ic_comparison/default_1e7_temp_cap/run"): "ART 2.0 SFE100 Bad Caps",
         full_dir("shangrila/hui/sfe_10"): "NBm SFE10",
         full_dir("shangrila/hui/sfe_50"): "NBm SFE50",
         full_dir("shangrila/hui/sfe_100"): "NBm SFE100",
         full_dir("shangrila/hui/sfe_200"): "NBm SFE200",
         full_dir("stampede2/production/sfe100"): "T&L SFE100",
         full_dir("stampede2/production/first_sfe_100_1e7_temp_cap"): "T&L SFE100 Bad Caps",
         full_dir("stampede2/ic_timing_tests/original_50_128"): "T&L - Original",
         full_dir("stampede2/ic_timing_tests/trim_12_128"): "T&L - 4x Trim 128",
         full_dir("stampede2/ic_timing_tests/trim_12_256"): "T&L - 4x Trim 256",
         full_dir("stampede2/ic_timing_tests/trim_25_256"): "T&L - 2x Trim 256"}
    
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
    out_dir = directory / "out"
    
    all_snapshots = [file.name for file in out_dir.iterdir()
                     if file.is_file()
                     and str(file.name).endswith(".art")
                     and str(file.name).startswith("continuous_a")]
    last_snapshots.append(sorted(all_snapshots)[-1])
    
earliest_last_snapshot = sorted(last_snapshots)[0]
common_scale = filename_to_scale_factor(earliest_last_snapshot) + 0.001  
# include fudge factor for scale comparisons (so 0.1801 and 0.1802 match)

# set up the dictionaries where we will store the datasets and halos
common_ds = dict()
last_ds = dict()
common_halos = dict()
last_halos = dict()
# and one to ensure each has the same color
colors = dict()

for idx, directory in enumerate(sys.argv[1:]):
    directory = Path(directory)
    colors[names[directory]] = bpl.color_cycle[idx]

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
                            if filename_to_scale_factor(file.name) <= common_scale]
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

def cimf(data_obj, include_initial_bound):
    """
    Make the cluster initial mass function. 
    
    :param data_obj: Sphere object representing a galaxy
    :param include_initial_bound: whether to incorporate the initial bound
                                  fraction of clusters, or just get the 
                                  distribution of initial particle masses
    :returns: Two lists. The first is f_i * M, representing the initial
              bound mass, of M if include_initial_bound=False. This will be 
              binned values suitable to plot. The second is dN/dLogM for each 
              of the bins in the first list.
    """
    star_initial_mass = data_obj[('STAR', 'INITIAL_MASS')].to("Msun").value
    if include_initial_bound:
        star_initial_bound = get_initial_bound_fraction(data_obj)
        star_initial_mass *= star_initial_bound
    
    # create bins with spacing of 0.16 dex
    bin_width = 0.16  # dex
    m_boundaries_log = np.arange(3-0.5*bin_width, 7, 0.16)
    m_centers_log = [np.mean([m_boundaries_log[idx], m_boundaries_log[idx+1]])
                     for idx in range(len(m_boundaries_log)-1)]
    
    m_boundaries = 10**m_boundaries_log
    m_centers = 10**np.array(m_centers_log)
    
    # then make the histogram showing how many there are per bin
    hist, edges = np.histogram(star_initial_mass, bins=m_boundaries)
    assert np.array_equiv(m_boundaries, edges)
    
    # We have dN, make it per dLogM
    hist = np.array(hist) / (bin_width * np.log(10))
    
    return m_centers, hist

def plot_cimf(ds_dict, halos_dict, plot_name_suffix):
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
    """
    if plot_name_suffix not in ["last", "common"]:
        raise ValueError("bad plot_name_suffix")
    
    fig, ax = bpl.subplots(figsize=[9, 7])

    for idx, name in enumerate(halos_dict):
        c = colors[name]
        for halo in halos_dict[name]:
            center = [halo.quantities["particle_position_x"],
                      halo.quantities["particle_position_y"],
                      halo.quantities["particle_position_z"]]
            
            sphere = ds_dict[name].sphere(center=center, radius=(100, "kpc"))
            mass_plot_with_bound, dn_dlogM_with_bound = cimf(sphere, True)
            mass_plot_without_bound, dn_dlogM_without_bound = cimf(sphere, False)
            
            # make the label only for the biggest halo
            if halo.quantities["rank"] == 1:
                # and include the redshift if it's different for each sim
                if plot_name_suffix == "last":
                    label = f"{name}: z = {1/ds_dict[name].scale_factor - 1:.1f}"
                else:
                    label = name
            else:
                label = None

            ax.plot(mass_plot_with_bound, dn_dlogM_with_bound, c=c, ls="-", label=label)
            ax.plot(mass_plot_without_bound, dn_dlogM_without_bound, c=c, ls=":")
            
    # plot the guiding lines
    plot_power_law(ax, -2, 1E6, 3E6, 1E4)
    plot_power_law(ax, -3, 1E6, 3E6, 1E4)

    ax.legend(loc=1, fontsize=14)
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.set_limits(1E3, 1E7, 10, 1E7)
    
    # if there is a common redshift, annotate it
    if plot_name_suffix == "common":
        ax.easy_add_text(f"z = {1/common_scale - 1:.1f}", "upper left")

    ax.add_labels("$f_i$M [$M_\odot$]", "dN/dlogM")

    fig.savefig(f"./comparison_plots/cimf_{plot_name_suffix}.png")

# then actually call this function to build the plots
plot_cimf(common_ds, common_halos, "common")
plot_cimf(last_ds, last_halos, "last")
