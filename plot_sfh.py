"""
plots_sfh.py

Creates a plot showing the comparative SFH of galaxies in an output

The parameters passed to this script must be directories with the 
"""
import sys
import os
from collections import defaultdict

import numpy as np
from yt.extensions.astro_analysis.halo_analysis.api import HaloCatalog
from astropy import cosmology
from astropy import units as u

import abundance_matching
um = abundance_matching.UniverseMachine()

import betterplotlib as bpl
bpl.set_style()
bpl.color_cycle.append("0.5")

import yt
yt.funcs.mylog.setLevel(50)  # ignore yt's output

# I have to hardcode some labels to make this easier, parsing them won't work
# nearly as well
def full_dir(partial_path):
    base_dir = "/u/home/gillenb/art_runs/runs/"
    return os.path.join(base_dir, partial_path)

names = {full_dir("shangrila/old_ic_comparison/test_default/run"): "ART 2.0 Not Final SFE100",
         full_dir("shangrila/old_ic_comparison/default/run"): "ART 2.0 Final SFE100",
         full_dir("shangrila/hui/sfe_100"): "NBm SFE100",
         full_dir("stampede2/ic_timing_tests/original_50_128"): "T&L - Original",
         full_dir("stampede2/ic_timing_tests/trim_12_128"): "T&L - 4x Trim 128",
         full_dir("stampede2/ic_timing_tests/trim_12_256"): "T&L - 4x Trim 256",
         full_dir("stampede2/ic_timing_tests/trim_25_256"): "T&L - 2x Trim 256"}
    
# make dictionary to store the resulting datasets
all_ds = dict()
all_halos = defaultdict(list)
for directory in sys.argv[1:]:
    out_dir = os.path.join(directory, "out")
    halos_dir = os.path.join(directory, "halos")

    last_out = sorted([file for file in  os.listdir(out_dir)
                       if file.endswith(".art")
                       and file.startswith("continuous_a")])[-1]

    last_halo = last_out.replace("continuous_", "halos_").replace(".art", ".0.bin")


    ds = yt.load(os.path.join(out_dir, last_out))
    halos_ds = yt.load(os.path.join(halos_dir, last_halo))

    all_ds[names[directory]] = ds

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
        all_halos[names[directory]].append(halo)


# =============================================================================
#         
# Set up cosmology
# 
# =============================================================================
temp_ds = all_ds[list(all_ds.keys())[0]]
# Initialize the cosmology object, used to put a redshift scale on the plots
h = temp_ds.artio_parameters["hubble"][0] 
H_0 = h * 100 * u.km / (u.Mpc * u.second)
omega_matter = temp_ds.artio_parameters["OmegaM"][0]
cosmo = cosmology.FlatLambdaCDM(H0=H_0, Om0=omega_matter, Tcmb0=2.725)
# and functions to transfer between redshift and age
def age_to_z(age):
    try:
        return cosmology.z_at_value(cosmo.age, age)
    except cosmology.funcs.CosmologyError:  # happens at very high z
        return 1000
def z_to_age_Gyr(z):
    # this is needed since the twin axis function always works with a scalar
    return z_to_age(z).to("Gyr").value
def z_to_age(z):
    return cosmo.age(z).to("Gyr")

# =============================================================================
#         
# Function that actually does this calculation
# 
# =============================================================================
def sfh(data_obj):
    masses = data_obj[("STAR", "INITIAL_MASS")].in_units("msun")
    creation_times = data_obj[("STAR", "creation_time")]
    
    # have the bins start with the last time, so that the final bin is not
    # strange due to only being incomplete
    dt = 100 * yt.units.Myr
    max_bin = data_obj.ds.current_time.to("Myr")
    bin_edges = [max_bin]
    while bin_edges[-1] > 0:
        bin_edges.append(bin_edges[-1] - dt)
    # then reverse it, so the biggest bin is at the end
    bin_edges = bin_edges[::-1]

    bin_centers = []
    sfr_values = []
    # go through each bin, seeing which stars were born there
    for left_idx in range(len(bin_edges) - 1):
        # get the age range
        min_age = bin_edges[left_idx]
        max_age = bin_edges[left_idx + 1]
        
        # then get the stars in that range
        age_more_idx = np.where(creation_times >= min_age)
        age_less_idx = np.where(creation_times <  max_age)
        age_good_idx = np.intersect1d(age_more_idx, age_less_idx)
        
        # and sum their masses
        this_formed_mass = np.sum(masses[age_good_idx])
        
        sfr_values.append(this_formed_mass / dt)
        bin_centers.append((min_age + max_age) / 2.0)
        
    return yt.YTArray(bin_centers), yt.YTArray(sfr_values)

# =============================================================================
#         
# Then plot!
# 
# =============================================================================
label_redshifts = [10, 5, 3, 2, 1, 0.5, 0.3, 0.2, 0.1]

fig, ax = bpl.subplots()
# store data about times
max_time = 0
for idx, name in enumerate(all_halos):
    c = bpl.color_cycle[idx]
    for halo in all_halos[name]:
        center = [halo.quantities["particle_position_x"],
                  halo.quantities["particle_position_y"],
                  halo.quantities["particle_position_z"]]
        sphere = all_ds[name].sphere(center=center, radius=(30, "kpc"))
        times, sfh_values = sfh(sphere)

        plot_times = times.to("Gyr").value
        plot_sfh = sfh_values.to("msun/yr").value
        dt = plot_times[1] - plot_times[0]
        if halo.quantities["rank"] == 1:
            label = name
        else:
            label = None
        ax.errorbar(plot_times, plot_sfh, xerr=0.5 * dt, markersize=8,
                    c=c, label=label)
        ax.plot(plot_times, plot_sfh, lw=1.0, c=c)

        # figure out the max time to use for the plot limit
        if max(plot_times) > max_time:
            max_time = max(plot_times)
        
# compare to Milky Way prediction
zs, sfhs, hi_lim, lo_lim = um.get_sfh("halo", 0, 1E12)
ages = [z_to_age(z).to("Gyr").value for z in zs]
ax.fill_between(x=ages, y1=lo_lim, y2=hi_lim, alpha=0.4, lw=0,
                color="0.3", label="MW-like")

ax.legend(loc=2)
ax.set_yscale("log")
ax.set_limits(0, 1.05*max_time, y_min=1E-1)
ax.add_labels("Time [Gyr]", "SFR  [$M_\odot$/yr]")

# then add the redshift axis. The process of selecting the labels raises
# warnings, so we can ignore that
# with warnings.catch_warnings():
#     warnings.simplefilter('ignore', UserWarning)
ax.twin_axis("x", label_redshifts, "Redshift", 
             new_to_old_func=z_to_age_Gyr)

fig.savefig("./comparison_plots/sfh_comparison.png")


# =============================================================================
#         
# Compare the cumulative mass - a very similar calculation
# 
# =============================================================================
def create_cumulative_mass(data_obj):
    """
    Create the cumulative stellar mass of the stars within some region.

    This will create many timesteps, then for each timestep record the mass
    that formed earlier than this time.

    :param data_obj: region of the simulation that stars will be selected from.
                     Can be something like a sphere, or even all_data()
    """
    masses = data_obj[("STAR", "INITIAL_MASS")].in_units("msun")
    creation_times = data_obj[("STAR", "creation_time")]

    sort_idxs = np.argsort(creation_times.to("Myr").value)

    times = creation_times[sort_idxs]
    mass_in_order = masses[sort_idxs]
    mass_cumulative = np.cumsum(mass_in_order)
    
    return times.in_units("Gyr"), yt.YTArray(mass_cumulative)

# =============================================================================
#         
# Then plot!
# 
# =============================================================================
label_redshifts = [10, 5, 3, 2, 1, 0.5, 0.3, 0.2, 0.1]

fig, ax = bpl.subplots()
# store data about times
max_time = 0
for idx, name in enumerate(all_halos):
    c = bpl.color_cycle[idx]
    for halo in all_halos[name]:
        center = [halo.quantities["particle_position_x"],
                  halo.quantities["particle_position_y"],
                  halo.quantities["particle_position_z"]]
        sphere = all_ds[name].sphere(center=center, radius=(30, "kpc"))
        times, cumulative_mass = create_cumulative_mass(sphere)

        plot_times = times.to("Gyr").value
        cumulative_mass = cumulative_mass.to("msun").value
        dt = plot_times[1] - plot_times[0]
        if halo.quantities["rank"] == 1:
            label = name
        else:
            label = None
        ax.plot(plot_times, cumulative_mass, c=c, label=label)

        # figure out the max time to use for the plot limit
        if max(plot_times) > max_time:
            max_time = max(plot_times)

ax.legend()
ax.set_yscale("log")
ax.set_limits(0, 1.05*max_time, y_min=1E6)
ax.add_labels("Time [Gyr]", "Stellar Mass  [$M_\odot$]")

# then add the redshift axis. The process of selecting the labels raises
# warnings, so we can ignore that
# with warnings.catch_warnings():
#     warnings.simplefilter('ignore', UserWarning)
ax.twin_axis("x", label_redshifts, "Redshift", 
             new_to_old_func=z_to_age_Gyr)

fig.savefig("./comparison_plots/mass_growth_comparison.png")


