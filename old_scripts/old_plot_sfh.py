"""
plots_sfh.py

Creates plots describing the star formation history of galaxies in an output.

Takes 2 required parameters.
1 - Location of the simulation output. Can be relative to the working directory
    where this code was called from, or an absolute path.
2 - Location of the simulation halo finder output. Can be relative to the 
    working directory where this code was called form, or an absolute path.
"""

import sys
import os
import warnings

import yt
from yt.extensions.astro_analysis.halo_analysis.api import HaloCatalog
import numpy as np
from astropy import cosmology
from astropy import units as u
from astropy.utils.exceptions import AstropyWarning
import colorcet
from matplotlib import colors
from matplotlib import cm
from matplotlib import patches

import abundance_matching
import betterplotlib as bpl

import utils

yt.funcs.mylog.setLevel(50)  # ignore yt's output
bpl.presentation_style()
bpl.presentation_style()  # for some reason this needs to be there twice

# initialize the abundance matching
um = abundance_matching.UniverseMachine()

# Check the parameters length
if len(sys.argv) > 3:  # one is the script name
    raise ValueError("Only two parameters can be passed to this script.")

# Load in the simulation and halo catalogs
ds_loc = os.path.abspath(sys.argv[1])
hc_loc = os.path.abspath(sys.argv[2])
scale_factor = ds_loc[-10:-4]
ds = yt.load(ds_loc)
halos_ds = yt.load(hc_loc)

# create halo catalog object
hc = HaloCatalog(halos_ds=halos_ds, data_ds=ds)
hc.add_callback("sphere", factor=0.5)
hc.create(save_halos=True, save_catalog=False)

# get the location of where to put the plots
sim_dir = os.path.dirname(ds_loc) + os.sep
plots_dir = sim_dir.replace("/out/", "/plots/")

# Initialize the cosmology object, used to put a redshift scale on the plots
H_0 = ds.artio_parameters["hubble"][0] * 100 * u.km / (u.Mpc * u.second)
omega_matter = ds.artio_parameters["OmegaM"][0]
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


# Do some parsing of the halo catalogs
# We get the indices that sort it. The reversing there makes the biggest halos
# first, like we want.
halo_masses = yt.YTArray([halo.quantities["particle_mass"] for halo in hc.halo_list])
rank_idxs = np.argsort(halo_masses)[::-1]
# Add this info to each halo object, and put the halos into a new sorted list,
# with the highest mass (lowest rank) halos first). We only keep ones with an
# appreciable amount of stellar mass.
halos = []
for rank, idx in enumerate(rank_idxs, start=1):
    halo = hc.halo_list[idx]
    halo.quantities["rank"] = rank
    # sometimes the halos can be too small for yt to handle. This seems to only
    # happen at extremely high redshift, before stars have formed, so we can
    # ignore it
    if halo.data_object is None:
        vr = halo.quantities["virial_radius"].to("pc")
        print("NO DATA OBJECT, VIRIAL RADIUS = {:.1f}".format(vr))
        continue
    sm = halo.data_object[("STAR", "INITIAL_MASS")].in_units("msun").sum()
    halo.quantities["stellar mass"] = sm
    if sm > 1e6 * yt.units.msun:
        halos.append(halo)

# Convenience function for annotating plots
def exponent_format(value, digits):
    exponent = np.floor(np.log10(value))
    constant = value / 10 ** exponent
    constant_str = "{number:.{digits}f}".format(number=constant, digits=digits)
    return constant_str + r"$\times 10^{" + str(int(exponent)) + r"}$"


# redshifts to label on plots later
label_redshifts = [10, 5, 3, 2, 1, 0.5, 0.3, 0.2, 0.1]

# =============================================================================
#
# Star Formation Histories
#
# =============================================================================
# Make plots showing the star formation histories of all galaxies with
# significant stellar mass, compared to abundance matching predictions
# first determine how big the figure should be, since we are including all
# galaxies in one figure
figsize = None
n_cols = 1
n_rows = 1
n_subplots = len(halos)
while figsize is None:
    if n_subplots <= n_rows * n_cols:
        figsize = [10 * n_cols, 7 * n_rows]
        break
    if n_rows > n_cols:
        n_cols += 1
    else:
        n_rows += 1

fig, axs = bpl.subplots(figsize=figsize, ncols=n_cols, nrows=n_rows, squeeze=False)
# squeeze being false there means axs is always an array, even if it's a single
# axis. We have to do this so flatten always works
axs = axs.flatten()

for ax, halo in zip(axs, halos):
    times, sfh_values = utils.sfh(halo.data_object)
    plot_times = times.to("Gyr").value
    plot_sfh = sfh_values.to("msun/yr").value
    dt = plot_times[1] - plot_times[0]
    ax.errorbar(plot_times, plot_sfh, xerr=0.5 * dt, markersize=8, c=bpl.almost_black)
    ax.plot(plot_times, plot_sfh, lw=1.0, c=bpl.almost_black)

    halo_mass = halo.quantities["particle_mass"].to("msun").value
    stellar_mass = halo.quantities["stellar mass"].to("msun").value
    # add this information
    msun = "$M_\odot$"
    ax.easy_add_text(
        "$M_*$ = "
        + exponent_format(stellar_mass, 2)
        + msun
        + "\n"
        + "$M_h$ = "
        + exponent_format(halo_mass, 2)
        + msun
        + "\n"
        + "z = {:.2f}".format(ds.current_redshift),
        "upper left",
    )

    # start figuring out the upper limit of the plot, which we only want to be
    # based on the data and matched abundance matching predictions
    upper_limit = 2 * plot_sfh.max()
    # then compare to abundance matching
    for name, mass in zip(["stellar", "halo"], [stellar_mass, halo_mass]):
        zs, sfhs, hi_lim, lo_lim = um.get_sfh(name, ds.current_redshift, mass)
        ages = [z_to_age(z).to("Gyr").value for z in zs]
        ax.fill_between(
            x=ages, y1=lo_lim, y2=hi_lim, alpha=0.4, label=name + " matched"
        )

        # see if we should change the max value
        if len(hi_lim) > 0:  # sometimes predictions don't exist?
            upper_limit = max(upper_limit, 2 * max(hi_lim))

    # compare to Milky Way prediction
    zs, sfhs, hi_lim, lo_lim = um.get_sfh("halo", 0, 1e12)
    ages = [z_to_age(z).to("Gyr").value for z in zs]
    ax.fill_between(x=ages, y1=lo_lim, y2=hi_lim, alpha=0.4, label="MW-like")

    # I want to set the lower limit to be the lowest non-zero sf
    lower_limit = 0.5 * min(plot_sfh[np.where(plot_sfh > 0)])

    ax.legend(loc=4)
    ax.set_yscale("log")
    ax.set_limits(0, max(plot_times) * 1.01, lower_limit, upper_limit)
    ax.add_labels("Time [Gyr]", "SFR  [$M_\odot$/yr]")

    # then add the redshift axis. The process of selecting the labels raises
    # warnings, so we can ignore that
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        ax.twin_axis("x", label_redshifts, "Redshift", new_to_old_func=z_to_age_Gyr)

# turn off any unused axes
for idx in range(len(halos), len(axs)):
    axs[idx].set_axis_off()
# only save the plot if there are actually galaxies.
if len(halos) > 0:
    fig.savefig(plots_dir + "sfh_a{}.png".format(scale_factor))

# =============================================================================
#
# Cumulative Stellar Mass histories
#
# =============================================================================
# Have the cumulative stellar mass for all galaxies on one plot
# I experimented with how to color the points. I used to use a colormap, but
# eventually went with colorcet's list of maximally different colors. I kept the
# old code as a reference in case I want to switch back
# cmap = cmocean.cm.thermal
# c_norm = colors.LogNorm(vmin=1E8, vmax=2E12)
# cmap = colorcet.cm.colorwheel
# cmap = colorcet.cm.rainbow
# c_norm = colors.Normalize(0, len(halos))
# scalar_map = cm.ScalarMappable(norm=c_norm, cmap=cmap)
# scalar_map.set_array([])

fig, ax = bpl.subplots()

for idx, halo in enumerate(halos):
    timesteps, mass = utils.create_cumulative_mass(halo.data_object)
    #     color = scalar_map.to_rgba(halo.quantities["particle_mass"].to("msun").value)
    #     color = scalar_map.to_rgba(idx)
    color = colorcet.glasbey_dark[idx]
    ax.plot(timesteps.to("Gyr").value, mass.to("msun").value, c=color)

ax.set_yscale("log")
ax.set_limits(0, ds.current_time.to("Gyr").value * 1.01)
ax.add_labels("Time [Gyr]", "Stellar Mass [$M_\odot$]")
# then add the redshift axis. The process of selecting the labels raises
# warnings, so we can ignore that
with warnings.catch_warnings():
    warnings.simplefilter("ignore", UserWarning)
    ax.twin_axis("x", label_redshifts, "Redshift", new_to_old_func=z_to_age_Gyr)
# cb = fig.colorbar(scalar_map)
# cb.set_label("Halo Mass [$M_\odot$] at z={:.1f}".format(ds.current_redshift))
# only save the plot if there are actually galaxies.
if len(halos) > 0:
    fig.savefig(plots_dir + "stellar_mass_a{}.png".format(scale_factor))

# =============================================================================
#
# Stellar Mass - Halo Mass relationship
#
# =============================================================================
# Plot all galaxies on the SMHM
fig, ax = bpl.subplots()
sms = [halo.quantities["stellar mass"].to("msun").value for halo in halos]
hms = [halo.quantities["particle_mass"].to("msun").value for halo in halos]

scatter = ax.scatter(hms, sms, c=bpl.almost_black, s=100)

# overplot SMHM relation from observations
z_smhm, masses, smhm, hi_lim, lo_lim = um.get_smhm(ds.current_redshift, "All")

# make the plot labels
um_label = "Universe Machine\nSMHM: z={:.1f}".format(abs(z_smhm))
sim_label = "Simulation: z={:.1f}".format(ds.current_redshift)

s_mass = masses * smhm
sm_up = masses * hi_lim
sm_down = masses * lo_lim

um_color = bpl.color_cycle[3]
line = ax.plot(masses, s_mass, c=um_color)[0]
ax.fill_between(x=masses, y1=sm_down, y2=sm_up, color=um_color, lw=0, alpha=0.1)
# two ways of having the shading in the legend. Either create a patch and put
# it in (uncommented option) or just manually put a shaded region in the
# background. I found that the second option behaved diferrently for pdf and png
# (don't know why!) so I went with the first option.
patch = patches.Patch(color=um_color, lw=0, alpha=0.1)
# ax.fill_between(x=[1.2724E9, 2.077E9], y1=[3.1E10, 3.1E10], y2=[7.5E10, 7.5E10],
#                 color=um_color, lw=0, alpha=0.1)
ax.legend(
    handles=((line, patch), scatter), labels=(um_label, sim_label), loc=2, frameon=False
)

ax.set_xscale("log")
ax.set_yscale("log")
ax.set_limits(1e9, 2e12, 1e5, 2e11)
ax.add_labels("Halo Mass [$M_\odot$]", "Stellar Mass [$M_\odot$]")

# save this even if there aren't galaxies, so we can have a sentinel file
fig.savefig(plots_dir + "smhm_a{}.png".format(scale_factor))
