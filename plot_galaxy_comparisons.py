"""
plot_galaxy_comparisons.py

Creates many plots showing time series comparisons of all the properties in the
galaxy summary files.

This script takes the list of all summary files.
"""

import sys
from pathlib import Path
from collections import defaultdict
import numpy as np
from astropy import units as u
import betterplotlib as bpl

from utils import load_galaxies

bpl.set_style()

sentinel_file = Path(sys.argv[1])
plot_dir = sentinel_file.parent

# =============================================================================
#
# functions to read the summary files
#
# =============================================================================
def check_for_rank(line, halo_dict):
    if line.startswith("Rank ") and line.strip().endswith(" halo:"):
        halo_dict["rank"] = int(line.split()[1])


def check_for_virial_mass(line, halo_dict):
    if line.startswith("Virial Mass: "):
        quantity = u.Quantity(line.split()[-2], line.split()[-1])
        halo_dict["virial_mass"] = quantity


def check_for_virial_radius(line, halo_dict):
    if line.startswith("Virial Radius: "):
        # have code units here, so need to get those, and throw awway the comma
        value = line.split()[-4]
        unit = line.split()[-3].replace(",", "")
        quantity = u.Quantity(value, unit)
        halo_dict["virial_radius"] = quantity


def check_for_position(line, halo_dict):
    if line.startswith("Position"):
        dimension = line.split()[1].replace(":", "")
        # have code units here, so need to get those, and throw awway the comma
        value = line.split()[-4]
        unit = line.split()[-3].replace(",", "")
        quantity = u.Quantity(value, unit)
        halo_dict[f"position_{dimension}"] = quantity


def check_for_stellar_mass(line, halo_dict):
    start = "Stellar Mass within "
    if line.startswith(start):
        radius_line = line.replace(start, "")
        radius_line = radius_line.replace("the ", "")
        radius_line = radius_line.replace(" ", "_")
        radius = radius_line.split(":")[0]
        quantity = u.Quantity(line.split()[-2], line.split()[-1])
        halo_dict[f"stellar_mass_{radius}"] = quantity


gas_types = ["Total", "HI", "HII", "H2", "He1", "HeII", "HeIII", "Metals"]


def check_for_gas_masses(line, halo_dict):
    for gas_type in gas_types:
        if line.startswith(gas_type + ": "):
            quantity = u.Quantity(line.split()[-2], line.split()[-1])
            halo_dict[f"gas_mass_{gas_type}"] = quantity


def check_for_metallicities(line, halo_dict):
    if " -> log(Z/Z_sun) = " in line:
        met_type = line.split("=")[0].strip()
        met_type = met_type.replace(" ", "_")
        met_value = u.Quantity(float(line.split()[-1]), u.dimensionless_unscaled)
        halo_dict[f"metallicity_{met_type}"] = met_value


def read_summary_file(file_loc, n_to_read):
    halos = dict()
    # set up the placeholder dict which will be filled later
    halo_dict = None
    with open(file_loc, "r") as summary:
        for line in summary:
            # first see if we find the separator. If so, we reset our dictionary
            if set(line.strip()) == {"="}:
                # if this isn't the first separator in the file, add the
                # existing halo to our list
                if halo_dict is not None:
                    halos[halo_dict["rank"]] = halo_dict

                    # if we have enough, quit
                    if len(halos) == n_to_read:
                        return halos
                # then reset it
                halo_dict = dict()
            # Then I have a bunch of separate functions that look for
            # items in each line. If they don't find what they're looking for,
            # they leave the dict unmodified
            check_for_rank(line, halo_dict)
            check_for_virial_mass(line, halo_dict)
            check_for_virial_radius(line, halo_dict)
            check_for_position(line, halo_dict)
            check_for_stellar_mass(line, halo_dict)
            check_for_gas_masses(line, halo_dict)
            check_for_metallicities(line, halo_dict)


# =============================================================================
#
# Then use these to parse the summaries
#
# =============================================================================
def get_n_halos_to_plot(summary_path):
    if "stampede2/production" in summary_path:
        return 2
    else:
        return 1


parsed_summaries = dict()
for summary_path in sys.argv[2:]:
    n_to_read = get_n_halos_to_plot(summary_path)

    parsed_summary = read_summary_file(summary_path, n_to_read)
    # some of these are None if there are no halos.
    if parsed_summary is not None:
        parsed_summaries[summary_path] = parsed_summary

# =============================================================================
#
# Then sort them based on their simulation of origin
#
# =============================================================================
binned_summaries = defaultdict(list)
for summary_path in parsed_summaries:
    for full_dir, short_name in load_galaxies.names.items():
        if summary_path.startswith(str(full_dir / "checks")):
            binned_summaries[short_name].append(summary_path)
# then sort them to get the summaries in order
for key in binned_summaries:
    binned_summaries[key] = sorted(binned_summaries[key])

# and have a function to get the scale factor from the filename
def get_scale_factor(summary_path):
    scale = summary_path.split("_")[-1]
    scale = scale.replace(".txt", "")
    scale = scale.replace("a", "")
    return float(scale)


# =============================================================================
#
# Then go through and plot everything!
#
# =============================================================================
# define some helper functions first
def plot_quantities(quantity, unit, ax):
    for idx, sim_name in enumerate(binned_summaries):
        color = load_galaxies.colors[sim_name]
        summaries = binned_summaries[sim_name]
        # the number of halos to plot will not change throughout the history
        # of the simulation
        n_to_plot = get_n_halos_to_plot(summaries[0])
        for rank in range(1, n_to_plot + 1):
            scale_factors = []
            quantities = []
            for summary in summaries:
                info = parsed_summaries[summary]
                # then get the quanties of interest
                try:  # check that the output has the desired quantity
                    quantities.append(info[rank][quantity].to(unit).value)
                    scale_factors.append(get_scale_factor(summary))
                except KeyError:  # output does not have it
                    continue
            # only add the label for rank 1, so I don't duplicate
            if rank == 1:
                label = sim_name
            else:
                label = None

            # only plot if at least one output has the quantity desired
            if len(scale_factors) > 0:
                ax.plot(scale_factors, quantities, label=label, c=color)


def plot_two_quantities(quantity_x, unit_x, quantity_y, unit_y, ax):
    for idx, sim_name in enumerate(binned_summaries):
        color = load_galaxies.colors[sim_name]
        summaries = binned_summaries[sim_name]
        # the number of halos to plot will not change throughout the history
        # of the simulation
        n_to_plot = get_n_halos_to_plot(summaries[0])
        for rank in range(1, n_to_plot + 1):
            quantities_x = []
            quantities_y = []
            for summary in summaries:
                info = parsed_summaries[summary]
                # then get the quanties of interest
                try:  # check that the output has the desired quantity
                    quantities_x.append(info[rank][quantity_x].to(unit_x).value)
                    quantities_y.append(info[rank][quantity_y].to(unit_y).value)
                except KeyError:  # output does not have it
                    continue
            # only add the label for rank 1, so I don't duplicate
            if rank == 1:
                label = sim_name
            else:
                label = None

            # only plot if at least one output has the quantities desired
            if len(quantities_x) > 0:
                ax.plot(quantities_x, quantities_y, label=label, c=color)


# -----------------------------------------------------------------------------
# halo mass plot
# -----------------------------------------------------------------------------
fig, ax = bpl.subplots()
plot_quantities("virial_mass", u.Msun, ax)
ax.set_yscale("log")
ax.legend(fontsize=10)
ax.add_labels("Scale Factor", "Virial Mass [M$_\odot$] ")
fig.savefig(plot_dir / "galaxy_comparison_virial_mass.png")

# -----------------------------------------------------------------------------
# stellar mass plot
# -----------------------------------------------------------------------------
fig, ax = bpl.subplots()
plot_quantities("stellar_mass_30_kpc", u.Msun, ax)
ax.set_yscale("log")
ax.legend(fontsize=10)
ax.add_labels("Scale Factor", "Stellar Mass [M$_\odot$] within 30 kpc")
ax.set_limits(y_min=2e7)
fig.savefig(plot_dir / "galaxy_comparison_stellar_mass.png")

# -----------------------------------------------------------------------------
# mass-metallicity plot
# -----------------------------------------------------------------------------
fig, ax = bpl.subplots()
plot_two_quantities(
    "stellar_mass_30_kpc", u.Msun, "metallicity_all_stars", u.dimensionless_unscaled, ax
)
# add Kirby 2013 line. Technically Kirby measures [Fe/H], I'll have to redo the
# summaries to calculate that instead. For now assuem log(Z/Z_sun) = [Fe/H]
masses = np.logspace(3, 11, 5)
fe_h_kirby = -1.69 + 0.3 * np.log10(masses / 1e6)
ax.plot(masses, fe_h_kirby, c=bpl.almost_black, ls=":", label="Kirby+2013 (z=0)")
ax.fill_between(masses, fe_h_kirby - 0.17, fe_h_kirby + 0.17, color="0.97", zorder=0)

ax.set_xscale("log")
ax.legend(fontsize=10, frameon=False)
ax.add_labels(
    "Stellar Mass [M$_\odot$] within 30 kpc",
    "Mean Stellar Metallicity [log(Z/$Z_\odot$)]",
)
ax.set_limits(1e5, 3e10)
fig.savefig(plot_dir / "galaxy_comparison_mass_metallicity.png")

# -----------------------------------------------------------------------------
# gas masses plots
# -----------------------------------------------------------------------------
for gas_type in gas_types:
    fig, ax = bpl.subplots()
    plot_quantities(f"gas_mass_{gas_type}", u.Msun, ax)
    ax.set_yscale("log")
    ax.legend(fontsize=10)
    ax.add_labels("Scale Factor", f"{gas_type} Gas Mass [M$_\odot$] within 30 kpc")
    fig.savefig(plot_dir / f"galaxy_comparison_gas_mass_{gas_type}.png")

# =============================================================================
#
# Touch the sentinel file
#
# =============================================================================
sentinel_file.touch()
