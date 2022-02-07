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
from matplotlib import pyplot as plt
import betterplotlib as bpl

from utils import load_galaxies, plot_utils, run_attributes

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
        # some versions of the summary printed both code units and kpc, while the
        # new version is only kpc. Need to do a check for that here.
        if "code" in line:
            value = line.split()[-4]
            unit = line.split()[3].replace(",", "")
        else:
            value = line.split()[-2]
            unit = line.split()[-1]
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
            quantity = u.Quantity(line.split()[1], line.split()[2])
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
for checks_dir in sys.argv[2:]:
    for summary_path in [
        str(f)
        for f in Path(checks_dir).iterdir()
        if f.name.startswith("galaxy_summaries")
    ]:
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
    for full_dir in run_attributes.names.keys():
        if summary_path.startswith(str(full_dir / "checks")):
            binned_summaries[full_dir].append(summary_path)
# then sort them to get the summaries in order
for key in binned_summaries:
    binned_summaries[key] = sorted(binned_summaries[key])

# and have a function to get the scale factor from the filename
def get_scale_factor(summary_path):
    scale = summary_path.split("_")[-1]
    scale = scale.replace(".txt", "")
    scale = scale.replace("a", "")
    return float(scale)


sim_groups = load_galaxies.get_plot_names_dirs(binned_summaries.keys())

# =============================================================================
#
# helper functions for plotting
#
# =============================================================================
def get_quantities(summaries, rank, quantity, unit):
    quantities = []
    scale_factors = []
    for summary in summaries:
        info = parsed_summaries[summary]
        # then get the quanties of interest
        try:  # check that the output has the desired quantity
            quantities.append(info[rank][quantity].to(unit).value)
            scale_factors.append(get_scale_factor(summary))
        except KeyError:  # output does not have it
            continue
    return quantities, scale_factors


def plot_quantities(quantity, unit, plot_name, ax):
    for idx, sim_dir in enumerate(binned_summaries):
        if plot_name not in run_attributes.names[sim_dir]:
            continue
        color = run_attributes.colors[sim_dir]
        summaries = binned_summaries[sim_dir]
        # the number of halos to plot will not change throughout the history
        # of the simulation
        n_to_plot = get_n_halos_to_plot(summaries[0])
        for rank in range(1, n_to_plot + 1):
            quantities, scale_factors = get_quantities(summaries, rank, quantity, unit)
            # only add the label for rank 1, so I don't duplicate
            if rank == 1:
                label = run_attributes.names[sim_dir][plot_name]
            else:
                label = None

            # only plot if at least one output has the quantity desired
            if len(scale_factors) > 0:
                ax.plot(scale_factors, quantities, label=label, c=color)


def plot_two_quantities(quantity_x, unit_x, quantity_y, unit_y, plot_name, ax):
    for idx, sim_dir in enumerate(binned_summaries):
        if plot_name not in run_attributes.names[sim_dir]:
            continue
        color = run_attributes.colors[sim_dir]
        summaries = binned_summaries[sim_dir]
        # the number of halos to plot will not change throughout the history
        # of the simulation
        n_to_plot = get_n_halos_to_plot(summaries[0])
        for rank in range(1, n_to_plot + 1):
            quantities_x, _ = get_quantities(summaries, rank, quantity_x, unit_x)
            quantities_y, _ = get_quantities(summaries, rank, quantity_y, unit_y)

            # only add the label for rank 1, so I don't duplicate
            if rank == 1:
                label = run_attributes.names[sim_dir][plot_name]
            else:
                label = None

            # only plot if at least one output has the quantities desired
            if len(quantities_x) > 0 and len(quantities_x) == len(quantities_y):
                ax.plot(quantities_x, quantities_y, label=label, c=color)


def plot_molecular_gas_fraction(plot_name, ax):
    for idx, sim_dir in enumerate(binned_summaries):
        if plot_name not in run_attributes.names[sim_dir]:
            continue
        color = run_attributes.colors[sim_dir]
        summaries = binned_summaries[sim_dir]
        # the number of halos to plot will not change throughout the history
        # of the simulation
        n_to_plot = get_n_halos_to_plot(summaries[0])
        for rank in range(1, n_to_plot + 1):
            gas_h2, scale_factors = get_quantities(
                summaries, rank, "gas_mass_H2", u.Msun
            )
            gas_hI, scale_factors = get_quantities(
                summaries, rank, "gas_mass_HI", u.Msun
            )

            ratio = [h2 / (h2 + hI) for h2, hI in zip(gas_h2, gas_hI)]
            # only add the label for rank 1, so I don't duplicate
            if rank == 1:
                label = run_attributes.names[sim_dir][plot_name]
            else:
                label = None

            # only plot if at least one output has the quantity desired
            if len(scale_factors) > 0:
                ax.plot(scale_factors, ratio, label=label, c=color)


# =============================================================================
#
# Then go through and plot everything!
#
# =============================================================================
for group in sim_groups:
    # -----------------------------------------------------------------------------
    # halo mass plot
    # -----------------------------------------------------------------------------
    fig, ax = bpl.subplots()
    plot_quantities("virial_mass", u.Msun, group, ax)
    ax.set_yscale("log")
    plot_utils.add_legend(ax, fontsize=10)
    ax.add_labels("Scale Factor", "Virial Mass [M$_\odot$] ")
    fig.savefig(plot_dir / f"galaxy_comparison_{group}_virial_mass.pdf")
    # then remove figure for memory purposes
    plt.close(fig)

    # -----------------------------------------------------------------------------
    # stellar mass plot
    # -----------------------------------------------------------------------------
    fig, ax = bpl.subplots()
    plot_quantities("stellar_mass_30_kpc", u.Msun, group, ax)
    ax.set_yscale("log")
    plot_utils.add_legend(ax, fontsize=10)
    ax.add_labels("Scale Factor", "Stellar Mass [M$_\odot$] within 30 kpc")
    ax.set_limits(y_min=2e7)
    fig.savefig(plot_dir / f"galaxy_comparison_{group}_stellar_mass.pdf")
    plt.close(fig)

    # -----------------------------------------------------------------------------
    # mass-metallicity plot
    # -----------------------------------------------------------------------------
    fig, ax = bpl.subplots()
    plot_two_quantities(
        "stellar_mass_30_kpc",
        u.Msun,
        "metallicity_all_stars",
        u.dimensionless_unscaled,
        group,
        ax,
    )
    # add Kirby 2013 line. Technically Kirby measures [Fe/H], I'll have to redo the
    # summaries to calculate that instead. For now assuem log(Z/Z_sun) = [Fe/H]
    masses = np.logspace(3, 11, 5)
    fe_h_kirby = -1.69 + 0.3 * np.log10(masses / 1e6)
    ax.plot(masses, fe_h_kirby, c=bpl.almost_black, ls=":", label="Kirby+2013 (z=0)")
    ax.fill_between(
        masses, fe_h_kirby - 0.17, fe_h_kirby + 0.17, color="0.97", zorder=0
    )

    ax.set_xscale("log")
    plot_utils.add_legend(ax, fontsize=10, frameon=False)
    ax.add_labels(
        "Stellar Mass [M$_\odot$] within 30 kpc",
        "Mean Stellar Metallicity [log(Z/$Z_\odot$)]",
    )
    ax.set_limits(1e5, 3e10)
    fig.savefig(plot_dir / f"galaxy_comparison_{group}_mass_metallicity.pdf")
    plt.close(fig)

    # -----------------------------------------------------------------------------
    # gas masses plots
    # -----------------------------------------------------------------------------
    for gas_type in gas_types:
        fig, ax = bpl.subplots()
        plot_quantities(f"gas_mass_{gas_type}", u.Msun, group, ax)
        ax.set_yscale("log")
        plot_utils.add_legend(ax, fontsize=10)
        ax.add_labels(
            "Scale Factor", f"{gas_type} Gas Mass [M$_\odot$] within" + " R$_{vir}$"
        )
        if "H2" in gas_type:
            ax.set_limits(y_min=1e5)

        fig.savefig(plot_dir / f"galaxy_comparison_{group}_gas_mass_{gas_type}.pdf")
        plt.close(fig)

    # -----------------------------------------------------------------------------
    # neutral fraction plot
    # -----------------------------------------------------------------------------
    fig, ax = bpl.subplots()
    plot_molecular_gas_fraction(group, ax)
    ax.set_yscale("log")
    plot_utils.add_legend(ax, fontsize=10)
    ax.add_labels("Scale Factor", f"H2 / (HI + H2) within" + " R$_{vir}$")
    ax.set_limits(y_min=1e-4, y_max=1)
    fig.savefig(plot_dir / f"galaxy_comparison_{group}_molecular_fraction.pdf")
    plt.close(fig)

# =============================================================================
#
# Touch the sentinel file
#
# =============================================================================
sentinel_file.touch()
