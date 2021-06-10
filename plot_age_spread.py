"""
plot_age_spread.py

Creates a plot showing the comparative cluster age spreads of clusters in 
different outputs. This code is very similar to plot_cimf.py, some of it was
copied. 
"""
import sys
from pathlib import Path

import numpy as np
import yt
from yt.extensions.astro_analysis.halo_analysis.api import HaloCatalog

import betterplotlib as bpl

bpl.set_style()
yt.funcs.mylog.setLevel(50)  # ignore yt's output

from plot_utils import names, colors

sentinel = Path(sys.argv[1]).resolve()
plot_dir = sentinel.parent


def filename_to_scale_factor(filename):
    return float(filename[-10:-4])


def get_ds_and_halos(ds_path):
    """ get the dataset and corresponding halo file """
    halo_path = ds_path.replace("out/continuous_", "halos/halos_").replace(
        ".art", ".0.bin"
    )

    ds = yt.load(ds_path)
    halos_ds = yt.load(halo_path)

    # if we are the old IC set, we have one galaxy, otherwise two
    # check what kind of particles are present
    if ("N-BODY_0", "MASS") in ds.derived_field_list:
        n_galaxies_each = 2
    else:
        n_galaxies_each = 1

    # create halo catalog object
    hc = HaloCatalog(halos_ds=halos_ds, data_ds=ds)
    hc.create(save_halos=True, save_catalog=False)

    # Do some parsing of the halo catalogs
    # We get the indices that sort it. The reversing there makes the biggest
    # halos first, like we want.
    halo_masses = yt.YTArray(
        [halo.quantities["particle_mass"] for halo in hc.halo_list]
    )
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
for directory in sys.argv[2:]:
    directory = Path(directory)
    if directory not in names:
        print(f"Skipping {directory}")
        continue

    out_dir = directory / "out"

    all_snapshots = [
        file.name
        for file in out_dir.iterdir()
        if file.is_file()
        and str(file.name).endswith(".art")
        and str(file.name).startswith("continuous_a")
    ]
    last_snapshots.append(sorted(all_snapshots)[-1])

earliest_last_snapshot = sorted(last_snapshots)[0]
common_scale = filename_to_scale_factor(earliest_last_snapshot) + 0.001
# include fudge factor for scale comparisons (so 0.1801 and 0.1802 match)

# set up the dictionaries where we will store the datasets and halos
common_ds = dict()
last_ds = dict()
common_halos = dict()
last_halos = dict()

for directory in sys.argv[2:]:
    directory = Path(directory)
    if directory not in names:
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

    ds_last, halos_last = get_ds_and_halos(str(last_snapshot))
    last_ds[names[directory]] = ds_last
    last_halos[names[directory]] = halos_last

    # then the last one that's in common with the other simulations
    all_common_snapshots = [
        file
        for file in all_snapshots
        if filename_to_scale_factor(file.name) <= common_scale
    ]
    # if there are no snapshots early enough for this, don't add them
    if len(all_common_snapshots) > 0:
        last_common_snapshot = sorted(all_common_snapshots)[-1]

        # then get the dataset and halo objects
        ds_common, halos_common = get_ds_and_halos(str(last_common_snapshot))
        common_ds[names[directory]] = ds_common
        common_halos[names[directory]] = halos_common

# Then functions to calculate the age spreads
def time_units(ds, array):
    return ds._handle.tphys_from_tcode_array(array) * yt.units.year


def duration(region, idxs):
    end_time = time_units(region.ds, region[("STAR", "TERMINATION_TIME")][idxs])
    return end_time - region[("STAR", "creation_time")][idxs]


def ave_time(region, idxs):
    art_units_ave_age = region[("STAR", "AVERAGE_AGE")][idxs]
    art_units_birth = region[("STAR", "BIRTH_TIME")][idxs]

    ave_time = time_units(region.ds, art_units_birth + art_units_ave_age)
    return ave_time - region[("STAR", "creation_time")][idxs]


def age_spread(region, idxs):
    initial_mass = region[("STAR", "initial_mass")][idxs]
    age_spread = region[("STAR", "AGE_SPREAD")][idxs] * region.ds.arr(1, "code_mass**2")
    birth_time = region[("STAR", "BIRTH_TIME")][idxs]
    creation_time = region[("STAR", "creation_time")][idxs]

    time = time_units(region.ds, (initial_mass ** 2 / age_spread) + birth_time)

    return time - creation_time


def time_cumulative_hist(data_obj, time_func, idxs):
    """
    Make the cluster initial mass function.

    :param data_obj: Sphere object representing a galaxy
    :param time_func: the function to be used to calculate the specific
                      implementation of the age spread
    :returns: Two lists. The first is the list of age spreads, in sorted order.
              The second is the fraction of clusters with age spreads less or
              equal to this value.
    """
    spreads = np.sort(time_func(data_obj, idxs).to("Myr").value)

    ranks = np.linspace(1 / len(spreads), 1, len(spreads))

    return spreads, ranks


def plot_age_growth(ds_dict, halos_dict, plot_name_suffix):
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

    fig, axs = bpl.subplots(figsize=[13, 17], ncols=2, nrows=3)

    funcs = [duration, ave_time, age_spread]
    names = ["Duration", "Average Age", "Age Spread"]

    for ax_row, func, name in zip(axs, funcs, names):
        # add the labels here
        ax_row[0].add_labels(
            name + " [Myr]", "Cumulative Fraction", "M < $10^5 M_\odot$"
        )
        ax_row[1].add_labels(
            name + " [Myr]", "Cumulative Fraction", "M > $10^5 M_\odot$"
        )
        ax_row[0].set_limits(0, 5, 0, 1)
        ax_row[1].set_limits(0, 8, 0, 1)
        for idx, name in enumerate(halos_dict):
            c = colors[name]
            for halo in halos_dict[name]:
                center = [
                    halo.quantities["particle_position_x"],
                    halo.quantities["particle_position_y"],
                    halo.quantities["particle_position_z"],
                ]

                sphere = ds_dict[name].sphere(center=center, radius=(100, "kpc"))

                # we split by low and high mass clusters, and some quality cuts
                good_idx = np.where(sphere[("STAR", "AGE_SPREAD")] > 0)
                cluster_masses = sphere[("STAR", "initial_mass")]
                cluster_cut = 1e5 * yt.units.Msun
                high_mask = cluster_masses > cluster_cut
                idx_lo = np.where(cluster_masses <= cluster_cut)
                idx_hi = np.where(cluster_masses > cluster_cut)
                idx_lo = np.intersect1d(idx_lo, good_idx)
                idx_hi = np.intersect1d(idx_hi, good_idx)

                # make the label only for the biggest halo
                if halo.quantities["rank"] == 1:
                    # and include the redshift if it's different for each sim
                    if plot_name_suffix == "last":
                        label = f"{name}: z = {1/ds_dict[name].scale_factor - 1:.1f}"
                    else:
                        label = name
                else:
                    label = None

                if len(idx_lo) > 0:
                    ages_lo, fractions_lo = time_cumulative_hist(sphere, func, idx_lo)
                    ax_row[0].plot(ages_lo, fractions_lo, c=c, lW=2, label=label)
                if len(idx_hi) > 0:
                    ages_hi, fractions_hi = time_cumulative_hist(sphere, func, idx_hi)
                    ax_row[1].plot(ages_hi, fractions_hi, c=c, lw=2, label=label)

    # put a legend in the top left panel
    axs[0][0].legend(loc=4, fontsize=14)

    # if there is a common redshift, annotate it
    if plot_name_suffix == "common":
        axs[0][0].easy_add_text(f"z = {1/common_scale - 1:.1f}", "upper left")

    fig.savefig(plot_dir / f"age_spread_{plot_name_suffix}.png")


# then actually call this function to build the plots
plot_age_growth(last_ds, last_halos, "last")
plot_age_growth(common_ds, common_halos, "common")

sentinel.touch()
