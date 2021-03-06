"""
plot_single_halo_full_nbody.py

Creates a projection plot of a single DM halo, where all species are used in
the projection.

Takes 4 required parameters.
1 - Location of the simulation output. Can be relative to the working directory
    where this code was called from, or an absolute path.
2 - Location of the halo finder outputs from ROCKSTAR. This file should be the
    *.0.bin file. Can be relative or absolute.
3 - Rank - which halo to plot. 1 is for the most massive, 2 for the second most
    massive, etc.
4 - Whether to split the plots by DM species (pass "split") or to show all 
    species on one plot (pass "full")
"""
from utils.nbody_projection_all_species import nbody_projection_all_species
from utils.nbody_projection_split_species import nbody_projection_split_species

import sys
import os

import yt
from yt.extensions.astro_analysis.halo_analysis.api import HaloCatalog
import numpy as np

import betterplotlib as bpl
bpl.set_style()

yt.funcs.mylog.setLevel(50)  # ignore yt's output

if len(sys.argv) != 5:
    raise ValueError("4 arguments required")
# the rank argument must be an integer greater than 0
rank = int(sys.argv[3])
if rank < 1 or rank != float(sys.argv[3]):
    raise ValueError("Rank must be an integer greater than zero.")
# the last argument must either be split or full
if sys.argv[4].lower() == "split":
    split = True
elif sys.argv[4].lower() == "full":
    split = False
else:
    raise ValueError("Last parameter not recognized.")

ds_loc = os.path.abspath(sys.argv[1])
scale_factor = ds_loc[-10:-4]
ds = yt.load(ds_loc)
a = ds.scale_factor
z = 1.0 / a - 1.0

# get the location of where to write the plot
sim_dir = os.path.dirname(ds_loc) + os.sep
plots_dir = sim_dir.replace("/out/", "/plots/")

if split:
    plot_name = plots_dir + "n_body_split_halo_rank_{}_a{}.png".format(rank, scale_factor)
else:
    plot_name = plots_dir + "n_body_halo_rank_{}_a{}.png".format(rank, scale_factor)

# =========================================================================
#         
# Read in halo catalogs
# 
# =========================================================================
halo_file = os.path.abspath(sys.argv[2])
ds_halos = yt.load(halo_file)

# Then create the halo catalogs
hc = HaloCatalog(halos_ds=ds_halos, data_ds=ds, output_dir="./")
# Restrict to things about LMC mass and above
hc.add_filter('quantity_value', 'particle_mass', '>', 3E10, 'Msun')
hc.create(save_catalog=False)

# check that we have any halos at all. If not, we can exit. This can happen
# for early outputs where nothing has collapsed yet.
halo_masses = yt.YTArray([item["particle_mass"] for item in hc.catalog])
if len(halo_masses) < rank:
    fig, ax = bpl.subplots(figsize=[8.0, 8.0])
    ax.easy_add_text("Rank {} does not exist at z={:.2f}".format(rank, z), 
                     "center left")
    ax.set_axis_off()
    fig.savefig(plot_name, dpi=400)
    exit()

# We get the indices that sort it. The reversing there makes the biggest halos
# first, like we want.
rank_idxs = np.argsort(halo_masses)[::-1]

# then we pick the one of the right rank to use. Rank 1 is the first item, so 
# we subtract one in the indexing
halo = hc.catalog[rank_idxs[rank-1]]

# then get the center to use
center = [halo["particle_position_x"],
          halo["particle_position_y"],
          halo["particle_position_z"]]

# then the plot is easy to call
plot_size = ds.quan(600, "kpccm").to("kpc")

if split:
    nbody_projection_split_species(ds, center, plot_size, "kpc", plot_name)
else:
    nbody_projection_all_species(ds, center, plot_size, "kpc", plot_name)
