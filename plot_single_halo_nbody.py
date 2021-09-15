"""
plot_single_halo_full_nbody.py

Creates a projection plot of a single DM halo, where all species are used in
the projection.

Takes 4 required parameters.
1 - Location of the simulation output. Can be relative to the working directory
    where this code was called from, or an absolute path.
2 - Rank - which halo to plot. 1 is for the most massive, 2 for the second most
    massive, etc.
3 - Whether to split the plots by DM species (pass "split") or to show all
    species on one plot (pass "full")
"""
from utils.nbody_projection_all_species import nbody_projection_all_species
from utils.nbody_projection_split_species import nbody_projection_split_species
from utils import load_galaxies

import sys
from pathlib import Path

import yt

import betterplotlib as bpl

bpl.set_style()

yt.funcs.mylog.setLevel(50)  # ignore yt's output

if len(sys.argv) != 4:
    raise ValueError("3 arguments required")
# the rank argument must be an integer greater than 0
rank = int(sys.argv[2])
if rank < 1 or rank != float(sys.argv[2]):
    raise ValueError("Rank must be an integer greater than zero.")
# the last argument must either be split or full
if sys.argv[3].lower() == "split":
    split = True
elif sys.argv[3].lower() == "full":
    split = False
else:
    raise ValueError("Last parameter not recognized.")

ds_loc = Path(sys.argv[1]).resolve()
sim = load_galaxies.Simulation(ds_loc, sphere_radius_kpc=None, n_galaxies=rank)
ad = sim.ds.all_data()

# get the location of where to write the plot
plots_dir = ds_loc.parent.parent / "plots"
scale_save = round(sim.scale_factor, 4)

if split:
    plot_name = plots_dir / "n_body_split_halo_rank_{}_a{}.png".format(rank, scale_save)
else:
    plot_name = plots_dir / "n_body_halo_rank_{}_a{}.png".format(rank, scale_save)

# =========================================================================
#
# Then start plotting
#
# =========================================================================
# find the galaxy of interest
for gal in sim.galaxies:
    if rank == gal.rank:
        break
else:  # no break
    # check that we have any halos at all. If not, we can exit. This can happen
    # for early outputs where nothing has collapsed yet.
    fig, ax = bpl.subplots(figsize=[8.0, 8.0])
    ax.easy_add_text(
        "Rank {} does not exist at z={:.2f}".format(rank, z), "center left"
    )
    ax.set_axis_off()
    fig.savefig(plot_name, dpi=400)
    exit()

# then get the center to use
center = gal.center.to("kpc")

# then the plot is easy to call
plot_size = sim.ds.quan(600, "kpccm").to("kpc")

if split:
    nbody_projection_split_species(sim.ds, center, plot_size, "kpc", plot_name)
else:
    nbody_projection_all_species(sim.ds, center, plot_size, "kpc", plot_name)
