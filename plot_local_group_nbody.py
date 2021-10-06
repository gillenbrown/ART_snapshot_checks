"""
plot_local_group_nbody.py

Creates a projection plot containing the top two halos.

Takes 4 required parameters.
1 - Location of the simulation output. Can be relative to the working directory
    where this code was called from, or an absolute path.
2- Whether to split the plots by DM species (pass "split") or to show all
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

if len(sys.argv) != 3:
    raise ValueError("2 arguments required")
# the last argument must either be split or full
if sys.argv[2].lower() == "split":
    split = True
elif sys.argv[2].lower() == "full":
    split = False
else:
    raise ValueError("Last parameter not recognized.")

ds_loc = Path(sys.argv[1]).resolve()
sim = load_galaxies.Simulation(ds_loc, sphere_radius_kpc=None, n_galaxies=2)
scale_factor = ds_loc.stem.split("_")[-1].replace("a", "")

# get the location of where to write the plot
plots_dir = ds_loc.parent.parent / "plots"

if split:
    plot_name = plots_dir / "n_body_split_local_group_a{}.png".format(scale_factor)
else:
    plot_name = plots_dir / "n_body_local_group_a{}.png".format(scale_factor)

# =========================================================================
#
# get the galaxies of interest
#
# =========================================================================
if len(sim.galaxies) < 2:
    fig, ax = bpl.subplots(figsize=[8.0, 8.0])
    ax.easy_add_text(
        "The Local Group does not exist at z={:.2f}".format(z), "center left"
    )
    ax.set_axis_off()
    fig.savefig(plot_name, dpi=400)
    exit()

# then we pick the one of the right rank to use
for gal in sim.galaxies:
    if gal.rank == 1:
        m31 = gal
    elif gal.rank == 2:
        mw = gal

# define the center to be the middle between the two halos. I don't use a numpy
# function for this so I can keep units.
center = [
    0.5 * (m31.center[0] + mw.center[0]),
    0.5 * (m31.center[1] + mw.center[1]),
    0.5 * (m31.center[2] + mw.center[2]),
]

# then the plot is easy to call
plot_size = sim.ds.quan(2, "Mpccm").to("Mpc")

if split:
    nbody_projection_split_species(sim.ds, center, plot_size, "Mpc", plot_name)
else:
    plot = nbody_projection_all_species(sim.ds, center, plot_size, "Mpc", None)

    # annotate the halos
    text_args = {"color": "k", "va": "center", "ha": "center"}
    for gal in [mw, m31]:
        plot.annotate_text(
            text=gal.rank, pos=gal.center.to("Mpc").value, text_args=text_args
        )

    plot.save(plot_name, mpl_kwargs={"dpi": 400})
