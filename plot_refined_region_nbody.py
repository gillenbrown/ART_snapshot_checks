"""
plot_refined_region_nbody.py

Creates a projection plot of the refined region.

Takes 2 required parameters.
1 - Location of the simulation output. Can be relative to the working directory
    where this code was called from, or an absolute path.
2 - Whether to split the plots by DM species (pass "split") or to show all
    species on one plot (pass "full")
"""
from utils.nbody_projection_all_species import nbody_projection_all_species
from utils.nbody_projection_split_species import nbody_projection_split_species
from utils import load_galaxies

import sys
from pathlib import Path

import yt
import numpy as np

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
sim = load_galaxies.Simulation(ds_loc, sphere_radius_kpc=None, n_galaxies=10)
ad = sim.ds.all_data()

# get the location of where to write the plot
plots_dir = ds_loc.parent.parent / "plots"
scale_save = ds_loc.stem.split("_")[-1].replace("a", "")

# =========================================================================
#
# Then start plotting
#
# =========================================================================
is_zoom = ("N-BODY_0", "POSITION_X") in sim.ds.field_list
if is_zoom:
    center = [
        np.median(ad[("N-BODY_0", "POSITION_X")]),
        np.median(ad[("N-BODY_0", "POSITION_Y")]),
        np.median(ad[("N-BODY_0", "POSITION_Z")]),
    ]
else:
    center = sim.ds.domain_center

# determine how big to make the plot window
box_length = sim.ds.domain_width[0]
max_length = sim.ds.quan(10.0, "Mpccm")
plot_size = min(max_length, box_length).to("Mpc")

if split:
    plot_name = plots_dir / "n_body_split_refined_a{}.png".format(scale_save)
    nbody_projection_split_species(sim.ds, center, plot_size, "Mpc", plot_name)
    # Here I don't do annotation of the halos, since it's too hard to do all
    # the pixel conversions and whatnot
else:
    plot_name = plots_dir / "n_body_refined_a{}.png".format(scale_save)
    plot = nbody_projection_all_species(sim.ds, center, plot_size, "Mpc", None)

    # annotate the halos
    text_args = {"color": "k", "va": "center", "ha": "center"}
    for halo in sim.galaxies:
        center = halo.center.to("Mpc").value
        rank = halo.rank
        plot.annotate_text(text=rank, pos=center, text_args=text_args)

    plot.save(plot_name, mpl_kwargs={"dpi": 400})
