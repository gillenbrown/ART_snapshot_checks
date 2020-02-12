"""
plot_refined_region_nbody.py

Creates a projection plot of the refined region.

Takes 2 required parameters.
1 - Location of the simulation output. Can be relative to the working directory
    where this code was called from, or an absolute path.
2 - Location of the halo finder outputs from ROCKSTAR. This file should be the
    *.0.bin file. Can be relative or absolute.
3 - Whether to split the plots by DM species (pass "split") or to show all 
    species on one plot (pass "full")
"""
from utils.nbody_projection_all_species import nbody_projection_all_species
from utils.nbody_projection_split_species import nbody_projection_split_species

import sys
import os

import yt
from yt.extensions.astro_analysis.halo_analysis.api import HaloCatalog
import numpy as np

yt.funcs.mylog.setLevel(50)  # ignore yt's output

if len(sys.argv) != 4:
    raise ValueError("3 arguments required")
# the last argument must either be split or full
if sys.argv[3].lower() == "split":
    split = True
elif sys.argv[3].lower() == "full":
    split = False
else:
    raise ValueError("Last parameter not recognized.")

ds_loc = os.path.abspath(sys.argv[1])
scale_factor = ds_loc[-10:-4]
ds = yt.load(ds_loc)
ad = ds.all_data()
a = ds.scale_factor
z = 1.0 / a - 1.0

# get the location of where to write the plot
sim_dir = os.path.dirname(ds_loc) + os.sep
plots_dir = sim_dir.replace("/out/", "/plots/")

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

# We get the indices that sort it. The reversing there makes the biggest halos
# first, like we want.
rank_idxs = np.argsort(halo_masses)[::-1]

# Add this info to each halo object, and put the halos into a new sorted list,
# with the highest mass (lowest rank) halos first)
halos = []
for rank, idx in enumerate(rank_idxs, start=1):
    halo = hc.catalog[idx]
    halo["rank"] = rank
    halos.append(halo)

# Then start plotting
is_zoom = ('N-BODY_0', 'POSITION_X') in ds.field_list
if is_zoom:
    center = [np.median(ad[('N-BODY_0', 'POSITION_X')]),
              np.median(ad[('N-BODY_0', 'POSITION_Y')]),
              np.median(ad[('N-BODY_0', 'POSITION_Z')])]
else:
    center = ds.domain_center

# determine how big to make the plot window
box_length = ds.domain_width[0]
max_length = ds.quan(10.0, "Mpccm")
plot_size = min(max_length, box_length).to("Mpc")

if split:
    plot_name = plots_dir + "n_body_split_refined_a{}.png".format(scale_factor)
    nbody_projection_split_species(ds, center, plot_size, "Mpc", plot_name)
    # Here I don't do annotation of the halos, since it's too hard to do all
    # the pixel conversions and whatnot
else:
    plot_name = plots_dir + "n_body_refined_a{}.png".format(scale_factor)
    plot = nbody_projection_all_species(ds, center, plot_size, "Mpc", None)

    # annotate the halos
    text_args={"color":"k", "va":"center", "ha":"center"}
    for halo in halos:
        coord = [halo["particle_position_x"], 
                 halo["particle_position_y"], 
                 halo["particle_position_z"]]
        rank = halo["rank"]
        if rank < 6:
            plot.annotate_text(text=rank, pos=coord, text_args=text_args)

    plot.save(plot_name, mpl_kwargs={"dpi": 400})