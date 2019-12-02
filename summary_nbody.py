"""
global_properties.py

Reports some global properties of the simulation, including grid structure, 
number and mass of dark matter particles, and simple halo information. Will
print this information to console and write an output file.

Takes 2 required and 1 optional parameter.
1 - Location of the simulation output. Can be relative to the working directory
    where this code was called from, or an absolute path.
2 - Location of the halo finder outputs from ROCKSTAR. This file should be the
    *.0.bin file. Can be relative or absolute.
3 - Optional argument. Must be "clobber" if included. This will bypass the 
    check whether we overwrite a previously existing output file.
4 - Optional argument. Must be "silent" if included. Will print info and 
    write it to the file if this is not included. Will only write to file
    if this is included.
"""

import sys
import os

import yt
from yt.extensions.astro_analysis.halo_analysis.api import HaloCatalog
import numpy as np

import betterplotlib as bpl
from matplotlib.patches import Circle

yt.funcs.mylog.setLevel(50)  # ignore yt's output
bpl.presentation_style()
bpl.presentation_style()  # for some reason this needs to be there twice

# Check that the third argument is correct
if len(sys.argv) == 4 and sys.argv[3] not in ["clobber", "silent"]:
    raise ValueError("Argument 3 (if used), must be 'clobber' or 'silent'")
if "silent" in sys.argv:
    silent = True
else:
    silent = False

def print_and_write(info, file_obj):
    if file_obj is not None:
        file_obj.write(info + "\n")
    if not silent:
        print(info)

ds_loc = os.path.abspath(sys.argv[1])
scale_factor = ds_loc[-10:-4]
ds = yt.load(ds_loc)
ad = ds.all_data()

# get the location of where to write the file.
sim_dir = os.path.dirname(ds_loc) + os.sep
file_dir = sim_dir.replace("/out/", "/checks/")
file_path = file_dir + "summary_nbody_a" + scale_factor + ".txt"
plots_dir = sim_dir.replace("/out/", "/plots/")

print_and_write("Output being written to:", None)
print_and_write(file_path, None)

# see if there is an existing file here that we don't want to replace.
if "clobber" not in sys.argv:
    if os.path.isfile(file_path):
        good = False
        while not good:
            choice = input("This output file already exists. "
                           "Overwrite? (y/n): ")
            if choice.lower() in ["y", "n"]:
                good = True
            else:
                print("Choose 'y' or 'n'")

        if choice.lower() == "n":
            # don't overwrite, so the code ends
            print("File won't be overwritten, so nothing will be done.")
            exit()
        # Don't need to do anything here inside this block if the user does 
        # want to overwrite, since the file created later will automatically 
        # overwrite any existing file.
# open the file
out_file = open(file_path, "w")

# =============================================================================
#         
# Basic information
# 
# =============================================================================
print_and_write("", None)  # no write since this is just for console formatting
print_and_write("Simulation location:", out_file)
print_and_write(ds_loc, out_file)

a = ds.scale_factor
z = 1.0 / a - 1.0
print_and_write("\na = {:.4f}".format(a), out_file)
print_and_write("z = {:.4f}".format(z), out_file)

# =============================================================================
#         
# Grid Structure
# 
# =============================================================================
box_size = ds.domain_width.to("Mpccm")[0]
full_grid_idxs = ad[('index', 'grid_level')].value
grid_levels, num_in_grid = np.unique(full_grid_idxs, return_counts=True)
total_cells = np.sum(num_in_grid)
cell_sizes = np.unique(ad["index", "dx"]).to("kpccm")[::-1]
# ^ np.unique returns the unique values in sorted order. We want the
# largest cells to correspond to the smallest level, so we reverse it
min_cell_size = np.min(cell_sizes)
max_cell_size = np.max(cell_sizes)
base_grid_size = np.log2(box_size / max_cell_size)

print_and_write("\nGrid Structure\n==============", out_file)
print_and_write("box size: {:.3f}".format(box_size), out_file)
print_and_write("base grid size: {:.3f}".format(base_grid_size), out_file)
print_and_write("total cells: {:,}".format(total_cells), out_file)
grid_out_str = "{:<5d}    {:>10,}    {:<8.3f}    {:<8.3f}"
grid_header_str = "{:<5}    {:>10}    {:<14}    {:<20}"
print_and_write(grid_header_str.format("Level", "Num Cells", "Cell Size",
                                        "Expected Cell Size"), out_file)
for level in grid_levels:
    level = int(level)
    num_cells = num_in_grid[level]
    cell_size = cell_sizes[level]
    expected_cell_size = box_size.to("kpccm") / (2**(level + base_grid_size))
    print_and_write(grid_out_str.format(level, num_cells, cell_size, 
                                        expected_cell_size), out_file)

# =============================================================================
#         
# Particles
# 
# =============================================================================
all_masses = ad[('N-BODY', 'MASS')].to("Msun")
masses, num_particles = np.unique(all_masses, return_counts=True)
# Do not reverse, since in ART species 0 is the lightest
# masses = masses[::-1] 
# num_particles = num_particles[::-1]
total_particles = np.sum(num_particles)

print_and_write("\nParticles\n=========", out_file)
print_and_write("total particles: {:,}".format(total_particles), out_file)
part_out_str = "{:<8d}    {:>15,}    {:<.2e}    {:<.8f}"
part_out_str_top = "{:<8d}    {:>15,}    {:<.2e}    ----------"
part_header_str = "{:<8}    {:>15}    {:<13}    {:<25}"
print_and_write(part_header_str.format("Species", "Num Particles", "Mass", 
                                        "Mass Ratio to Next"), out_file)

# print info about the species, but also add a dictionary that can convert
# from a particle's mass to its species
mass_to_species = dict()
for species in range(len(masses)):
    species = int(species)
    dm_mass = masses[species]
    dm_num = num_particles[species]

    mass_to_species["{:.0f}".format(dm_mass)] = species
    
    if species == len(masses) - 1:
        print_and_write(part_out_str_top.format(species, dm_num, dm_mass),
                        out_file)
    else:
        ratio = (masses[species + 1] / dm_mass).value
        print_and_write(part_out_str.format(species, dm_num, dm_mass, ratio),
                        out_file)

# =========================================================================
#         
# Plots
# 
# =========================================================================
normals = {"x": [1, 0, 0], 
           "y": [0, 1, 0], 
           "z": [0, 0, 1]}
for direction in normals:
    grid_plot_name = plots_dir + "grid_idxs_{}_{}.png".format(direction,
                                                              scale_factor)
    n_body_plot_name = plots_dir + "n_body_{}_{}.png".format(direction,
                                                             scale_factor)

    n_body_field = ("deposit", "N-BODY_density")
    grid_level_field = ('index', 'grid_level')

    grid_plot = yt.SlicePlot(ds, normal=normals[direction], 
                             fields=grid_level_field, width=(15, "Mpccm"))
    grid_plot.set_log(grid_level_field, False)
    grid_plot.set_cmap(grid_level_field, "Pastel1")
    grid_plot.set_zlim(grid_level_field, -0.5, 8.5)
    grid_plot.save(grid_plot_name)

    n_body_plot = yt.SlicePlot(ds, normal=normals[direction], 
                               fields=n_body_field, width=(15, "Mpccm"))
    n_body_plot.save(n_body_plot_name)

# =========================================================================
#         
# Halo analysis - printing app
# 
# =========================================================================
halo_file = os.path.abspath(sys.argv[2])
ds_halos = yt.load(halo_file)

# Then create the halo catalogs
hc = HaloCatalog(halos_ds=ds_halos, data_ds=ds, output_dir="./")
# Restrict to things about LMC mass and above
hc.add_filter('quantity_value', 'particle_mass', '>', 3E10, 'Msun')
hc.create(save_catalog=False)

# make better names for the quantities in the halo catalog
quantity_names = {"virial_radius": "Virial Radius", 
                  "particle_position_x": "Position X", 
                  "particle_position_y": "Position Y", 
                  "particle_position_z": "Position Z",
                  "particle_mass": "Virial Mass"}
# and decide what order we want to print them in
ordered_quantities = ["particle_mass", "virial_radius", "particle_position_x",
                      "particle_position_y", "particle_position_z"]

# check that we have any halos at all. If not, we can exit. This can happen
# for early outputs where nothing has collapsed yet.
halo_masses = yt.YTArray([item["particle_mass"] for item in hc.catalog])
if len(halo_masses) == 0:
    print_and_write("No halos at this redshift", out_file)
    out_file.close()
    exit()

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

# get the N-body particle locations. These are different in the old and new
# simulations, so we have to check 
species_x = dict()
species_y = dict()
species_z = dict()
if ('N-BODY_0', 'POSITION_X') in ds.field_list:
    check_contamination = True
    idx = 0
    while True:
        try:
            species_x[idx] = ad[('N-BODY_'.format(idx), 'POSITION_X')].to("Mpc").value
            species_y[idx] = ad[('N-BODY_'.format(idx), 'POSITION_Y')].to("Mpc").value
            species_z[idx] = ad[('N-BODY_'.format(idx), 'POSITION_Z')].to("Mpc").value
        except yt.utilities.exceptions.YTFieldNotFound:
            break
else:  # old sims
    check_contamination = False
    species_x[0] = ad[('N-BODY', 'POSITION_X')].to("Mpc").value
    species_y[0] = ad[('N-BODY', 'POSITION_Y')].to("Mpc").value
    species_z[0] = ad[('N-BODY', 'POSITION_Z')].to("Mpc").value

# define some helper functions that will be used for contamination calculations
def get_center(halo, with_units):
    x = halo["particle_position_x"]
    y = halo["particle_position_y"]
    z = halo["particle_position_z"]
    if not with_units:
        x = x.to("Mpc").value
        y = y.to("Mpc").value
        z = z.to("Mpc").value
    return (x, y, z)

def distance(x_0, y_0, z_0, x_1, y_1, z_1):
    return np.sqrt((x_1 - x_0)**2 + (y_1 - y_0)**2 + (z_1 - z_0)**2)

def mass_fractions(sphere):
    # Prints the fraction of mass inside a sphere that comes from different
    # species of N-body particles. This is useful for checking contamination
    masses = sphere[('N-BODY', 'MASS')]
    unique_masses, num_particles = np.unique(masses, return_counts=True)
    total_mass = masses.sum()
    for m, n_m in zip(unique_masses, num_particles):
        this_m_tot = m * n_m  # total mass in this species of particle
        frac = (this_m_tot / total_mass).value
        s = mass_to_species["{:.0f}".format(m.to("Msun"))]
        print_and_write("{}: N = {:>10,} {:>6.2f}%".format(s, n_m, frac*100), out_file)

# Print information about all the halos present
for halo in halos:
    print_and_write("\n==================================\n", out_file)
    print_and_write("Rank {} halo:".format(halo["rank"]), out_file)
    # First print the important quantities
    for quantity in ordered_quantities:   
        q_name = quantity_names[quantity]
        if q_name == "Virial Radius" or "Position" in q_name:
            value = halo[quantity].to("kpc")
            print_and_write("{}: {:<7.3f}".format(q_name, value), 
                            out_file)
        elif "Mass" in q_name:
            value = halo[quantity].to("Msun")
            print_and_write("{}: {:<2.3e}".format(q_name, value), 
                            out_file)

    # then print information about contamination, if we need to
    if check_contamination:
        print_and_write("\nClosest particle of each low-res DM species", out_file)
         # First we calculate the closest particle of each unrefined DM species
        x, y, z = get_center(halo, with_units=False)
        for idx in species_x:
            distances = distance(x, y, z, species_x[idx], species_y[idx], species_z[idx])
            print_and_write("{}: {:.0f} kpc".format(idx, np.min(distances_1)*1000), out_file)

        virial_radius = halo["virial_radius"]
        center = get_center(halo, with_units=True)
        virial_sphere = ds.sphere(center=center, radius=virial_radius)
        mpc_sphere = ds.sphere(center=center, radius=1*yt.units.Mpc)
        print_and_write("\nFraction of the DM mass from different species", out_file)
        print_and_write("Within the virial radius", out_file)
        mass_fractions(virial_sphere)
        
        print_and_write("\nWithin 1 Mpc", out_file)
        mass_fractions(mpc_sphere)
print_and_write("\n==================================\n", out_file)
        
# Then print the separation of the two biggest halos
if len(halos) >= 2:
    halo_1 = halos[0]
    halo_2 = halos[1]
    dx = halo_1["particle_position_x"] - halo_2["particle_position_x"]
    dy = halo_1["particle_position_y"] - halo_2["particle_position_y"]
    dz = halo_1["particle_position_z"] - halo_2["particle_position_z"]
    dist = np.sqrt(dx**2 + dy**2 + dz**2).to("kpc")
    print_and_write("\nSeparation of two largest halos: {:.2f}".format(dist), 
                    out_file)

# Then do a plot of the halos
halos_plot_name = plots_dir + "halos_{}.png".format(scale_factor)
def add_virial_radii(hc, axis_1, axis_2, ax):
    for halo in halos:
        coord_1 = halo["particle_position_{}".format(axis_1)].to("Mpc").value
        coord_2 = halo["particle_position_{}".format(axis_2)].to("Mpc").value
        radius = halo["virial_radius"].to("Mpc").value
        rank = halo["rank"]
        c = Circle((coord_1, coord_2), radius, 
                   fill=False, clip_on=True, ls="--")
        ax.add_artist(c)
        ax.add_text(coord_1, coord_2, rank, ha="center", va="center", 
                    fontsize=10, color="w")
        ax.scatter([coord_1], [coord_2], s=20, c=bpl.color_cycle[1])


fig, axs = bpl.subplots(ncols=3, figsize=[15, 5])
ax_xy, ax_xz, ax_yz = axs.flatten()

ax_xy.scatter(species_x[0][::10000], species_y[0][::10000], s=10)
add_virial_radii(hc, "x", "y", ax_xy)
ax_xy.add_labels("X [Mpc]", "Y [Mpc]")
ax_xy.equal_scale()

ax_xz.scatter(species_x[0][::10000], species_z[0][::10000], s=10)
add_virial_radii(hc, "x", "z", ax_xz)
ax_xz.add_labels("X [Mpc]", "Z [Mpc]")
ax_xz.equal_scale()

ax_yz.scatter(species_y[0][::10000], species_z[0][::10000], s=10)
add_virial_radii(hc, "y", "z", ax_yz)
ax_yz.add_labels("Y [Mpc]", "Z [Mpc]")
ax_yz.equal_scale()

fig.savefig(halos_plot_name)

print_and_write("\nPlots will be saved to:", out_file)
print_and_write(grid_plot_name, out_file)
print_and_write(n_body_plot_name, out_file)
print_and_write(halos_plot_name, out_file)

out_file.close()
