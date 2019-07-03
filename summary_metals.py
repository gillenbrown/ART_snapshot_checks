"""
summmary_metals.py

Reports some properties related to the metals in the simulation

Takes 1 required and 2 optional parameters.
1 - Location of the simulation output. Can be relative to the working directory
    where this code was called from, or an absolute path.
2 - Optional argument. Must be "clobber" if included. This will bypass the 
    check whether we overwrite a previously existing output file.
3 - Optional argument. Must be "silent" if included. Will print info and 
    write it to the file if this is not included. Will only write to file
    if this is included.
"""

import sys
import os

import yt
import numpy as np

import betterplotlib as bpl

yt.funcs.mylog.setLevel(50)  # ignore yt's output
bpl.presentation_style()
bpl.presentation_style()  # for some reason this needs to be there twice

# Check that the third argument is correct
if len(sys.argv) == 3 and sys.argv[2] not in ["clobber", "silent"]:
    raise ValueError("Argument 2 (if used), must be 'clobber' or 'silent'")
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
file_path = file_dir + "summary_metals_a" + scale_factor + ".txt"
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
# Gas cells
# 
# =============================================================================
# See some statistics about the metal properties of the gas at various levels

if ('artio', 'HVAR_METAL_DENSITY_Mg') in ds.field_list:
    elements = ["II", "Ia", "AGB", "C", "N", "O", "Mg", "S", "Ca", "Fe"]
elif ('artio', 'HVAR_METAL_DENSITY_C') in ds.field_list:
    elements = ["II", "Ia", "AGB", "C", "N", "O", "Fe"]
else:
    elements = ["II", "Ia"]

full_grid_levels = ad[('index', 'grid_level')].value
grid_levels = np.unique(full_grid_levels)
cell_sizes = np.unique(ad["index", "dx"]).to("pc")[::-1]
# ^ np.unique returns the unique values in sorted order. We want the
# largest cells to correspond to the smallest level, so we reverse it

# Get all the cells at various levels
level_idxs = {level:np.where(full_grid_levels == level)[0] 
              for level in grid_levels}

# Get the necessary data beforehand to reduce data accessing needs
densities = ad[('gas', 'density')]
volumes = ad[('gas', 'cell_volume')]
metal_densities = {elt:ad[('artio', 'HVAR_METAL_DENSITY_{}'.format(elt))]
                   for elt in elements}

# yt doesn't know about the units for my new fields
code_density = ds.mass_unit / (ds.length_unit)**3
for elt in metal_densities:
    if elt not in ["II", "Ia"]:
        metal_densities[elt] *= code_density
    metal_densities[elt] = metal_densities[elt].to("g/cm**3")

# Then go level by level to get the properties of the gas at that level
header_str = "{:<10s}\t{:>10s}\t{:>10s}\t{:>10s}\t{:>10s}"
row_str = "{:<10s}\t{:>10.3E}\t{:>10.3E}\t{:>10.3E}\t{:>10.3E}"
for level, cell_size in zip(level_idxs, cell_sizes):
    print_and_write("level={:.0f}, cell size={:.2f}".format(level, cell_size), 
                    out_file)
    print_and_write(header_str.format("Element", "Minimum Z", "Median Z", 
                                      "Mean Z", "Maximum Z"), out_file)
    
    # get the info at this level
    idxs = level_idxs[level]  # indices of cells at this level
    density_level = densities[idxs]  # densities of cell on this level
    volume_level = volumes[idxs]  # volumes of cell on this level
    total_mass_level = np.sum(density_level * volume_level) 
    # ^ The total gas mass in cells on this level
    
    # then go element by element
    for element in elements:
        densities_elt_level = metal_densities[element][idxs]
        z_elt = (densities_elt_level / density_level).value
        total_elt_mass_level = np.sum(densities_elt_level * volume_level)
        true_mean = (total_elt_mass_level / total_mass_level).value

        print_and_write(row_str.format(element, np.min(z_elt), 
                                       np.median(z_elt), true_mean, 
                                       np.max(z_elt)), out_file)
    print_and_write("", out_file)  # for spacing
# =============================================================================
#         
# Stars
# 
# =============================================================================
# We do the equivalent thing here, although it's different because stars have
# metallicity and mass, not densities
# We also have to modify the elements
star_elements = [elt.replace("II", "SNII").replace("Ia", "SNIa")
                 for elt in elements]

stellar_masses = ad[("STAR", "MASS")]
total_mass = np.sum(stellar_masses)

if len(stellar_masses) == 0:
    print_and_write("No stars at this redshift", out_file)
else:
    print_and_write("Stars", out_file)
    print_and_write(header_str.format("Element", "Minimum Z", "Median Z", 
                                    "Mean Z", "Maximum Z"), out_file)
    for element in star_elements:
        z_elt = ad[("STAR", "METALLICITY_{}".format(element))].value
        total_elt_mass = np.sum(z_elt * stellar_masses)
        true_mean = (total_elt_mass / total_mass).value

        print_and_write(row_str.format(element, np.min(z_elt), 
                                       np.median(z_elt), true_mean, 
                                       np.max(z_elt)), out_file)
print_and_write("", out_file)  # for spacing

# =========================================================================
#         
# Cell Masses
# 
# =========================================================================
# Checking the mass of cells can help debug refinement, so make sure that 
# the Lagrangian refinement is actually working correctly
refine_top_header = "Percentiles of cell gas mass distribution, units of {}"
refine_header_str = "{:<10s}\t{:>10d}\t{:>10d}\t{:>10d}\t{:>10d}\t{:>10d}"
refine_row_str = "{:<10.0f}" + 5 * "\t{:>10.3E}"

percentiles = [0, 25, 50, 75, 100]

for unit in ["code_mass", "Msun"]:
    # write header info
    print_and_write(refine_top_header.format(unit), out_file)
    print_and_write(refine_header_str.format("Level", *percentiles), out_file)

    gas_mass = ad[('gas', 'cell_mass')].to(unit).value
    for level in level_idxs:
        idxs = level_idxs[level]  # indices of cells at this level
        mass_percentiles = np.percentile(gas_mass[idxs], percentiles)
        print_and_write(refine_row_str.format(level, *mass_percentiles), 
                        out_file)
    print_and_write("", out_file)  # for spacing

# =========================================================================
#         
# Plots
# 
# =========================================================================
normals = {"x": [1, 0, 0], 
           "y": [0, 1, 0], 
           "z": [0, 0, 1]}
for direction in normals:
    gas_plot_name = plots_dir + "gas_density_{}_{}.png".format(direction,
                                                               scale_factor)

    gas_plot = yt.SlicePlot(ds, normal=normals[direction], 
                            fields=("gas", "density"), width=(15, "Mpccm"))
    gas_plot.save(gas_plot_name)

print_and_write("\nPlots will be saved to:", out_file)
print_and_write(gas_plot_name, out_file)

out_file.close()
