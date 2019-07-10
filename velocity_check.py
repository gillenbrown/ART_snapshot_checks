"""
velocity_check.py

Reports the velocities of gas cells, N-boody particles, and star particles, to 
help with debugging

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

import numpy as np
import yt

import utils as u 

yt.funcs.mylog.setLevel(50)  # ignore yt's output

# Check that the third argument is correct
if len(sys.argv) == 3 and sys.argv[2] not in ["clobber", "silent"]:
    raise ValueError("Argument 2 (if used), must be 'clobber' or 'silent'")
if "silent" in sys.argv:
    silent = True
else:
    silent = False

# get the dataset info
ds_loc = os.path.abspath(sys.argv[1])
scale_factor = ds_loc[-10:-4]
ds = yt.load(ds_loc)
ad = ds.all_data()

# get the location of where to write the file.
sim_dir = os.path.dirname(ds_loc) + os.sep
file_dir = sim_dir.replace("/out/", "/checks/")
file_path = file_dir + "velocities_a" + scale_factor + ".txt"
plots_dir = sim_dir.replace("/out/", "/plots/")

# inform the user of where the output is going
u.print_and_write("Output being written to:", None, silent)
u.print_and_write(file_path, None, silent)

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

# define a shorthand function that handles the parameters to print_and_write,
# so we don't have to duplicate that each time
def out(info):
    return u.print_and_write(info, out_file, silent)

ds = yt.load(sys.argv[1])
ad = ds.all_data()
levels_gas = ad[('index', 'grid_level')].value
grid_levels, num_in_grid = np.unique(levels_gas, return_counts=True)
cell_sizes = np.unique(ad["index", "dx"]).to("pc").value[::-1]
# ^ np.unique returns the unique values in sorted order. We want the
# largest cells to correspond to the smallest level, so we reverse it

# get the gas velocity
velocity_gas = ad[('gas', 'velocity_magnitude')].to("km/s").value

# do this for stars and DM too
velocity_star = ad[('STAR', 'particle_velocity_magnitude')].to("km/s").value
velocity_dm   = ad[('N-BODY', 'particle_velocity_magnitude')].to("km/s").value

# We'll need to restrict this to a few particles. We wnat to find the cell a 
# given particle is in, which is computationally expensive, so we'll restrict
# to only the fastest velocities
n_vels = 100  # how many cells/DM particles/star particles to save
# Sort the velocities, then get the biggest ones (that will be at the end).
idx_fast_dm   = np.argsort(velocity_dm)[-n_vels::]
idx_fast_star = np.argsort(velocity_star)[-n_vels::]
# then get only those velocities
velocity_dm = velocity_dm[idx_fast_dm]
velocity_star = velocity_star[idx_fast_star]

# We want to figure out what cells the particles are in, which we can do if we
# get their locations. Again only get the fastest ones
position_star = ad[('STAR', 'particle_position')][idx_fast_star]
position_dm   = ad[('N-BODY', 'particle_position')][idx_fast_dm]

levels_star = ds.find_field_values_at_points([('index', 'grid_level')], 
                                             position_star).value
levels_dm   = ds.find_field_values_at_points([('index', 'grid_level')], 
                                             position_dm).value

# print the max velocities in each level
# We have to have this ugly code to handle what happens when there are no stars
# or DM on a given level. This is all for the string that gets printed
header_str = "{:<10s}" + 5 * "{:>10s}"
level_str = "{:<10.0f}"
not_empty = "{:>10.2f}"
time = "{:>10.2E}"
empty = "{:>10s}".format("---")
row_str = level_str + 4 * not_empty + time
row_str_no_star = level_str + 2 * not_empty + empty + not_empty + time
row_str_no_dm = level_str + empty + 3 * not_empty + time
row_str_no_both = level_str + empty + not_empty + empty + not_empty + time

out("\nThis shows the highest velocity present in the following components at " 
    "each level.\nAll velocities are in km/s, cell size in pc, dt in years.")
out(header_str.format("Level", "DM", "Gas", "Stars", "Cell Size", "dt"))

for level, cell_size in zip(grid_levels, cell_sizes):
    idx_gas = np.where(levels_gas == level)
    idx_star = np.where(levels_star == level)
    idx_dm = np.where(levels_dm == level)

    # get the maximum velocity at this level, if we have particles here
    vel_max_gas = np.max(velocity_gas[idx_gas])
    if len(idx_star[0]) > 0:  # we do have stars at this level
        vel_max_star = np.max(velocity_star[idx_star])
    if len(idx_dm[0]) > 0:  # DM at this level
        vel_max_dm = np.max(velocity_dm[idx_dm])

    # then decide what to print
    if len(idx_star[0]) > 0:  # stars 
        if len(idx_dm[0]) > 0:  # stars and DM 
            vel_max_all = max([vel_max_gas, vel_max_dm, vel_max_star])
            dt = cell_size * yt.units.pc / (vel_max_all * yt.units.km / yt.units.s)
            dt = dt.to("yr").value

            out_str = row_str.format(level, vel_max_dm, vel_max_gas, vel_max_star, cell_size, dt)
        else:  # stars, no DM
            vel_max_all = max([vel_max_gas, vel_max_star])
            dt = cell_size * yt.units.pc / (vel_max_all * yt.units.km / yt.units.s)
            dt = dt.to("yr").value
            out_str = row_str_no_dm.format(level, vel_max_gas, vel_max_star, cell_size, dt)
    else:  # no stars
        if len(idx_dm[0]) > 0:  # no stars, but DM 
            vel_max_all = max([vel_max_gas, vel_max_dm])
            dt = cell_size * yt.units.pc / (vel_max_all * yt.units.km / yt.units.s)
            dt = dt.to("yr").value
            out_str = row_str_no_star.format(level, vel_max_dm, vel_max_gas, cell_size, dt)
        else:  # no stars, no dm
            dt = cell_size * yt.units.pc / (vel_max_gas * yt.units.km / yt.units.s)
            dt = dt.to("yr").value
            out_str = row_str_no_both.format(level, vel_max_gas, cell_size, dt)
        
    out(out_str)

# print the information of the highest-velocity cells/particles
out("\nHere are the cells/particles with the highest velocities.")
n_each = 10  # how many cells/DM particles/star particles to print
# Sort the velocities, then get the biggest ones (that will be at the end).
# The indexing is weird. We want to go backwards, (biggest first), so the step
# must be -1. Since we're going backwards, we end once we count back n_each, so 
# that goes in the end part of the slice. We have to subtract one, since the 
# final value isn't included in the slice
idx_sort_gas  = np.argsort(velocity_gas)[:-n_each-1:-1]
idx_sort_dm   = np.argsort(velocity_dm)[:-n_each-1:-1]
idx_sort_star = np.argsort(velocity_star)[:-n_each-1:-1]

header_str = "{:>20s}\t" + 2 * "{:<10s}" +3 * "{:>20s}"
row_str = "{:>20.10f}\t" + "{:<10.0f}" + "{:<10.3E}" + 3 * "{:>20.10f}"

for name in ["DM", "Gas", "Stars"]:
    # get the info for each of the components
    if name == "Gas":
        idxs = idx_sort_gas
        velocities = velocity_gas
        levels = levels_gas
        masses = ad[('gas', 'cell_mass')].to("Msun").value
        # haven't gotten the positions for the gas yet, so do that now. This 
        # is a bit ugly, but is useful to make the later code easier
        pos_gas_x = ad[('gas', 'x')].to("code_length").value
        pos_gas_y = ad[('gas', 'y')].to("code_length").value
        pos_gas_z = ad[('gas', 'z')].to("code_length").value
        positions = np.stack([pos_gas_x, pos_gas_y, pos_gas_z], axis=1)
        positions = ds.arr(positions, "code_length")
    elif name == "DM":
        idxs = idx_sort_dm
        velocities = velocity_dm
        levels = levels_dm
        masses = ad[('N-BODY', 'particle_mass')][idx_fast_dm].to("Msun").value
        positions = position_dm
    elif name == "Stars":
        idxs = idx_sort_star
        velocities = velocity_star
        levels = levels_star
        masses = ad[('STAR', 'MASS')][idx_fast_star].to("Msun").value
        positions = position_star

    out("\n{} Highest Velocities".format(name))
    out(header_str.format("Velocity [km/s]", "Level", "Mass", "Location X [code]", 
                          "Location Y [code]", "Location Z [code]"))
    # go through each of the highest velocity cells and print their information
    for idx in idxs:
        x, y, z = positions[idx].to("code_length").value
        out(row_str.format(velocities[idx], levels[idx], masses[idx], x, y, z))


