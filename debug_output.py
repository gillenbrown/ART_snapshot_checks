"""
debug_output.py

Reports some global properties of the simulation that are useful for looking
at the nitty gritty of the simulations. 

Takes 2 required and 1 optional parameter.
1 - Location of the simulation output. Can be relative to the working directory
    where this code was called from, or an absolute path.
2 - Optional argument. Must be "clobber" if included. This will bypass the 
    check whether we overwrite a previously existing output file.
3 - Optional argument. Must be "silent" if included. Will print info and 
    write it to the file if this is not included. Will only write to file
    if this is included.
"""

import sys, os, gc

import yt
import numpy as np

yt.funcs.mylog.setLevel(50)  # ignore yt's output

# Check that the third argument is correct
if len(sys.argv) == 4 and sys.argv[3] not in ["clobber", "silent"]:
    raise ValueError("Argument 3 (if used), must be 'clobber' or 'silent'")
if "silent" in sys.argv:
    silent = True
else:
    silent = False


def print_and_write(info, file_obj):
    if file_obj is not None:
        file_obj.write(str(info) + "\n")
    if not silent:
        print(info)


ds_loc = os.path.abspath(sys.argv[1])
scale_factor = ds_loc[-10:-4]
ds = yt.load(ds_loc)
ad = ds.all_data()

# get the location of where to write the file.
sim_dir = os.path.dirname(ds_loc) + os.sep
file_dir = sim_dir.replace("/out/", "/checks/")
file_path = file_dir + "debug_a" + scale_factor + ".txt"
plots_dir = sim_dir.replace("/out/", "/plots/")

print_and_write("Output being written to:", None)
print_and_write(file_path, None)

# see if there is an existing file here that we don't want to replace.
if "clobber" not in sys.argv:
    if os.path.isfile(file_path):
        good = False
        while not good:
            choice = input("This output file already exists. " "Overwrite? (y/n): ")
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

# define a shorthand function that handles the parameters to out,
# so we don't have to duplicate that each time
def out(info):
    return print_and_write(info, out_file)


# =============================================================================
#
# Basic information
#
# =============================================================================
print_and_write("", None)  # no write since this is just for console formatting
out("Simulation location:")
out(ds_loc)

a = ds.scale_factor
z = 1.0 / a - 1.0
out("\na = {:.4f}".format(a))
out("z = {:.4f}".format(z))

# =============================================================================
#
# Grid Structure
#
# =============================================================================
box_size = ds.domain_width.to("Mpccm")[0]
levels_gas = ad[("index", "grid_level")].value
grid_levels, num_in_grid = np.unique(levels_gas, return_counts=True)
total_cells = np.sum(num_in_grid)
cell_sizes = np.unique(ad["index", "dx"]).to("pc")[::-1]
# ^ np.unique returns the unique values in sorted order. We want the
# largest cells to correspond to the smallest level, so we reverse it
min_cell_size = np.min(cell_sizes)
max_cell_size = np.max(cell_sizes)
base_grid_size = np.log2(box_size / max_cell_size)

out("\nGrid Structure\n==============")
out("box size: {:.3f}".format(box_size))
out("base grid size: {:.3f}".format(base_grid_size))
out("total cells: {:,}".format(total_cells))
grid_out_str = "{:<5d}    {:>13,}    {:>16.3f}    {:>16}    {:>16,}"
grid_header_str = "{:<5}    {:>13}    {:>16}    {:>16}    {:>16}"
out(
    grid_header_str.format(
        "Level", "Num Cells", "Cell Size [pc]", "Time Refinement", "Num Timesteps"
    )
)
for level in grid_levels:
    level = int(level)
    num_cells = num_in_grid[level]
    cell_size = cell_sizes[level]
    refinement_factor = 1
    for factor in ds.artio_parameters["time_refinement_factor"][: level + 1]:
        refinement_factor *= factor
    num_timesteps = num_cells * refinement_factor
    out(
        grid_out_str.format(
            level, num_cells, cell_size.to("pc").value, refinement_factor, num_timesteps
        )
    )

# =============================================================================
#
# Particles
#
# =============================================================================
masses, num_particles = np.unique(ad[("N-BODY", "MASS")].to("Msun"), return_counts=True)
# Do not reverse, since in ART species 0 is the lightest
total_particles = np.sum(num_particles)

out("\nParticles\n=========")
out("total particles: {:,}".format(total_particles))
part_out_str = "{:<8d}    {:>15,}    {:<12e}    {:<12e}    {:<.8f}"
part_out_str_top = "{:<8d}    {:>15,}    {:<12e}    {:<12e}    ----------"
part_header_str = "{:<8s}    {:>15s}    {:<12s}    {:<12s}    {:<25s}"
out(
    part_header_str.format(
        "Species", "Num Particles", "Mass [Code]", "Mass [Msun]", "Mass Ratio to Next"
    )
)

# print info about the species, but also add a dictionary that can convert
# from a particle's mass to its species
mass_to_species = dict()
for species in range(len(masses)):
    species = int(species)
    dm_mass = masses[species]
    m_sol = dm_mass.to("Msun").value
    m_code = dm_mass.to("code_mass").value
    dm_num = num_particles[species]

    mass_to_species["{:.0f}".format(dm_mass)] = species

    if species == len(masses) - 1:
        out(part_out_str_top.format(species, dm_num, m_code, m_sol))
    else:
        ratio = (masses[species + 1] / dm_mass).value
        out(part_out_str.format(species, dm_num, m_code, m_sol, ratio))


# =============================================================================
#
# Gas cells
#
# =============================================================================
# See some statistics about the metal properties of the gas at various levels
# first check if we actually have hydro cells. If not we exit
if ("gas", "density") not in ds.derived_field_list:
    out_file.close()
    exit()

out("\nElement Abundances\n==================")

if ("artio", "HVAR_METAL_DENSITY_Mg") in ds.field_list:
    elements = ["II", "Ia", "AGB", "C", "N", "O", "Mg", "S", "Ca", "Fe"]
elif ("artio", "HVAR_METAL_DENSITY_Fe") in ds.field_list:
    elements = ["II", "Ia", "AGB", "C", "N", "O", "Fe"]
elif ("artio", "HVAR_METAL_DENSITY_FE") in ds.field_list:
    elements = ["II", "Ia", "AGB", "C", "N", "O", "FE"]
else:
    elements = ["II", "Ia"]

# Get all the cells at various levels
level_idxs = {level: np.where(levels_gas == level)[0] for level in grid_levels}

# Get the necessary data beforehand to reduce data accessing needs
densities = ad[("gas", "density")]
volumes = ad[("gas", "cell_volume")]
metal_densities = {
    elt: ad[("artio", "HVAR_METAL_DENSITY_{}".format(elt))] for elt in elements
}

# yt doesn't know about the units for my new fields
code_density = ds.mass_unit / (ds.length_unit) ** 3
for elt in metal_densities:
    if elt not in ["II", "Ia"]:
        metal_densities[elt] *= code_density
    metal_densities[elt] = metal_densities[elt].to("g/cm**3")

# Then go level by level to get the properties of the gas at that level
header_str = "{:<10s}{:>12s}{:>12s}{:>12s}{:>12s}"
row_str = "{:<10s}{:>12.3E}{:>12.3E}{:>12.3E}{:>12.3E}"
for level, cell_size, n_cells in zip(level_idxs, cell_sizes, num_in_grid):
    out(
        "level={:.0f}, cell size={:.2f}, "
        "number of cells={:,.0f}".format(level, cell_size, n_cells)
    )
    out(header_str.format("Element", "Minimum Z", "Median Z", "Mean Z", "Maximum Z"))

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

        out(
            row_str.format(
                element, np.min(z_elt), np.median(z_elt), true_mean, np.max(z_elt)
            )
        )
    out("")  # for spacing
# delete unneeded variables
del densities, volumes, metal_densities
gc.collect()

# =============================================================================
#
# Stars
#
# =============================================================================
# We do the equivalent thing here, although it's different because stars have
# metallicity and mass, not densities
# We also have to modify the elements
star_elements = [elt.replace("II", "SNII").replace("Ia", "SNIa") for elt in elements]

stellar_masses = ad[("STAR", "MASS")]
total_mass = np.sum(stellar_masses)

if len(stellar_masses) == 0:
    out("No stars at this redshift")
else:
    out("Stars")
    out(header_str.format("Element", "Minimum Z", "Median Z", "Mean Z", "Maximum Z"))
    for element in star_elements:
        z_elt = ad[("STAR", "METALLICITY_{}".format(element))].value
        total_elt_mass = np.sum(z_elt * stellar_masses)
        true_mean = (total_elt_mass / total_mass).value

        out(
            row_str.format(
                element, np.min(z_elt), np.median(z_elt), true_mean, np.max(z_elt)
            )
        )
        del z_elt
out("")  # for spacing
# delete unnecessary variables
del star_elements, stellar_masses
gc.collect()

# =========================================================================
#
# Cell Masses
#
# =========================================================================
out("\nCell Masses\n===========")
# Checking the mass of cells can help debug refinement, so make sure that
# the Lagrangian refinement is actually working correctly
refine_top_header = "Percentiles of cell gas mass distribution, units of {}"
percentiles = [0, 0.1, 1, 25, 50, 75, 99, 99.9, 100]
refine_header_str = "{:<10s}" + "{:>13s}" + len(percentiles) * "{:>12.1f}"
refine_row_str = "{:<10.0f}" + "{:>13,.0f}" + len(percentiles) * "{:>12.3E}"

for unit in ["code_mass", "Msun"]:
    # write header info
    out(refine_top_header.format(unit))
    out(refine_header_str.format("Level", "Number", *percentiles))

    gas_mass = ad[("gas", "cell_mass")].to(unit).value
    for level in level_idxs:
        idxs = level_idxs[level]  # indices of cells at this level
        n_cell = len(level_idxs[level])
        mass_percentiles = np.percentile(gas_mass[idxs], percentiles)
        out(refine_row_str.format(level, n_cell, *mass_percentiles))
    del gas_mass
    out("")  # for spacing
gc.collect()

# =========================================================================
#
# Checking feedback caps
#
# =========================================================================
# Check the number of cells containing young stars with temperatures above
# certain limits, to check how the feedback caps are behaving.
out("\nTemperatures around young stars\n===============================")
out("This shows the fraction of stars with temperatures above the given limit")
star_ages = ad[("STAR", "age")]
star_positions = ad[("STAR", "particle_position")]
for max_age in [15, 40]:
    out(f"\nMax age = {max_age} Myr")
    # Get the young stars (defined as less than 15 Myr). This value is chosen to
    young_idxs = np.where(star_ages < 40 * yt.units.Myr)[0]
    if len(young_idxs) == 0:
        out("No stars!")
    else:
        out("Temperature    Fraction")
        # cap the number at 1000
        if len(young_idxs) > 1000:
            np.random.shuffle(young_idxs)  # get a random 1000
            young_idxs = young_idxs[:1000]
        # then get the temperature at the positions of all these young stars
        young_positions = star_positions[young_idxs]
        young_temps = ds.find_field_values_at_points(
            [("gas", "temperature")], young_positions
        )
        # then see how many stars have temperatures above certain values
        def find_fraction_above(temps, temp_max):
            num_above = np.sum(temps > temp_max)
            return 100 * num_above / len(temps)

        temp_caps = [1e6, 1e7, 2e7, 4e7, 4.3e7, 5e7, 1e8, 1e9]
        for cap in temp_caps:
            out(f"{cap:11.2e}    {find_fraction_above(young_temps, cap):.2f}")
        del young_positions, young_temps, young_idxs
del star_ages, star_positions
gc.collect()

# =========================================================================
#
# Velocities
#
# =========================================================================
out("\nVelocities\n==========")
# get the gas velocity
vx_gas = np.abs(ad[("gas", "velocity_x")].to("km/s").value)
vy_gas = np.abs(ad[("gas", "velocity_y")].to("km/s").value)
vz_gas = np.abs(ad[("gas", "velocity_z")].to("km/s").value)
# the velocity ART uses for the bulk gas velocith is the maximum of the x, y,
# and z components. The next function call does this for each cell.
velocity_bulk_gas = np.amax([vx_gas, vy_gas, vz_gas], axis=0)
# and sound speed
try:
    gamma = ad[("gas", "gamma")]
except:  # field not found in the entropy based scheme. May want to double check this.
    gamma = 5 / 3
sound_speed_gas = (
    np.sqrt(
        gamma
        * (gamma - 1.0)
        * ad[("artio", "HVAR_INTERNAL_ENERGY")]
        / ad[("gas", "density")]
    )
    .to("km/s")
    .value
)
v_tot_gas = velocity_bulk_gas + sound_speed_gas

# do this for stars and DM too
vx_star = np.abs(ad[("STAR", "particle_velocity_x")].to("km/s").value)
vy_star = np.abs(ad[("STAR", "particle_velocity_y")].to("km/s").value)
vz_star = np.abs(ad[("STAR", "particle_velocity_z")].to("km/s").value)
velocity_star = np.amax([vx_star, vy_star, vz_star], axis=0)

vx_dm = np.abs(ad[("N-BODY", "particle_velocity_x")].to("km/s").value)
vy_dm = np.abs(ad[("N-BODY", "particle_velocity_y")].to("km/s").value)
vz_dm = np.abs(ad[("N-BODY", "particle_velocity_z")].to("km/s").value)
velocity_dm = np.amax([vx_dm, vy_dm, vz_dm], axis=0)

# delete full arrays that aren't needed
del vx_gas, vy_gas, vz_gas
del vx_star, vy_star, vz_star
del vx_dm, vy_dm, vz_dm
gc.collect()

# We'll need to restrict this to a few particles. We wnat to find the cell a
# given particle is in, which is computationally expensive, so we'll restrict
# to only the fastest velocities
n_vels = 100  # how many cells/DM particles/star particles to save
# Sort the velocities, then get the biggest ones (that will be at the end).
idx_fast_dm = np.argsort(velocity_dm)[-n_vels::]
idx_fast_star = np.argsort(velocity_star)[-n_vels::]
# then get only those velocities
velocity_dm = velocity_dm[idx_fast_dm]
velocity_star = velocity_star[idx_fast_star]

# We want to figure out what cells the particles are in, which we can do if we
# get their locations. Again only get the fastest ones
position_star = ad[("STAR", "particle_position")][idx_fast_star]
position_dm = ad[("N-BODY", "particle_position")][idx_fast_dm]

levels_star = ds.find_field_values_at_points(
    [("index", "grid_level")], position_star
).value
levels_dm = ds.find_field_values_at_points([("index", "grid_level")], position_dm).value

# print the max velocities in each level
# We have to have this ugly code to handle what happens when there are no stars
# or DM on a given level. This is all for the string that gets printed
header_str = "{:<12s} {:>12s}" + 8 * "{:>12s}"
level_str = "{:<12.0f} {:>12,.0f}"
not_empty = "{:>12.2f}"
empty = "{:>12s}".format("---")
time = "{:>12.2E}"
row_str = level_str + 7 * not_empty + time
row_str_no_star = level_str + 5 * not_empty + empty + not_empty + time
row_str_no_dm = level_str + empty + 6 * not_empty + time
row_str_no_both = level_str + empty + 4 * not_empty + empty + not_empty + time

out(
    "- This shows the highest velocity present in the following components \n"
    "  at each level.\n"
    "- Velocities reported (other than sound speed) are the maximum of the \n"
    "  x, y, and z component of the velocity for the star, gas, or DM. \n"
    "  This is because ART does it this way and I want to be consistent.\n"
    "- All velocities are in km/s, cell size in pc, dt in years.\n"
    "- The gas cell with the highest v_tot = bulk + c_s is selected, then the\n"
    "  bulk motion and sound speed for that cell are reported.\n"
    "- dt is the equivalent to the timestep in ART, i.e. \n"
    "  0.5 * cell_size / v_tot, where v_tot is defined above.\n"
    "- Only the fastest {} DM and star particles in the entire simulation box\n"
    "  are selected, so some levels may have no DM or stars listed.\n"
    "".format(n_vels)
)
out(
    header_str.format(
        "Level",
        "num_cells",
        "DM",
        "Gas bulk",
        "Gas c_s",
        "Gas v_tot",
        "Fast %",
        "Stars",
        "Cell Size",
        "dt",
    )
)

for level, cell_size, n_cells in zip(grid_levels, cell_sizes, num_in_grid):
    idx_gas = np.where(levels_gas == level)
    idx_star = np.where(levels_star == level)
    idx_dm = np.where(levels_dm == level)

    # get the maximum total gas velocity, then report the components of that
    # that contribute to the total
    vel_max_gas_tot_idx = np.argmax(v_tot_gas[idx_gas])
    vel_max_gas_tot = v_tot_gas[idx_gas][vel_max_gas_tot_idx]
    vel_max_gas_bulk = velocity_bulk_gas[idx_gas][vel_max_gas_tot_idx]
    vel_max_gas_cs = sound_speed_gas[idx_gas][vel_max_gas_tot_idx]

    # Calculate the fraction of cells above 1000 km/s
    fast_mask = velocity_bulk_gas[idx_gas] > 975 * yt.units.km / yt.units.s
    fast_percent = 100 * np.sum(fast_mask) / len(idx_gas[0])

    # the timestep ART chooses only depends on the gas
    dt = 0.5 * cell_size / (vel_max_gas_tot * yt.units.km / yt.units.s)
    dt = dt.to("yr").value

    # get the maximum velocity at this level, if we have particles here
    if len(idx_star[0]) > 0:  # we do have stars at this level
        vel_max_star = np.max(velocity_star[idx_star])
    if len(idx_dm[0]) > 0:  # DM at this level
        vel_max_dm = np.max(velocity_dm[idx_dm])

    # then decide what to print
    if len(idx_star[0]) > 0:  # stars
        if len(idx_dm[0]) > 0:  # stars and DM
            out_str = row_str.format(
                level,
                n_cells,
                vel_max_dm,
                vel_max_gas_bulk,
                vel_max_gas_cs,
                vel_max_gas_tot,
                fast_percent,
                vel_max_star,
                cell_size.to("pc").value,
                dt,
            )
        else:  # stars, no DM
            out_str = row_str_no_dm.format(
                level,
                n_cells,
                vel_max_gas_bulk,
                vel_max_gas_cs,
                vel_max_gas_tot,
                fast_percent,
                vel_max_star,
                cell_size.to("pc").value,
                dt,
            )
    else:  # no stars
        if len(idx_dm[0]) > 0:  # no stars, but DM
            out_str = row_str_no_star.format(
                level,
                n_cells,
                vel_max_dm,
                vel_max_gas_bulk,
                vel_max_gas_cs,
                vel_max_gas_tot,
                fast_percent,
                cell_size.to("pc").value,
                dt,
            )
        else:  # no stars, no dm
            out_str = row_str_no_both.format(
                level,
                n_cells,
                vel_max_gas_bulk,
                vel_max_gas_cs,
                vel_max_gas_tot,
                fast_percent,
                cell_size.to("pc").value,
                dt,
            )

    out(out_str)

# print the information of the highest-velocity cells/particles
out("\nHere are the cells/particles with the highest velocities.")
n_each = 10  # how many cells/DM particles/star particles to print
# Sort the velocities, then get the biggest ones (that will be at the end).
# The indexing is weird. We want to go backwards, (biggest first), so the step
# must be -1. Since we're going backwards, we end once we count back n_each, so
# that goes in the end part of the slice. We have to subtract one, since the
# final value isn't included in the slice
idx_sort_gas_bulk = np.argsort(velocity_bulk_gas)[: -n_each - 1 : -1]
idx_sort_gas_c_s = np.argsort(sound_speed_gas)[: -n_each - 1 : -1]
idx_sort_gas_tot = np.argsort(v_tot_gas)[: -n_each - 1 : -1]
idx_sort_dm = np.argsort(velocity_dm)[: -n_each - 1 : -1]
idx_sort_star = np.argsort(velocity_star)[: -n_each - 1 : -1]

header_str = "{:>20s}\t" + 2 * "{:<10s}" + 3 * "{:>20s}"
row_str = "{:>20.10f}\t" + "{:<10.0f}" + "{:<10.3E}" + 3 * "{:>20.10f}"

for name in ["DM", "Gas bulk", "Gas c_s", "Gas total", "Stars"]:
    # get the info for each of the components
    if "Gas" in name:
        levels = levels_gas
        masses = ad[("gas", "cell_mass")].to("Msun").value
        # haven't gotten the positions for the gas yet, so do that now. This
        # is a bit ugly, but is useful to make the later code easier
        pos_gas_x = ad[("gas", "x")].to("code_length").value
        pos_gas_y = ad[("gas", "y")].to("code_length").value
        pos_gas_z = ad[("gas", "z")].to("code_length").value
        positions = np.stack([pos_gas_x, pos_gas_y, pos_gas_z], axis=1)
        positions = ds.arr(positions, "code_length")
    if name == "Gas bulk":
        idxs = idx_sort_gas_bulk
        velocities = velocity_bulk_gas
    elif name == "Gas c_s":
        idxs = idx_sort_gas_c_s
        velocities = sound_speed_gas
    elif name == "Gas total":
        idxs = idx_sort_gas_tot
        velocities = v_tot_gas
    elif name == "DM":
        idxs = idx_sort_dm
        velocities = velocity_dm
        levels = levels_dm
        masses = ad[("N-BODY", "particle_mass")][idx_fast_dm].to("Msun").value
        positions = position_dm
    elif name == "Stars":
        idxs = idx_sort_star
        velocities = velocity_star
        levels = levels_star
        masses = ad[("STAR", "MASS")][idx_fast_star].to("Msun").value
        positions = position_star

    out("\n{} Highest Velocities".format(name))
    out(
        header_str.format(
            "Velocity [km/s]",
            "Level",
            "Mass",
            "Location X [code]",
            "Location Y [code]",
            "Location Z [code]",
        )
    )
    # go through each of the highest velocity cells and print their information
    for idx in idxs:
        x, y, z = positions[idx].to("code_length").value
        out(row_str.format(velocities[idx], levels[idx], masses[idx], x, y, z))

out_file.close()
