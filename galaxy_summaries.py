"""
galaxy_summaries.py

Reports the properties of the galaxies in the simulation.

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
        file_obj.write(info + "\n")
    if not silent:
        print(info)

ds_loc = os.path.abspath(sys.argv[1])
scale_factor = ds_loc[-10:-4]
ds = yt.load(ds_loc)
ad = ds.all_data()

is_zoom = ('N-BODY_0', 'POSITION_X') in ds.field_list

# get the location of where to write the file.
sim_dir = os.path.dirname(ds_loc) + os.sep
file_dir = sim_dir.replace("/out/", "/checks/")
file_path = file_dir + "galaxy_summaries_a" + scale_factor + ".txt"
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

# define a shorthand function that handles the parameters to out,
# so we don't have to duplicate that each time
def out(info):
    return print_and_write(info, out_file)

# =========================================================================
#         
# Halo analysis
# 
# =========================================================================
halo_file = os.path.abspath(sys.argv[2])
ds_halos = yt.load(halo_file)

# Then create the halo catalogs
hc = HaloCatalog(halos_ds=ds_halos, data_ds=ds, output_dir="./")
# Restrict to things about LMC mass and above
hc.add_filter('quantity_value', 'particle_mass', '>', 1E9, 'Msun')
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
    out("No DM halos above 10^9 Msun at this redshift.")
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

# double check that the boundary conditions have been correctly respected. I'm not sure
# why this happens, but it seems not to be respected at times.

for halo in hc.catalog:
    for position_name in ["particle_position_x", "particle_position_y", "particle_position_z"]:
        while halo[position_name] > ds.domain_width[0]:
            # make sure subtraction doesn't happen in code units
            print(f"fixing halo {halo['rank']} position")
            halo[position_name] = halo[position_name].to("cm") - ds.domain_width[0].to("cm")

# get the N-body particle locations. These are different in the old and new
# simulations, so we have to check 
species_x = dict()
species_y = dict()
species_z = dict()
if is_zoom:
    idx = 0
    while True:
        try:
            species_x[idx] = ad[('N-BODY_{}'.format(idx), 'POSITION_X')].to("Mpc").value
            species_y[idx] = ad[('N-BODY_{}'.format(idx), 'POSITION_Y')].to("Mpc").value
            species_z[idx] = ad[('N-BODY_{}'.format(idx), 'POSITION_Z')].to("Mpc").value
            idx += 1
        except yt.utilities.exceptions.YTFieldNotFound:
            break
else:  # old sims
    species_x[0] = ad[('N-BODY', 'POSITION_X')].to("Mpc").value
    species_y[0] = ad[('N-BODY', 'POSITION_Y')].to("Mpc").value
    species_z[0] = ad[('N-BODY', 'POSITION_Z')].to("Mpc").value

# define some helper functions that will be used for contamination calculations
def get_center(halo, with_units):
    x_cen = halo["particle_position_x"]
    y_cen = halo["particle_position_y"]
    z_cen = halo["particle_position_z"]
    if not with_units:
        x_cen = x_cen.to("Mpc").value
        y_cen = y_cen.to("Mpc").value
        z_cen = z_cen.to("Mpc").value
    return (x_cen, y_cen, z_cen)

def distance(x_0, y_0, z_0, x_1, y_1, z_1):
    return np.sqrt((x_1 - x_0)**2 + (y_1 - y_0)**2 + (z_1 - z_0)**2)

def mass_fractions(sphere):
    # Prints the fraction of mass inside a sphere that comes from different
    # species of N-body particles. This is useful for checking contamination
    masses = sphere[('N-BODY', 'MASS')]
    unique_masses, num_particles = np.unique(masses, return_counts=True)
    total_mass = masses.sum()
    out("Total Mass: {:.3e}".format(total_mass.to("Msun")))
    out("{:<7s} {:>10s} {:>15s} {:>15s}".format("Species", "Number", "Mass", "Mass Fraction"))
    for m, n_m in zip(unique_masses, num_particles):
        this_m_tot = (m * n_m).to("Msun")  # total mass in this species of particle
        frac = (this_m_tot / total_mass).value
        # then get the species. We take advantage of the fact the the unique masses
        # are sorted, with the smallest mass (last species) first.
        s = len(unique_masses) - np.where(unique_masses == m)[0][0] - 1

        out("{:<7} {:>10,} {:>10.2e} {:>14.2f}%".format(s, n_m, this_m_tot, frac*100))

# =========================================================================
#         
# Then we go through and print information about the halos present
# 
# =========================================================================
n_halos = 10
for halo in halos[:n_halos]:
    out("\n==================================\n")
    out("Rank {} halo:".format(halo["rank"]))
    # First print the important quantities
    for quantity in ordered_quantities:   
        q_name = quantity_names[quantity]
        if q_name == "Virial Radius" or "Position" in q_name:
            value_kpc = halo[quantity].to("kpc").value
            value_code = (halo[quantity] / ds.length_unit).value
            out("{}: {:>7.3f} kpc, {:>7.3f} code".format(q_name, value_kpc, value_code))
        elif "Mass" in q_name:
            value = halo[quantity].to("Msun")
            out("{}: {:<2.3e}".format(q_name, value))

    # get some spheres that will be used in the calculations of galaxy 
    # properties
    virial_radius = halo["virial_radius"]
    center = get_center(halo, with_units=True)
    sphere_virial = ds.sphere(center=center, radius=virial_radius)
    sphere_mpc = ds.sphere(center=center, radius=1*yt.units.Mpc)
    sphere_30_kpc = ds.sphere(center=center, radius=30*yt.units.kpc)

    # then print information about contamination, if we need to
    if is_zoom:
        out("\nClosest particle of each low-res DM species")
         # First we calculate the closest particle of each unrefined DM species
        x_cen, y_cen, z_cen = get_center(halo, with_units=False)
        for idx in species_x:
            distances = distance(x_cen, y_cen, z_cen, 
                                 species_x[idx], species_y[idx], species_z[idx])
            out("{}: {:.0f} kpc".format(idx, np.min(distances)*1000))

        for sphere, name in zip([sphere_virial, sphere_mpc],
                                ["the virial radius", "1 Mpc"]):
            out("\nDM mass and fraction of total from different "
                "species within {}".format(name))
            mass_fractions(sphere)
    else:
        total_mass = sphere_virial[('N-BODY', 'MASS')].to("Msun").value.sum()
        out("\nTotal DM particle mass within the virial radius = {:.3e} Msun".format(total_mass))

    # print the stellar mass 
    stellar_mass_30kpc = np.sum(sphere_30_kpc[('STAR', 'MASS')].to("Msun").value)
    stellar_mass_virial = np.sum(sphere_virial[('STAR', 'MASS')].to("Msun").value)
    out("")
    out("Stellar Mass within 30 kpc: {:.2e} Msun".format(stellar_mass_30kpc))
    out("Stellar Mass within the virial radius: {:.2e} Msun".format(stellar_mass_virial))

    # Print the masses in various states
    cell_volumes = sphere_30_kpc[('index', 'cell_volume')]
    gas_mass = sphere_30_kpc[('gas', 'cell_mass')]
    gas_mass_HI = sphere_30_kpc[('gas', 'HI density')] * cell_volumes
    gas_mass_HII = sphere_30_kpc[('gas', 'HII density')] * cell_volumes
    gas_mass_H2 = 2 * sphere_30_kpc[('artio', 'RT_HVAR_H2')] * ds.arr(1, "code_mass/code_length**3") * cell_volumes
    gas_mass_HeI = 4 * sphere_30_kpc[('gas', 'HeI density')] * cell_volumes
    gas_mass_HeII = 4 * sphere_30_kpc[('gas', 'HeII density')] * cell_volumes
    gas_mass_HeIII = 4 * sphere_30_kpc[('gas', 'HeIII density')] * cell_volumes
    gas_mass_metals = sphere_30_kpc[('gas', 'metal_density')] * cell_volumes

    gas_mass_total = np.sum(gas_mass.to("Msun").value)
    gas_mass_HI_total = np.sum(gas_mass_HI.to("Msun").value)
    gas_mass_HII_total = np.sum(gas_mass_HII.to("Msun").value)
    gas_mass_H2_total = np.sum(gas_mass_H2.to("Msun").value)
    gas_mass_HeI_total = np.sum(gas_mass_HeI.to("Msun").value)
    gas_mass_HeII_total = np.sum(gas_mass_HeII.to("Msun").value)
    gas_mass_HeIII_total = np.sum(gas_mass_HeIII.to("Msun").value)
    gas_mass_metals_total = np.sum(gas_mass_metals.to("Msun").value)

    out("")
    out("Gas masses within 30 kpc:")
    out("Total:  {:.2e} Msun".format(gas_mass_total))
    out("HI:     {:.2e} Msun".format(gas_mass_HI_total))
    out("HII:    {:.2e} Msun".format(gas_mass_HII_total))
    out("H2:     {:.2e} Msun".format(gas_mass_H2_total))
    out("He1:    {:.2e} Msun".format(gas_mass_HeI_total))
    out("HeII:   {:.2e} Msun".format(gas_mass_HeII_total))
    out("HeIII:  {:.2e} Msun".format(gas_mass_HeIII_total))
    out("Metals: {:.2e} Msun".format(gas_mass_metals_total))


    # Metallicity
    metallicity = sphere_30_kpc[('STAR', 'METALLICITY_SNII')].value
    metallicity += sphere_30_kpc[('STAR', 'METALLICITY_SNIa')].value
    if ('STAR', 'METALLICITY_AGB') in ds.field_list:
        metallicity += sphere_30_kpc[('STAR', 'METALLICITY_AGB')].value
    # get the masses so we can mass-weight it
    masses = sphere_30_kpc[('STAR', 'MASS')].to("Msun").value
    metal_mass = masses * metallicity

    out("")
    out("Stellar metallicity (nans mean no stars)")
    out("Mass-weighted mean metallicity of:")
    for max_age in [np.inf, 40]:
        if np.isinf(max_age):
            descriptor = "all stars"
        else:
            descriptor = f"stars younger than {max_age} Myr"

        star_ages = sphere_30_kpc[('STAR', 'age')].to("Myr").value
        # get the mass-weighted metallicity
        if len(star_ages) == 0 or np.min(star_ages) > max_age:
            mean_metallicity = np.nan
        else:
            star_mask = star_ages < max_age
            mean_metallicity = np.sum(metal_mass[star_mask]) / np.sum(masses[star_mask])

        out(f"{descriptor:>25} = {mean_metallicity:.2e} -> log(Z/Z_sun) = {np.log10(mean_metallicity/0.02):.2f}")

        
out("\n==================================\n")
        
# Then print the separation of the two biggest halos
if len(halos) >= 2:
    halo_1 = halos[0]
    halo_2 = halos[1]
    dx = halo_1["particle_position_x"] - halo_2["particle_position_x"]
    dy = halo_1["particle_position_y"] - halo_2["particle_position_y"]
    dz = halo_1["particle_position_z"] - halo_2["particle_position_z"]
    dist = np.sqrt(dx**2 + dy**2 + dz**2).to("kpc")
    out("\nSeparation of two largest halos: {:.2f}".format(dist))


# Saving this for posterity - way to get multiple species into one field
# n_body_density_field = ("deposit", "High-Res_Dark_Matter_Density")
# # We want to set up a deposit N-body density field that has only the high res
# # particles. 
# def _n_density_high_res(field, data):
#     # We need to check for sims without the high res particles
#     try:
#         return data[("deposit", "N-BODY_0_density")] + \
#                data[("deposit", "N-BODY_1_density")]
#     except yt.utilities.exceptions.YTFieldNotFound:
#         return data[("deposit", "N-BODY_density")]
# ds.add_field(n_body_density_field, function=_n_density_high_res, 
#              units="g/cm**3", sampling_type="cell")

# then find the center, which will be the median of the high res particles, or
# the domain center, if there are no high res particles


out_file.close()
