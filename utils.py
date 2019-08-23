import numpy as np
import yt

def print_and_write(info, file_obj, silent):
    if file_obj is not None:
        file_obj.write(str(info) + "\n")
    if not silent:
        print(info)

def a_to_z(a):
    return (1.0 / a) - 1.0

def z_to_a(z):
    return 1.0 / (1.0 + z)

def get_stellar_births(data_obj, initial=False):
    """
    Get the birth times and initial mass for all star particles in a region.

    :param data_obj: region of the simulation that stars will be selected from.
                     Can be something like a sphere, or even all_data()
    :param initial: Whether or not to use the initial mass of star particles,
                    rather than their current masses (NOT cluster bound masses).
                    Using initial masses is better when calculating SFH, while
                    turning this off is better for the cumulative stellar mass
                    history, since it sums to the present stellar mass.
    :return: Two arrays: first the masses in solar masses, then the 
             birth times of the stars.
    """
    if initial:
        masses = data_obj[("STAR", "INITIAL_MASS")].in_units("msun")
    else:
        masses = data_obj[("STAR", "MASS")].in_units("msun")
    creation_times = data_obj[("STAR", "creation_time")]
    return masses, creation_times

def create_cumulative_mass(data_obj):
    """
    Create the cumulative stellar mass of the stars within some region.

    This will create many timesteps, then for each timestep record the mass
    that formed earlier than this time.

    :param data_obj: region of the simulation that stars will be selected from.
                     Can be something like a sphere, or even all_data()
    """
    masses, creation_times = get_stellar_births(data_obj, initial=True)
    
    n_steps = 1000
    # the last timestep should include the current time of the simulation. We
    # round the current time up to make sure that happens
    max_time = np.ceil(data_obj.ds.current_time)
    timesteps = np.linspace(0, max_time, n_steps)
    
    cumulative_mass = []
    for time in timesteps:
        # get all the stars that were born earlier than the current time
        idx = np.where(creation_times <= time)
        # then sum their masses
        cumulative_mass.append(np.sum(masses[idx]).in_units("msun"))
    
    return timesteps.in_units("Gyr"), yt.YTArray(cumulative_mass)

def sfh(data_obj):
    masses, creation_times = get_stellar_births(data_obj, initial=True)
    
    dt = 100 # Myr
    # add one extra bin, so that any current SF is captured, and not 
    # thrown away by arange not including the end point
    max_bin = data_obj.ds.current_time.to("Myr").value + dt
    bin_edges_age = np.arange(0, max_bin, dt) * yt.units.Myr
    timestep = bin_edges_age[1] - bin_edges_age[0]

    bin_centers = []
    sfr_values = []
    # go through each bin, seeing which stars were born there
    for left_idx in range(len(bin_edges_age) - 1):
        # get the age range
        min_age = bin_edges_age[left_idx]
        max_age = bin_edges_age[left_idx + 1]
        
        # then get the stars in that range
        age_more_idx = np.where(creation_times >= min_age)
        age_less_idx = np.where(creation_times <  max_age)
        age_good_idx = np.intersect1d(age_more_idx, age_less_idx)
        
        # and sum their masses
        this_formed_mass = np.sum(masses[age_good_idx])
        
        sfr_values.append(this_formed_mass / timestep)
        bin_centers.append((min_age + max_age) / 2.0)
        
    return yt.YTArray(bin_centers), yt.YTArray(sfr_values)