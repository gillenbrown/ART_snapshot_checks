from pathlib import Path
import numpy as np
from scipy import special
from astropy import cosmology, table
from astropy import units as u
import yt

from tqdm import tqdm

# ======================================================================================
#
# CIMF calculations
#
# ======================================================================================
# Then the functions to calculate the CIMF. Here we need to do some analysis
# of the bound fraction.
def f_bound(eps_int):
    # Li et al 2019: https://ui.adsabs.harvard.edu/abs/2019MNRAS.487..364L/abstract
    # equation 17
    alpha_star = 0.48
    f_sat = 0.94
    term_a = special.erf(np.sqrt(3 * eps_int / alpha_star))
    term_b = np.sqrt(12 * eps_int / (np.pi * alpha_star))
    term_c = np.exp(-3 * eps_int / alpha_star)
    return (term_a - (term_b * term_c)) * f_sat


def get_initial_bound_fraction(galaxy):
    star_initial_mass = galaxy[("STAR", "INITIAL_MASS")].to("Msun").value
    # the variable named INITIAL_BOUND_FRACTION is not the initial_bound fraction,
    # it's actually the accumulated mass nearby through the course of accretion, in
    # code masses. This is used to calculate the formation efficiency, which is then
    # used to get the bound fraction.
    star_accumulated_mass = galaxy[("STAR", "INITIAL_BOUND_FRACTION")].to("").value
    star_accumulated_mass *= galaxy.ds.mass_unit
    star_accumulated_mass = star_accumulated_mass.to("Msun").value

    eps_int = star_initial_mass / star_accumulated_mass

    return f_bound(eps_int)


# ======================================================================================
#
# tidal evolution to z=0
#
# ======================================================================================
omega_tid = 175 / yt.units.Gyr  # calibrated to match MW GC number
# dt = 100 * yt.units.Myr


def t_tidal(M):
    term_1 = 10 * yt.units.Gyr
    term_2 = (M / (2e5 * yt.units.Msun)) ** (2 / 3)
    term_3 = 100 / (omega_tid * yt.units.Gyr)
    total = term_1 * term_2 * term_3
    return total.to("Myr")


# def get_n_steps(t_now, t_z_0):
#     return int(np.ceil((t_z_0 - t_now).to("Myr") / dt))
#
#
# precalculated_disruption = dict()
#
#
# def precalculate_disruption(n_steps):
#     if n_steps in precalculated_disruption:
#         return
#     else:
#         precalculated_disruption[n_steps] = dict()
#
#     # otherwise, we need to build this
#     print(f"\n\n\nprecalculating\n{n_steps}\n\n\n")
#     for log_m in tqdm(np.arange(2, 8, 0.01)):
#         key_m = str(round(log_m, 2))
#
#         m = 10 ** log_m * yt.units.Msun
#         for _ in range(n_steps):
#             if m < 100 * yt.units.Msun:
#                 m = 0
#                 break
#             m *= np.exp(-dt / t_tidal(m))
#         precalculated_disruption[n_steps][key_m] = m
#
#
# def get_evolution_single_cluster(n_steps, log_m):
#     try:
#         key_m = str(round(log_m, 2))
#         return precalculated_disruption[n_steps][key_m]
#     except KeyError:  # only happens for initial M below 100, which will be fully
#         # disrupted by z=0
#         return 0
#
#
# def evolve_cluster_population_old(galaxy):
#     # get cosmology to get times
#     # Initialize the cosmology object, used to put a redshift scale on the plots
#     H_0 = galaxy.ds.artio_parameters["hubble"][0] * 100 * u.km / (u.Mpc * u.second)
#     omega_matter = galaxy.ds.artio_parameters["OmegaM"][0]
#     cosmo = cosmology.FlatLambdaCDM(H0=H_0, Om0=omega_matter, Tcmb0=2.725)
#
#     t_now = galaxy.ds.current_time
#     # need to convert the astropy units into yt units.
#     t_z_0 = cosmo.age(0).to("Gyr").value * yt.units.Gyr
#     n_steps = get_n_steps(t_now, t_z_0)
#
#     # precalculate the disruption. This won't do anything if it already exists.
#     precalculate_disruption(n_steps)
#
#     # then get star masses (include stellar evolution)
#     raw_mass = galaxy[("STAR", "MASS")].to("Msun")
#     star_initial_bound = get_initial_bound_fraction(galaxy)
#     tidal_bound_fraction = galaxy[("STAR", "BOUND_FRACTION")].value
#     cluster_masses = raw_mass * star_initial_bound * tidal_bound_fraction
#     # get the log here. This makes using it in precalculation easier, as this is
#     # vectorized. Set a minimum of 0.1 to avoid log of zero errors.
#     log_cluster_masses_msun = np.log10(np.maximum(0.1, cluster_masses.to("Msun").value))
#
#     evolved_masses = [
#         get_evolution_single_cluster(n_steps, log_m)
#         for log_m in tqdm(log_cluster_masses_msun)
#     ]
#     return np.array(evolved_masses)


def evolve_cluster_population(galaxy):
    # get cosmology to get times
    # Initialize the cosmology object, used to put a redshift scale on the plots
    H_0 = galaxy.ds.artio_parameters["hubble"][0] * 100 * u.km / (u.Mpc * u.second)
    omega_matter = galaxy.ds.artio_parameters["OmegaM"][0]
    cosmo = cosmology.FlatLambdaCDM(H0=H_0, Om0=omega_matter, Tcmb0=2.725)

    t_now = galaxy.ds.current_time
    # need to convert the astropy units into yt units.
    t_z_0 = cosmo.age(0).to("Gyr").value * yt.units.Gyr
    dt = t_z_0 - t_now

    # then get star masses (include stellar evolution)
    raw_mass = galaxy[("STAR", "MASS")].to("Msun")
    star_initial_bound = get_initial_bound_fraction(galaxy)
    tidal_bound_fraction = galaxy[("STAR", "BOUND_FRACTION")].value
    cluster_masses = raw_mass * star_initial_bound * tidal_bound_fraction

    timescales_initial = t_tidal(cluster_masses)
    evolved_masses = cluster_masses * (1 - 2 * dt / (3 * timescales_initial)) ** (3 / 2)
    # replace nans with zeros
    bad_idx = np.isnan(evolved_masses)
    evolved_masses[bad_idx] = 0

    return evolved_masses.to("Msun").value


# ======================================================================================
#
# CIMF itself
#
# ======================================================================================
def _cimf_base(masses, bin_width=0.16):
    # create bins with spacing of 0.16 dex.
    # I use the maximum mass as the second to last bin, so that it drops to
    # zero after that. But I do have to check for having zero clusters, which I think
    # only happens for the evolved to z=0 plot or one restricting the age range
    if len(masses) == 0:
        max_mass = 1e5
    else:
        max_mass = max(np.max(masses), 1e4)  # guard against no high mass clusters
    m_centers_log = np.arange(np.log10(max_mass) + bin_width, 2, -bin_width)[::-1]
    m_boundaries_log = np.arange(
        np.min(m_centers_log) - 0.5 * bin_width,
        np.max(m_centers_log) + 0.6 * bin_width,
        bin_width,
    )

    m_boundaries = 10 ** m_boundaries_log
    m_centers = 10 ** np.array(m_centers_log)

    # then make the histogram showing how many there are per bin
    hist, edges = np.histogram(masses, bins=m_boundaries)
    assert np.array_equiv(m_boundaries, edges)

    # We have dN, make it per dLogM
    hist = np.array(hist) / (bin_width * np.log(10))

    return m_centers, hist


def cimf(sim, mass_type, max_age_myr, max_z):
    """
    Make the cluster initial mass function.

    Notes on the different mass variables:
    ('STAR', 'INITIAL_MASS') - initial stellar mass of the star particle
    ('STAR', 'MASS') - current stellar mass of the star particle, accounting for
        stellar evolution
    ('STAR', 'INITIAL_BOUND_FRACTION') - NOT the actual initial_bound fraction. See
        the `get_initial_bound_fraction` function above for more on how to use this,
        but this variable is the accumulated mass near the cluster over the course
        of accretion. This is used to calculate formation efficiency, which is then
        used to get the actual initial_bound fraction
    ('STAR', 'BOUND_FRACTION') - This is the actual bound fraction at the current
        time, but NOT accounting for the proper initial_bound fraction

    :param sim: simulation object object
    :param mass_type: String encoding which mass to get here. The options are:
                      "initial" - just the initial stellar masses
                      "initial_bound" - initial masses including initial_bound
                                        fraction
                      "current" - current bound mass, accounting for the initial
                                  bound fraction, tidal disruption, and stellar
                                  death
    :param include_initial_bound: whether to incorporate the initial bound
                                  fraction of clusters, or just get the
                                  distribution of initial particle masses
    :param max_age_myr: The maximum age to restrict the plot to. Is infinity as the
                        default, which plots all stars.
    :param max_z: The maximum metallicity to restrict the plot to. Is infinity as the
                  default, which plots all stars.
    :returns: Two lists. The first is f_i * M, representing the initial
              bound mass, of M if include_initial_bound=False. This will be
              binned values suitable to plot. The second is dN/dLogM for each
              of the bins in the first list.
    """
    id_string = f"{mass_type}_{max_age_myr}_{max_z}"

    if id_string not in sim.precalculated:
        mass = []
        for galaxy in sim.galaxies:
            if mass_type == "initial":
                this_mass = galaxy[("STAR", "INITIAL_MASS")].to("Msun").value
            elif mass_type == "initial_bound":
                initial_mass = galaxy[("STAR", "INITIAL_MASS")].to("Msun").value
                star_initial_bound = get_initial_bound_fraction(galaxy)
                this_mass = initial_mass * star_initial_bound
            elif mass_type == "current":
                # use mass that accounts for stellar evolution
                raw_mass = galaxy[("STAR", "MASS")].to("Msun").value
                star_initial_bound = get_initial_bound_fraction(galaxy)
                tidal_bound_fraction = galaxy[("STAR", "BOUND_FRACTION")].value
                this_mass = raw_mass * star_initial_bound * tidal_bound_fraction
            elif mass_type == "evolved":
                this_mass = evolve_cluster_population(galaxy)
            else:
                raise ValueError("Mass not recognized")

            # then restrict the metallicity and age
            metallicity = galaxy[("STAR", "METALLICITY_SNII")].value
            metallicity += galaxy[("STAR", "METALLICITY_SNIa")].value
            if ("STAR", "METALLICITY_AGB") in sim.ds.field_list:
                metallicity += galaxy[("STAR", "METALLICITY_AGB")].value
            mask_z = metallicity < max_z
            mask_age = galaxy[("STAR", "age")] < max_age_myr * yt.units.Myr
            this_mass = this_mass[np.logical_and(mask_z, mask_age)]

            mass = np.concatenate([this_mass, mass])

        m_centers, hist = _cimf_base(mass)

        # Also normalize for the number of galaxies
        hist = hist / sim.n_galaxies

        sim.precalculated[id_string] = m_centers, hist
    return sim.precalculated[id_string]


def harric_gc_mass_function():
    # Note that I carefully checked this against Figure 3 from Hui paper 3. My
    # normalization is different than his. I think he forgot the ln(10) factor when
    # using dlogM. When I remove that factor, I get similar normalizations as him.
    cat_path = Path(__file__).resolve().parent.parent / "data" / "harris_gc_catalog.txt"
    harris_catalog = table.Table.read(str(cat_path), format="ascii")
    abs_mag = harris_catalog["M_V,t"]
    harris_catalog["m_to_l"] = 1.3 + 4.5 / (1 + np.exp(2 * abs_mag + 21.4))
    # convert abs_mag to luminosity
    # M - M_sun = -2.5 log (L / L_sun)
    # (M_sun - M) / 2.5 = log(L / L_sun)
    # L = L_sun * 10**((M - M_sun) / 2.5)
    M_sun = 4.81  # Willmer 2018 ApJS 236 47
    harris_catalog["luminosity"] = 10 ** ((M_sun - abs_mag) / 2.5)
    harris_catalog["mass"] = harris_catalog["luminosity"] * harris_catalog["m_to_l"]
    return _cimf_base(harris_catalog["mass"])
