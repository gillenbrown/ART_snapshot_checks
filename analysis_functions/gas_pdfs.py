import numpy as np
import yt

# ==============================================================================
#
# Gas PDFs
#
# ==============================================================================
def get_gas_density(region):
    n = region[("gas", "density")] / yt.units.mass_hydrogen
    return n.to("cm**(-3)").value


def get_gas_volume(region):
    return region[("gas", "cell_volume")].to("pc**3").value


def get_gas_mass(region):
    return region[("gas", "cell_mass")].to("Msun").value


def get_gas_temp(region):
    return region[("gas", "temperature")].to("K").value


def get_gas_p_over_k(region):
    p_over_k = region[("gas", "pressure")] / yt.units.kboltz
    return p_over_k.to("K/cm**3").value


def get_gas_turbulent_temp(region):
    # Taken from src/core/hydro_sgst.c
    gamma = 5 / 3
    wmu = 4.0 / (8.0 - 5.0 * 0.24)
    density = region[("artio", "HVAR_GAS_DENSITY")]
    # add units. energy density is
    # energy / length^3 = mass * velocity^2 / length^3
    turbulent_energy_density = (
        region[("artio", "HVAR_GAS_TURBULENT_ENERGY")]
        * region.ds.mass_unit
        * region.ds.velocity_unit ** 2
        / region.ds.length_unit ** 3
    )
    # then use the equation from ART. Also need to use unit conversions.
    temp = (
        (gamma - 1)
        * wmu
        * turbulent_energy_density
        * yt.units.mass_hydrogen  # technically is mb in ART, but that's basically m_H
        / (density * yt.units.kboltz)
    )
    return temp.to("K").value


def get_gas_turbulent_velocity_dispersion(region):
    # taken from sf_recipe.turbulent.c
    density = region[("artio", "HVAR_GAS_DENSITY")]
    # add units. energy density is
    # energy / length^3 = mass * velocity^2 / length^3
    turbulent_energy_density = (
        region[("artio", "HVAR_GAS_TURBULENT_ENERGY")]
        * region.ds.mass_unit
        * region.ds.velocity_unit ** 2
        / region.ds.length_unit ** 3
    )

    sigma_squared = 2 * turbulent_energy_density / density
    sigma = np.sqrt(sigma_squared).to("km/s")
    return sigma


def get_gas_metallicity_ii(region):
    z_density = region[("artio", "HVAR_METAL_DENSITY_II")]
    gas_density = region[("artio", "HVAR_GAS_DENSITY")]
    return (z_density / gas_density).to("").value


def get_gas_velocity(region):
    vx = region[("gas", "velocity_x")].to("km/s").value
    vy = region[("gas", "velocity_y")].to("km/s").value
    vz = region[("gas", "velocity_z")].to("km/s").value

    mass = get_gas_mass(region)

    mean_vx = np.average(vx, weights=mass)
    mean_vy = np.average(vy, weights=mass)
    mean_vz = np.average(vz, weights=mass)

    dvx = vx - mean_vx
    dvy = vy - mean_vy
    dvz = vz - mean_vz

    gas_velocity = np.sqrt(dvx ** 2 + dvy ** 2 + dvz ** 2)
    return gas_velocity


def get_gas_level(region):
    return region[("index", "grid_level")].to("").value


def cumulative_property(region, prop_func, weight_func):
    """
    Create the cumulative distribution of the given property.

    :param region: The yt region to extract the data from.
    :param prop_func: A function to extract a property from `region`. This function must
                      take the region as the only argument, and return dimensionless
                      values.
    :param weight_func: A function to extract the weights to be used in the cumuative
                        distribution. Passing None will provide no weighting, meaning
                        the histogram will show the number of cells. Passing volume
                        or mass will result in the cumulative mass fraction being
                        returned.
    :returns: Two arrays. One for the sorted values of properties, then the second with
              the cumulative weights that are less than or equal to this value. This
              is done so that these two lists can simply be plotted as x and y.
    """
    values = prop_func(region)

    if weight_func is None:
        weights = np.ones(values.shape)
    else:
        weights = weight_func(region)

    sort_idxs = np.argsort(values)
    values = values[sort_idxs]
    weights = weights[sort_idxs]
    cumulative_weight = np.cumsum(weights)
    cumulative_weight = cumulative_weight / cumulative_weight[-1]

    return values, cumulative_weight


def pdf_property(region, prop_func, weight_func, x_min, x_max, bin_size):
    """
    Create a pdf of gas properties.

    Specifically, this will return d(Weight) / d(log(property)).

    :param region: The yt region to extract the data from.
    :param prop_func: A function to extract a property from `region`. This function must
                      take the region as the only argument, and return dimensionless
                      values.
    :param weight_func: A function to extract the weights to be used in the cumuative
                        distribution. Passing None will provide no weighting, meaning
                        the histogram will show the number of cells per bin.
    :param x_min: Minimum bin center.
    :param x_max: Maximum bin center.
    :param bin_size: Bin size (in dex)
    :returns: Two arrays. One for the bin centers of property, the second the value of
              the pdf at that value. This is so they can simply be plotted as x and y.
    """
    boundaries_log = np.arange(x_min - 0.5 * bin_size, x_max + 0.5 * bin_size, bin_size)
    centers_log = [
        np.mean([boundaries_log[idx], boundaries_log[idx + 1]])
        for idx in range(len(boundaries_log) - 1)
    ]

    boundaries = 10 ** boundaries_log
    centers = 10 ** np.array(centers_log)

    values = prop_func(region)

    if weight_func is None:
        weights = np.ones(values.shape)
    else:
        weights = weight_func(region)

    dx = np.histogram(a=values, bins=boundaries, weights=weights)[0]
    # make dx per log rho
    return centers, dx / (bin_size * np.log(10))
