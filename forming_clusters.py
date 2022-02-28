"""
forming_clusters.py

Reports the properties of the clusters with ages less than 15 Myr in the current
snapshot.

Takes 2 required parameter.
1 - Location of the output file
2 - Location of the simulation output. Can be relative to the working directory
    where this code was called from, or an absolute path.

"""
from utils import load_galaxies as lg

import sys
from pathlib import Path

import yt
from analysis_functions import gas_pdfs

yt.funcs.mylog.setLevel(50)  # ignore yt's output

# =========================================================================
#
# Add certain fields to output so I can access them more cleanly
#
# =========================================================================
# I do this before loading the sim so that they're properly registered.
# I tried using lambdas to shorten this, but it didn't work for some reason
def _my_number_density(field, data):
    return gas_pdfs.get_gas_number_density(data)


def _h2_mass(field, data):
    return gas_pdfs.get_h2_mass(data)


def _h2_frac(field, data):
    return gas_pdfs.get_h2_frac(data)


def _my_temp(field, data):
    return gas_pdfs.get_gas_temp(data)


yt.add_field(
    name=("gas", "my_number_density"),
    function=_my_number_density,
    sampling_type="local",
    units="cm**(-3)",
)

yt.add_field(
    name=("gas", "H2_mass"),
    function=_h2_mass,
    sampling_type="local",
    units="g",
)

yt.add_field(
    name=("gas", "H2_frac"),
    function=_h2_frac,
    sampling_type="local",
    units="",
)

yt.add_field(
    name=("gas", "my_temperature"),
    function=_my_temp,
    sampling_type="local",
    units="K",
)

# =========================================================================
#
# load sims and prep output file
#
# =========================================================================
out_file = open(sys.argv[1], "w")
ds_loc = Path(sys.argv[2]).resolve()
sim = lg.Simulation(ds_loc, sphere_radius_kpc=30)

# define a shorthand function that handles the parameters to out,
# so we don't have to duplicate that each time
def out(info):
    out_file.write(info + "\n")


# print header for this output file
out("# This catalog contains columns with the following info: ")
out("# - Particle ID")
out("# - Cluster Age [Myr]")
out("# - number density of cluster cell [cm^(-3)]")
out("# - H2 mass of cluster cell [Msun]")
out("# - H2 fraction in cluster cell")
out("# - Temperature in cluster cell [K]")

# =========================================================================
#
# Analysis
#
# =========================================================================
for galaxy in sim.galaxies:
    # make sure the mask is created
    if galaxy.mask_done_forming is None:
        galaxy.make_finished_cluster_mask()

    # get positions of clusters
    x = galaxy.prop_all_clusters(("STAR", "POSITION_X"))[~galaxy.mask_done_forming]
    y = galaxy.prop_all_clusters(("STAR", "POSITION_Y"))[~galaxy.mask_done_forming]
    z = galaxy.prop_all_clusters(("STAR", "POSITION_Z"))[~galaxy.mask_done_forming]
    # then their ages and particle IDs
    ages = (
        galaxy.prop_all_clusters(("STAR", "age"))
        .to("Myr")
        .value[~galaxy.mask_done_forming]
    )
    ids = (
        galaxy.prop_all_clusters(("STAR", "PID"))
        .to("")
        .value[~galaxy.mask_done_forming]
    )

    # turn this into 3-tuples for position, which is what the function needs
    xyz = [(x[i], y[i], z[i]) for i in range(len(x))]

    # list the fields we request
    fields = [
        ("gas", "my_number_density"),
        ("gas", "H2_mass"),
        ("gas", "H2_frac"),
        ("gas", "my_temperature"),
    ]
    # then get these fields
    ns, h2_mass, h2_frac, temp = galaxy.ds.find_field_values_at_points(fields, xyz)
    # clean up units
    ns = ns.to("cm**(-3)").value
    h2_mass = h2_mass.to("Msun").value
    h2_frac = h2_frac.to("").value
    temp = temp.to("K").value

    # then print info
    for i in range(len(x)):
        out(
            f"{ids[i]:.0f} "
            f"{ages[i]:.2f} "
            f"{ns[i]:.2e} "
            f"{h2_mass[i]:.2e} "
            f"{h2_frac[i]:.2f} "
            f"{temp[i]:.2e}"
        )


out_file.close()