import sys
from pathlib import Path

from astropy import units as u
from astropy import constants as c
import numpy as np
from tqdm import tqdm

from matplotlib import colors, cm
import cmocean
import betterplotlib as bpl
bpl.set_style()

# Get the log directory specified by the user
log_dir = Path(sys.argv[1]).absolute()

# parse the output file first
# I used to use the stdout.000.log file, as it had the timestep info, but had the run
# info (which we discard) from only one MPI rank. Now, since we're investigating the
# CFL violations, we have to use the full output file. I need to copy this into the
# log directory with this name so the code can find it.
log_file = log_dir / "stdout.full.log"

# ======================================================================================
#
# Read in the file
#
# ======================================================================================
# first get the number of lines in the file, so we can tell the user how its going
with open(log_file, "r") as in_file:
    for i, l in enumerate(in_file):
        pass
total_lines = i + 1

def get_rid_of_beginning_of_line(line):
    line = line.strip()
    while True:
        try:
            int(line[0])
            line = line[1:]
        except ValueError:
            break

    # also get rid of the colon, space, and >
    line = line.replace(": > ", "")
    # get rid of the level indicators
    level = 0
    while line[0] in [".", ":"]:
        line = line[2:]
        level += 1

    return line, level

in_cfl = False
raw_violations = []

with open(log_file, "r") as in_file:
    for line in tqdm(in_file, total=total_lines):
        if "CFL CONDITION VIOLATED:" in line:
            in_cfl = True
            this_violation = dict()
        # sometimes with multiple processors we get contaminating timestep lines
        elif in_cfl and "timestep(" in line:
            continue

        elif in_cfl and "---------------------------------------------------------" in line:
            in_cfl = False
            raw_violations.append(this_violation)

        elif in_cfl:
            line, level = get_rid_of_beginning_of_line(line)
            this_violation["level"] = level

            if line == "courant cell information:" or line.startswith("CFL tolerance"):
                continue
            else:
                key, value = line.split("=")
                this_violation[key.strip()] = value.strip()

# ======================================================================================
#
# Then parse the information
#
# ======================================================================================
units_dict = {"Myr": u.Myr,
              "K": u.K,
              "cm^-3": u.cm**-3,
              "cm/s": u.cm / u.s}


def parse_simple_quan(value_string):
    value, unit = value_string.split()
    return float(value) * units_dict[unit]


def parse_pressure(pressure_string):
    value, unit_1, unit_2 = pressure_string.split()
    assert unit_1 == "ergs"
    assert unit_2 == "cm^-3"
    return float(value) * u.erg * u.cm ** -3


def parse_velocity(velocity_string):
    vx, vy, vz, unit = velocity_string.split()
    vx = float(vx) * units_dict[unit]
    vy = float(vy) * units_dict[unit]
    vz = float(vz) * units_dict[unit]
    # I calculate this with and without the absolute values.
    # the CFL reporting does not use the absolute value (which
    # is dumb), but the timestep (correctly) does
    return max(abs(vx), abs(vy), abs(vz))


violations = []
for item in raw_violations:
    new_item = item.copy()
    del new_item["v"]

    new_item["current dt"] = parse_simple_quan(item["current dt"])
    new_item["needed  dt"] = parse_simple_quan(item["needed  dt"])
    new_item["T"] = parse_simple_quan(item["T"])
    new_item["P"] = parse_pressure(item["P"])
    new_item["n"] = parse_simple_quan(item["n"])
    new_item["cs"] = parse_simple_quan(item["cs"])
    new_item["v_max"] = parse_velocity(item["v"])
    new_item["dt sound crossing"] = parse_simple_quan(item["dt sound crossing"])
    new_item["dt bulk velocity"] = parse_simple_quan(item["dt bulk velocity"])
    violations.append(new_item)

# validate the reported crossing times, which can be used to calculate the cell size
is_new_run = True # set this is the default, we'll check until its wrong. We can't just
# check it once because sometimes it is accurate, as all the velocities are positive,
# so ART's broken calculation still works
for item in violations:
    level_size_cs = item["cs"] * item["dt sound crossing"]
    level_size_bulk = item["v_max"] * item["dt bulk velocity"]
    item["cell size"] = level_size_cs
    # we'll save whether or not this agrees. In the new runs it should, but not the old.
    # The calculated v_max is always fine, just not the crossing time ART reports. In
    # addition, it uses a slightly different sound speec calculation that messes up
    # the comparison to the true dt
    if is_new_run:  # check if our assumption needs to change
        is_new_run = np.isclose(level_size_cs, level_size_bulk)

# then validate the required dt
for item in violations:
    item["total v"] = item["v_max"] + item["cs"]
    item["my needed dt"] = 0.7 * item["cell size"] / item["total v"]
    if is_new_run:
        assert np.isclose(item["my needed dt"], item["needed  dt"])


# Add cell mass. I currently don't do anything with this, but I keep it for later
for item in violations:
    item["cell mass"] =  (c.m_p * item["n"] * item["cell size"]**3).to("Msun")

# ======================================================================================
#
# plot
#
# ======================================================================================
mass_norm = colors.LogNorm(vmin=1, vmax=1E5)
mass_mappable = cm.ScalarMappable(norm=mass_norm, cmap=cmocean.cm.haline_r)
mass_colors = [mass_mappable.to_rgba(v["cell mass"].to("Msun").value) for v in violations]

fig, ax = bpl.subplots(figsize=[10, 6])
ax.scatter([v["cs"].to("km/s").value for v in violations],
           [v["v_max"].to("km/s").value for v in violations],
           c=mass_colors, alpha=1)

ax.add_labels("Sound Speed [km/s]", "Bulk Velocity [km/s]")
ax.equal_scale()

ax.fill_between([0.001, 1E6], y1=c.c.to("km/s").value, y2=1E6,
                color="0.8", zorder=0)
ax.fill_between([c.c.to("km/s").value, 1E6], y1=0, y2=1E6,
                color="0.8", zorder=0)
ax.axhline(1E3, ls="--", zorder=0)
ax.axvline(1E3, ls="--", zorder=0)
ax.plot([1E-3, 1E8], [1E-3, 1E8], c=bpl.almost_black, ls=":", lw=1, zorder=0)
ax.set_yscale("log")
ax.set_xscale("log")
ax.set_limits(10, 1E6, 10, 1E6)

# ax.add_text(x=0.2, y=1.01E3, text="Velocity Cap",
#             ha="left", va="bottom",
#             fontsize=16)
# ax.add_text(y=0.2, x=1.01E3, text="Velocity Cap",
#             ha="left", va="bottom", rotation=270,
#             fontsize=16)
ax.add_text(x=2000, y=4E5, text="Speed of Light!",
            ha="left", va="bottom",
            fontsize=16)

cb = fig.colorbar(mass_mappable, ax=ax)
cb.set_label("Cell Mass [$M_\odot$]")
fig.savefig(log_dir / "cfl_cell_speeds.png", bbox_inches="tight")
