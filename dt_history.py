from pathlib import Path
import sys

import numpy as np
from astropy import units as u
from matplotlib import colors
from matplotlib import cm
import colorcet as cc
from tqdm import tqdm

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
# To do this I'll go through the file. We start with a timestep. We'll get the 
# timesteps on all the levels. This ends with "Spatial Resolution" appears.
# Then the next thing is to figure out whether the timestep was successfull or 
# not. We'll either get "Error: could not complete timestep" or 
# "done with timestep". Then we know which timestep comes next.
# One complication: we don't know what the first timestep is. We'll have to go 
# through the file to start with to find that info.
# ======================================================================
# 
# Convenience functions
# 
# ======================================================================
def strip_mpi_rank_from_line(line):
    # get rid of the beginning thing. First we get rid of the MPI rank number
    line = line.strip()
    while True:
        try:
            int(line[0])
            line = line[1:]
        except ValueError:
            break

    # also get rid of the colon and space
    line = line[1:]
    return line.strip()

def find_first_timestep_number(log_file):
    """
    Find the timestep number the simulation started at.
    
    This is done by finding the first completed timestep and seeing
    what the completed timestep was.
    """
    with open(log_file, "r") as stdout:
        for line in stdout:
            line = strip_mpi_rank_from_line(line)
            if line.startswith("done with timestep"):
                return int(line.split()[3])
            
def is_global_timestep_line(line):
    return line.startswith("chose") and line.endswith("Myr as our next global time-step")

def is_level_timestep_line(line):
    return line.startswith("level") and \
           ", dt = " in line and \
           " Myr, time-ref = " in line and \
           ", global time-ref = " in line

def is_end_of_level_timesteps(line):
    return line.startswith("Spatial resolution: ")

def is_timestep_success_line(line):
    return line.startswith("done with timestep") or \
           line.strip().endswith("CFL CONDITION VIOLATED:") or \
           line.startswith("level=0 vFac(aexpv)/vFac(aexp0)")
           # This last one is only used for the intial condition

def is_end_of_cfl_violation_info(line):
    return "---------------------------------------------------------" in line

def get_global_timestep(line):
    # Returns a value in years, rather than Myr as is in the file
    if not is_global_timestep_line(line):
        raise ValueError("This is not a global timestep line!")
    # convert from Myr to yr
    return float(line.split()[1]) * 1E6

def get_level_timestep(line):
    # Returns a value in years, rather than Myr as is in the file
    if not is_level_timestep_line(line):
        raise ValueError("This is not a level timestep line!")
    
    # get the level and the dt
    split = line.split()
    level = int(split[1].replace(",", ""))
    
    # convert from Myr to yr
    dt = float(split[4]) * 1E6
    
    return level, dt

def did_timestep_succeed(line):
    if not is_timestep_success_line(line):
        raise ValueError("This is not a timestep success line!")
    
    # The first item returned will be the success or failure. If there was a success, it
    # will be the timestep number. If not, it will be the level of the CFL violation.
    if line.startswith("done with timestep"):
        return True
    elif line.startswith("level=0 vFac(aexpv)/vFac(aexp0)"):
        return True
    else:
        return False

        # idx_start = line.find("[") + 1  # add 1 to get to the number in brackets
        # idx_end = line.find("]")
        # return False, int(line[idx_start:idx_end])

def get_successful_timestep_number(line):
    if not is_timestep_success_line(line):
        raise ValueError("This is not a timestep success line!")

    if line.startswith("level=0 vFac(aexpv)/vFac(aexp0)"):
        return 0  # IC doesn't really count
    elif line.startswith("done with timestep"):
        return int(line.split()[3])
    else:
        raise ValueError("Line not recognized!")

def get_rid_of_level_indicators(line):
    line = line.replace("> ", "")
    # get rid of the level indicators
    level = 0
    while line[0] in [".", ":"]:
        line = line[2:]
        level += 1

    return line, level

# Functions for parsing the CFL violation things
units_dict = {"Myr": u.Myr,
              "K": u.K,
              "cm^-3": u.cm**-3,
              "cm/s": u.cm / u.s}

def parse_simple_quan(value_string):
    value, unit = value_string.split()
    return float(value) * units_dict[unit]

def parse_velocity(velocity_string):
    vx, vy, vz, unit = velocity_string.split()
    vx = float(vx) * units_dict[unit]
    vy = float(vy) * units_dict[unit]
    vz = float(vz) * units_dict[unit]
    # the CFL reporting does not use the absolute value (which
    # is dumb), but the timestep (correctly) does
    return max(abs(vx), abs(vy), abs(vz))

# ======================================================================
# 
# then the actual reading of the log file
# 
# ======================================================================

# We'll keep track of the following things:
timestep_numbers = []  # the actual number
# the level on which CFL violations happened in a given level. Will be -99 if there was
# no CFL violation in that timestep
timestep_successes = []
cfl_violation_levels = []
cfl_violation_reasons = []

level_dts = {l:[] for l in range(20)}

# count number of lines in the file
with open(log_file, "r") as in_file:
    for i, l in enumerate(in_file):
        pass
total_lines = i + 1

# first we need to know where we're starting
timestep_number = find_first_timestep_number(log_file)
# open the file
stdout = open(log_file, "r")

# Have state variables indicating where we are within the file
looking_for_global_timestep = True
looking_for_level_dt = False
looking_for_timestep_success = False
inside_cfl_info = False
# We start by looking for the dts used for the first timestep
for line in tqdm(stdout, total=total_lines):
    line = strip_mpi_rank_from_line(line)
    
    # we start by throwing out all lines that indicate work. We have to be careful
    # since the CFL violation lines also start with "> "
    if line.startswith("> ") and "timestep(" in line:
        continue

    # get the global timestep if we need to
    if looking_for_global_timestep:
        if is_global_timestep_line(line):
            # start a clean dictionary for this new timestep
            this_level_dt = dict()
            this_level_dt[0] = get_global_timestep(line)
            # then set the state variables
            looking_for_global_timestep = False
            looking_for_level_dt = True
            continue  # go to next line
        else:    
            # check for lines that shouldn't happen
            try:
                assert not is_level_timestep_line(line)
                assert not is_end_of_level_timesteps(line)
                # timestep success can happen, since each MPI processor indicates that
                # the timestep finished
                # assert not is_timestep_success_line(line)
                continue
            except AssertionError:
                raise RuntimeError(f"Expected global timestep line, found: \n{line}")
            
        
    # if we need to look for the level info, do that
    if looking_for_level_dt:
        # look for lines with actual timestep info
        if is_level_timestep_line(line):
            level, dt = get_level_timestep(line)
            this_level_dt[level] = dt
        # also look for the end of the timestep info
        elif is_end_of_level_timesteps(line):
            # store all of this info
            timestep_numbers.append(timestep_number)

            for level in level_dts:
                # If this level exists now, append the timestep used
                if level in this_level_dt:
                    level_dts[level].append(this_level_dt[level])
                # otherwise append zero
                else:
                    level_dts[level].append(0)    

            # then reset the state variables
            looking_for_level_dt = False
            looking_for_timestep_success = True
            
            # then go to the next line
            continue
            
        else:
            try:
                assert not is_global_timestep_line(line)
                assert not is_end_of_level_timesteps(line)
                assert not is_timestep_success_line(line)
                continue
            except AssertionError:
                # this shouldn't happen, something is wrong
                raise RuntimeError(f"Was expecting timestep info, did not find in line:\n{line}")
    

    # Figure out whether the timestep succeeded or not
    if looking_for_timestep_success:
        if is_timestep_success_line(line):
            # figure out what to do with the counters
            if did_timestep_succeed(line):
                this_timestep_number = get_successful_timestep_number(line)
                # the IC is defined as a success, but we don't do anything with it.
                if this_timestep_number != 0:
                    # check that this matches what we think the current timestep is
                    if not timestep_number == this_timestep_number:
                        raise ValueError(f"Something is wrong with the timestep numbers at {this_timestep_number}!")
                    # otherwise just increment it
                    timestep_number += 1
                    # and att this info to our lists
                    timestep_successes.append(True)
                    cfl_violation_levels.append(None)
                    cfl_violation_reasons.append("")

                # reset state variables
                looking_for_timestep_success = False
                looking_for_global_timestep = True
            else:  # CFL violation
                # set up CFL info
                this_violation = dict()
                # reset state variables, we'll find the info later
                looking_for_timestep_success = False
                inside_cfl_info = True
            continue  # go to next line
        else:  # was not this line
            # check for lines that shouldn't happen
            try:
                assert not is_global_timestep_line(line)
                assert not is_level_timestep_line(line)
                assert not is_end_of_level_timesteps(line)
                continue
            except AssertionError:
                raise RuntimeError(f"Expected timestep success line, found: \n{line}")

    # if we need to get the CFL information
    if inside_cfl_info:
        # first check if we're done, as this line breaks the parsers
        if is_end_of_cfl_violation_info(line):
            # First figure out what caused the violation. The simplest way to do this
            # would be to just see which crossing times are smaller than the current
            # timestep. But because the needed velocity is bulk + c_s, typically it
            # will be hard for one to get there all by itself, unless it's drastiaclly
            # higher than the other. I tried making this plot, and most of the time
            # neither got there by myself. So instead we just determine which one is
            # larger
            c_s = parse_simple_quan(this_violation["cs"])
            bulk = parse_velocity(this_violation["v"])
            # see if one is significantly bigger than the other
            threshold = 1.25
            if c_s > threshold * bulk:
                this_reason = "sound"
            elif bulk > threshold * c_s:
                this_reason = "bulk"
            else:
                this_reason = "sound bulk"

            # append the info we collected to the lists
            timestep_successes.append(False)
            cfl_violation_levels.append(this_violation["level"])
            cfl_violation_reasons.append(this_reason)

            # reset state variables
            inside_cfl_info = False
            looking_for_global_timestep = True
            continue # go to next timestep

        # otherwise parse what we have
        line, level = get_rid_of_level_indicators(line)
        this_violation["level"] = level

        if line == "courant cell information:" or line.startswith("CFL tolerance"):
            # dont' need this info
            continue
        else:
            key, value = line.split("=")
            this_violation[key.strip()] = value.strip()


# close the file
stdout.close()

# there may be one extra dt available compared to the CFL levels if
# the run died in the middle of a timestep (as it usually does). Check that
if len(timestep_numbers) == len(timestep_successes) + 1:
    # throw away the last timestep
    timestep_numbers = timestep_numbers[:-1]
    for level in level_dts:
        level_dts[level] = level_dts[level][:-1]
    
# then they should all be equal
assert len(timestep_numbers) == len(timestep_successes)
assert len(timestep_numbers) == len(cfl_violation_levels)
assert len(timestep_numbers) == len(cfl_violation_reasons)
for level in level_dts:
    assert len(level_dts[level]) == len(timestep_numbers)

# ======================================================================
# 
# plotting
# 
# ======================================================================
def get_max_level(level_dts):
    # iterate backwards, and find the first level with any data
    for level in range(max(level_dts.keys()), 0, -1):
        if max(level_dts[level]) > 0:
            return level

cmap = cc.m_CET_R1
# make boundaries for colormap. They will be centered on the levels in the simulation.
boundaries = [-0.5 + level for level in range(get_max_level(level_dts)+2)]
norm = colors.BoundaryNorm(boundaries, cmap.N)
mappable = cm.ScalarMappable(cmap=cmap, norm=norm)
mappable.set_array([])

def _plot_base(ax, xs, ys, successes, cfl_violation_level, cfl_violation_reason, level):
    """
    Underlying function to plot the timesteps as points, with lines connecting
    them. Points are filled if the timestep was successful, open if not. There will also
    be markers showing which level was responsible and which velocity (sound or bulk)
    was the biggest contributor
    """
    color = mappable.to_rgba(level)

    ax.plot(xs, ys, lw=1, c=color, zorder=1)
    
    # then plot the successes and failures as points separately
    good_idx = [idx for idx in range(len(cfl_violation_level))
                if successes[idx]]
    bad_idx =  [idx for idx in range(len(cfl_violation_level))
                if (not successes[idx]) and cfl_violation_level[idx] != level]
    cfl_idx =  [idx for idx in range(len(cfl_violation_level))
                if (not successes[idx]) and cfl_violation_level[idx] == level]
    cfl_bulk = [idx for idx in range(len(cfl_violation_level))
                if idx in cfl_idx and "bulk" in cfl_violation_reason[idx]]
    cfl_sound = [idx for idx in range(len(cfl_violation_level))
                 if idx in cfl_idx and "sound" in cfl_violation_reason[idx]]
    
    x_good = [xs[idx] for idx in good_idx]
    y_good = [ys[idx] for idx in good_idx]
    
    x_bad = [xs[idx] for idx in bad_idx]
    y_bad = [ys[idx] for idx in bad_idx]

    x_cfl = [xs[idx] for idx in cfl_idx]
    y_cfl = [ys[idx] for idx in cfl_idx]

    x_cfl_bulk = [xs[idx] for idx in cfl_bulk]
    y_cfl_bulk = [ys[idx] for idx in cfl_bulk]

    x_cfl_sound = [xs[idx] for idx in cfl_sound]
    y_cfl_sound = [ys[idx] for idx in cfl_sound]
    
    # good points are filled circles, bad ones filled with white.
    # Make a list of colors (that are all the same) so scatter doesn't get
    # confused by the rgb tuples
    colors_good = [color for _ in x_good]
    colors_bad = [color for _ in x_bad]
    colors_cfl = [color for _ in x_cfl]
    colors_cfl_bulk = [color for _ in x_cfl_bulk]
    colors_cfl_sound = [color for _ in x_cfl_sound]

    common_circle = {"alpha":1.0}
    # the CFL indicators should go behind the bad markers
    common_x = {"lw":2, "zorder": 2, "alpha":1.0}
    size_cfl = 150  # leave this separate since the X needs to be smaller than +
    ax.scatter(x_good, y_good, edgecolors=colors_good, c=colors_good,
               s=30, zorder=3, lw=1, **common_circle)
    ax.scatter(x_bad,  y_bad,  edgecolors=colors_bad,  c=bpl.light_grey,
               s=10, zorder=3, lw=1, **common_circle)
    ax.scatter(x_cfl,  y_cfl,  edgecolors=colors_cfl,  c="w",
               s=30, zorder=4, lw=2, **common_circle)
    # My mnemonic is "P" for plus and for pressure (which relates to c_s)
    ax.scatter(x_cfl_sound,  y_cfl_sound,  c=colors_cfl_sound, marker="+",
               s=size_cfl, **common_x)
    ax.scatter(x_cfl_bulk, y_cfl_bulk,   c=colors_cfl_bulk, marker="x",
               s=size_cfl*0.5, **common_x)

def plot_level(ax, timestep_successes, cfl_violation_level, cfl_violation_reasons,
               dts, level):
    # what we do here is go through the dts and break it into contiguous segments
    # where that level existed
    
    xs = list(range(len(dts)))
    
    # check that they're the same length
    assert len(xs) == len(dts)
    
    # we'll keep a list of the contiguous points without a break. 
    # when we hit a zero in dt (indicating that level isn't present),
    # we'll plot it and reset
    contiguous_xs = []
    contiguous_ys = []
    contiguous_success = []
    contiguous_cfl_level = []
    contiguous_cfl_reason = []
    
    for idx in range(len(xs)):
        y = dts[idx]
        if y > 0:
            contiguous_xs.append(xs[idx])
            contiguous_ys.append(y)
            contiguous_success.append(timestep_successes[idx])
            contiguous_cfl_level.append(cfl_violation_level[idx])
            contiguous_cfl_reason.append(cfl_violation_reasons[idx])
            
        else:
            if len(contiguous_xs) > 0:
                _plot_base(ax, contiguous_xs, contiguous_ys, contiguous_success,
                           contiguous_cfl_level, contiguous_cfl_reason, level)
                contiguous_xs = []
                contiguous_ys = []
                contiguous_success = []
                contiguous_cfl_level = []
                contiguous_cfl_reason = []
            
    # if we've reached the end, we may need to plot too
    if len(contiguous_xs) > 0:
        _plot_base(ax, contiguous_xs, contiguous_ys, contiguous_success,
                   contiguous_cfl_level, contiguous_cfl_reason, level)


# ======================================================================
#
# actual figure
#
# ======================================================================
fig, ax = bpl.subplots(figsize=[3 + len(timestep_numbers) / 8, 7])
ax.make_ax_dark()

for level in range(get_max_level(level_dts)+1):
    plot_level(ax, timestep_successes, cfl_violation_levels, cfl_violation_reasons,
               level_dts[level], level)

ax.set_yscale("log")

# label every 5 timesteps
x_labels = [x for x in timestep_numbers
            if x % 5 == 0]
plot_xs = list(range(len(timestep_numbers)))
for label in x_labels:
    # find the right x_value to place this at
    for idx in range(len(timestep_numbers)):
        if timestep_numbers[idx] == label and timestep_successes[idx]:
            x = plot_xs[idx]
            y = level_dts[0][idx]
            
            ax.add_text(x, y*1.2, label, va="bottom", ha="center", fontsize=12)
            ax.plot([x, x], [y, y*1.2], c=mappable.to_rgba(0), lw=1, zorder=1)
            
# Then add the colorbar
cbar = fig.colorbar(mappable, ax=ax, ticks=list(level_dts.keys()), pad=0)
cbar.set_label("Level")
cbar.ax.invert_yaxis()

ax.set_limits(-2, max(plot_xs)+2)
ax.add_labels("Timesteps Attempted This Run", "dt [years]")
fig.savefig(log_dir / "timestep_history.png", bbox_inches="tight")
