from pathlib import Path
import sys

from matplotlib import colors
from matplotlib import cm
import cmocean

import betterplotlib as bpl
bpl.set_style()

# Get the log directory specified by the user
log_dir = Path(sys.argv[1]).absolute()

# parse the output file first
# Always use the stdout_000 file, as it has the timestep info.
log_file = log_dir / "stdout.000.log"
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
def find_first_timestep_number(log_file):
    """
    Find the timestep number the simulation started at.
    
    This is done by finding the first completed timestep and seeing
    what the completed timestep was.
    """
    with open(log_file, "r") as stdout:
        for line in stdout:
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
           line == "Error: could not complete timestep, restarting from previous timestep" or \
           line.startswith("level=0 vFac(aexpv)/vFac(aexp0)")
           # This last one is only used for the intial condition

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

def get_timestep_that_succeeded(line):
    if not is_timestep_success_line(line):
        raise ValueError("This is not a timestep success line!")
    
    # CFL violations will give -1
    if line.startswith("done with timestep"):
        return int(line.split()[3])
    else:
        return -1

# ======================================================================
# 
# then the actual reading of the log file
# 
# ======================================================================

# We'll keep track of the following things:
timestep_numbers = []  # the actual number
timestep_successes = []  # whether the timestep succeeded or had CFL error
level_dts = {l:[] for l in range(20)}

# first we need to know where we're starting
timestep_number = find_first_timestep_number(log_file)

# open the file
stdout = open(log_file, "r")

# Have state variables indicating where we are within the file
looking_for_global_timestep = True
looking_for_level_dt = False
looking_for_timestep_success = False
# We start by looking for the dts used for the first timestep
for line in stdout:
    # get rid of spaces
    line = line.strip()
    
    # we start by throwing out all lines that indicate work
    if line.startswith("> "):
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
                assert not is_timestep_success_line(line)
                continue
            except AssertionError:
                raise RuntimeError(f"Expected timestep success line, found: \n{line}")
            
        
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
            # this shouldn't happen, something is wrong
            raise RuntimeError(f"Was expecting timestep info, did not find in line:\n{line}")
    

    # Figure out whether the timestep succeeded or not
    if looking_for_timestep_success:
        if is_timestep_success_line(line):
            # figure out what to do with the counters
            new_timestep = get_timestep_that_succeeded(line)
            if new_timestep < 0:
                # CFL violation
                timestep_successes.append(False)
            else:
                # check that this matches what we think the current timestep is
                if not timestep_number == new_timestep:
                    raise ValueError(f"Something is wrong with the timestep numbers at {new_timestep}!")
                # otherwise just increment it
                timestep_number += 1
                timestep_successes.append(True)
        
            # reset state variables
            looking_for_timestep_success = False
            looking_for_global_timestep = True
        else:  # was not this line
            # check for lines that shouldn't happen
            try:
                assert not is_global_timestep_line(line)
                assert not is_level_timestep_line(line)
                assert not is_end_of_level_timesteps(line)
                continue
            except AssertionError:
                raise RuntimeError(f"Expected timestep success line, found: \n{line}")
        
    

# close the file
stdout.close()

# there may be one extra dt available compared to the timestep success if 
# the run died in the middle of a timestep (as it usually does). Check that
if len(timestep_numbers) == len(timestep_successes) + 1:
    # throw away the last timestep
    timestep_numbers = timestep_numbers[:-1]
    for level in level_dts:
        level_dts[level] = level_dts[level][:-1]
    
# then they should all be equal
assert len(timestep_numbers) == len(timestep_successes)
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

def _plot_base(ax, xs, ys, successes, color):
    """
    Underlying function to plot the timesteps as points, with lines connecting
    them. Points are filled if the timestep was successful, open if not.
    """
    ax.plot(xs, ys, lw=1, c=color, zorder=1)
    
    # then plot the successes and failures as points separately
    good_idx = [idx for idx in range(len(successes)) 
                if successes[idx] is True]
    bad_idx =  [idx for idx in range(len(successes)) 
                if successes[idx] is False]
    
    x_good = [xs[idx] for idx in good_idx]
    y_good = [ys[idx] for idx in good_idx]
    
    x_bad = [xs[idx] for idx in bad_idx]
    y_bad = [ys[idx] for idx in bad_idx]
    
    # good points are filled circles, bad ones filled with white
    ax.scatter(x_good, y_good, s=20, edgecolors=[color], c=[color], lw=1, zorder=2)
    ax.scatter(x_bad,  y_bad,  s=20, edgecolors=[color], c="w",     lw=1, zorder=2)

def plot_level(ax, timestep_successes, dts, color):
    # what we do here is go through the values and break it into a list
    # of nonzero values, to be plotted separately, without connection
    
    xs = list(range(len(timestep_successes)))
    
    # check that they're the same length
    assert len(xs) == len(dts)
    
    # we'll keep a list of the contiguous points without a break. 
    # when we hit a zero in dt (indicating that level isn't present),
    # we'll plot it and reset
    contiguous_xs = []
    contiguous_ys = []
    contiguous_successes = []
    
    for idx in range(len(xs)):
        y = dts[idx]
        if y > 0:
            contiguous_xs.append(xs[idx])
            contiguous_ys.append(y)
            contiguous_successes.append(timestep_successes[idx])
            
        else:
            if len(contiguous_xs) > 0:
                _plot_base(ax, contiguous_xs, contiguous_ys, contiguous_successes, color=color)
                contiguous_xs = []
                contiguous_ys = []
                contiguous_successes = []
            
    # if we've reached the end, we may need to plot too
    if len(contiguous_xs) > 0:
        _plot_base(ax, contiguous_xs, contiguous_ys, contiguous_successes, color=color) 

cmap = cmocean.cm.thermal
# make boundaries for colormap. They will be centered on the levels in the simulation.
boundaries = [-0.5 + level for level in range(get_max_level(level_dts)+2)]
norm = colors.BoundaryNorm(boundaries, cmap.N*0.98)
mappable = cm.ScalarMappable(cmap=cmap, norm=norm)
mappable.set_array([])

fig, ax = bpl.subplots(figsize=[2 + len(timestep_successes) / 8, 7])
ax.make_ax_dark()

for level in range(get_max_level(level_dts)+1):
    plot_level(ax, timestep_successes, level_dts[level], mappable.to_rgba(level))

ax.set_yscale("log")

# label every 5 timesteps
x_labels = [x for x in timestep_numbers
            if x % 5 == 0]
plot_xs = list(range(len(timestep_successes)))
for label in x_labels:
    # find the right x_value to place this at
    for idx in range(len(timestep_successes)):
        if timestep_numbers[idx] == label and timestep_successes[idx]:
            x = plot_xs[idx]
            y = level_dts[0][idx]
            
            ax.add_text(x, y*1.2, label, va="bottom", ha="center", fontsize=12)
            ax.plot([x, x], [y, y*1.2], c=bpl.almost_black, lw=1, zorder=0)
            
# Then add the colorbar
cbar = fig.colorbar(mappable, ax=ax, ticks=list(level_dts.keys()))
cbar.set_label("Level")

ax.set_limits(-2, max(plot_xs)+2)
ax.add_labels("Timesteps Attempted This Run", "dt [years]")
fig.savefig(log_dir / "timestep_history.png", bbox_inches="tight")