"""
plot_cimf_evolution.py

Plot the evolution of the CIMF for one simulation.

This takes two parameters:
- the location to save the plot to. The name of this plot should tell the directory
  where the sim lives
"""
import sys
from pathlib import Path

import numpy as np
import yt
import betterplotlib as bpl

from utils import load_galaxies, plot_utils, run_attributes
from analysis_functions import cimf

bpl.set_style()
yt.funcs.mylog.setLevel(50)  # ignore yt's output

plot = Path(sys.argv[1])
# then get the simulation directory. This is based on reverse engineering the function
# in the makefile to get the plot name.
sim_dir = Path(
    plot.name.replace("cimf_zevolution_", "").replace("__", "/").replace(".pdf", "")
)

# ======================================================================================
#
# redshift evolution of a single run
#
# ======================================================================================
# This is pretty simple. We load the simulation outputs at a range of redshift and
# plot the observable mass function at that redshift.
fig, ax = bpl.subplots()
# colors = ["#E996CA", "#9E729F", "#534D74", "#082949"]
colors = [
    run_attributes.h(*hsv)
    for hsv in [
        (0.95, 0.20, 0.95),
        (0.85, 0.30, 0.85),
        (0.75, 0.40, 0.72),
        (0.65, 0.50, 0.60),
        (0.65, 0.00, 0.15),  # consistent with bpl.almost_black
    ]
]

# first figure out which redshifts to show. I'll show z=0, the last output, and 3 others
sim_last = load_galaxies.get_simulations_last([sim_dir])
if len(sim_last) == 0:  # when running on shangrila, Stampede2 sims are not there
    fig.savefig(plot)
    exit()

assert len(sim_last) == 1  # check that there aren't 2 for some reason
sim_last = sim_last[0]
last_z = sim_last.z
# then figure out other zs to show
if last_z < 2:
    zs = [6, 5, 4, 3, last_z]
elif last_z < 3:
    zs = [7, 6, 5, 4, last_z]
elif last_z < 4:
    zs = [7, 6, 5, last_z]
else:
    zs = [last_z]

max_yvalue = 0
for z, c in zip(zs, colors):
    sim = load_galaxies.get_simulations_same_scale([sim_dir], z)
    if len(sim) == 0:  # no sim at this redshift found
        continue
    assert len(sim) == 1
    plot_masses, dn_dlogM = cimf.cimf(sim[0], "current", np.inf, np.inf)

    # determine whether to just show integers in the legend
    if np.isclose(int(z), z, atol=0.05):
        label = f"z={z:.0f}"
    else:
        label = f"z={z:.1f}"

    plot_utils.plot_line_with_cut(
        ax,
        plot_masses,
        dn_dlogM,
        cut_x=sim_last.unreliable_mass,
        c=c,
        label=label,
        ls2="--",
        zorder=5,
    )

    if np.max(dn_dlogM) > max_yvalue:
        max_yvalue = np.max(dn_dlogM)

# then show analytical disruption to z=0
plot_masses, dn_dlogM = cimf.cimf(sim_last, "evolved", np.inf, np.inf)
ax.plot(plot_masses, dn_dlogM, c=bpl.almost_black, ls=":", label=f"z=0", zorder=100)

# and the observed GC mass function
m_mw_gc, y_mw_gc = cimf.harric_gc_mass_function()
# ax.plot(m_mw_gc, y_mw_gc, c="0.5", ls=":", label="MW GCs", zorder=0)
ax.fill_between(x=m_mw_gc, y1=0, y2=y_mw_gc, color="0.5", alpha=0.5, zorder=1)

# format axis
plot_utils.add_legend(ax, loc=1, fontsize=18)
ax.set_yscale("log")
ax.set_xscale("log")
# put 1 cluster as the y lower limit of the plot
ax.set_limits(
    1e3, 1e7, 0.8 / (0.16 * np.log(10)), 10 ** (np.ceil(np.log10(max_yvalue)))
)
ax.add_labels("$f_b$M [$M_\odot$]", "dN/dlogM")
fig.savefig(plot)
