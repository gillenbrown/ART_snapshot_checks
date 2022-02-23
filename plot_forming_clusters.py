"""
plot_forming_clusters.py

Creates several plots examining the properties of clusters as they form. This parses
the scripts previously generated

This script takes the list of all directories containing the checks
"""

import sys
from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt
from astropy import table
import betterplotlib as bpl

from utils import load_galaxies, plot_utils, run_attributes
from analysis_functions import age_spreads

bpl.set_style()

sentinel_file = Path(sys.argv[1])
plot_dir = sentinel_file.parent

# ======================================================================================
#
# read the files
#
# ======================================================================================
def read_file(file_loc):
    # check that this table actually has data
    with open(file_loc, "r") as in_file:
        good = False
        for line in in_file:
            if not line.startswith("#"):
                good = True
                break
    if not good:
        return None
    # If the table does have data, read it with an astropy table
    cat = table.Table.read(file_loc, format="ascii")
    # columns need to be renamed
    cat["col1"].name = "ID"
    cat["col2"].name = "Age_Myr"
    cat["col3"].name = "number_density"
    cat["col4"].name = "H2_mass"
    cat["col5"].name = "H2_frac"
    cat["col6"].name = "temp"

    return cat


tables = dict()
for checks_dir in sys.argv[2:]:
    this_cats = []
    for summary_path in [
        str(f)
        for f in Path(checks_dir).iterdir()
        if f.name.startswith("forming_clusters")
    ]:
        this_cat = read_file(summary_path)
        if this_cat is not None:
            this_cats.append(this_cat)

    # then stack all tables together
    if len(this_cats) > 0:
        tables[Path(checks_dir).parent] = table.vstack(this_cats)

# ======================================================================================
#
# Add info to the table based on the last output
#
# ======================================================================================
for run_dir, this_cat in tables.items():
    # add dummy columns
    this_cat["matched"] = False
    this_cat["duration"] = -99.9
    this_cat["M_initial"] = -99.9

    # then load the simulation and get the relevant data
    sim = load_galaxies.get_simulations_last([run_dir])[0]
    ids = sim.func_all_galaxies(lambda g: g[("STAR", "PID")].to("").value)
    durations = sim.func_all_galaxies(lambda g: age_spreads.duration(g).to("Myr").value)
    m_initial = sim.func_all_galaxies(
        lambda g: g[("STAR", "INITIAL_MASS")].to("Msun").value
    )

    # then go through and match the clusters from runtime to those in the output files
    # To do this, we rely on the IDs in both the table and the sim output being sorted.
    # Then we can simply go through and match without double iteration
    this_cat.sort("ID")
    assert np.array_equal(np.sort(this_cat["ID"]), np.array(this_cat["ID"]))
    # now sort the simulation data
    sim_sort_idx = np.argsort(ids)
    ids = ids[sim_sort_idx]
    durations = durations[sim_sort_idx]
    m_initial = m_initial[sim_sort_idx]
    assert np.array_equal(np.sort(ids), ids)

    # Then go through and work on the matching
    idx_table = 0
    idx_sim = 0
    while idx_table < len(this_cat) and idx_sim < len(ids):
        id_tab = this_cat["ID"][idx_table]
        id_sim = ids[idx_sim]
        # check if we have a match
        if id_tab == id_sim:
            row = this_cat[idx_table]
            row["matched"] = True
            row["duration"] = durations[idx_sim]
            row["M_initial"] = m_initial[idx_sim]
            idx_table += 1
            idx_sim += 1
        # If we got here, we did not match. Increase the counter of the index that
        # points to a smaller value. This way we keep the indices roughly together and
        # avoid missing a match
        elif id_tab > id_sim:
            idx_sim += 1
        else:
            idx_table += 1

    # add derived properties: whether the cluster failed formation, and whether or not
    # we recorded its properties while it was still forming.
    this_cat["failed"] = this_cat["duration"] > 14
    this_cat["still_forming"] = this_cat["Age_Myr"] < this_cat["duration"]
    # restrict to only clusters that matched and are still forming
    this_cat = this_cat[this_cat["matched"]]
    this_cat = this_cat[this_cat["still_forming"]]

    # Delete unnecessary columns
    del this_cat["ID"]
    del this_cat["matched"]
    del this_cat["still_forming"]
    del this_cat["duration"]


# ======================================================================================
#
# Functions for analysis
#
# ======================================================================================
class YAxisProp(object):
    def __init__(self, table_name, label, y_min, y_max, log):
        self.table_name = table_name
        self.label = label
        self.y_min = y_min
        self.y_max = y_max
        if log:
            self.y_scale = "log"
        else:
            self.y_scale = "linear"

    def plot_shaded_region_age_vs_prop(self, ax, data_table, color, idxs):
        plot_utils.shaded_region(
            ax,
            xs=data_table["Age_Myr"][idxs],
            ys=data_table[self.table_name][idxs],
            color=color,
            p_lo=10,
            p_hi=90,
            dx=1,
            cut_x=np.inf,
            log_x=False,
            label=self.label,
        )
        self.format(ax)

    def format(self, ax):
        ax.add_labels("Cluster Age [Myr]", self.label)
        ax.set_yscale(self.y_scale)
        ax.set_limits(0, 15, self.y_min, self.y_max)


y_attributes = [
    YAxisProp("number_density", "Number Density", 1, 1e5, True),
    YAxisProp("H2_mass", "H$_2$ Mass [M$_\odot$]", 1, 1e5, True),
    YAxisProp("H2_frac", "Fraction of Mass in H$_2$", 0, 1, False),
    YAxisProp("temp", "Temperature [K]", 1, 1e8, True),
]


def plot_single_sim_massive_split_failure(run_dir):
    this_table = tables[run_dir]
    for prop in y_attributes:
        fig, axs = bpl.subplots(figsize=[14, 7], ncols=2)

        idxs_left = (this_table["M_initial"] > 1e5) & ~this_table["failed"]
        idxs_right = (this_table["M_initial"] > 1e5) & this_table["failed"]

        prop.plot_shaded_region_age_vs_prop(
            axs[0], this_table, run_attributes.colors[run_dir], idxs_left
        )
        prop.plot_shaded_region_age_vs_prop(
            axs[1], this_table, run_attributes.colors[run_dir], idxs_right
        )

        axs[0].easy_add_text("M > $10^5 M_\odot$\nSuccessfull", "upper right")
        axs[1].easy_add_text("M > $10^5 M_\odot$\nFailures", "upper right")

        sim_name = plot_utils.get_sim_dirname(run_dir)
        plot_name = f"forming_clusters_age_{prop.table_name}_{sim_name}.pdf"
        fig.savefig(plot_dir / plot_name)
        # close figures to save memory
        plt.close(fig)


# ======================================================================================
#
# Do analysis
#
# ======================================================================================
for run_dir in tables.keys():
    plot_single_sim_massive_split_failure(run_dir)

# figure out which groupings are needed when plotting
# sim_groups = load_galaxies.get_plot_names_dirs(tables.keys())
# for group in sim_groups:
#     pass

sentinel_file.touch()
