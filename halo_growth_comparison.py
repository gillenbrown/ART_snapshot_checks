import sys
from pathlib import Path

import numpy as np
import yt
from astropy import cosmology
from astropy import units as u
from scipy import interpolate

import matplotlib.pyplot as plt
from matplotlib import gridspec
import betterplotlib as bpl

from utils import plot_utils

bpl.set_style()
yt.funcs.mylog.setLevel(50)  # ignore yt's output

plot_name = Path(sys.argv[1])
# ==============================================================================
#
# We need the time, so we need to set up cosmology info
#
# ==============================================================================
# I looked at the MUSIC config file to get the cosmology information
cosmo = {
    "R&J": cosmology.FlatLambdaCDM(H0=68, Om0=0.310, Tcmb0=2.725),
    "T&L": cosmology.FlatLambdaCDM(H0=71, Om0=0.266, Tcmb0=2.725),
}

# define some helper functions to handle ages
def a_to_z(a):
    if a <= 0:
        return np.inf
    else:
        return (1.0 / a) - 1


def z_to_age_gyr(z, ic):
    return cosmo[ic].age(z).to("Gyr").value


def age_gyr_to_z(age, ic):
    age *= u.Gyr
    if age >= 0.95 * cosmo[ic].age(0):
        return -1 * (
            (age / cosmo[ic].age(0)).value - 1
        )  # dumb guess that will be used in minimize
    elif age <= 0:
        return np.inf
    elif age < 1 * u.Myr:
        return 608 / age.to(u.Myr).value  # scale to real values
    return cosmology.z_at_value(cosmo[ic].age, age)


# ==============================================================================
#
# Class to handle the halo catalogs
#
# ==============================================================================
class Halo(object):
    def __init__(self, mergers_file_loc, growth_file_loc, name):
        self.name = name
        self.ic = self.name.split()[0]

        self.h = cosmo[self.ic].h

        self._init_mergers(mergers_file_loc)
        self._init_growth(growth_file_loc)

        self.mergers_z = [a_to_z(a) for a in self.mergers_scale]
        self.growth_z = [a_to_z(a) for a in self.growth_scale]

        self.mergers_age = [z_to_age_gyr(z, self.ic) for z in self.mergers_z]
        self.growth_age = [z_to_age_gyr(z, self.ic) for z in self.growth_z]

        # then create an interpolation object, which will be used for calculating the ratios
        # of the growth history
        if len(self.growth_age) == 0:
            raise ValueError("bad halo growth")
        elif len(self.growth_age) == 1:
            self.growth_age_interp = None
        else:
            self.growth_age_interp = interpolate.interp1d(
                x=self.growth_age, y=self.growth_mass, kind="linear"
            )

    def _init_mergers(self, mergers_file_loc):
        self.mergers_scale = []
        self.mergers_sat_mass = []
        self.mergers_central_mass = []
        self.mergers_mass_ratio = []
        self.mergers_scale_final = []

        with open(mergers_file_loc, "r") as in_file:
            for line in in_file:
                if line.startswith("#"):
                    continue
                # otherwise, columns are scale factor, satellite mass, central mass,
                # then ids
                split_line = line.split()
                a, m_sat, m_cen, id_sat, id_cen = split_line[0:5]
                a_f, m_sat_f, m_cen_f, id_sat_f, id_cen_f = split_line[-5:]
                self.mergers_scale.append(float(a))
                self.mergers_sat_mass.append(float(m_sat) / self.h)
                self.mergers_central_mass.append(float(m_cen) / self.h)
                self.mergers_mass_ratio.append(float(m_cen) / float(m_sat))
                self.mergers_scale_final.append(float(a_f))

    def _init_growth(self, growth_file_loc):
        self.growth_scale = []
        self.growth_mass = []
        with open(growth_file_loc, "r") as in_file:
            for line in in_file:
                if line.startswith("#"):
                    continue
                # otherwise, the data is scale factor, halo mass, then id
                a, m, id = line.split()
                # put h back in
                m = float(m) / self.h
                a = float(a)

                self.growth_scale.append(a)
                self.growth_mass.append(m)


# ==============================================================================
#
# Then get the two halos in each of these simulations
#
# ==============================================================================
class Simulation(object):
    def __init__(self, output_dir):
        self.dirname = Path(output_dir).resolve() / "checks"
        mergers_1 = self.dirname / "mergers_1.txt"
        mergers_2 = self.dirname / "mergers_2.txt"

        growth_1 = self.dirname / "growth_1.txt"
        growth_2 = self.dirname / "growth_2.txt"

        try:
            name_base = plot_utils.names[Path(output_dir)]
        except KeyError:
            raise ValueError("Sim name not found in plot_utils.py", output_dir)

        self.g1 = Halo(mergers_1, growth_1, name_base)
        self.g2 = Halo(mergers_2, growth_2, name_base)


# fmt: off
tl_nbody = Simulation("/u/home/gillenb/art_runs/runs/pleiades/nbody/new_ic_trim_25mpc/root_08/run/outputs/vel_offset/")
tl_sfe001_hn20 = Simulation("/u/home/gillenb/art_runs/runs/stampede2/production/tl_sfe001_hn20/run/")
tl_sfe010_hn20 = Simulation("/u/home/gillenb/art_runs/runs/stampede2/production/tl_sfe010_hn20/run/")
tl_sfe100_hn00 = Simulation("/u/home/gillenb/art_runs/runs/stampede2/production/tl_sfe100_hn00/run/")
tl_sfe100_hn05 = Simulation("/u/home/gillenb/art_runs/runs/stampede2/production/tl_sfe100_hn05/run/")
tl_sfe100_hn20 = Simulation("/u/home/gillenb/art_runs/runs/stampede2/production/tl_sfe100_hn20/run/")
rj_nbody_92 = Simulation("/u/home/gillenb/art_runs/runs/stampede2/rj_nbody/original_92.48mpc_level07/run/")
rj_nbody_46 = Simulation("/u/home/gillenb/art_runs/runs/stampede2/rj_nbody/hybrid_46.24mpc_level08/run/")
rj_nbody_23 = Simulation("/u/home/gillenb/art_runs/runs/stampede2/rj_nbody/hybrid_23.12mpc_level08/run/")
# fmt: on

# ==============================================================================
#
# Then we plot! First a bunch of helper functions
#
# ==============================================================================
def format_main_ax(ax, ic):
    ax.add_labels("", "Mass [$M_\odot$]")
    ax.set_limits(0, 14, 1e9, 3e12)
    ax.set_xticklabels([])
    ax.tick_params(axis="x", direction="in")

    ax.twin_axis(
        "x",
        [0, 0.5, 1, 2, 3, 5],
        label="Redshift",
        old_to_new_func=lambda age: age_gyr_to_z(age, ic),
    )
    ax.set_yscale("log")
    ax.yaxis.set_ticks_position("both")

    ax.yaxis.set_ticks([1e10, 1e11, 1e12])
    ax.yaxis.set_ticklabels(["$10^{10}$", "$10^{11}$", "$10^{12}$"])


def format_diff_ax(ax):
    ax.add_labels("Cosmic Time [Gyr]", "$\Delta$ M [dex]")
    ax.set_limits(0, 14, -0.5, 0.5)
    ax.yaxis.set_ticks_position("both")


def plot_mass_ratio(ax, halo1, halo2, color):
    age_start = max(halo1.growth_age[-1], halo2.growth_age[-1])
    age_end = min(halo1.growth_age[0], halo2.growth_age[0])

    ages = np.linspace(age_start, age_end, 1000)

    if halo1.growth_age_interp is not None and halo2.growth_age_interp is not None:
        ratios = [halo2.growth_age_interp(a) / halo1.growth_age_interp(a) for a in ages]
    else:
        ages = [age_start]
        ratios = [halo2.growth_mass[0] / halo1.growth_mass[0]]

    diffs_dex = [np.log10(r) for r in ratios]

    ax.plot(ages, diffs_dex, lw=3, c=color)


def plot_merger_marker(ax, x, y, color, marker, zorder=1):
    msize = 300
    if marker in ["+", "x", "o", "v"]:
        linewidth = 3
        if marker == "+":
            msize *= np.sqrt(2)
        if marker in ["x", "+"]:
            facecolor = color
        else:
            facecolor = "none"
    else:
        linewidth = 0
        facecolor = color
        color = bpl.almost_black

        if marker == "s":
            msize /= 3
        elif marker == "P":
            msize /= 1.5
        elif marker == ".":
            msize *= 2
        elif marker == ">":
            msize /= 1.5

    ax.scatter(
        x,
        y,
        color=color,
        marker=marker,
        facecolor=facecolor,
        linewidth=linewidth,
        alpha=1.0,
        s=msize,
        zorder=zorder,
    )


def plot_merger(ax, halo, color, marker, ic):
    ax.plot(
        halo.growth_age,
        halo.growth_mass,
        lw=3,
        label=halo.name,
        color=color,
        zorder=100,
    )

    for a_i, a_f, m_sat, m_ratio in zip(
        halo.mergers_scale,
        halo.mergers_scale_final,
        halo.mergers_sat_mass,
        halo.mergers_mass_ratio,
    ):
        if merger_ratio_cut(m_ratio, a_i) and m_sat > 1e9:
            age = z_to_age_gyr(a_to_z(a_i), ic)
            plot_merger_marker(ax, [age], [m_sat], color, marker)


def merger_ratio_cut(m_ratio, a):
    return m_ratio < 3 and a > 0.1


# ==============================================================================
#
# Then the full plot itself
#
# ==============================================================================
fig = plt.figure(figsize=[16, 17])
gs = gridspec.GridSpec(
    nrows=5, ncols=2, hspace=0, wspace=0.3, height_ratios=[4, 1, 1.5, 4, 1]
)  # include dummy spacer
ax_l = fig.add_subplot(gs[0, 0], projection="bpl")
ax_t = fig.add_subplot(gs[0, 1], projection="bpl")
ax_dl = fig.add_subplot(gs[1, 0], projection="bpl")
ax_dt = fig.add_subplot(gs[1, 1], projection="bpl")
ax_r = fig.add_subplot(gs[3, 0], projection="bpl")
ax_j = fig.add_subplot(gs[3, 1], projection="bpl")
ax_dr = fig.add_subplot(gs[4, 0], projection="bpl")
ax_dj = fig.add_subplot(gs[4, 1], projection="bpl")

c1, m1 = bpl.color_cycle[0], "."  # "o"
c2, m2 = bpl.color_cycle[1], "s"  # "v"
c3, m3 = bpl.color_cycle[2], "P"  # "+"
c4, m4 = bpl.color_cycle[3], ">"  # "x"
c5, m5 = bpl.color_cycle[4], "v"
c6, m6 = bpl.color_cycle[5], "x"

print(
    "\nThe most massive galaxy at z=0 is not the most massive galaxy at all times,\n"
    "so in this plot I manually adjust which galaxy goes on which axis.\n"
    "Make sure this is correct.\n"
)

plot_merger(ax_t, tl_nbody.g1, c1, m1, "T&L")
plot_merger(ax_t, tl_sfe001_hn20.g2, c2, m2, "T&L")
plot_merger(ax_t, tl_sfe010_hn20.g2, c3, m3, "T&L")
plot_merger(ax_t, tl_sfe100_hn00.g2, c4, m4, "T&L")
plot_merger(ax_t, tl_sfe100_hn05.g2, c5, m5, "T&L")
plot_merger(ax_t, tl_sfe100_hn20.g2, c6, m6, "T&L")

plot_mass_ratio(ax_dt, tl_nbody.g1, tl_nbody.g1, c1)
plot_mass_ratio(ax_dt, tl_nbody.g1, tl_sfe001_hn20.g2, c2)
plot_mass_ratio(ax_dt, tl_nbody.g1, tl_sfe010_hn20.g2, c3)
plot_mass_ratio(ax_dt, tl_nbody.g1, tl_sfe100_hn00.g2, c4)
plot_mass_ratio(ax_dt, tl_nbody.g1, tl_sfe100_hn05.g2, c5)
plot_mass_ratio(ax_dt, tl_nbody.g1, tl_sfe100_hn20.g2, c6)

plot_merger(ax_l, tl_nbody.g2, c1, m1, "T&L")
plot_merger(ax_l, tl_sfe001_hn20.g1, c2, m2, "T&L")
plot_merger(ax_l, tl_sfe010_hn20.g1, c3, m3, "T&L")
plot_merger(ax_l, tl_sfe100_hn00.g1, c4, m4, "T&L")
plot_merger(ax_l, tl_sfe100_hn05.g1, c5, m5, "T&L")
plot_merger(ax_l, tl_sfe100_hn20.g1, c6, m6, "T&L")

plot_mass_ratio(ax_dl, tl_nbody.g2, tl_nbody.g2, c1)
plot_mass_ratio(ax_dl, tl_nbody.g2, tl_sfe001_hn20.g1, c2)
plot_mass_ratio(ax_dl, tl_nbody.g2, tl_sfe010_hn20.g1, c3)
plot_mass_ratio(ax_dl, tl_nbody.g2, tl_sfe100_hn00.g1, c4)
plot_mass_ratio(ax_dl, tl_nbody.g2, tl_sfe100_hn05.g1, c5)
plot_mass_ratio(ax_dl, tl_nbody.g2, tl_sfe100_hn20.g1, c6)

plot_merger(ax_r, rj_nbody_92.g1, c1, m1, "R&J")
plot_merger(ax_r, rj_nbody_46.g1, c2, m2, "R&J")
plot_merger(ax_r, rj_nbody_23.g1, c3, m3, "R&J")

plot_mass_ratio(ax_dr, rj_nbody_92.g1, rj_nbody_92.g1, c1)
plot_mass_ratio(ax_dr, rj_nbody_92.g1, rj_nbody_46.g1, c2)
plot_mass_ratio(ax_dr, rj_nbody_92.g1, rj_nbody_23.g1, c3)

plot_merger(ax_j, rj_nbody_92.g2, c1, m1, "R&J")
plot_merger(ax_j, rj_nbody_46.g2, c2, m2, "R&J")
plot_merger(ax_j, rj_nbody_23.g2, c3, m3, "R&J")

plot_mass_ratio(ax_dj, rj_nbody_92.g2, rj_nbody_92.g2, c1)
plot_mass_ratio(ax_dj, rj_nbody_92.g2, rj_nbody_46.g2, c2)
plot_mass_ratio(ax_dj, rj_nbody_92.g2, rj_nbody_23.g2, c3)

format_main_ax(ax_l, "T&L")
format_main_ax(ax_t, "T&L")
format_diff_ax(ax_dl)
format_diff_ax(ax_dt)
format_main_ax(ax_r, "R&J")
format_main_ax(ax_j, "R&J")
format_diff_ax(ax_dr)
format_diff_ax(ax_dj)

ax_t.easy_add_text("Andromeda Analog - Thelma", "upper left")
ax_l.easy_add_text("Milky Way Analog - Louise", "upper left")
ax_r.easy_add_text("Andromeda Analog - Romeo", "upper left")
ax_j.easy_add_text("Milky Way Analog - Juliet", "upper left")

ax_l.legend(frameon=False, loc=4)
ax_r.legend(frameon=False, loc=4)
# Commented out code below shows how to make markers be next to the legend. I'll ignore
# that for now, since the colors should be enough
# # I want both the line color and markers in the legend.
# # I'll accomplish this with a bit of a hack
# # I'll increase the spacing between the legend entry and labels,
# # then manually plot points
# ax_l.legend(frameon=False, loc=4, handletextpad=2.0)
# x_fake = 5.5
# # plot_merger_marker(ax_l, [x_fake], [1.162E10], c1, m1)# 1.16 to 1.165
# plot_merger_marker(ax_l, [x_fake], [6.24e9], c1, m1)
# plot_merger_marker(ax_l, [x_fake], [3.49e9], c2, m2)  # 3.48 to 3.5
# plot_merger_marker(ax_l, [x_fake], [1.97e9], c3, m3)

fig.savefig(plot_name, bbox_inches="tight")
