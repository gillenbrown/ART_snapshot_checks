import sys
from pathlib import Path

import numpy as np
import yt
from astropy import cosmology
from astropy import units as u

import betterplotlib as bpl

bpl.set_style()
yt.funcs.mylog.setLevel(50)  # ignore yt's output

out_plot = Path(sys.argv[1]).resolve()

# ==============================================================================
#
# Set up cosmology to be able to have the second time axis
#
# ==============================================================================
class Cosmology(cosmology.FlatLambdaCDM):
    def __init__(self, ds_loc):

        ds = yt.load(ds_loc)

        H0 = ds.artio_parameters["hubble"][0] * 100
        omega_m = ds.artio_parameters["OmegaM"][0]

        super().__init__(H0, omega_m, Tcmb0=2.725)

    # helper functions
    @staticmethod
    def a_to_z(a):
        if a <= 0:
            return np.inf
        else:
            return (1.0 / a) - 1

    def z_to_age_gyr(self, z):
        return self.age(z).to("Gyr").value

    def age_gyr_to_z(self, age):
        age *= u.Gyr
        if age >= 0.95 * self.age(0):
            # dumb guess that will be used in minimize
            return -1 * ((age / self.age(0)).value - 1)
        elif age <= 0:
            return np.inf
        elif age < 1 * u.Myr:
            return 608 / age.to(u.Myr).value  # scale to real values
        return cosmology.z_at_value(self.age, age)


dir_tl = Path(
    "/u/home/gillenb/art_runs/runs/pleiades/nbody/new_ic_trim_25mpc/root_08/run/outputs/vel_offset/"
)
dir_rj = Path(
    "/u/home/gillenb/art_runs/runs/stampede2/rj_nbody/hybrid_23.12mpc_level08/run/"
)
dir_hui = Path("/u/home/gillenb/art_runs/runs/shangrila/hui/dm_only")

cosmo_tl = Cosmology(str(dir_tl / "out" / "continuous_a1.0007.art"))
cosmo_rj = Cosmology(str(dir_rj / "out" / "continuous_a1.0018.art"))
cosmo_hui = Cosmology(str(dir_hui / "out" / "continuous_a1.0005.art"))

# ==============================================================================
#
# Code to read in halos
#
# ==============================================================================
class Halo(object):
    def __init__(self, mergers_file_loc, growth_file_loc, name, cosmo, color):
        self.color = color
        self.cosmo = cosmo
        self.h = cosmo.h
        self._init_mergers(mergers_file_loc)
        self._init_growth(growth_file_loc)

        self.mergers_z = [self.cosmo.a_to_z(a) for a in self.mergers_scale]
        self.growth_z = [self.cosmo.a_to_z(a) for a in self.growth_scale]

        self.mergers_age = [self.cosmo.z_to_age_gyr(z) for z in self.mergers_z]
        self.growth_age = [self.cosmo.z_to_age_gyr(z) for z in self.growth_z]

        self.name = name

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
                (
                    a_final,
                    m_sat_final,
                    m_cen_final,
                    id_sat_final,
                    id_cen_final,
                ) = split_line[-5:]
                self.mergers_scale.append(float(a))
                self.mergers_sat_mass.append(float(m_sat) / self.h)
                self.mergers_central_mass.append(float(m_cen) / self.h)
                self.mergers_mass_ratio.append(float(m_cen) / float(m_sat))
                self.mergers_scale_final.append(float(a_final))

    def _init_growth(self, growth_file_loc):
        self.growth_scale = []
        self.growth_mass = []
        with open(growth_file_loc, "r") as in_file:
            for line in in_file:
                if line.startswith("#"):
                    continue
                # otherwise, the data is scale factor, halo mass, then id
                a, m, id = line.split()
                self.growth_scale.append(float(a))
                self.growth_mass.append(float(m) / self.h)


class Simulation(object):
    def __init__(self, output_dir, name_1, name_2, color_1, color_2, cosmo):
        self.cosmo = cosmo

        mergers_1 = output_dir / "checks" / "mergers_1.txt"
        mergers_2 = output_dir / "checks" / "mergers_2.txt"

        growth_1 = output_dir / "checks" / "growth_1.txt"
        growth_2 = output_dir / "checks" / "growth_2.txt"

        self.halo_1 = Halo(mergers_1, growth_1, name_1, self.cosmo, color_1)
        self.halo_2 = Halo(mergers_2, growth_2, name_2, self.cosmo, color_2)


# ==============================================================================
#
# do the actual read in
#
# ==============================================================================
tl = Simulation(
    dir_tl, "Thelma", "Louise", bpl.color_cycle[0], bpl.color_cycle[1], cosmo_tl
)
rj = Simulation(
    dir_rj, "Romeo", "Juliet", bpl.color_cycle[3], bpl.color_cycle[4], cosmo_rj
)
hui = Simulation(dir_hui, "Isolated MW", "", bpl.color_cycle[2], "", cosmo_hui)


# ==============================================================================
#
# Then plot
#
# ==============================================================================
msize = 100


def plot_merger(ax, halo, zorder, marker):
    if marker == "+":
        msize_plot = msize * np.sqrt(2)
    elif marker == "o":
        msize_plot = 0.7 * msize
    else:
        msize_plot = msize

    ax.plot(
        halo.growth_age,
        halo.growth_mass,
        lw=2,
        label=halo.name,
        color=halo.color,
        zorder=zorder + 100,
    )

    for a_i, a_f, m_sat, m_ratio in zip(
        halo.mergers_scale,
        halo.mergers_scale_final,
        halo.mergers_sat_mass,
        halo.mergers_mass_ratio,
    ):
        if m_ratio < 4 and m_sat > 1e9:
            age = halo.cosmo.z_to_age_gyr(halo.cosmo.a_to_z(a_i))
            ax.scatter(
                [age],
                [m_sat],
                c=halo.color,
                marker=marker,
                linewidth=1.5,
                alpha=1.0,
                s=msize_plot,
                zorder=60 + 1 / m_ratio,
            )


fig, ax = bpl.subplots(figsize=[8, 7])

plot_merger(ax, tl.halo_1, 1, "+")
plot_merger(ax, tl.halo_2, 1, "x")
plot_merger(ax, rj.halo_1, 1, "^")
plot_merger(ax, rj.halo_2, 1, "v")
plot_merger(ax, hui.halo_1, 1, "o")

ax.add_labels("Cosmic Time [Gyr]", "Halo Mass [$M_\odot$]")
ax.set_limits(0, 14, 1e10, 2e12)
ax.twin_axis(
    "x", [0, 0.5, 1, 2, 3, 5], label="Redshift", old_to_new_func=cosmo_tl.age_gyr_to_z
)
ax.set_yscale("log")
ax.yaxis.set_ticks_position("both")

# I want both the line color and markers in the legend.
# I'll accomplish this with a bit of a hack
# I'll increase the spacing between the legend entry and labels,
# then manually plot points
ax.legend(frameon=False, loc=4, handletextpad=1.5)
x_marker = 9.8  # 9.0-9.9
ax.scatter(
    x_marker,
    6.48e10,
    c=tl.halo_1.color,
    marker="+",
    linewidth=1.5,
    alpha=1.0,
    s=msize * np.sqrt(2),
)
ax.scatter(
    x_marker,
    4.515e10,  # 4.51-4.53
    c=tl.halo_2.color,
    marker="x",
    linewidth=1.5,
    alpha=1.0,
    s=msize * np.sqrt(2),
)
ax.scatter(
    x_marker,
    3.15e10,  # 3.10-3.20
    c=rj.halo_1.color,
    marker="^",
    linewidth=1.5,
    alpha=1.0,
    s=msize,
)
ax.scatter(
    x_marker,
    2.195e10,  # 2.19-2.20
    c=rj.halo_2.color,
    marker="v",
    linewidth=1.5,
    alpha=1.0,
    s=msize,
)
ax.scatter(
    x_marker,
    1.535e10,  # 1.53-1.54
    c=hui.halo_1.color,
    marker="o",
    linewidth=1.5,
    alpha=1.0,
    s=msize * 0.7,
)

# I also want to reconfigure the x labels
x_values = [0, 2, 4, 6, 8, 10, 12, 14]
x_labels = [str(i) for i in x_values]
ax.set_xticks(x_values)
ax.set_xticklabels(x_labels)

fig.savefig(out_plot)
