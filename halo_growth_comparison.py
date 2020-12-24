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

bpl.set_style()
yt.funcs.mylog.setLevel(50)  # ignore yt's output

plot_name = Path(sys.argv[1])
# ==============================================================================
# 
# We need the time, so we need to load the dataset to get the cosmology info
#
# ==============================================================================
out_dir = Path("/u/home/gillenb/art_runs/runs/stampede2/production/sfe100/out/")
# first validate that this exists
if not out_dir.exists():
    plot_name.touch()
    exit()
# get a simulation from this path
found_sims = []
for item in out_dir.iterdir():
    if item.name.endswith(".art"):
        found_sims.append(item)
# if we didn't find anything, exit
if len(found_sims) == 0:
    plot_name.touch()
    exit()
# get the earliest one so there are less stars and it loads quicker
ds = yt.load(str(sorted(found_sims)[0]))
H0 = ds.artio_parameters["hubble"][0]*100
omega_m = ds.artio_parameters["OmegaM"][0]
cosmo = cosmology.FlatLambdaCDM(H0, omega_m, Tcmb0=2.725)

# define some helper functions to handle ages
def a_to_z(a):
    if a <= 0:
        return np.inf
    else:
        return (1.0 / a) - 1
    
def z_to_age_gyr(z):
    return cosmo.age(z).to("Gyr").value

def age_gyr_to_z(age):
    age *= u.Gyr
    if age >= 0.95 * cosmo.age(0):
        return -1 * ((age / cosmo.age(0)).value - 1)  # dumb guess that will be used in minimize
    elif age <= 0:
        return np.inf
    elif age < 1 * u.Myr:
        return 608 / age.to(u.Myr).value  # scale to real values
    return cosmology.z_at_value(cosmo.age, age)

# ==============================================================================
# 
# Class to handle the halo catalogs
#
# ==============================================================================
class Halo(object):
    def __init__(self, mergers_file_loc, growth_file_loc, name):
        self.h = 0.710000  
        self._init_mergers(mergers_file_loc)
        self._init_growth(growth_file_loc)
        
        self.mergers_z = [a_to_z(a) for a in self.mergers_scale]
        self.growth_z = [a_to_z(a) for a in self.growth_scale]
        
        self.mergers_age = [z_to_age_gyr(z) for z in self.mergers_z]
        self.growth_age = [z_to_age_gyr(z) for z in self.growth_z]
        
        self.name = name
        
        # then create an interpolation object, which will be used for calculating the ratios
        # of the growth history
        if len(self.growth_age) == 0:
            raise ValueError("bad halo growth")
        elif len(self.growth_age) == 1:
            self.growth_age_interp = None
        else:
            self.growth_age_interp = interpolate.interp1d(x=self.growth_age, y=self.growth_mass, kind="linear")

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
                a_final, m_sat_final, m_cen_final, id_sat_final, id_cen_final = split_line[-5:]
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
    def __init__(self, output_dir, name_base):
        self.dirname = Path(output_dir).resolve()
        mergers_thelma = self.dirname / "mergers_1.txt"
        mergers_louise = self.dirname / "mergers_2.txt"
        
        growth_thelma = self.dirname / "growth_1.txt"
        growth_louise = self.dirname / "growth_2.txt"
        
        self.thelma = Halo(mergers_thelma, growth_thelma, name_base)
        self.louise = Halo(mergers_louise, growth_louise, name_base)

l1_128 = Simulation("/u/home/gillenb/art_runs/runs/pleiades/nbody/new_ic_50mpc/root_07/run/outputs/vel_offset/checks/", 
                    "Collisionless L1-128")
l2_256 = Simulation("/u/home/gillenb/art_runs/runs/pleiades/nbody/new_ic_trim_25mpc/root_08/run/outputs/vel_offset/checks/", 
                    "Collisionless L2-256")
l2_256_sfe100 = Simulation("/u/home/gillenb/art_runs/runs/stampede2/production/sfe100/checks/", 
                           "Hydro SFE 100 L2-256")

# ==============================================================================
# 
# Then we plot! First a bunch of helper functions
#
# ==============================================================================
def format_main_ax(ax):
    ax.add_labels("", "Mass [$M_\odot$]")
    ax.set_limits(0, 14, 1E9, 3E12)
    ax.set_xticklabels([])
    ax.tick_params(axis="x", direction="in")
    
    ax.twin_axis("x", [0, 0.5, 1, 2, 3, 5], label="Redshift", old_to_new_func=age_gyr_to_z)
    ax.set_yscale("log")
    ax.yaxis.set_ticks_position("both")
    
    ax.yaxis.set_ticks([1E10, 1E11, 1E12])
    ax.yaxis.set_ticklabels(["$10^{10}$", "$10^{11}$", "$10^{12}$"])
    
def format_diff_ax(ax):
    ax.add_labels("Cosmic Time [Gyr]", "$\Delta$ M [dex]")
    ax.set_limits(0, 14, -0.5, 0.5)
    ax.yaxis.set_ticks_position("both")

def plot_mass_ratio(ax, halo1, halo2, color):
    age_start = max(halo1.growth_age[-1], halo2.growth_age[-1])
    age_end   = min(halo1.growth_age[0],  halo2.growth_age[0])
    
    ages = np.linspace(age_start, age_end, 1000)
    
    if halo1.growth_age_interp is not None and halo2.growth_age_interp is not None:
        ratios = [halo2.growth_age_interp(a) / halo1.growth_age_interp(a)
                  for a in ages]
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
            facecolor=color
        else:
            facecolor="none"
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
        
        
    ax.scatter(x, y, color=color, marker=marker, facecolor=facecolor,
               linewidth=linewidth, alpha=1.0, s=msize, zorder=zorder)

def plot_merger(ax, halo, zorder, color, marker):
    ax.plot(halo.growth_age, halo.growth_mass, lw=3, 
            label=halo.name, color=color, zorder=zorder+100)

    for a_i, a_f, m_sat, m_ratio in zip(halo.mergers_scale, halo.mergers_scale_final, halo.mergers_sat_mass, halo.mergers_mass_ratio):
        if merger_ratio_cut(m_ratio, a_i) and m_sat > 1E9:
            age = z_to_age_gyr(a_to_z(a_i))
            plot_merger_marker(ax, [age], [m_sat], color, marker, zorder=zorder)

def merger_ratio_cut(m_ratio, a):
    return m_ratio < 3 and a > 0.1

# ==============================================================================
# 
# Then the full plot itself
#
# ==============================================================================
fig = plt.figure(figsize=[15, 8])
gs = gridspec.GridSpec(nrows=2, ncols=2, hspace=0, wspace=0.3, height_ratios=[4, 1])
ax_l = fig.add_subplot(gs[0,0], projection="bpl")
ax_t = fig.add_subplot(gs[0,1], projection="bpl")
ax_dl = fig.add_subplot(gs[1,0], projection="bpl")
ax_dt = fig.add_subplot(gs[1,1], projection="bpl")

c1, m1 = bpl.color_cycle[0], "." #"o"
c2, m2 = bpl.color_cycle[1], "s" #"v"
c3, m3 = bpl.color_cycle[3], "P" #"+"
c4, m4 = bpl.color_cycle[4], ">" #"x"

plot_merger(ax_t, l1_128.thelma, 2, c1, m1)
plot_merger(ax_t, l2_256.thelma, 1, c2, m2)
plot_merger(ax_t, l2_256_sfe100.thelma, 1, c3, m3)

plot_mass_ratio(ax_dt, l2_256.thelma, l1_128.thelma, c1)
plot_mass_ratio(ax_dt, l2_256.thelma, l2_256.thelma, c2)
plot_mass_ratio(ax_dt, l2_256.thelma, l2_256_sfe100.thelma, c3)

plot_merger(ax_l, l1_128.louise, 2, c1, m1)
plot_merger(ax_l, l2_256.louise, 1, c2, m2)
plot_merger(ax_l, l2_256_sfe100.louise, 1, c3, m3)

plot_mass_ratio(ax_dl, l2_256.louise, l1_128.louise, c1)
plot_mass_ratio(ax_dl, l2_256.louise, l2_256.louise, c2)
plot_mass_ratio(ax_dl, l2_256.louise, l2_256_sfe100.louise, c3)

format_main_ax(ax_l)
format_main_ax(ax_t)
format_diff_ax(ax_dl)
format_diff_ax(ax_dt)
    
ax_t.easy_add_text("Andromeda Analog", "upper left")
ax_l.easy_add_text("Milky Way Analog", "upper left")

# I want both the line color and markers in the legend. 
# I'll accomplish this with a bit of a hack
# I'll increase the spacing between the legend entry and labels,
# then manually plot points
ax_l.legend(frameon=False, loc=4, handletextpad=2.0)
x_fake = 5.5
# plot_merger_marker(ax_l, [x_fake], [1.162E10], c1, m1)# 1.16 to 1.165
plot_merger_marker(ax_l, [x_fake], [6.24E9], c1, m1) 
plot_merger_marker(ax_l, [x_fake], [3.49E9], c2, m2) # 3.48 to 3.5
plot_merger_marker(ax_l, [x_fake], [1.97E9], c3, m3)
    
fig.savefig(plot_name, bbox_inches="tight")