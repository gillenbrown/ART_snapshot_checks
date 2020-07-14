# https://www.gnu.org/software/make/manual/html_node/Automatic-Variables.html
# My notation is that no paths end in a slash
# ------------------------------------------------------------------------------
#
#  Flag to figure out what machine we're on
# 
# ------------------------------------------------------------------------------
# This tells us which directories to look for simulations in, which will 
# be used later. We'll use the hostname of the machine
hostname = $(shell hostname)
# findstring returns the matching part of the string. If it's not empty when
# we try to find the shangrila hostname, we know we're on shangrila
ifneq (,$(findstring shangrila,$(hostname)))
	machine = shangrila
endif
ifneq (,$(findstring ldan,$(hostname)))
	machine = lou
endif
ifneq (,$(findstring lfe,$(hostname)))
	machine = lou
endif
ifneq (,$(findstring gl-login,$(hostname)))
	machine = great_lakes
endif
ifneq (,$(findstring stampede2,$(hostname)))
	machine = stampede2
endif

# ------------------------------------------------------------------------------
#
#  name of python executable - depends on the machine
#
# ------------------------------------------------------------------------------
ifeq ($(machine),shangrila)
	python=python
endif
ifeq ($(machine),lou)
	python=python
endif
ifeq ($(machine),great_lakes)
	python=python
endif
ifeq ($(machine),stampede2)
	python=python3
endif

# ------------------------------------------------------------------------------
#
#  Code locations - depend on the machine
# 
# ------------------------------------------------------------------------------
ifeq ($(machine),shangrila)
	tree_config_script = /u/home/gillenb/code/not_mine/rockstar/scripts/gen_merger_cfg.pl
	tree_dir = /u/home/gillenb/code/not_mine/consistent-trees
	halo_finding_script = ./run_rockstar.sh
	timing_script = /u/home/gillenb/code/mine/cart/utils/scripts/parse_timing2.pl
endif
ifeq ($(machine),lou)
	tree_config_script = /u/gbrown12/code/rockstar/scripts/gen_merger_cfg.pl
	tree_dir = /u/gbrown12/code/consistent-trees
	halo_finding_script = ./run_rockstar_ldan.sh
	timing_script = /u/gbrown12/code/parse_timing2.pl
endif
ifeq ($(machine),great_lakes)
	tree_config_script = /home/gillenb/code/rockstar/scripts/gen_merger_cfg.pl
	tree_dir = /home/gillenb/code/consistent-trees
	halo_finding_script = ./run_rockstar_gl.sh
	timing_script = /home/gillenb/code/art_cluster/utils/scripts/parse_timing2.pl
endif
ifeq ($(machine),stampede2)
	tree_config_script = $(HOME)/code/rockstar-galaxies/scripts/gen_merger_cfg.pl
	tree_dir = $(HOME)/code/consistent-trees
	halo_finding_script = ./run_rockstar_stampede2.sh
	timing_script = $(HOME)/code/cart/utils/scripts/parse_timing2.pl
endif

# ------------------------------------------------------------------------------
#
#  Code locations that are relative to this file
# 
# ------------------------------------------------------------------------------
debug_script = ./debug_output.py
galaxies_script = ./galaxy_summaries.py
nbody_single_halo_plots_script = ./plot_single_halo_nbody.py
nbody_refined_plot_script = ./plot_refined_region_nbody.py
nbody_local_group_plot_script = ./plot_local_group_nbody.py
nbody_full_plot_script = ./utils/nbody_projection_all_species.py
nbody_split_plot_script = ./utils/nbody_projection_split_species.py
halo_growth_comp_script = ./halo_growth_comparison.py
sfh_plots_script = ./plot_sfh.py
cimf_plots_script = ./plot_cimf.py
dt_history_script = ./dt_history.py
prep_for_merger_tree_script = ./merger_prep.py
read_tree_dir = ./read_tree
read_tree_exe = $(read_tree_dir)/halo_history
read_tree_src = $(read_tree_dir)/halo_history.c
comparison_plots_dir = ./comparison_plots

# ------------------------------------------------------------------------------
#
#  Simulation outputs to run this on - depend on the machine
# 
# ------------------------------------------------------------------------------
ifeq ($(machine),shangrila)
	runs_home = /u/home/gillenb/art_runs/runs
	sim_dirs_nbody = 
	sim_dirs_hydro = $(runs_home)/shangrila/old_ic_comparison/default_1e7_temp_cap/run
	#$(runs_home)/shangrila/old_ic_comparison/default/run \
	                 #$(runs_home)/shangrila/hui/sfe_10 \
	                 #$(runs_home)/shangrila/hui/sfe_50 \
	                 #$(runs_home)/shangrila/hui/sfe_100 \
	                 #$(runs_home)/shangrila/hui/sfe_200 \
	                 #$(runs_home)/stampede2/production/sfe100
endif

ifeq ($(machine),lou)
	runs_home = /u/gbrown12/art_runs/runs
	sim_dirs_nbody = $(runs_home)/nbody/new_ic_trim_12.5mpc/root_05/run/outputs/fix_vel \
	                 $(runs_home)/nbody/new_ic_trim_12.5mpc/root_06/run/outputs/fix_vel \
	                 $(runs_home)/nbody/new_ic_trim_12.5mpc/root_07/run/outputs/fix_vel \
	                 $(runs_home)/nbody/new_ic_trim_12.5mpc/root_08/run/outputs/fix_vel \
	                 $(runs_home)/nbody/new_ic_trim_12.5mpc/root_05/run/outputs/vel_offset \
	                 $(runs_home)/nbody/new_ic_trim_12.5mpc/root_06/run/outputs/vel_offset \
	                 $(runs_home)/nbody/new_ic_trim_12.5mpc/root_07/run/outputs/vel_offset \
	                 $(runs_home)/nbody/new_ic_trim_12.5mpc/root_08/run/outputs/vel_offset \
	                 $(runs_home)/nbody/new_ic_trim_25mpc/root_06/run/outputs/vel_offset \
	                 $(runs_home)/nbody/new_ic_trim_25mpc/root_07/run/outputs/vel_offset \
	                 $(runs_home)/nbody/new_ic_trim_25mpc/root_08/run/outputs/vel_offset \
	                 $(runs_home)/nbody/new_ic_50mpc/root_07/run/outputs/fix_vel
	sim_dirs_hydro = $(runs_home)/hydro/new_ic_50mpc/root_07/run/outputs/first \
	                 $(runs_home)/hydro/new_ic_trim_12.5mpc/root_05/run/outputs/first \
	                 $(runs_home)/hydro/new_ic_trim_12.5mpc/root_06/run/outputs/first \
	                 $(runs_home)/hydro/new_ic_trim_12.5mpc/root_07/run/outputs/first \
	                 $(runs_home)/hydro/new_ic_trim_12.5mpc/root_08/run/outputs/first \
	                 $(runs_home)/hydro/old_ic/default/run/outputs/run_1 \
	                 $(runs_home)/hydro/old_ic/no_hn/run/outputs/run_1
	                 #$(runs_home)/hydro/old_ic/continuous/run/outputs/run_1 \
	                 #$(runs_home)/hydro/old_ic/continuous_no_hn/run/outputs/run_1 \
	                 #$(runs_home)/hydro/old_ic/no_changes/run/outputs/run_1 \
	                 #$(runs_home)/hydro/old_ic/no_alpha/run/outputs/run_1
endif

ifeq ($(machine),great_lakes)                
	runs_home = /home/gillenb/art_runs/runs
	sim_dirs_nbody = 
	sim_dirs_hydro = $(runs_home)/hydro_test/discrete_128/run/outputs/first \
	                 $(runs_home)/hydro_test/discrete_256/run/outputs/first \
	                 $(runs_home)/hydro_test/scaling_discrete_256/run/outputs/1_cores \
	                 $(runs_home)/hydro_test/scaling_discrete_256/run/outputs/2_cores \
	                 $(runs_home)/hydro_test/scaling_discrete_256/run/outputs/4_cores \
	                 $(runs_home)/hydro_test/scaling_discrete_256/run/outputs/6_cores \
	                 $(runs_home)/hydro_test/scaling_discrete_256/run/outputs/8_cores 
endif

ifeq ($(machine),stampede2)
	runs_home = $(SCRATCH)/art_runs/runs
	sim_dirs_nbody =
	sim_dirs_hydro = $(runs_home)/production/sfe_100/run
endif

# combine the N-Body and Hydro into one big list
sim_dirs = $(sim_dirs_nbody) $(sim_dirs_hydro)
# ------------------------------------------------------------------------------
#
#  Directories for each simulation, where outputs will be stored
# 
# ------------------------------------------------------------------------------
sim_out_dirs = $(foreach dir,$(sim_dirs),$(dir)/out)
sim_out_dirs_hydro = $(foreach dir,$(sim_dirs_hydro),$(dir)/out)
sim_rockstar_halos_dirs = $(foreach dir,$(sim_dirs),$(dir)/rockstar_halos)
sim_checks_dirs = $(foreach dir,$(sim_dirs),$(dir)/checks)
sim_plots_dirs = $(foreach dir,$(sim_dirs),$(dir)/plots)
my_directories = $(sim_checks_dirs) $(sim_rockstar_halos_dirs) $(sim_plots_dirs) $(comparison_plots_dir)

# ------------------------------------------------------------------------------
#
#  List of all simulation outputs and their corresponding halo catalogs
# 
# ------------------------------------------------------------------------------
dir_to_sims = $(wildcard $(1)/*_a*.art)
snapshots = $(foreach dir,$(sim_out_dirs),$(call dir_to_sims,$(dir)))
snapshots_hydro = $(foreach dir,$(sim_out_dirs_hydro),$(wildcard $(dir)/*_a*.art))
# Parse the snapshot names into halo catalogs. These will all be in their
# own special directory
sim_to_halo_dir = $(subst .art,,$(subst out/continuous,rockstar_halos/halos,$(1)))
halos_catalogs = $(foreach snapshot,$(snapshots),$(call sim_to_halo_dir,$(snapshot)/halos_0.0.bin))
halos_dirs = $(foreach cat,$(halos_catalogs),$(dir $(cat)))
# then some functions to go backwards to determine the simulation output
halo_catalog_to_sim = $(subst /halos_0.0.bin,.art,$(subst rockstar_halos/halos,out/continuous,$(1)))
# add these directories to what we need to create
my_directories += $(halos_dirs)

# ------------------------------------------------------------------------------
#
#  Debug output files
# 
# ------------------------------------------------------------------------------
sim_to_debug = $(subst .art,.txt,$(subst out/continuous,checks/debug, $(1)))
debug_to_sim = $(subst .txt,.art,$(subst checks/debug,out/continuous, $(1)))
debugs = $(foreach snapshot,$(snapshots),$(call sim_to_debug,$(snapshot)))

# ------------------------------------------------------------------------------
#
#  Galaxies output files
# 
# ------------------------------------------------------------------------------
sim_to_galaxies = $(subst .art,.txt,$(subst out/continuous,checks/galaxy_summaries, $(1)))
galaxies_to_sim = $(subst .txt,.art,$(subst checks/galaxy_summaries,out/continuous, $(1)))
galaxies_to_halo = $(subst .txt,/halos_0.0.bin,$(subst checks/galaxy_summaries,rockstar_halos/halos, $(1)))
galaxies = $(foreach snapshot,$(snapshots_hydro),$(call sim_to_galaxies,$(snapshot)))

# ------------------------------------------------------------------------------
#
#  comparison plots that have all sims on one plot. Only on shangrila!
# 
# ------------------------------------------------------------------------------
ifeq ($(machine),shangrila)
	sfh_plots = $(comparison_plots_dir)/mass_growth_comparison.png \
				$(comparison_plots_dir)/sfh_comparison.png
	cimf_plots = $(comparison_plots_dir)/cimf_common.png \
				 $(comparison_plots_dir)/cimf_last.png
	halo_growth_plot = $(comparison_plots_dir)/halo_growth.png
else
	sfh_plots =
	cimf_plots =
	halo_growth_plot =
endif

# ------------------------------------------------------------------------------
#
#  Nbody plots
# 
# ------------------------------------------------------------------------------
sim_to_halo_1_full = $(subst .art,.png,$(subst out/continuous,plots/n_body_halo_rank_1, $(1)))
sim_to_halo_2_full = $(subst .art,.png,$(subst out/continuous,plots/n_body_halo_rank_2, $(1)))
sim_to_halo_1_split = $(subst .art,.png,$(subst out/continuous,plots/n_body_split_halo_rank_1, $(1)))
sim_to_halo_2_split = $(subst .art,.png,$(subst out/continuous,plots/n_body_split_halo_rank_2, $(1)))
sim_to_refined_full  = $(subst .art,.png,$(subst out/continuous,plots/n_body_refined, $(1)))
sim_to_refined_split = $(subst .art,.png,$(subst out/continuous,plots/n_body_split_refined, $(1)))
sim_to_local_group_full = $(subst .art,.png,$(subst out/continuous,plots/n_body_local_group, $(1)))
sim_to_local_group_split = $(subst .art,.png,$(subst out/continuous,plots/n_body_split_local_group, $(1)))

halo_1_full_to_sim = $(subst .png,.art,$(subst plots/n_body_halo_rank_1,out/continuous, $(1)))
halo_2_full_to_sim = $(subst .png,.art,$(subst plots/n_body_halo_rank_2,out/continuous, $(1)))
halo_1_split_to_sim = $(subst .png,.art,$(subst plots/n_body_split_halo_rank_1,out/continuous, $(1)))
halo_2_split_to_sim = $(subst .png,.art,$(subst plots/n_body_split_halo_rank_2,out/continuous, $(1)))
refined_full_to_sim = $(subst .png,.art,$(subst plots/n_body_refined,out/continuous, $(1)))
refined_split_to_sim = $(subst .png,.art,$(subst plots/n_body_split_refined,out/continuous, $(1)))
local_group_full_to_sim = $(subst .png,.art,$(subst plots/n_body_local_group,out/continuous, $(1)))
local_group_split_to_sim = $(subst .png,.art,$(subst plots/n_body_split_local_group,out/continuous, $(1)))

halo_1_full_to_halo = $(subst .png,/halos_0.0.bin,$(subst plots/n_body_halo_rank_1,rockstar_halos/halos, $(1)))
halo_2_full_to_halo = $(subst .png,/halos_0.0.bin,$(subst plots/n_body_halo_rank_2,rockstar_halos/halos, $(1)))
halo_1_split_to_halo = $(subst .png,/halos_0.0.bin,$(subst plots/n_body_split_halo_rank_1,rockstar_halos/halos, $(1)))
halo_2_split_to_halo = $(subst .png,/halos_0.0.bin,$(subst plots/n_body_split_halo_rank_2,rockstar_halos/halos, $(1)))
refined_full_to_halo = $(subst .png,/halos_0.0.bin,$(subst plots/n_body_refined,rockstar_halos/halos, $(1)))
refined_split_to_halo = $(subst .png,/halos_0.0.bin,$(subst plots/n_body_split_refined,rockstar_halos/halos, $(1)))
local_group_full_to_halo = $(subst .png,/halos_0.0.bin,$(subst plots/n_body_local_group,rockstar_halos/halos, $(1)))
local_group_split_to_halo = $(subst .png,/halos_0.0.bin,$(subst plots/n_body_split_local_group,rockstar_halos/halos, $(1)))

halo_1_full_plots = $(foreach snapshot,$(snapshots),$(call sim_to_halo_1_full,$(snapshot)))
halo_2_full_plots = $(foreach snapshot,$(snapshots),$(call sim_to_halo_2_full,$(snapshot)))
halo_1_split_plots = $(foreach snapshot,$(snapshots),$(call sim_to_halo_1_split,$(snapshot)))
halo_2_split_plots = $(foreach snapshot,$(snapshots),$(call sim_to_halo_2_split,$(snapshot)))
refined_full_plots = $(foreach snapshot,$(snapshots),$(call sim_to_refined_full,$(snapshot)))
refined_split_plots = $(foreach snapshot,$(snapshots),$(call sim_to_refined_split,$(snapshot)))
local_group_full_plots = $(foreach snapshot,$(snapshots),$(call sim_to_local_group_full,$(snapshot)))
local_group_split_plots = $(foreach snapshot,$(snapshots),$(call sim_to_local_group_split,$(snapshot)))

nbody_plots = $(halo_1_full_plots) $(halo_2_full_plots) $(halo_1_split_plots) $(halo_2_split_plots) $(refined_full_plots) $(refined_split_plots) $(local_group_full_plots) $(local_group_split_plots)

# ------------------------------------------------------------------------------
#
#  Consistent trees
# 
# ------------------------------------------------------------------------------
# we need to get the correct rockstar config files
rockstar_cfgs = $(foreach d,$(sim_rockstar_halos_dirs),$(d)/rockstar.cfg)
# This is a bit tricky here, but we basically just use this to find all the 
# subdirectories in the same directory as the rockstar config file
rockstar_cfg_to_halos = $(foreach d,$(shell find $(dir $(1)) -mindepth 1 -maxdepth 1 -type d),$(d)/halos_0.0.bin)
# ^ Here I use the .bin files, even though I really need the out_0.list files.
# this is because I require the .bin files elsewhere, so use those as the 
# proxies for the halo creation. 

# then we use these to make the merger tree config files
tree_cfgs = $(foreach dir,$(sim_rockstar_halos_dirs),$(dir)/outputs/merger_tree.cfg)
tree_cfg_to_rockstar_cfg = $(subst outputs/merger_tree.cfg,rockstar.cfg,$(1))

# then the actual halo trees
trees = $(foreach dir,$(sim_rockstar_halos_dirs),$(dir)/trees/tree_0_0_0.dat)
tree_to_tree_cfg = $(subst trees/tree_0_0_0.dat,outputs/merger_tree.cfg,$(1))

# ------------------------------------------------------------------------------
#
#  Parsing merger trees
# 
# ------------------------------------------------------------------------------
# Here I again use a sentinel file, since there will be multiple files created
# by the one C file:
merger_sentinels = $(foreach dir,$(sim_checks_dirs),$(dir)/merger_sentinel.txt)
# need to get the trees that the sim should parse
merger_to_tree = $(subst checks/merger_sentinel.txt,rockstar_halos/trees/tree_0_0_0.dat,$(1))

# ------------------------------------------------------------------------------
#
#  timing output and dt history
# 
# ------------------------------------------------------------------------------
# Here we assume all log files are located as follows. Each `log` directory
# has several subdirectories holding the logs for each run, starting with
# `log_`. Any other directories will be ignored. We have to do this weird
# thing we we get all items in each directory, then get the parent of each of
# those, to ensure that what we got is really a directory. 
log_dirs = $(foreach dir,$(sim_dirs),$(foreach item,$(wildcard $(dir)/log/log_*/*),$(dir $(item))))
# sorting this removes any duplicates. Order doesn't matter to me anyway
timing_dirs = $(sort $(log_dirs))
# the output of $(dir ...) keeps the last slash, so don't include it here
timings = $(foreach t_dir,$(timing_dirs),$(t_dir)timing_debug.txt)
dt_history_plots = $(foreach t_dir,$(timing_dirs),$(t_dir)timestep_history.png)

# ------------------------------------------------------------------------------
#
#  Movies
# 
# ------------------------------------------------------------------------------
movies_all = $(foreach dir,$(sim_dirs),$(dir)/plots/$(1).mp4)
movies_hydro = $(foreach dir,$(sim_dirs_hydro),$(dir)/plots/$(1).mp4)

# Here 1 is the full location of the movie
movie_to_out_dir = $(subst plots/,out,$(dir $(1)))
movie_to_sims = $(call dir_to_sims,$(call movie_to_out_dir,$(1)))
# in this next function argument 1 is the full location where the movie will 
# be saved, and 2 is the function that turns a simulation into a plot file 
# (these are all defined above)
movie_to_plots = $(foreach sim,$(call movie_to_sims,$(1)),$(call $(2),$(sim)))
movie_to_metal_summary = $(call sim_to_summary_metal,$(call dir_to_sims,$(subst plots/$(1).mp4,out,$(2))))
movie_to_plot_dir = $(subst /$(1).mp4,,$(2))

# ------------------------------------------------------------------------------
#
#  Rules
# 
# ------------------------------------------------------------------------------
movies = $(call movies_all,n_body_refined) $(call movies_all,n_body_split_refined) $(call movies_all,n_body_local_group) $(call movies_all,n_body_split_local_group)
all: $(my_directories) $(halo_growth_plot) # $(timings) $(dt_history_plots) $(sfh_plots) $(cimf_plots) $(debugs) $(galaxies) $(movies) $(halo_growth_plot)

.PHONY: clean
clean:
	rm -r $(sim_checks_dirs) $(sim_human_halos_dirs) $(sim_plots_dirs)

.PHONY: clean_all
clean_all:
	rm -r $(my_directories)  # includes rockstar halos

# Make directories if they don't exist
.PHONY: dirs
dirs: $(my_directories)

$(my_directories):
	mkdir $@

# have phony target for movies
.PHONY: movies
movies: $(movies)

# Phony target to allow us to make only the halos, for example on 
# a PBS script where we don't want to waste computation on the serial aspects
# that come later
.PHONY: halos
halos: $(halos_catalogs)

# Do the halo finding. We do this on individual snapshots, rather than running
# it on all of the halos at once, as it allows us to better manage the output
# files and feed them to the merger tree code
.SECONDEXPANSION:
$(halos_catalogs): %: $$(call halo_catalog_to_sim, %)
	$(halo_finding_script) $(call halo_catalog_to_sim, $@) $(dir $@)

# Make the debug files
.SECONDEXPANSION:
$(debugs): %: $(debug_script)
	$(python) $(debug_script) $(call debug_to_sim, $@) clobber silent

# Make the summary files
.SECONDEXPANSION:
$(galaxies): %: $$(call galaxies_to_halo, %) $(galaxies_script)
	$(python) $(galaxies_script) $(call galaxies_to_sim, $@) $(call galaxies_to_halo, $@) clobber silent

# Make the CIMF plots. We could use &: instead of : to indicate a grouped
# target, but that required Make 4.3 or higher
$(sfh_plots): $(snapshots_hydro) $(halos_catalogs)
	$(python) $(sfh_plots_script) $(sim_dirs_hydro)

# Make the CIMF plots. We could use &: instead of : to indicate a grouped
# target, but that required Make 4.3 or higher
$(cimf_plots): $(snapshots_hydro) $(halos_catalogs)
	$(python) $(cimf_plots_script) $(sim_dirs_hydro)

# Make the individual nbody plots - several examples of very similar things here
.SECONDEXPANSION:
$(halo_1_full_plots): %: $$(call halo_1_full_to_halo, %) $(nbody_single_halo_plots_script) $(nbody_full_plot_script)
	$(python) $(nbody_single_halo_plots_script) $(call halo_1_full_to_sim, $@) $(call halo_1_full_to_halo, $@) 1 full

.SECONDEXPANSION:
$(halo_2_full_plots): %: $$(call halo_2_full_to_halo, %) $(nbody_single_halo_plots_script) $(nbody_full_plot_script)
	$(python) $(nbody_single_halo_plots_script) $(call halo_2_full_to_sim, $@) $(call halo_2_full_to_halo, $@) 2 full

.SECONDEXPANSION:
$(halo_1_split_plots): %: $$(call halo_1_split_to_halo, %) $(nbody_single_halo_plots_script) $(nbody_split_plot_script)
	$(python) $(nbody_single_halo_plots_script) $(call halo_1_split_to_sim, $@) $(call halo_1_split_to_halo, $@) 1 split

.SECONDEXPANSION:
$(halo_2_split_plots): %: $$(call halo_2_split_to_halo, %) $(nbody_single_halo_plots_script) $(nbody_split_plot_script)
	$(python) $(nbody_single_halo_plots_script) $(call halo_2_split_to_sim, $@) $(call halo_2_split_to_halo, $@) 2 split

.SECONDEXPANSION:
$(refined_full_plots): %: $$(call refined_full_to_halo, %) $(nbody_refined_plot_script) $(nbody_full_plot_script)
	$(python) $(nbody_refined_plot_script) $(call refined_full_to_sim, $@) $(call refined_full_to_halo, $@) full

.SECONDEXPANSION:
$(refined_split_plots): %: $$(call refined_split_to_halo, %) $(nbody_refined_plot_script) $(nbody_split_plot_script)
	$(python) $(nbody_refined_plot_script) $(call refined_split_to_sim, $@) $(call refined_split_to_halo, $@) split

.SECONDEXPANSION:
$(local_group_full_plots): %: $$(call local_group_full_to_halo, %) $(nbody_local_group_plot_script) $(nbody_full_plot_script)
	$(python) $(nbody_local_group_plot_script) $(call local_group_full_to_sim, $@) $(call local_group_full_to_halo, $@) full

.SECONDEXPANSION:
$(local_group_split_plots): %: $$(call local_group_split_to_halo, %) $(nbody_local_group_plot_script) $(nbody_split_plot_script)
	$(python) $(nbody_local_group_plot_script) $(call local_group_split_to_sim, $@) $(call local_group_split_to_halo, $@) split


# Make the movies. We can't parametrize this all since there is no third expansion 
.SECONDEXPANSION:
$(call movies_all,n_body_refined): %: $$(call movie_to_plots,%,sim_to_refined_full)
	ffmpeg -framerate 2 -pattern_type glob -i '$(call movie_to_plot_dir,n_body_refined,$@)/n_body_refined*.png' -c:v h264 -pix_fmt yuv420p -y $@

SECONDEXPANSION:
$(call movies_all,n_body_split_refined): %: $$(call movie_to_plots,%,sim_to_refined_split)
	ffmpeg -framerate 2 -pattern_type glob -i '$(call movie_to_plot_dir,n_body_split_refined,$@)/n_body_split_refined*.png' -c:v h264 -pix_fmt yuv420p -y $@

.SECONDEXPANSION:
$(call movies_all,n_body_local_group): %: $$(call movie_to_plots,%,sim_to_local_group_full)
	ffmpeg -framerate 2 -pattern_type glob -i '$(call movie_to_plot_dir,n_body_local_group,$@)/n_body_local_group*.png' -c:v h264 -pix_fmt yuv420p -y $@

SECONDEXPANSION:
$(call movies_all,n_body_split_local_group): %: $$(call movie_to_plots,%,sim_to_local_group_split)
	ffmpeg -framerate 2 -pattern_type glob -i '$(call movie_to_plot_dir,n_body_split_local_group,$@)/n_body_split_local_group*.png' -c:v h264 -pix_fmt yuv420p -y $@

# Get ready to make the merger trees by moving the files correctly. The 
# rockstar config files will be used as the sentinel here
.SECONDEXPANSION:
$(rockstar_cfgs): %: $$(call rockstar_cfg_to_halos,%) $(prep_for_merger_tree_script)
	$(python) $(prep_for_merger_tree_script) $@

# Make the consistent trees config files
.SECONDEXPANSION:
$(tree_cfgs): %: $$(call tree_cfg_to_rockstar_cfg,%)
	perl $(tree_config_script) $(call tree_cfg_to_rockstar_cfg,$@)

# Then build the merger trees
.SECONDEXPANSION:
$(trees): %: $$(call tree_to_tree_cfg,%)
	cd $(tree_dir) && perl do_merger_tree.pl $<

# Build the merger tree reading code
$(read_tree_exe): $(read_tree_src)
	cd $(read_tree_dir) && make

# Build the accretion history output files
.SECONDEXPANSION:
$(merger_sentinels): %: $$(call merger_to_tree,%) $(read_tree_exe)
	$(read_tree_exe) $(call merger_to_tree,$@) && touch $@

# Then use those to make the halo comparison plot
$(halo_growth_plot): $(merger_sentinels) $(halo_growth_comp_script)
	$(python)  $(halo_growth_comp_script)

# timing output
.PHONY: timing
timing: $(timings)

# Here I have all timing outputs dependent on all snapshots. This isn't right.
# But it's really hard to get the actual directory, since the timing scripts
# have many subdirectories with unknown names. This is techincaly wrong, but 
# costs essentially no time, so I keep it. 
$(timings): $(snapshots)
	$(timing_script) $(dir $@) > $@ 

$(dt_history_plots): $(snapshots) $(dt_history_script)
	$(python) $(dt_history_script) $(dir $@)

