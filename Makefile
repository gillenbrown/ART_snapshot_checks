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

# ------------------------------------------------------------------------------
#
#  Code locations - depend on the machine
# 
# ------------------------------------------------------------------------------
ifeq ($(machine),shangrila)
	tree_config_script = /u/home/gillenb/code/not_mine/rockstar/scripts/gen_merger_cfg.pl
	tree_dir = /u/home/gillenb/code/not_mine/consistent-trees
	halo_finding_script = ./run_rockstar.sh
	timing_script = /u/home/gillenb/code/mine/art_cluster/utils/scripts/parse_timing2.pl
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

# ------------------------------------------------------------------------------
#
#  Code locations that are relative to this file
# 
# ------------------------------------------------------------------------------
halo_finding_py_file = ./halo_finding_rockstar.py
rename_script = ./rename_halos.py
summary_nbody_script = ./summary_nbody.py
summary_metal_script = ./summary_metals.py
summary_vel_script = ./summary_velocity.py
sfh_plots_script = ./plots_sfh.py
read_tree_dir = ./read_tree
read_tree_exe = $(read_tree_dir)/halo_history
read_tree_src = $(read_tree_dir)/halo_history.c

# ------------------------------------------------------------------------------
#
#  Simulation outputs to run this on - depend on the machine
# 
# ------------------------------------------------------------------------------
ifeq ($(machine),shangrila)
	runs_home = /u/home/gillenb/art_runs/runs
	sim_dirs_nbody = $(runs_home)/shangrila/test_new_ic/outputs 
#                   $(runs_home)/shangrila/nbody/run/outputs/tl 
#                   $(runs_home)/shangrila/nbody/run/outputs/rj 
#                   $(runs_home)/shangrila/nbody/run/outputs/br_no_refine_1 
#                   $(runs_home)/shangrila/nbody/run/outputs/br_no_refine_2
	sim_dirs_hydro = $(runs_home)/shangrila/discrete_old/run/outputs \
	                 $(runs_home)/shangrila/test_music_256/run/outputs \
	                 $(runs_home)/shangrila/test_music_128/run/outputs \
	                 $(runs_home)/great_lakes/hydro_test/discrete_128/run/outputs/first \
	                 $(runs_home)/great_lakes/hydro_test/discrete_256/run/outputs/first 
#                   $(runs_home)/pleiades/hydro/intel_broadwell_discrete/run/outputs/third 
#                   $(runs_home)/pleiades/hydro/intel_broadwell/run/outputs/alpha_restrict 
#                   $(runs_home)/pleiades/hydro/intel_broadwell/run/outputs/no_virial
#                   $(runs_home)/shangrila/test_all_elts/run
#                   $(runs_home)/shangrila/test_mine/run
#                   $(runs_home)/shangrila/NBm_10SFE_tidal_writeout/run
endif

ifeq ($(machine),lou)                
	runs_home = /u/gbrown12/art_runs/runs
	sim_dirs_nbody = $(runs_home)/nbody/trim_ic_05/run/outputs/first \
	                 $(runs_home)/nbody/trim_ic_06/run/outputs/first \
	                 $(runs_home)/nbody/trim_ic_07/run/outputs/first \
	                 $(runs_home)/nbody/trim_ic_08/run/outputs/first \
	                 $(runs_home)/nbody/no_trim/run/outputs/first
#                   $(runs_home)/nbody/trim_ic/run/outputs/first/ 
#                   $(runs_home)/nbody/trim_ic/run/outputs/second/ 
#                   $(runs_home)/nbody/intel/run/outputs/tl_production 
#                   $(runs_home)/nbody/intel/run/outputs/br_production 
#                   $(runs_home)/nbody/intel/run/outputs/rj_production 
#                   $(runs_home)/nbody/intel/run/outputs/change_core 
#                   $(runs_home)/nbody/intel/run/outputs/br_1.1.28_pleiades_no_refine 
#                   $(runs_home)/nbody/intel/run/outputs/br_2.1.28_pleiades_no_refine 
#                   $(runs_home)/nbody/intel/run/outputs/br_4.1.28_pleiades_no_refine 
#                   $(runs_home)/nbody/intel/run/outputs/br_8.1.28_electra_no_refine 
#                   $(runs_home)/nbody/intel/run/outputs/br_8.1.28_pleiades_no_refine 
#                   $(runs_home)/nbody/intel/run/outputs/br_8.2.14_pleiades_no_refine 
	sim_dirs_hydro = #$(runs_home)/hydro/intel_broadwell/run/outputs/no_virial \
                   #$(runs_home)/hydro/intel_broadwell/run/outputs/alpha_restrict \
                   #$(runs_home)/hydro/intel_broadwell_discrete/run/outputs/third \
#                   $(runs_home)/hydro/timing_test/current_code_no_elts/run/outputs/first \
#                   $(runs_home)/hydro/timing_test/current_code_with_elts/run/outputs/first \
#                   $(runs_home)/hydro/timing_test/current_code_with_elts_continuous/run/outputs/first \
#                   $(runs_home)/hydro/timing_test/old_code/run/outputs/first
#                   $(runs_home)/hydro/intel_broadwell/run/outputs/tl_first 
#                   $(runs_home)/hydro/intel_broadwell/run/outputs/tl_first_restart 
#                   $(runs_home)/hydro/intel_broadwell/run/outputs/tl_second 
#                   $(runs_home)/hydro/intel_broadwell/run/outputs/tl_debug 
#                   $(runs_home)/hydro/intel_broadwell/run/outputs/detail 
#                   $(runs_home)/hydro/intel_broadwell/run/outputs/detail_cfl_restart 
#                   $(runs_home)/hydro/intel_broadwell/run/outputs/refinement_delete 
#                   $(runs_home)/hydro/intel_broadwell/run/outputs/detail_low_res 
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
sim_human_halos_dirs = $(foreach dir,$(sim_dirs),$(dir)/halos)
sim_checks_dirs = $(foreach dir,$(sim_dirs),$(dir)/checks)
sim_plots_dirs = $(foreach dir,$(sim_dirs),$(dir)/plots)
my_directories = $(sim_checks_dirs) $(sim_human_halos_dirs) $(sim_rockstar_halos_dirs) $(sim_plots_dirs)

# ------------------------------------------------------------------------------
#
#  List of all simulation outputs and their corresponding halo catalogs
# 
# ------------------------------------------------------------------------------
dir_to_sims = $(wildcard $(1)/*_a*.art)
snapshots = $(foreach dir,$(sim_out_dirs),$(call dir_to_sims,$(dir)))
snapshots_hydro = $(foreach dir,$(sim_out_dirs_hydro),$(wildcard $(dir)/*_a*.art))
# Parse the snapshot names into halo catalogs 
# replace the directory and suffix
sim_to_halo = $(subst .art,.0.bin,$(subst out/continuous,halos/halos,$(1)))
halos_catalogs = $(foreach snapshot,$(snapshots),$(call sim_to_halo,$(snapshot)))

# ------------------------------------------------------------------------------
#
#  Rockstar sentinel files. 
# 
# ------------------------------------------------------------------------------
# A sentinel file just shows that rockstar has run for a given set of 
# simulation outputs. I do it this way since we run ROCKSTAR on all outputs
# at once, and the sentinel file is an easier way to show that dependency
rockstar_sentinels = $(foreach dir,$(sim_rockstar_halos_dirs),$(dir)/sentinel.txt)
# Function to get all simulations that are prerequisites to this sentinel
sentinel_to_sims = $(wildcard $(subst rockstar_halos/sentinel.txt,out/,$(1))*_a*.art)
sentinel_to_out_dir = $(subst rockstar_halos/sentinel.txt,out/,$(1))
sentinel_to_rh_dir = $(subst rockstar_halos/sentinel.txt,rockstar_halos/,$(1))

# Then some complex functions to find the rockstar sentinel file for a given 
# halo catalog. This is hard and ugly since we have to mess around with
# wildcards, which require $(patsubst), which requires words, not a path.
empty :=
space := $(empty) $(empty)
path_to_words = $(subst /,$(space),$(1))
words_to_path = /$(subst $(space),/,$(1))
halo_words_to_rockstar_words = $(patsubst halos,rockstar_halos,$(1))
rockstar_words_to_sentinel_words = $(patsubst %.0.bin,sentinel.txt,$(1))
halo_words_to_sentinel_words = $(call halo_words_to_rockstar_words,$(call rockstar_words_to_sentinel_words,$(1)))
halo_to_sentinel = $(call words_to_path,$(call halo_words_to_sentinel_words,$(call path_to_words,$(1))))

# ------------------------------------------------------------------------------
#
#  Summary files - Nbody
# 
# ------------------------------------------------------------------------------
sim_to_summary_nbody = $(subst .art,.txt,$(subst out/continuous,checks/summary_nbody, $(1)))
summary_nbody_to_sim = $(subst .txt,.art,$(subst checks/summary_nbody,out/continuous, $(1)))
summary_nbody_to_halo = $(subst .txt,.0.bin,$(subst checks/summary_nbody,halos/halos, $(1)))
summaries_nbody = $(foreach snapshot,$(snapshots),$(call sim_to_summary_nbody,$(snapshot)))

# ------------------------------------------------------------------------------
#
#  Nbody movies
# 
# ------------------------------------------------------------------------------
movies_nbody = $(foreach dir,$(sim_dirs),$(dir)/plots/n_body.mp4)
movie_nbody_to_summary_nbody = $(call sim_to_summary_nbody,$(call dir_to_sims,$(subst plots/n_body.mp4,out,$(1))))
movie_nbody_to_plot_dir = $(subst /n_body.mp4,,$(1))

movies_grid = $(foreach dir,$(sim_dirs),$(dir)/plots/grid.mp4)
movie_grid_to_summary_nbody = $(call sim_to_summary_nbody,$(call dir_to_sims,$(subst plots/grid.mp4,out,$(1))))
movie_grid_to_plot_dir = $(subst /grid.mp4,,$(1))

# ------------------------------------------------------------------------------
#
#  Consistent trees
# 
# ------------------------------------------------------------------------------
# first are the config files
tree_cfgs = $(foreach dir,$(sim_rockstar_halos_dirs),$(dir)/outputs/merger_tree.cfg)
tree_cfg_to_rockstar_cfg = $(subst outputs/merger_tree.cfg,rockstar.cfg,$(1))
tree_cfg_to_sentinel = $(subst outputs/merger_tree.cfg,sentinel.txt,$(1))

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
#  Summary files - metals
# 
# ------------------------------------------------------------------------------
sim_to_summary_metal = $(subst .art,.txt,$(subst out/continuous,checks/summary_metals, $(1)))
summary_metal_to_sim = $(subst .txt,.art,$(subst checks/summary_metals,out/continuous, $(1)))
summaries_metal = $(foreach snapshot,$(snapshots_hydro),$(call sim_to_summary_metal,$(snapshot)))

# ------------------------------------------------------------------------------
#
#  Summary files - velocity
# 
# ------------------------------------------------------------------------------
sim_to_summary_vel = $(subst .art,.txt,$(subst out/continuous,checks/summary_velocity, $(1)))
summary_vel_to_sim = $(subst .txt,.art,$(subst checks/summary_velocity,out/continuous, $(1)))
summaries_vel = $(foreach snapshot,$(snapshots_hydro),$(call sim_to_summary_vel,$(snapshot)))

# ------------------------------------------------------------------------------
#
#  Plots - SFH
# 
# ------------------------------------------------------------------------------
# Here the SMHM relation plot is used as the sentinel
sim_to_smhm = $(subst .art,.png,$(subst out/continuous,plots/smhm, $(1)))
smhm_to_sim = $(subst .png,.art,$(subst plots/smhm,out/continuous, $(1)))
smhm_plots = $(foreach snapshot,$(snapshots_hydro),$(call sim_to_smhm,$(snapshot)))

# ------------------------------------------------------------------------------
#
#  timing output
# 
# ------------------------------------------------------------------------------
timings = $(foreach dir,$(sim_dirs),$(dir)/checks/timing_debug.txt)
timing_to_sims = $(call dir_to_sims,$(subst checks/timing_debug.txt,out,$(1)))
timing_to_dir = $(subst checks/timing_debug.txt,log,$(1))

# ------------------------------------------------------------------------------
#
#  Rules
# 
# ------------------------------------------------------------------------------
all: $(my_directories) $(summaries_nbody) $(movies_nbody) $(movies_grid) $(summaries_metal) $(summaries_vel) $(merger_sentinels) $(smhm_plots) $(timings)

.PHONY: clean
clean:
	rm -r $(my_directories)

# Make directories if they don't exist
.PHONY: dirs
dirs: $(my_directories)

$(my_directories):
	mkdir $@

# Rule to make the rockstar sentinel files
# Phony target to allow us to make only the rockstar halos, for example on 
# a PBS script where we don't want to waste computation on the serial aspects
# that come later
.PHONY: halos
halos: $(rockstar_sentinels)
# We run the script with parameters to the output directory and rockstar halo 
# directory
.SECONDEXPANSION:
$(rockstar_sentinels): %: $$(call sentinel_to_sims, %) 
	$(halo_finding_script) $(call sentinel_to_out_dir, $@) $(call sentinel_to_rh_dir, $@)

# Rule to rename the halo catalogs into something more user-friendly
.SECONDEXPANSION:
$(halos_catalogs): %: $(rename_script) $$(call halo_to_sentinel,%)
	python $(rename_script) $@

# Make the summary files for N-body
.SECONDEXPANSION:
$(summaries_nbody): %: $$(call summary_nbody_to_halo, %) $(summary_nbody_script)
	python $(summary_nbody_script) $(call summary_nbody_to_sim, $@) $(call summary_nbody_to_halo, $@) clobber silent

# Make the moves from the N-body plots
.SECONDEXPANSION:
$(movies_nbody): %: $$(call movie_nbody_to_summary_nbody, %)
	ffmpeg -framerate 2 -pattern_type glob -i '$(call movie_nbody_to_plot_dir,$@)/n_body_*.png' -c:v h264 -pix_fmt yuv420p -y $@

.SECONDEXPANSION:
$(movies_grid): %: $$(call movie_grid_to_summary_nbody, %)
	ffmpeg -framerate 2 -pattern_type glob -i '$(call movie_grid_to_plot_dir,$@)/grid_idxs_*.png' -c:v h264 -pix_fmt yuv420p -y $@

# Make the consistent trees config files
.SECONDEXPANSION:
$(tree_cfgs): %: $$(call tree_cfg_to_sentinel,%)
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

# Make the summary files for metals
.SECONDEXPANSION:
$(summaries_metal): %: $$(call summary_metal_to_sim,%) $(summary_metal_script)
	python $(summary_metal_script) $(call summary_metal_to_sim, $@) clobber silent

# Make the summary files for velocities
.SECONDEXPANSION:
$(summaries_vel): %: $$(call summary_vel_to_sim,%) $(summary_vel_script)
	python $(summary_vel_script) $(call summary_vel_to_sim, $@) clobber silent

# SFH plots
.SECONDEXPANSION:
$(smhm_plots): %: $$(call smhm_to_sim,%) $(sfh_plots_script)
	python $(sfh_plots_script) $(call smhm_to_sim, $@) $(call sim_to_halo, $(call smhm_to_sim, $@))

# timing output
.PHONY: timing
timing: $(timings)

.SECONDEXPANSION:
$(timings): %: $$(call timing_to_sims, %)
	$(timing_script) $(call timing_to_dir, $@) > $@

