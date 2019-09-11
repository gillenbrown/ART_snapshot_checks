# ------------------------------------------------------------------------------
#
#  Flag to figure out what machine we're on
# 
# ------------------------------------------------------------------------------
# This tells us which directories to look for simulations in, which will 
# be used later. We'll use the path of the current directory to figure this out.
# My username is different on the two machines, so that part of the path will
# be unique
username_shangrila = gillenb
username_lou = gbrown12
current_dir = $(shell pwd)
# findstring returns the matching part of the string. If it's not empty when
# we try to find the shangrila username, we know we're on shangrila
ifneq (,$(findstring $(username_shangrila),$(current_dir)))
	machine = shangrila
endif
ifneq (,$(findstring $(username_lou),$(current_dir)))
	machine = lou
endif

# ------------------------------------------------------------------------------
#
#  Code locations - depend on the machine
# 
# ------------------------------------------------------------------------------
ifeq ($(machine),shangrila)
	tree_config_script = /u/home/gillenb/code/not_mine/rockstar/scripts/gen_merger_cfg.pl
	tree_dir = /u/home/gillenb/code/not_mine/consistent-trees/
	halo_finding_script = ./run_rockstar.sh
endif
ifeq ($(machine),lou)
	tree_config_script = /u/gbrown12/yt-conda/src/rockstar/scripts/gen_merger_cfg.pl
	tree_dir = /u/gbrown12/code/consistent-trees/
	halo_finding_script = ./run_rockstar_ldan.sh
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
	runs_home = /u/home/gillenb/art_runs/runs/
	sim_dirs_nbody = $(runs_home)shangrila/nbody/run/outputs/tl 
#                   $(runs_home)shangrila/nbody/run/outputs/rj 
#                   $(runs_home)shangrila/nbody/run/outputs/br_no_refine_1 
#                   $(runs_home)shangrila/nbody/run/outputs/br_no_refine_2
	sim_dirs_hydro = $(runs_home)shangrila/test_discrete/run/outputs \
                   $(runs_home)pleiades/hydro/intel_broadwell_discrete/run/outputs/second \
                   $(runs_home)pleiades/hydro/intel_broadwell/run/outputs/alpha_restrict \
                   $(runs_home)pleiades/hydro/intel_broadwell/run/outputs/no_virial
#                   $(runs_home)shangrila/test_all_elts/run
#                   $(runs_home)shangrila/test_mine/run
#                   $(runs_home)shangrila/NBm_10SFE_tidal_writeout/run
#                   $(runs_home)shangrila/test_mine_music/run/outputs 
endif

ifeq ($(machine),lou)                
	runs_home = /u/gbrown12/art_runs/runs/
	sim_dirs_nbody = $(runs_home)nbody/intel/run/outputs/tl_production 
#                   $(runs_home)nbody/intel/run/outputs/br_production 
#                   $(runs_home)nbody/intel/run/outputs/rj_production 
#                   $(runs_home)nbody/intel/run/outputs/change_core 
#                   $(runs_home)nbody/intel/run/outputs/br_1.1.28_pleiades_no_refine 
#                   $(runs_home)nbody/intel/run/outputs/br_2.1.28_pleiades_no_refine 
#                   $(runs_home)nbody/intel/run/outputs/br_4.1.28_pleiades_no_refine 
#                   $(runs_home)nbody/intel/run/outputs/br_8.1.28_electra_no_refine 
#                   $(runs_home)nbody/intel/run/outputs/br_8.1.28_pleiades_no_refine 
#                   $(runs_home)nbody/intel/run/outputs/br_8.2.14_pleiades_no_refine 
	sim_dirs_hydro = $(runs_home)hydro/intel_broadwell/run/outputs/no_virial \
                   $(runs_home)hydro/intel_broadwell/run/outputs/alpha_restrict \
                   $(runs_home)hydro/intel_broadwell_discrete/run/outputs/second
#                   $(runs_home)hydro/intel_broadwell_debug_timestep/run/outputs/detail_dt
#                   $(runs_home)hydro/intel_broadwell/run/outputs/tl_first 
#                   $(runs_home)hydro/intel_broadwell/run/outputs/tl_first_restart 
#                   $(runs_home)hydro/intel_broadwell/run/outputs/tl_second 
#                   $(runs_home)hydro/intel_broadwell/run/outputs/tl_debug 
#                   $(runs_home)hydro/intel_broadwell/run/outputs/detail 
#                   $(runs_home)hydro/intel_broadwell/run/outputs/detail_cfl_restart 
#                   $(runs_home)hydro/intel_broadwell/run/outputs/refinement_delete 
#                   $(runs_home)hydro/intel_broadwell_no_elts/run/outputs/no_elts 
#                   $(runs_home)hydro/pgi_broadwell/run/outputs/detail 
#                   $(runs_home)hydro/pgi_broadwell/run/outputs/tl_first 
#                   $(runs_home)hydro/pgi_broadwell/run/outputs/tl_second 
#                   $(runs_home)hydro/intel_skylake/tl_first 
#                   $(runs_home)hydro/intel_broadwell/run/outputs/detail_low_res 
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
snapshots = $(foreach dir,$(sim_out_dirs),$(wildcard $(dir)/*_a*.art))
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
#  Rules
# 
# ------------------------------------------------------------------------------
all: $(my_directories) $(summaries_nbody) $(summaries_metal) $(summaries_vel) $(merger_sentinels) $(smhm_plots)

.PHONY: clean
clean:
	rm -r $(my_directories)

# Make directories if they don't exist
$(my_directories):
	mkdir $@

# Rule to make the rockstar sentinel files
# Phony target to allow us to make only the rockstar halos, for example on 
# a PBS script where we don't want to waste computation on the serial aspects
# that come later
.PHONY: halos
halos: $(rockstar_sentinels)
# We run the script with parameters to the output directory and rockstar halo 
# directory, plus the remove keyword to replace previous halo catalogs
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

