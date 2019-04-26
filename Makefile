# ------------------------------------------------------------------------------
#
#  Code locations - depend on the machine
# 
# ------------------------------------------------------------------------------
tree_config_script_lou = /u/gbrown12/yt-conda/src/rockstar/scripts/gen_merger_cfg.pl
tree_config_script_shangrila = /u/home/gillenb/code/not_mine/rockstar/scripts/gen_merger_cfg.pl
tree_dir_lou = /u/gbrown12/code/consistent-trees/
tree_dir_shangrila = /u/home/gillenb/code/not_mine/consistent-trees/

tree_config_script = $(tree_config_script_shangrila)
tree_dir = $(tree_dir_shangrila)

# ------------------------------------------------------------------------------
#
#  Code locations that are relative to this file
# 
# ------------------------------------------------------------------------------
halo_finding_script = ./run_rockstar.sh
rename_script = ./rename_halos.py
summary_script = ./global_properties.py
read_tree_dir = ./read_tree
read_tree_exe = $(read_tree_dir)/halo_history
read_tree_src = $(read_tree_dir)/halo_history.c

# ------------------------------------------------------------------------------
#
#  Simulation outputs to run this on - depend on the machine
# 
# ------------------------------------------------------------------------------
runs_home_shangrila = /u/home/gillenb/art_runs/runs/
sim_dirs_shangrila = $(runs_home_shangrila)shangrila/test_mine/run \
                     $(runs_home_shangrila)shangrila/nbody/run/outputs/rj \
         			 $(runs_home_shangrila)shangrila/nbody/run/outputs/tl \
       			     $(runs_home_shangrila)shangrila/nbody/run/outputs/br_no_refine_1 \
         			 $(runs_home_shangrila)pleiades/nbody/intel/br_1.1.28_pleiades_no_refine \
        		     $(runs_home_shangrila)pleiades/nbody/intel/br_2.1.28_pleiades_no_refine \
         			 $(runs_home_shangrila)pleiades/nbody/intel/br_4.1.28_pleiades_no_refine \
         			 $(runs_home_shangrila)pleiades/nbody/intel/br_8.1.28_pleiades_no_refine \
         			 $(runs_home_shangrila)pleiades/nbody/intel/br_8.1.28_electra_no_refine \
         			 $(runs_home_shangrila)pleiades/nbody/intel/br_8.2.14_pleiades_no_refine \
         			 $(runs_home_shangrila)pleiades/nbody/intel/br_production \
         			 $(runs_home_shangrila)pleiades/nbody/intel/rj_production \
         			 $(runs_home_shangrila)pleiades/nbody/intel/tl_production 
runs_home_lou = /u/gbrown12/art_runs/runs/nbody/intel/run/outputs/
sim_dirs_lou = $(runs_home_lou)br_production \
          	   $(runs_home_lou)tl_production \
         	   $(runs_home_lou)rj_production 

runs_home = $(runs_home_shangrila)
sim_dirs = $(sim_dirs_shangrila)

# ------------------------------------------------------------------------------
#
#  Directories for each simulation, where outputs will be stored
# 
# ------------------------------------------------------------------------------
sim_out_dirs = $(foreach dir,$(sim_dirs),$(dir)/out)
sim_rockstar_halos_dirs = $(foreach dir,$(sim_dirs),$(dir)/rockstar_halos)
sim_human_halos_dirs = $(foreach dir,$(sim_dirs),$(dir)/halos)
sim_checks_dirs = $(foreach dir,$(sim_dirs),$(dir)/checks)
sim_plots_dirs = $(foreach dir,$(sim_dirs),$(dir)/plots)
all_directories = $(sim_checks_dirs) $(sim_human_halos_dirs) $(sim_rockstar_halos_dirs) $(sim_plots_dirs)

# ------------------------------------------------------------------------------
#
#  List of all simulation outputs and their corresponding halo catalogs
# 
# ------------------------------------------------------------------------------
snapshots = $(foreach dir,$(sim_out_dirs),$(wildcard $(dir)/*_a*.art))
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
#  Summary files
# 
# ------------------------------------------------------------------------------
sim_to_summary = $(subst .art,.txt,$(subst out/continuous,checks/summary, $(1)))
summary_to_sim = $(subst .txt,.art,$(subst checks/summary,out/continuous, $(1)))
summary_to_halo = $(subst .txt,.0.bin,$(subst checks/summary,halos/halos, $(1)))
summaries = $(foreach snapshot,$(snapshots),$(call sim_to_summary,$(snapshot)))

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
#  Rules
# 
# ------------------------------------------------------------------------------
all: $(all_directories) $(halos_catalogs) $(summaries) $(merger_sentinels)

# Make directories if they don't exist
$(all_directories):
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
	$(halo_finding_script) $(call sentinel_to_out_dir, $@) $(call sentinel_to_rh_dir, $@) remove

# Rule to rename the halo catalogs into something more user-friendly
.SECONDEXPANSION:
$(halos_catalogs): %: $(rename_script) $$(call halo_to_sentinel,%)
	python $(rename_script) $@

# Make the summary files
.SECONDEXPANSION:
$(summaries): %: $$(call summary_to_halo, %) $(summary_script)
	python $(summary_script) $(call summary_to_sim, $@) $(call summary_to_halo, $@) clobber silent

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
