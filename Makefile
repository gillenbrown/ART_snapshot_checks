# ------------------------------------------------------------------------------
#
#  Code locations
# 
# ------------------------------------------------------------------------------
halo_finding_script = /u/home/gillenb/code/mine/halo_scripts/run_rockstar.sh
rename_script = /u/home/gillenb/code/mine/ART_snapshot_checks/rename_halos.py
summary_script = /u/home/gillenb/code/mine/ART_snapshot_checks/global_properties.py
tree_config_script = /u/home/gillenb/code/not_mine/rockstar/scripts/gen_merger_cfg.pl

# ------------------------------------------------------------------------------
#
#  Simulation outputs to run this on
# 
# ------------------------------------------------------------------------------
runs_home = /u/home/gillenb/art_runs/runs/
sim_dirs = $(runs_home)shangrila/test_mine/run \
           $(runs_home)shangrila/nbody/run/outputs/rj \
           $(runs_home)shangrila/nbody/run/outputs/tl \
           $(runs_home)shangrila/nbody/run/outputs/br_no_refine_1 \
           $(runs_home)pleiades/nbody/intel/br_1.1.28_pleiades_no_refine \
           $(runs_home)pleiades/nbody/intel/br_2.1.28_pleiades_no_refine \
           $(runs_home)pleiades/nbody/intel/br_4.1.28_pleiades_no_refine \
           $(runs_home)pleiades/nbody/intel/br_8.1.28_pleiades_no_refine \
           $(runs_home)pleiades/nbody/intel/br_8.1.28_electra_no_refine \
           $(runs_home)pleiades/nbody/intel/br_8.2.14_pleiades_no_refine \
           $(runs_home)pleiades/nbody/intel/br_production \
           $(runs_home)pleiades/nbody/intel/rj_production \
           $(runs_home)pleiades/nbody/intel/tl_production 

# ------------------------------------------------------------------------------
#
#  Directories for each simulation, where outputs will be stored
# 
# ------------------------------------------------------------------------------
sim_out_dirs = $(foreach dir, $(sim_dirs), $(dir)/out/)
sim_rockstar_halos_dirs = $(foreach dir, $(sim_dirs), $(dir)/rockstar_halos/)
sim_human_halos_dirs = $(foreach dir, $(sim_dirs), $(dir)/halos/)
sim_checks_dirs = $(foreach dir, $(sim_dirs), $(dir)/checks/)
sim_plots_dirs = $(foreach dir, $(sim_dirs), $(dir)/plots/)
all_directories = $(sim_checks_dirs) $(sim_human_halos_dirs) $(sim_rockstar_halos_dirs) $(sim_plots_dirs)

# ------------------------------------------------------------------------------
#
#  List of all simulation outputs and their corresponding halo catalogs
# 
# ------------------------------------------------------------------------------
snapshots = $(foreach dir, $(sim_out_dirs), $(wildcard $(dir)*_a*.art))
# Parse the snapshot names into halo catalogs 
# replace the directory and suffix
sim_to_halo = $(subst .art,.0.bin, $(subst out/continuous,halos/halos, $(1)))
halos_catalogs = $(foreach snapshot, $(snapshots), $(call sim_to_halo, $(snapshot)))

# ------------------------------------------------------------------------------
#
#  Rockstar sentinel files. 
# 
# ------------------------------------------------------------------------------
# A sentinel file just shows that rockstar has run for a given set of 
# simulation outputs. I do it this way since we run ROCKSTAR on all outputs
# at once, and the sentinel file is an easier way to show that dependency
rockstar_sentinels = $(foreach dir, $(sim_rockstar_halos_dirs), $(dir)sentinel.txt)
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
sim_to_summary = $(subst .art,.txt, $(subst out/continuous,checks/summary, $(1)))
summary_to_sim = $(subst .txt,.art, $(subst checks/summary,out/continuous, $(1)))
summary_to_halo = $(subst .txt,.0.bin, $(subst checks/summary,halos/halos, $(1)))
summaries = $(foreach snapshot, $(snapshots), $(call sim_to_summary, $(snapshot)))

# ------------------------------------------------------------------------------
#
#  Consistent trees
# 
# ------------------------------------------------------------------------------
tree_cfgs = $(foreach dir, $(sim_rockstar_halos_dirs), $(dir)/outputs/merger_tree.cfg)
tree_cfg_to_rockstar_cfg = $(subst outputs/merger_tree.cfg,rockstar.cfg,$(1))
tree_cfg_to_sentinel = $(subst outputs/merger_tree.cfg,sentinel.txt,$(1))

# ------------------------------------------------------------------------------
#
#  Rules
# 
# ------------------------------------------------------------------------------
all: $(all_directories) $(summaries)

# Make directories if they don't exist
$(all_directories):
	mkdir $@

# Rule to make the rockstar sentinel files
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
