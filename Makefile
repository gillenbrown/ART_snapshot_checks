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
ifneq (,$(findstring stampede2,$(hostname)))
	machine = stampede2
endif

# ------------------------------------------------------------------------------
#
#  Code locations - depend on the machine
# 
# ------------------------------------------------------------------------------
ifeq ($(machine),shangrila)
	tree_config_script = /u/home/gillenb/code/not_mine/rockstar/scripts/gen_merger_cfg.pl
	tree_dir = /u/home/gillenb/code/not_mine/consistent-trees
	timing_script = /u/home/gillenb/code/mine/cart/utils/scripts/parse_timing2.pl
endif
ifeq ($(machine),stampede2)
	tree_config_script = $(HOME)/code/rockstar-galaxies/scripts/gen_merger_cfg.pl
	tree_dir = $(HOME)/code/consistent-trees
	timing_script = $(HOME)/code/cart/utils/scripts/parse_timing2.pl
endif

# ------------------------------------------------------------------------------
#
#  Code locations that are relative to this file
# 
# ------------------------------------------------------------------------------
halo_finding_script = ./run_rockstar.py
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
age_plot_script = ./plot_age_spread.py
galaxies_comparison_script = ./plot_galaxy_comparisons.py
gal_readin_script = ./utils/load_galaxies.py
plot_utils_script = ./utils/plot_utils.py
dt_history_script = ./dt_history.py
cfl_script = ./cfl_violations.py
tidal_consolidation_script = ./tidal_consolidation.py
merger_prep_script = ./merger_prep.py
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
	sim_dirs_nbody = $(runs_home)/stampede2/rj_nbody/original_92.48mpc_level07/run \
	                 $(runs_home)/stampede2/rj_nbody/hybrid_46.24mpc_level08/run \
	                 $(runs_home)/stampede2/rj_nbody/hybrid_23.12mpc_level08/run 
	sim_dirs_hydro = $(runs_home)/shangrila/hui/sfe_10 \
	                 $(runs_home)/shangrila/hui/sfe_50 \
	                 $(runs_home)/shangrila/hui/sfe_100 \
 	                 $(runs_home)/shangrila/hui/sfe_200 \
	                 $(runs_home)/stampede2/production/tl_sfe001_hn20/run \
	                 $(runs_home)/stampede2/production/tl_sfe010_hn20/run \
	                 $(runs_home)/stampede2/production/tl_sfe100_hn20/run \
	                 $(runs_home)/stampede2/production/tl_sfe100_hn05/run \
	                 $(runs_home)/stampede2/production/tl_sfe100_hn00/run \
	                 $(runs_home)/stampede2/production/tl_sfe100_hn00_fboost3/run \
 	                 $(runs_home)/stampede2/production/rj_sfe010_hn20/run \
 	                 $(runs_home)/stampede2/production/rj_sfe100_hn20/run \
 	                 $(runs_home)/stampede2/old_ic_comparison_production_analog/continuous_hn00_novirial/run \
 	                 $(runs_home)/stampede2/old_ic_comparison_production_analog/continuoushui_hn00_novirial/run \
 	                 $(runs_home)/stampede2/old_ic_comparison_production_analog/continuouspopmcluster_hn00_novirial/run \
 	                 $(runs_home)/stampede2/old_ic_comparison_production_analog/continuoussnr_hn00_novirial/run \
 	                 $(runs_home)/stampede2/old_ic_comparison_production_analog/discrete_hn00_novirial/run \
 	                 $(runs_home)/stampede2/old_ic_comparison_production_analog/discrete_hn00_virial10/run \
 	                 $(runs_home)/stampede2/old_ic_comparison_production_analog/discrete_hn00_virial10_19/run \
 	                 $(runs_home)/stampede2/old_ic_comparison_production_analog/discrete_hn00_virial10_19_advect/run \
 	                 $(runs_home)/stampede2/old_ic_comparison_production_analog/discrete_hn00_virial10_19_old/run \
 	                 $(runs_home)/stampede2/old_ic_comparison_production_analog/discrete_hn00_virial10_advect/run \
 	                 $(runs_home)/stampede2/old_ic_comparison_production_analog/discrete_hn00_virial10_advect_nostars/run \
 	                 $(runs_home)/stampede2/old_ic_comparison_production_analog/discrete_hn00_virial10_elements/run \
 	                 $(runs_home)/stampede2/old_ic_comparison_production_analog/discrete_hn00_virial10_entropy/run \
 	                 $(runs_home)/stampede2/old_ic_comparison_production_analog/discrete_hn00_virial10_entropy_fboost1/run \
 	                 $(runs_home)/stampede2/old_ic_comparison_production_analog/discrete_hn00_virial10_entropy_fboost3/run \
 	                 $(runs_home)/stampede2/old_ic_comparison_production_analog/discrete_hn00_virial10_entropy_hybridagediff/run \
 	                 $(runs_home)/stampede2/old_ic_comparison_production_analog/discrete_hn00_virial10_entropy_molvadim/run \
 	                 $(runs_home)/stampede2/old_ic_comparison_production_analog/discrete_hn00_virial10_entropy_molvadim_fboost1/run \
 	                 $(runs_home)/stampede2/old_ic_comparison_production_analog/discrete_hn00_virial10_entropy_noagediff/run \
 	                 $(runs_home)/stampede2/old_ic_comparison_production_analog/discrete_hn00_virial10_entropy_nosync/run \
 	                 $(runs_home)/stampede2/old_ic_comparison_production_analog/discrete_hn00_virial10_fboost3/run \
 	                 $(runs_home)/stampede2/old_ic_comparison_production_analog/discrete_hn00_virial10_noadvoradia_nostars/run \
 	                 $(runs_home)/stampede2/old_ic_comparison_production_analog/discrete_hn00_virial10_nostars/run \
 	                 $(runs_home)/stampede2/old_ic_comparison_production_analog/discrete_hn00_virial10_noturb_adi/run \
 	                 $(runs_home)/stampede2/old_ic_comparison_production_analog/discrete_hn00_virial10_noturb_adv/run \
 	                 $(runs_home)/stampede2/old_ic_comparison_production_analog/discrete_hn20_virial10/run 
endif
ifeq ($(machine),stampede2)
	runs_home = $(SCRATCH)/art_runs/runs
	sim_dirs_nbody =
	sim_dirs_hydro = $(runs_home)/production/tl_sfe001_hn20/run \
	                 $(runs_home)/production/tl_sfe010_hn20/run \
	                 $(runs_home)/production/tl_sfe100_hn00/run \
	                 $(runs_home)/production/tl_sfe100_hn00_fboost3/run \
	                 $(runs_home)/production/tl_sfe100_hn05/run \
	                 $(runs_home)/production/tl_sfe100_hn20/run \
	                 $(runs_home)/production/rj_sfe010_hn20/run \
	                 $(runs_home)/production/rj_sfe100_hn20/run 
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
# tidal directories are only needed for production runs
prod_sim_dirs = $(filter $(runs_home)/stampede2/production%,$(sim_dirs))
sim_tidal_dirs = $(foreach dir,$(prod_sim_dirs),$(dir)/tidal)
my_directories = $(sim_checks_dirs) $(sim_human_halos_dirs) $(sim_rockstar_halos_dirs) $(sim_plots_dirs) $(sim_tidal_dirs) $(comparison_plots_dir)

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
#sim_to_halo = $(subst .art,.0.bin,$(subst out/continuous,halos/halos,$(1)))
#halos_catalogs = $(foreach snapshot,$(snapshots),$(call sim_to_halo,$(snapshot)))

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
sentinel_to_halos_dir = $(subst rockstar_halos/sentinel.txt,halos/,$(1))

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
galaxies_to_halo = $(subst .txt,.0.bin,$(subst checks/galaxy_summaries,halos/halos, $(1)))
galaxies = $(foreach snapshot,$(snapshots),$(call sim_to_galaxies,$(snapshot)))

# ------------------------------------------------------------------------------
#
#  comparison plots that have all sims on one plot. Only on shangrila!
# 
# ------------------------------------------------------------------------------
ifeq ($(machine),shangrila)
	sfh_sentinel = $(comparison_plots_dir)/sfh_sentinel.txt
	cimf_sentinel = $(comparison_plots_dir)/cimf_sentinel.txt
	halo_growth_plot = $(comparison_plots_dir)/halo_growth.png
	age_spread_sentinel = $(comparison_plots_dir)/age_spread_sentinel.txt
	# one sentinel file for the script that makes a bunch of plots
	galaxy_comparison_sentinel = $(comparison_plots_dir)/comparison_sentinel.txt
else
	sfh_sentinel =
	cimf_sentinel =
	halo_growth_plot =
	age_spread_sentinel =
	galaxy_comparison_sentinel =
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
# first are the config files. We need these any any machine, since they'll need
# to be copied over to be used before creating the merger trees
tree_cfgs = $(foreach dir,$(sim_rockstar_halos_dirs),$(dir)/outputs/merger_tree.cfg)
tree_cfg_to_rockstar_cfg = $(subst outputs/merger_tree.cfg,rockstar.cfg,$(1))
tree_cfg_to_sentinel = $(subst outputs/merger_tree.cfg,sentinel.txt,$(1))
# I also have some sentinel files for the prepping of the out.list files
tree_cfg_to_prep_sentinel = $(subst outputs/merger_tree.cfg,sentinel_prep.txt,$(1))
merger_prep_sentinels = $(foreach dir,$(sim_rockstar_halos_dirs),$(dir)/sentinel_prep.txt)
merger_prep_sentinel_to_halo_sentinel = $(subst sentinel_prep.txt,sentinel.txt,$(1))

# we only do the actual merger tree creation on shangrila
ifeq ($(machine),shangrila)
	# then the actual halo trees
	trees = $(foreach dir,$(sim_rockstar_halos_dirs),$(dir)/trees/tree_0_0_0.dat)
	tree_to_tree_cfg = $(subst trees/tree_0_0_0.dat,outputs/merger_tree.cfg,$(1))

	# ----------------------------------------------------------------------------
	#
	#  Parsing merger trees
	# 
	# ----------------------------------------------------------------------------
	# Here I again use a sentinel file, since there will be multiple files created
	# by the one C file:
	merger_sentinels = $(foreach dir,$(sim_checks_dirs),$(dir)/merger_sentinel.txt)
	# need to get the trees that the sim should parse
	merger_to_tree = $(subst checks/merger_sentinel.txt,rockstar_halos/trees/tree_0_0_0.dat,$(1))
else
	trees =
	merger_sentinels = 
endif

# ------------------------------------------------------------------------------
#
#  timing output and dt history
# 
# ------------------------------------------------------------------------------
# Here we assume all log files are located as follows. Each `log` directory
# has several subdirectories holding the runtime outputs for each run, starting
# with `runtime_`. The log directories are inside that directory. Any other
# directories will be ignored. We have to do this weird thing we we get all
# items in each directory, then get the parent of each of those, to ensure that
# what we got is really a directory.
log_dirs = $(foreach dir,$(sim_dirs),$(foreach item,$(wildcard $(dir)/log/runtime_*/*),$(dir $(item))log))
# sorting this removes any duplicates. Order doesn't matter to me anyway
timing_dirs = $(sort $(log_dirs))
timings = $(foreach t_dir,$(timing_dirs),$(t_dir)/timing_debug.txt)
dt_history_plots = $(foreach t_dir,$(timing_dirs),$(t_dir)/timestep_history.png)
cfl_plots = $(foreach t_dir,$(timing_dirs),$(t_dir)/cfl_cell_speeds.png)

# ------------------------------------------------------------------------------
#
#  tidal file consolidation
# 
# ------------------------------------------------------------------------------
# Take the tidal ouput files from each run and consolidate them together into one
# big file for each star.
# only do this on shangrila
ifeq ($(machine),shangrila)
	tidal_sentinels = $(foreach dir,$(sim_tidal_dirs),$(dir)/sentinel.txt)
else
	tidal_sentinels =
endif


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

# only make the movies on shangrila
ifeq ($(machine),shangrila)
	movies = $(call movies_all,n_body_refined) $(call movies_all,n_body_split_refined) $(call movies_all,n_body_local_group) $(call movies_all,n_body_split_local_group)
else
	movies =
endif

# ------------------------------------------------------------------------------
#
#  Rules
# 
# ------------------------------------------------------------------------------
all: $(my_directories) $(timings) $(sfh_sentinel) $(cimf_sentinel) $(age_spread_sentinel) $(halo_growth_plot) $(galaxy_comparison_sentinel) $(debugs) $(movies)
# $(dt_history_plots) $(cfl_plots)

.PHONY: clean
clean:
	rm -r $(sim_checks_dirs) $(sim_plots_dirs)

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

# Rule to make the rockstar sentinel files
# Phony target to allow us to make only the rockstar halos, for example on 
# a PBS script where we don't want to waste computation on the serial aspects
# that come later
.PHONY: halos
halos: $(rockstar_sentinels)

.PHONY: tidal
tidal: $(tidal_sentinels)

# Throughout I use Make's static pattern rules to parameterize over my different
# outputs
# We run the script with parameters to the output directory and rockstar halo 
# directory
.SECONDEXPANSION:
$(rockstar_sentinels): %: $$(call sentinel_to_sims, %) 
	python $(halo_finding_script) $(call sentinel_to_out_dir, $@) $(call sentinel_to_rh_dir, $@) $(call sentinel_to_halos_dir, $@) $(machine)

# Make the debug files
$(debugs): %: $(debug_script)
	python $(debug_script) $(call debug_to_sim, $@) clobber silent

# Make the summary files
$(galaxies): %: $$(call galaxies_to_halo, %) $(galaxies_script)
	python $(galaxies_script) $(call galaxies_to_sim, $@) $(call galaxies_to_halo, $@) clobber silent

# Make the CIMF plots. We could use &: instead of : to indicate a grouped
# target, but that required Make 4.3 or higher
$(sfh_sentinel): $(sfh_plots_script) $(gal_readin_script) $(plot_utils_script) $(snapshots_hydro) $(rockstar_sentinels)
	python $(sfh_plots_script) $(sfh_sentinel) $(sim_dirs_hydro)

# Make the CIMF plots. We could use &: instead of : to indicate a grouped
# target, but that required Make 4.3 or higher
$(cimf_sentinel): $(cimf_plots_script) $(gal_readin_script) $(plot_utils_script) $(snapshots_hydro) $(rockstar_sentinels)
	python $(cimf_plots_script) $(cimf_sentinel) $(sim_dirs_hydro)

# Make the age spread plots. 
$(age_spread_sentinel): $(age_plot_script) $(gal_readin_script) $(plot_utils_script) $(snapshots_hydro) $(rockstar_sentinels)
	python $(age_plot_script) $(age_spread_sentinel) $(sim_dirs_hydro)

# and the galaxy comparison plots
$(galaxy_comparison_sentinel): $(galaxies) $(galaxies_comparison_script) $(gal_readin_script) $(plot_utils_script)
	python $(galaxies_comparison_script) $(galaxy_comparison_sentinel) $(galaxies)

# Make the individual nbody plots - several examples of very similar things here
$(halo_1_full_plots): %: $(nbody_single_halo_plots_script) $(nbody_full_plot_script)
	python $(nbody_single_halo_plots_script) $(call halo_1_full_to_sim, $@) 1 full

$(halo_2_full_plots): %: $(nbody_single_halo_plots_script) $(nbody_full_plot_script)
	python $(nbody_single_halo_plots_script) $(call halo_2_full_to_sim, $@) 2 full

$(halo_1_split_plots): %: $(nbody_single_halo_plots_script) $(nbody_split_plot_script)
	python $(nbody_single_halo_plots_script) $(call halo_1_split_to_sim, $@) 1 split

$(halo_2_split_plots): %: $(nbody_single_halo_plots_script) $(nbody_split_plot_script)
	python $(nbody_single_halo_plots_script) $(call halo_2_split_to_sim, $@) 2 split

$(refined_full_plots): %: $(nbody_refined_plot_script) $(nbody_full_plot_script)
	python $(nbody_refined_plot_script) $(call refined_full_to_sim, $@) full

$(refined_split_plots): %: $(nbody_refined_plot_script) $(nbody_split_plot_script)
	python $(nbody_refined_plot_script) $(call refined_split_to_sim, $@) split

$(local_group_full_plots): %: $(nbody_local_group_plot_script) $(nbody_full_plot_script)
	python $(nbody_local_group_plot_script) $(call local_group_full_to_sim, $@) full

$(local_group_split_plots): %: $(nbody_local_group_plot_script) $(nbody_split_plot_script)
	python $(nbody_local_group_plot_script) $(call local_group_split_to_sim, $@) split


# Make the movies. We can't parametrize this all since there is no third expansion 
$(call movies_all,n_body_refined): %: $$(call movie_to_plots,%,sim_to_refined_full)
	ffmpeg -framerate 2 -pattern_type glob -i '$(call movie_to_plot_dir,n_body_refined,$@)/n_body_refined*.png' -c:v h264 -pix_fmt yuv420p -y $@

$(call movies_all,n_body_split_refined): %: $$(call movie_to_plots,%,sim_to_refined_split)
	ffmpeg -framerate 2 -pattern_type glob -i '$(call movie_to_plot_dir,n_body_split_refined,$@)/n_body_split_refined*.png' -c:v h264 -pix_fmt yuv420p -y $@

$(call movies_all,n_body_local_group): %: $$(call movie_to_plots,%,sim_to_local_group_full)
	ffmpeg -framerate 2 -pattern_type glob -i '$(call movie_to_plot_dir,n_body_local_group,$@)/n_body_local_group*.png' -c:v h264 -pix_fmt yuv420p -y $@

$(call movies_all,n_body_split_local_group): %: $$(call movie_to_plots,%,sim_to_local_group_split)
	ffmpeg -framerate 2 -pattern_type glob -i '$(call movie_to_plot_dir,n_body_split_local_group,$@)/n_body_split_local_group*.png' -c:v h264 -pix_fmt yuv420p -y $@

# and format the out.list files
$(merger_prep_sentinels): %: $$(call merger_prep_sentinel_to_halo_sentinel,%) $(merger_prep_script)
	python $(merger_prep_script) $@

# Make the consistent trees config files
$(tree_cfgs): %: $$(call tree_cfg_to_sentinel,%) $$(call tree_cfg_to_prep_sentinel,%)
	perl $(tree_config_script) $(call tree_cfg_to_rockstar_cfg,$@)

# Then build the merger trees
$(trees): %: $$(call tree_to_tree_cfg,%)
	cd $(tree_dir) && perl do_merger_tree.pl $<

# Build the merger tree reading code
$(read_tree_exe): $(read_tree_src)
	cd $(read_tree_dir) && make

# Build the accretion history output files
$(merger_sentinels): %: $$(call merger_to_tree,%) $(read_tree_exe)
	$(read_tree_exe) $(call merger_to_tree,$@) && touch $@

# Then use those to make the halo comparison plot
$(halo_growth_plot): $(merger_sentinels) $(halo_growth_comp_script) $(gal_readin_script)
	python $(halo_growth_comp_script) $@

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
	python $(dt_history_script) $(dir $@)

$(cfl_plots): $(snapshots) $(cfl_script)
	python $(cfl_script) $(dir $@)

# tidal output consolidation
$(tidal_sentinels): $(snapshots) $(tidal_consolidation_script)
	python $(tidal_consolidation_script) $@