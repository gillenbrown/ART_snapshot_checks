import os

makefile = open("makefile_generated", "w")

current_dir = os.path.dirname(os.path.abspath(__file__))

# https://www.gnu.org/software/make/manual/html_node/Automatic-Variables.html
this_target = "$@"

# ------------------------------------------------------------------------------
#  functions to write to the makefile
# ------------------------------------------------------------------------------
def write(text=""):
    makefile.write(text)
    makefile.write("\n")


def write_list_one_per_line(name, items):
    # check if there's a single value
    if len(items) == 1:
        write("{} = {}".format(name, items[0]))
        write()
        return
    # if there are multiple we need to escape it with a slash
    write("{} = {} \\".format(name, items[0]))
    # have to pad the beginning to make it look nice
    empty = " " * (len(name) + 3)  # 3 for the =
    for i in range(1, len(items)):
        makefile.write(empty)
        makefile.write(items[i])
        if i != len(items) - 1:
            makefile.write(" \\")
        makefile.write("\n")

    write()


def write_rule(target, prerequisite, recipe):
    if prerequisite is None:
        prerequisite = ""
    write("{}: {}".format(target, prerequisite))
    write("\t{}".format(recipe))
    write()


def list_to_str(items):
    return " ".join(items)


# ------------------------------------------------------------------------------
#
# functions to navigate the directories
#
# ------------------------------------------------------------------------------
def join(home_dir=current_dir, *args):
    return os.path.join(home_dir, *args)


def base_dir_to_out_dir(base_dir):
    return join(base_dir, "out")


def base_dir_to_rockstar_dir(base_dir):
    return join(base_dir, "rockstar_halos")


def base_dir_to_halos_dir(base_dir):
    return join(base_dir, "halos")


def base_dir_to_checks_dir(base_dir):
    return join(base_dir, "checks")


def base_dir_to_plots_dir(base_dir):
    return join(base_dir, "plots")


# ------------------------------------------------------------------------------
#  go back and forth through various directories. These are ugly but useful
# ------------------------------------------------------------------------------
def base_to_sims(base_dir):
    out_dir = base_dir_to_out_dir(base_dir)
    return [
        join(out_dir, item) for item in os.listdir(out_dir) if item.endswith(".art")
    ]


def sim_to_halo(sim_path):
    halo_path = sim_path.replace("out/continuous", "halos/halos")
    return halo_path.replace(".art", ".0.bin")


def halo_to_sentinel(halo_path):
    return halo_path.replace("halos/halos", "rockstar_halos/sentinel.txt")


# ------------------------------------------------------------------------------
#
#  Flag to figure out what machine we're on
#
# ------------------------------------------------------------------------------
# This tells us which directories to look for simulations in, which will
# be used later. We'll use the hostname of the machine
hostname = os.uname()[1]
if "shangrila" in hostname:
    machine = "shangrila"
elif "ldan" in hostname or "lfe" in hostname:
    machine = "lou"
elif "gl-login" in hostname:
    machine = "great_lakes"
elif "gillenb" in hostname:
    machine = "macbook"
else:
    raise ValueError("hostname not recognized")

# ------------------------------------------------------------------------------
#
#  Code locations - depend on the machine, will be used later
#
# ------------------------------------------------------------------------------
if machine == "shangrila":
    tree_config_script = (
        "/u/home/gillenb/code/not_mine/rockstar/scripts/gen_merger_cfg.pl"
    )
    tree_dir = "/u/home/gillenb/code/not_mine/consistent-trees"
    halo_finding_script = "./run_rockstar.sh"
    timing_script = (
        "/u/home/gillenb/code/mine/art_cluster/utils/scripts/parse_timing2.pl"
    )
elif machine == "lou":
    tree_config_script = "/u/gbrown12/code/rockstar/scripts/gen_merger_cfg.pl"
    tree_dir = "/u/gbrown12/code/consistent-trees"
    halo_finding_script = "./run_rockstar_ldan.sh"
    timing_script = "/u/gbrown12/code/parse_timing2.pl"
elif machine == "great lakes":
    tree_config_script = "/home/gillenb/code/rockstar/scripts/gen_merger_cfg.pl"
    tree_dir = "/home/gillenb/code/consistent-trees"
    halo_finding_script = "./run_rockstar_gl.sh"
    timing_script = "/home/gillenb/code/art_cluster/utils/scripts/parse_timing2.pl"
elif machine == "macbook":
    tree_config_script = "/home/gillenb/code/rockstar/scripts/gen_merger_cfg.pl"
    tree_dir = "/home/gillenb/code/consistent-trees"
    halo_finding_script = "./run_rockstar_gl.sh"
    timing_script = "/home/gillenb/code/art_cluster/utils/scripts/parse_timing2.pl"
else:
    raise ValueError("hostname not recognized")

# ------------------------------------------------------------------------------
#
#  Simulation outputs to run this on - depend on the machine
#
# ------------------------------------------------------------------------------
if machine == "shangrila":
    runs_home = "/u/home/gillenb/art_runs/runs"
    sim_dirs_nbody = []
    sim_dirs_hydro = [
        join(runs_home, "shangrila/hui/sfe_100"),
        join(runs_home, "shangrila/old_ic_comparison/default_test/run"),
        join(runs_home, "stampede2/ic_timing_tests/original_50_128"),
        join(runs_home, "stampede2/ic_timing_tests/trim_12_256"),
        join(runs_home, "stampede2/ic_timing_tests/trim_25_256"),
    ]

elif machine == "lou":
    runs_home = "/u/gbrown12/art_runs/runs"
    sim_dirs_nbody = [
        join(runs_home, "nbody/new_ic_trim_12.5mpc/root_05/run/outputs/fix_vel"),
        join(runs_home, "nbody/new_ic_trim_12.5mpc/root_06/run/outputs/fix_vel"),
        join(runs_home, "nbody/new_ic_trim_12.5mpc/root_07/run/outputs/fix_vel"),
        join(runs_home, "nbody/new_ic_trim_12.5mpc/root_08/run/outputs/fix_vel"),
        join(runs_home, "nbody/new_ic_trim_12.5mpc/root_05/run/outputs/vel_offset"),
        join(runs_home, "nbody/new_ic_trim_12.5mpc/root_06/run/outputs/vel_offset"),
        join(runs_home, "nbody/new_ic_trim_12.5mpc/root_07/run/outputs/vel_offset"),
        join(runs_home, "nbody/new_ic_trim_12.5mpc/root_08/run/outputs/vel_offset"),
        join(runs_home, "nbody/new_ic_trim_25mpc/root_06/run/outputs/vel_offset"),
        join(runs_home, "nbody/new_ic_trim_25mpc/root_07/run/outputs/vel_offset"),
        join(runs_home, "nbody/new_ic_trim_25mpc/root_08/run/outputs/vel_offset"),
        join(runs_home, "nbody/new_ic_50mpc/root_07/run/outputs/fix_vel"),
    ]
    sim_dirs_hydro = [
        join(runs_home, "hydro/new_ic_50mpc/root_07/run/outputs/first"),
        join(runs_home, "hydro/new_ic_trim_12.5mpc/root_05/run/outputs/first"),
        join(runs_home, "hydro/new_ic_trim_12.5mpc/root_06/run/outputs/first"),
        join(runs_home, "hydro/new_ic_trim_12.5mpc/root_07/run/outputs/first"),
        join(runs_home, "hydro/new_ic_trim_12.5mpc/root_08/run/outputs/first"),
        join(runs_home, "hydro/old_ic/default/run/outputs/run_1"),
        join(runs_home, "hydro/old_ic/no_hn/run/outputs/run_1"),
    ]

elif machine == "great lakes":
    runs_home = "/home/gillenb/art_runs/runs"
    sim_dirs_nbody = []
    sim_dirs_hydro = [
        join(runs_home, "hydro_test/discrete_128/run/outputs/first"),
        join(runs_home, "hydro_test/discrete_256/run/outputs/first"),
        join(runs_home, "hydro_test/scaling_discrete_256/run/outputs/1_cores"),
        join(runs_home, "hydro_test/scaling_discrete_256/run/outputs/2_cores"),
        join(runs_home, "hydro_test/scaling_discrete_256/run/outputs/4_cores"),
        join(runs_home, "hydro_test/scaling_discrete_256/run/outputs/6_cores"),
        join(runs_home, "hydro_test/scaling_discrete_256/run/outputs/8_cores"),
    ]

elif machine == "macbook":
    sim_dirs_nbody = []
    sim_dirs_hydro = [
        "/Users/gillenb/google_drive/research/agb/agb_wrapper/testing/stdout_tests/art_dataset"
    ]
else:
    raise ValueError("hostname not recognized")

# combine the N-Body and Hydro into one big list
sim_dirs = sim_dirs_nbody + sim_dirs_hydro

# ------------------------------------------------------------------------------
#
#  Directories for each simulation, where outputs will be stored
#
# ------------------------------------------------------------------------------
sim_out_dirs = [base_dir_to_out_dir(d) for d in sim_dirs]
sim_out_dirs_hydro = [base_dir_to_out_dir(d) for d in sim_dirs_hydro]
sim_rockstar_halos_dirs = [base_dir_to_rockstar_dir(d) for d in sim_dirs]
sim_human_halos_dirs = [base_dir_to_halos_dir(d) for d in sim_dirs]
sim_checks_dirs = [base_dir_to_checks_dir(d) for d in sim_dirs]
sim_plots_dirs = [base_dir_to_plots_dir(d) for d in sim_dirs]
comparison_plots_dir = join(current_dir, "comparison_plots")
most_directories = (
    sim_human_halos_dirs + sim_checks_dirs + sim_plots_dirs + [comparison_plots_dir]
)

# targets to create all these directories
write_list_one_per_line("most_directories", most_directories)
write_list_one_per_line("rockstar_directories", sim_rockstar_halos_dirs)

write("all_directories = $(most_directories) $(rockstar_directories)")
write()

# have a phony name to create the directories\
write(".PHONY: dirs")
write("dirs: $(all_directories)")
write()
write_rule("$(all_directories)", None, "mkdir $@")

# ------------------------------------------------------------------------------
#
#  The clean target removes these directories
#
# ------------------------------------------------------------------------------
write(".PHONY: clean")
write_rule("clean", None, "rm -r $(most_directories)")

write(".PHONY: clean_all")
write_rule("clean_all", None, "rm -r $(all_directories)")

# ------------------------------------------------------------------------------
#
#  Halo Catalogs
#
# ------------------------------------------------------------------------------
halo_finding_py_file = join("halo_finding_rockstar.py")
rename_script = join("rename_halos.py")

# Here I just use one of them as the indicator, since I don't know how many
# .bin files there are per simulation

# Rule to make the rockstar sentinel files
# Phony target to allow us to make only the rockstar halos, for example on
# a PBS script where we don't want to waste computation on the serial aspects
# that come later
write(".PHONY: halos")
write("halos: $(rockstar_sentinels)")
write()

# # Rule to rename the halo catalogs into something more user-friendly
# .SECONDEXPANSION:
# $(halos_catalogs): %: $(rename_script) $$(call halo_to_sentinel,%)
#     python $(rename_script) $@


# first make the rule for the sentinels
for base_dir in sim_dirs:
    sentinel = join(base_dir_to_rockstar_dir(base_dir), "sentinel.txt")
    prereq = list_to_str(base_to_sims(base_dir))

    sentinel_rule = (
        f"{halo_finding_script} "
        f"{base_dir_to_out_dir(base_dir)} "
        f"{base_dir_to_rockstar_dir(base_dir)}"
    )
    write_rule(sentinel, prereq, sentinel_rule)

# then the rule to turn those into human readable halo catalogs
for base_dir in sim_dirs:
    sentinel = join(base_dir_to_rockstar_dir(base_dir), "sentinel.txt")
    for sim in base_to_sims(base_dir):
        hc = sim_to_halo(sim)

        prereqs = f"{rename_script} {sentinel}"
        rule = f"python {rename_script} {this_target}"

        write_rule(hc, prereqs, rule)


# # ------------------------------------------------------------------------------
# #
# #  Rockstar sentinel files.
# #
# # ------------------------------------------------------------------------------

# # ------------------------------------------------------------------------------
# #
# #  Debug output files
# #
# # ------------------------------------------------------------------------------
# sim_to_debug = $(subst .art,.txt,$(subst out/continuous,checks/debug, $(1)))
# debug_to_sim = $(subst .txt,.art,$(subst checks/debug,out/continuous, $(1)))
# debugs = $(foreach snapshot,$(snapshots),$(call sim_to_debug,$(snapshot)))

# # ------------------------------------------------------------------------------
# #
# #  Galaxies output files
# #
# # ------------------------------------------------------------------------------
# sim_to_galaxies = $(subst .art,.txt,$(subst out/continuous,checks/galaxies, $(1)))
# galaxies_to_sim = $(subst .txt,.art,$(subst checks/galaxies,out/continuous, $(1)))
# galaxies_to_halo = $(subst .txt,.0.bin,$(subst checks/galaxies,halos/halos, $(1)))
# galaxies = $(foreach snapshot,$(snapshots),$(call sim_to_galaxies,$(snapshot)))

# # ------------------------------------------------------------------------------
# #
# #  SFH plots
# #
# # ------------------------------------------------------------------------------
# sfh_plots = $(comparison_plots_dir)/mass_growth_comparison.png \
#             $(comparison_plots_dir)/sfh_comparison.png

# # ------------------------------------------------------------------------------
# #
# #  Nbody plots
# #
# # ------------------------------------------------------------------------------
# sim_to_halo_1_full = $(subst .art,.png,$(subst out/continuous,plots/n_body_halo_rank_1, $(1)))
# sim_to_halo_2_full = $(subst .art,.png,$(subst out/continuous,plots/n_body_halo_rank_2, $(1)))
# sim_to_halo_1_split = $(subst .art,.png,$(subst out/continuous,plots/n_body_split_halo_rank_1, $(1)))
# sim_to_halo_2_split = $(subst .art,.png,$(subst out/continuous,plots/n_body_split_halo_rank_2, $(1)))
# sim_to_refined_full  = $(subst .art,.png,$(subst out/continuous,plots/n_body_refined, $(1)))
# sim_to_refined_spilt = $(subst .art,.png,$(subst out/continuous,plots/n_body_split_refined, $(1)))
# sim_to_local_group_full = $(subst .art,.png,$(subst out/continuous,plots/n_body_local_group, $(1)))
# sim_to_local_group_spilt = $(subst .art,.png,$(subst out/continuous,plots/n_body_split_local_group, $(1)))

# halo_1_full_to_sim = $(subst .png,.art,$(subst plots/n_body_halo_rank_1,out/continuous, $(1)))
# halo_2_full_to_sim = $(subst .png,.art,$(subst plots/n_body_halo_rank_2,out/continuous, $(1)))
# halo_1_split_to_sim = $(subst .png,.art,$(subst plots/n_body_split_halo_rank_1,out/continuous, $(1)))
# halo_2_split_to_sim = $(subst .png,.art,$(subst plots/n_body_split_halo_rank_2,out/continuous, $(1)))
# refined_full_to_sim = $(subst .png,.art,$(subst plots/n_body_refined,out/continuous, $(1)))
# refined_split_to_sim = $(subst .png,.art,$(subst plots/n_body_split_refined,out/continuous, $(1)))
# local_group_full_to_sim = $(subst .png,.art,$(subst plots/n_body_local_group,out/continuous, $(1)))
# local_group_split_to_sim = $(subst .png,.art,$(subst plots/n_body_split_local_group,out/continuous, $(1)))

# halo_1_full_to_halo = $(subst .png,.0.bin,$(subst plots/n_body_halo_rank_1,halos/halos, $(1)))
# halo_2_full_to_halo = $(subst .png,.0.bin,$(subst plots/n_body_halo_rank_2,halos/halos, $(1)))
# halo_1_split_to_halo = $(subst .png,.0.bin,$(subst plots/n_body_split_halo_rank_1,halos/halos, $(1)))
# halo_2_split_to_halo = $(subst .png,.0.bin,$(subst plots/n_body_split_halo_rank_2,halos/halos, $(1)))
# refined_full_to_halo = $(subst .png,.0.bin,$(subst plots/n_body_refined,halos/halos, $(1)))
# refined_split_to_halo = $(subst .png,.0.bin,$(subst plots/n_body_split_refined,halos/halos, $(1)))
# local_group_full_to_halo = $(subst .png,.0.bin,$(subst plots/n_body_local_group,halos/halos, $(1)))
# local_group_split_to_halo = $(subst .png,.0.bin,$(subst plots/n_body_split_local_group,halos/halos, $(1)))

# halo_1_full_plots = $(foreach snapshot,$(snapshots),$(call sim_to_halo_1_full,$(snapshot)))
# halo_2_full_plots = $(foreach snapshot,$(snapshots),$(call sim_to_halo_2_full,$(snapshot)))
# halo_1_split_plots = $(foreach snapshot,$(snapshots),$(call sim_to_halo_1_split,$(snapshot)))
# halo_2_split_plots = $(foreach snapshot,$(snapshots),$(call sim_to_halo_2_split,$(snapshot)))
# refined_full_plots = $(foreach snapshot,$(snapshots),$(call sim_to_refined_full,$(snapshot)))
# refined_split_plots = $(foreach snapshot,$(snapshots),$(call sim_to_refined_spilt,$(snapshot)))
# local_group_full_plots = $(foreach snapshot,$(snapshots),$(call sim_to_local_group_full,$(snapshot)))
# local_group_split_plots = $(foreach snapshot,$(snapshots),$(call sim_to_local_group_spilt,$(snapshot)))

# nbody_plots = $(halo_1_full_plots) $(halo_2_full_plots) $(halo_1_split_plots) $(halo_2_split_plots) $(refined_full_plots) $(refined_split_plots) $(local_group_full_plots) $(local_group_split_plots)

# # ------------------------------------------------------------------------------
# #
# #  Consistent trees
# #
# # ------------------------------------------------------------------------------
# # first are the config files
# tree_cfgs = $(foreach dir,$(sim_rockstar_halos_dirs),$(dir)/outputs/merger_tree.cfg)
# tree_cfg_to_rockstar_cfg = $(subst outputs/merger_tree.cfg,rockstar.cfg,$(1))
# tree_cfg_to_sentinel = $(subst outputs/merger_tree.cfg,sentinel.txt,$(1))

# # then the actual halo trees
# trees = $(foreach dir,$(sim_rockstar_halos_dirs),$(dir)/trees/tree_0_0_0.dat)
# tree_to_tree_cfg = $(subst trees/tree_0_0_0.dat,outputs/merger_tree.cfg,$(1))

# # ------------------------------------------------------------------------------
# #
# #  Parsing merger trees
# #
# # ------------------------------------------------------------------------------
# # Here I again use a sentinel file, since there will be multiple files created
# # by the one C file:
# merger_sentinels = $(foreach dir,$(sim_checks_dirs),$(dir)/merger_sentinel.txt)
# # need to get the trees that the sim should parse
# merger_to_tree = $(subst checks/merger_sentinel.txt,rockstar_halos/trees/tree_0_0_0.dat,$(1))

# # ------------------------------------------------------------------------------
# #
# #  timing output
# #
# # ------------------------------------------------------------------------------
# timings = $(foreach dir,$(sim_dirs),$(dir)/checks/timing_debug.txt)
# timing_to_sims = $(call dir_to_sims,$(subst checks/timing_debug.txt,out,$(1)))
# timing_to_dir = $(subst checks/timing_debug.txt,log,$(1))

# # ------------------------------------------------------------------------------
# #
# #  Movies
# #
# # ------------------------------------------------------------------------------
# movies_all = $(foreach dir,$(sim_dirs),$(dir)/plots/$(1).mp4)
# movies_hydro = $(foreach dir,$(sim_dirs_hydro),$(dir)/plots/$(1).mp4)

# # Here 1 is the full location of the movie
# movie_to_out_dir = $(subst plots/,out,$(dir $(1)))
# movie_to_sims = $(call dir_to_sims,$(call movie_to_out_dir,$(1)))
# # in this next function argument 1 is the full location where the movie will
# # be saved, and 2 is the function that turns a simulation into a plot file
# # (these are all defined above)
# movie_to_plots = $(foreach sim,$(call movie_to_sims,$(1)),$(call $(2),$(sim)))
# movie_to_metal_summary = $(call sim_to_summary_metal,$(call dir_to_sims,$(subst plots/$(1).mp4,out,$(2))))
# movie_to_plot_dir = $(subst /$(1).mp4,,$(2))

# ------------------------------------------------------------------------------
#
#  Rules
#
# ------------------------------------------------------------------------------


# movies = $(call movies_all,n_body_refined) $(call movies_all,n_body_split_refined) $(call movies_all,n_body_local_group) $(call movies_all,n_body_split_local_group) $(call movies_hydro,gas_density) $(call movies_hydro,gas_velocity_x) $(call movies_hydro,gas_velocity_y) $(call movies_hydro,gas_velocity_z)
# all: $(my_directories) $(debugs) $(sfh_plots) $(movies) $(smhm_plots) $(timings) $(merger_sentinels)


# # have phony target for movies
# .PHONY: movies
# movies: $(movies)

# # Make the debug files
# .SECONDEXPANSION:
# $(debugs): %: $(debug_script)
#     python $(debug_script) $(call debug_to_sim, $@) clobber silent

# # Make the summary files
# .SECONDEXPANSION:
# $(galaxies): %: $$(call galaxies_to_halo, %) $(galaxies_script)
#     python $(galaxies_script) $(call galaxies_to_sim, $@) $(call galaxies_to_halo, $@) clobber silent

# # Make the SFH plots. That & indicates a grouped target, telling make that the recipe
# # produces all the targets.
# $(sfh_plots) &: $(snapshots_hydro)
#     python $(sfh_plots_script) $(sim_dirs_hydro)

# # Make the individual nbody plots - several examples of very similar things here
# .SECONDEXPANSION:
# $(halo_1_full_plots): %: $$(call halo_1_full_to_halo, %) $(nbody_single_halo_plots_script) $(nbody_full_plot_script)
#     python $(nbody_single_halo_plots_script) $(call halo_1_full_to_sim, $@) $(call halo_1_full_to_halo, $@) 1 full

# .SECONDEXPANSION:
# $(halo_2_full_plots): %: $$(call halo_2_full_to_halo, %) $(nbody_single_halo_plots_script) $(nbody_full_plot_script)
#     python $(nbody_single_halo_plots_script) $(call halo_2_full_to_sim, $@) $(call halo_2_full_to_halo, $@) 2 full

# .SECONDEXPANSION:
# $(halo_1_split_plots): %: $$(call halo_1_split_to_halo, %) $(nbody_single_halo_plots_script) $(nbody_split_plot_script)
#     python $(nbody_single_halo_plots_script) $(call halo_1_split_to_sim, $@) $(call halo_1_split_to_halo, $@) 1 split

# .SECONDEXPANSION:
# $(halo_2_split_plots): %: $$(call halo_2_split_to_halo, %) $(nbody_single_halo_plots_script) $(nbody_split_plot_script)
#     python $(nbody_single_halo_plots_script) $(call halo_2_split_to_sim, $@) $(call halo_2_split_to_halo, $@) 2 split

# .SECONDEXPANSION:
# $(refined_full_plots): %: $$(call refined_full_to_halo, %) $(nbody_refined_plot_script) $(nbody_full_plot_script)
#     python $(nbody_refined_plot_script) $(call refined_full_to_sim, $@) $(call refined_full_to_halo, $@) full

# .SECONDEXPANSION:
# $(refined_split_plots): %: $$(call refined_split_to_halo, %) $(nbody_refined_plot_script) $(nbody_split_plot_script)
#     python $(nbody_refined_plot_script) $(call refined_split_to_sim, $@) $(call refined_split_to_halo, $@) split

# .SECONDEXPANSION:
# $(local_group_full_plots): %: $$(call local_group_full_to_halo, %) $(nbody_local_group_plot_script) $(nbody_full_plot_script)
#     python $(nbody_local_group_plot_script) $(call local_group_full_to_sim, $@) $(call local_group_full_to_halo, $@) full

# .SECONDEXPANSION:
# $(local_group_split_plots): %: $$(call local_group_split_to_halo, %) $(nbody_local_group_plot_script) $(nbody_split_plot_script)
#     python $(nbody_local_group_plot_script) $(call local_group_split_to_sim, $@) $(call local_group_split_to_halo, $@) split

# # parameters:
# # 1 - movie
# # 2 - sim to plot name
# define movie_set_template =
#     $(1):
#         $(info $@)
#         $(info $<)
#         $(info "===========")
#         ffmpeg -framerate 2 -pattern_type glob -i '$(call movie_to_plot_dir,$(2),$@)/$(2)*.png' -c:v h264 -pix_fmt yuv420p -y $@
# endef
# # $$(call movie_to_plots,%,$(2))
# $(foreach movie,$(call movies_all,n_body_refined),$(eval $(call movie_set_template,$(movie),n_body_refined)))

# # # Make the movies. We can't parametrize this all since there is no third expansion
# # .SECONDEXPANSION:
# # $(call movies_all,n_body_refined): %: $$(call movie_to_plots,$@,sim_to_refined_full)
# #   $(info $@)
# #   $(info $<)
# # #     $(info $(call movie_to_plots,$@,sim_to_refined_full))
# #   $(info "===========")
# # #     ffmpeg -framerate 2 -pattern_type glob -i '$(call movie_to_plot_dir,n_body_refined,$@)/n_body_refined*.png' -c:v h264 -pix_fmt yuv420p -y $@

# # SECONDEXPANSION:
# # $(call movies_all,n_body_split_refined): %: $$(call movie_to_plots,$@,sim_to_refined_split)
# #   $(info $@)
# #   $(info $<)
# # #     $(info $(call movie_to_plots,$@,sim_to_refined_split))
# #   $(info "===========")
# #   ffmpeg -framerate 2 -pattern_type glob -i '$(call movie_to_plot_dir,n_body_split_refined,$@)/n_body_split_refined*.png' -c:v h264 -pix_fmt yuv420p -y $@

# # .SECONDEXPANSION:
# # $(call movies_all,n_body_local_group): %: $$(call movie_to_plots,%,sim_to_local_group_full)
# #   $(info $@)
# #   $(info $<)
# #   $(info "===========")
# # #     ffmpeg -framerate 2 -pattern_type glob -i '$(call movie_to_plot_dir,n_body_local_group,$@)/n_body_local_group*.png' -c:v h264 -pix_fmt yuv420p -y $@

# # SECONDEXPANSION:
# # $(call movies_all,n_body_split_local_group): %: $$(call movie_to_plots,%,sim_to_refined_split)
# #   ffmpeg -framerate 2 -pattern_type glob -i '$(call movie_to_plot_dir,n_body_split_local_group,$@)/n_body_split_local_group*.png' -c:v h264 -pix_fmt yuv420p -y $@

# # .SECONDEXPANSION:
# # $(call movies_hydro,gas_density): %: $$(call movie_to_metal_summary,gas_density,%)
# #   ffmpeg -framerate 2 -pattern_type glob -i '$(call movie_to_plot_dir,gas_density,$@)/gas_density_*.png' -c:v h264 -pix_fmt yuv420p -y $@

# # .SECONDEXPANSION:
# # $(call movies_hydro,gas_velocity_x): %: $$(call movie_to_metal_summary,gas_velocity_x,%)
# #   ffmpeg -framerate 2 -pattern_type glob -i '$(call movie_to_plot_dir,gas_velocity_x,$@)/gas_velocity_x_*.png' -c:v h264 -pix_fmt yuv420p -y $@

# # SECONDEXPANSION:
# # $(call movies_hydro,gas_velocity_y): %: $$(call movie_to_metal_summary,gas_velocity_y,%)
# #   ffmpeg -framerate 2 -pattern_type glob -i '$(call movie_to_plot_dir,gas_velocity_y,$@)/gas_velocity_y_*.png' -c:v h264 -pix_fmt yuv420p -y $@

# # SECONDEXPANSION:
# # $(call movies_hydro,gas_velocity_z): %: $$(call movie_to_metal_summary,gas_velocity_z,%)
# #   ffmpeg -framerate 2 -pattern_type glob -i '$(call movie_to_plot_dir,gas_velocity_z,$@)/gas_velocity_z_*.png' -c:v h264 -pix_fmt yuv420p -y $@

# # Make the consistent trees config files
# .SECONDEXPANSION:
# $(tree_cfgs): %: $$(call tree_cfg_to_sentinel,%)
#     perl $(tree_config_script) $(call tree_cfg_to_rockstar_cfg,$@)

# # Then build the merger trees
# .SECONDEXPANSION:
# $(trees): %: $$(call tree_to_tree_cfg,%)
#     cd $(tree_dir) && perl do_merger_tree.pl $<

# # Build the merger tree reading code
# $(read_tree_exe): $(read_tree_src)
#     cd $(read_tree_dir) && make

# # Build the accretion history output files
# .SECONDEXPANSION:
# $(merger_sentinels): %: $$(call merger_to_tree,%) $(read_tree_exe)
#     $(read_tree_exe) $(call merger_to_tree,$@) && touch $@

# # SFH plots
# .SECONDEXPANSION:
# $(smhm_plots): %: $$(call smhm_to_sim,%) $(sfh_plots_script)
#     python $(sfh_plots_script) $(call smhm_to_sim, $@) $(call sim_to_halo, $(call smhm_to_sim, $@))

# # timing output
# .PHONY: timing
# timing: $(timings)

# .SECONDEXPANSION:
# $(timings): %: $$(call timing_to_sims, %)
#     $(timing_script) $(call timing_to_dir, $@) > $@ || echo "Timing failed"

makefile.close()


# ------------------------------------------------------------------------------

#  Code locations that are relative to this file

# ------------------------------------------------------------------------------


debug_script = join("debug_output.py")
galaxies_script = join("galaxy_summaries.py")
nbody_single_halo_plots_script = join("plot_single_halo_nbody.py")
nbody_refined_plot_script = join("plot_refined_region_nbody.py")
nbody_local_group_plot_script = join("plot_local_group_nbody.py")
nbody_full_plot_script = join("utils/nbody_projection_all_species.py")
nbody_split_plot_script = join("utils/nbody_projection_split_species.py")
sfh_plots_script = join("plot_sfh.py")
read_tree_dir = join("read_tree")
read_tree_exe = join(read_tree_dir, "halo_history")
read_tree_src = join(read_tree_dir, "halo_history.c")


snapshots = [base_to_sims(d) for d in sim_dirs]
snapshots_hydro = [base_to_sims(d) for d in sim_dirs_hydro]
