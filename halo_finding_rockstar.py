import sys
import os
import pathlib
import shutil
import yt
from yt.extensions.astro_analysis.halo_finding.rockstar.api import RockstarHaloFinder
from yt.data_objects.particle_filters import add_particle_filter
yt.funcs.mylog.setLevel(50)  # ignore yt's output
yt.enable_parallelism()

# format of sys.argv:
# idx 0: name of script
# ids 1: directory of simulation outputs to run halos on
# idx 2: directory to store halo finding outputs in
# idx 3: must be `remove`, `silent` or not present

# do error checking on root only.
if yt.is_root():
    if len(sys.argv) < 3:
        raise ValueError("Please provide the proper command line argument.")
    elif len(sys.argv) == 4:
        if sys.argv[3] != "remove":
            raise ValueError("If present, the third command line "
                             "argument can only be `remove`")
    elif len(sys.argv) > 4:
        raise ValueError("Extra command line arguments not recognized.")
    # check the order. The input will be checked elsewhere, so we only have to 
    # check remove
    if sys.argv[2] == "remove":
        raise ValueError("The third command line argument must be `remove`, "
                         "not the second.")

# turn the directories the user passes into the absolute path
sim_dir = os.path.abspath(sys.argv[1])
out_dir = os.path.abspath(sys.argv[2])
if not sim_dir.endswith(os.sep):
    sim_dir += os.sep
if not out_dir.endswith(os.sep):
    out_dir += os.sep

# check to see if there is a currently existing halo catalog already here
# again, do this only on the root.
if yt.is_root():
    if os.path.exists(out_dir):
        if len(sys.argv) == 3:  # remove not present
            raise ValueError("Existing halo catalogs already present.\n"
                             "            Use keyword `remove` to overwrite.")
        else:  # user does want to remove
            shutil.rmtree(out_dir)
    else:
        os.mkdir(out_dir)

if yt.is_root():
    print("Reading simulations from: {}".format(sim_dir))
    print("Writing halo catalogs to: {}".format(out_dir))

ts = yt.load(sim_dir + 'continuous_a?.????.art')

# check what kind of particles are present
if ('N-BODY_0', 'MASS') in ts[0].derived_field_list:
    particle_type = "N-BODY_0"
else:
    particle_type = "N-BODY"

# Lou analysis machines have Skylake nodes with 20 cores. We have one master
# process too.
rh = RockstarHaloFinder(ts, num_readers=5, num_writers=14, outbase=out_dir,
                        particle_type=particle_type)
rh.run()

# update the sentinel file
pathlib.Path(out_dir + "sentinel.txt").touch()
