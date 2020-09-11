import sys
import os
import pathlib
import yt
from yt.extensions.astro_analysis.halo_finding.rockstar.api import RockstarHaloFinder
yt.funcs.mylog.setLevel(0)  # ignore yt's output
yt.enable_parallelism()

# format of sys.argv:
# idx 0: name of script
# ids 1: directory of simulation outputs to run halos on
# idx 2: directory to store halo finding outputs in

# print function that always flushes
import functools
print = functools.partial(print, flush=True)

# do error checking on root only.
if yt.is_root():
    if len(sys.argv) != 4:
        raise ValueError("Please provide the proper command line argument.")

# turn the directories the user passes into the absolute path
sim_dir = os.path.abspath(sys.argv[1])
out_dir = os.path.abspath(sys.argv[2])
cores_to_use = int(sys.argv[3])
if not sim_dir.endswith(os.sep):
    sim_dir += os.sep
if not out_dir.endswith(os.sep):
    out_dir += os.sep

# check to see if there is a currently existing halo catalog already here
# to restart from. 
if os.path.exists(out_dir + "restart.cfg"):
    restart = True
else:
    restart = False

if yt.is_root():
    print("Reading simulations from: {}".format(sim_dir))
    print("Writing halo catalogs to: {}".format(out_dir))

ts = yt.load(sim_dir + 'continuous_a?.????.art')

# check what kind of particles are present
if ('N-BODY_0', 'MASS') in ts[0].derived_field_list:
    particle_type = "N-BODY_0"
else:
    particle_type = "N-BODY"

# determine how many readers and writers to use, depending on what machine we're on. We
# have one master process, plus the number of readers and writers. Reading is quicker
# then analysis, so we typically have fewer readers
processing_cores = cores_to_use - 1
if processing_cores < 5:
    readers = 1
elif processing_cores > 30:
    readers = 16
else:
    raise NotImplementedError("Choose writers for 5 < cores < 30")
writers = processing_cores - readers
if yt.is_root():
    print("Rockstar will be running with:")
    print(f"\t- {cores_to_use} total cores")
    print(f"\t- 1 master process")
    print(f"\t- {readers} readers")
    print(f"\t- {writers} writers")
rh = RockstarHaloFinder(ts, num_readers=readers, num_writers=writers, outbase=out_dir,
                        particle_type=particle_type)
rh.run(restart=restart)

# update the sentinel file
if yt.is_root():
    pathlib.Path(out_dir + "sentinel.txt").touch()
