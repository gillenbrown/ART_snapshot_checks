import sys
from pathlib import Path
import yt
from yt.extensions.astro_analysis.halo_finding.rockstar.api import RockstarHaloFinder
yt.funcs.mylog.setLevel(50)  # ignore yt's output - can disable to debug
yt.enable_parallelism()

# format of sys.argv:
# idx 0: name of script
# ids 1: location of simulation output to run halos on
# idx 2: directory to store halo finding outputs in

# do error checking on root only.
if yt.is_root():
    if len(sys.argv) != 3:
        raise ValueError("Please provide the proper command line argument.")

# turn the directories the user passes into the absolute path
sim_loc = str(Path(sys.argv[1]).absolute())
out_dir = str(Path(sys.argv[2]).absolute())

if yt.is_root():
    print("Reading simulation from: {}".format(sim_loc))
    print("Writing halo catalogs to: {}".format(out_dir))

ds = yt.load(sim_loc)

# check what kind of particles are present
if ('N-BODY_0', 'MASS') in ds.derived_field_list:
    particle_type = "N-BODY_0"
else:
    particle_type = "N-BODY"

# The number of readers and writers is relatively arbitrary in terms of 
# performance, I've found (although without any rigorous testing). The important
# this is the number of files created, which is proportional to num_writers.
# on systems with limited file counts this should be a lower number. 
rh = RockstarHaloFinder(ds, num_readers=1, num_writers=2, 
                        outbase=out_dir, particle_type=particle_type)
rh.run()
