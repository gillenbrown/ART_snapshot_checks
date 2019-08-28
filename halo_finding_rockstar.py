import sys
import os
import pathlib
import shutil
import time
import yt
from yt.extensions.astro_analysis.halo_finding.rockstar.api import RockstarHaloFinder
from yt.data_objects.particle_filters import add_particle_filter
yt.funcs.mylog.setLevel(0)  # ignore yt's output
yt.enable_parallelism()

# format of sys.argv:
# idx 0: name of script
# ids 1: directory of simulation outputs to run halos on
# idx 2: directory to store halo finding outputs in
# idx 3: must be `silent` or not present

# print function that always flushes
import functools
print = functools.partial(print, flush=True)

# do error checking on root only.
if yt.is_root():
    if len(sys.argv) < 3:
        raise ValueError("Please provide the proper command line argument.")
    elif len(sys.argv) == 4:
        if sys.argv[3] != "silent":
            raise ValueError("If present, the third command line "
                             "argument can only be `silent`")
    elif len(sys.argv) > 4:
        raise ValueError("Extra command line arguments not recognized.")

# turn the directories the user passes into the absolute path
sim_dir = os.path.abspath(sys.argv[1])
out_dir = os.path.abspath(sys.argv[2])
if not sim_dir.endswith(os.sep):
    sim_dir += os.sep
if not out_dir.endswith(os.sep):
    out_dir += os.sep

# check to see if there is a currently existing halo catalog already here
# to restart from. 
if os.path.exists(out_dir + "datasets.txt"):
    restart = True
else:
    restart = False

# rockstar has a bug. When restarting, it overwrites the existing datasets.txt 
# file, but doesn't include the original datasets that were in that file. This 
# is only a bug with the datasets.txt file. The outputs are correctly numbered,
# but the datasets.txt file doesn't reflect that. To fix this, we will copy the
# datasets.txt file to a new location so that we can save the info until the 
# end of this script, when we'll parse the two datasets.txt
if restart and yt.is_root():
    shutil.copyfile(out_dir + "datasets.txt", out_dir + "old_datasets.txt")


if yt.is_root():
    print("Reading simulations from: {}".format(sim_dir))
    print("Writing halo catalogs to: {}".format(out_dir))

time.sleep(20)# wait to make sure things are done

ts = yt.load(sim_dir + 'continuous_a?.????.art')

# check what kind of particles are present
if ('N-BODY_0', 'MASS') in ts[0].derived_field_list:
    particle_type = "N-BODY_0"
else:
    particle_type = "N-BODY"

# Lou analysis machines have Skylake nodes with 20 cores. We have one master
# process too.
rh = RockstarHaloFinder(ts, num_readers=9, num_writers=10, outbase=out_dir,
                        particle_type=particle_type)
rh.run(restart=restart)

time.sleep(20)  # wait for all to be done

# then we can fix the dataset.txt nonsense
if restart and yt.is_root():
    # first move the new datasets.txt to a place where we won't overwrite it
    shutil.move(out_dir + "datasets.txt", out_dir + "new_datasets.txt")

    # the original file's numbering is fine. We don't need to mess with that.
    # But we do need to get the number of the last file there, so we'll read 
    # that file in.
    with open(out_dir + "datasets.txt", "w") as final_datasets:
        with open(out_dir + "old_datasets.txt", "r") as old_datasets:
            for line in old_datasets:
                # put all the lines in the new file, including the header
                final_datasets.write(line)
                if not line.startswith("#"):
                    # get the index of this line. They are in order, so the
                    # last one processed will be the highest index
                    last_index = int(line.split()[-1])
                
        # then we can go through the new file and replace the indices properly
        with open(out_dir + "new_datasets.txt", "r") as new_datasets:
            for line in new_datasets:
                # we don't need to copy the header for the second file
                if line.startswith("#"):
                    continue  

                new_idx = int(line.split()[-1])
                # make the proper index. The first in this file will currently 
                # have index 0, but needs last_index + 1.
                correct_idx = new_idx + last_index + 1 
                # then we can replace that in the string. Tabs are used to 
                # separate the index from the sim location, so we can use
                # that to avoid accidentally replacing part of the scale
                # factor in the sim name
                line = line.replace("\t{}".format(str(new_idx)), 
                                    "\t{}".format(str(correct_idx)))
                # put the line with the updated index in the new file
                final_datasets.write(line)

    # and finally delete the temporary files
    os.remove(out_dir + "old_datasets.txt")
    os.remove(out_dir + "new_datasets.txt")

# update the sentinel file
if yt.is_root():
    pathlib.Path(out_dir + "sentinel.txt").touch()
