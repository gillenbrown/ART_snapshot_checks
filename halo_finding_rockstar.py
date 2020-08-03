"""
halo_finding_rockstar.py - Run the Rockstar halo finder on a set of outputs

This is a rather complex script now. It basically works as its own makefile. It
looks for all the existing simulations, then sees which ones still need to have
a halo file created for them. It then does the halo processing for all that 
need it. This is done in chunks to avoid too much memory with the production
runs, which are huge.
"""
import sys
import shutil
import re
from pathlib import Path
import gc
import yt
from yt.extensions.astro_analysis.halo_finding.rockstar.api import RockstarHaloFinder
yt.funcs.mylog.setLevel(0)  # ignore yt's output
yt.enable_parallelism()

# print function that always flushes
import functools
print = functools.partial(print, flush=True)

# yt.is_root has issues, I'll make my own. This runs at the beginning of the 
# script so any future changes to the rank numbers (which I've seen happen
# somehow!) don't affect anything
is_root = yt.is_root()

# format of sys.argv:
# idx 0: name of script
# ids 1: directory of simulation outputs to run halos on
# idx 2: directory to store rockstar halo finding outputs in
# idx 3: directory to store processed halo finding outputs in
# idx 4: how many cores to use for Rockstar

# do error checking on root only.
if is_root:
    print(len(sys.argv))
    if len(sys.argv) != 5:
        raise ValueError("Please provide the proper command line argument.")

# turn the directories the user passes into the absolute path
sim_dir = Path(sys.argv[1]).resolve()
rockstar_dir = Path(sys.argv[2]).resolve()
halos_dir = Path(sys.argv[3]).resolve()
cores_to_use = int(sys.argv[4])

if is_root:
    print("Reading simulations from: {}".format(sim_dir))
    print("Writing intermediate halo catalogs to: {}".format(rockstar_dir))
    print("Writing final halo catalogs to: {}".format(halos_dir))

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
if is_root:
    print("Rockstar will be running with:")
    print(f"\t- {cores_to_use} total cores")
    print(f"\t- 1 master process")
    print(f"\t- {readers} readers")
    print(f"\t- {writers} writers")

# check to see if there is a currently existing config file already here
# to restart from. 
if (rockstar_dir / "restart.cfg").exists():
    restart = True
else:
    restart = False

# ==============================================================================
# 
# determining which snapshots need to have Rockstar run on them
# 
# ==============================================================================
# first get all the simulations
all_sims = []
for item in sorted(sim_dir.iterdir()):
    if item.name.startswith("continuous_a") and item.name.endswith(".art"):
        all_sims.append(item)
# This should already be sorted, but be sure
all_sims = sorted(all_sims)

# then find the first one that does not have an existing halo catalog. This is
# what we will start with
def sim_to_halo(sim_path, idx=0, extension="bin"):
    halo_name = Path(sim_path).name.replace("continuous", "halos")
    halo_name = halo_name.replace(".art", f".{idx}.{extension}")
    return halos_dir / halo_name

def sim_to_out(sim_path):
    out_name = Path(sim_path).name.replace("continuous", "out")
    out_name = out_name.replace(".art", ".list") 
    return halos_dir / out_name

needed_sims = []
for sim in all_sims:
    halo = sim_to_halo(sim)
    # if we have found previous catalogs that need making, any further ones
    # will need to be remade
    if len(needed_sims) > 0:
        if halo.exists() and is_root:
            halo.unlink()
    # if there is no halo catalog, we'll need to make it
    if not halo.exists():
        needed_sims.append(sim)

# ==============================================================================
# 
# splitting simulations into chunks
# 
# ==============================================================================
# I do not use the simple method of loading all as one time series object. This
# turns out to be too memory intensive for use on stampede2 with the very large
# output files. Instead, I'll create a list of several at a time and do it in 
# chunks. There needs to be overlap between the chunks so that Rockstar can 
# link them together.
chunk_size = 2
chunks = [[]]
for sim in sorted(needed_sims):
    # then add the simulation to the last chunk. Convert to a string so
    # yt knows how to read it
    chunks[-1].append(str(sim))
    # if the last chunk is full, add another, and put this sim at the front of
    # that chunk. Don't do this if we're at the last simulation.
    if len(chunks[-1]) == chunk_size and sim != sorted(needed_sims)[-1]:
        chunks.append([])
        chunks[-1].append(str(sim))

if is_root:
    print("Running the halo finder on these chunks:")
    print("------")
    for c in chunks:
        for s in c:
            print(s)
        print("------")

# ==============================================================================
# 
# Running the analysis
# 
# ==============================================================================
# then do the halo finding on the chunks, one at a time
for chunk in chunks:
    # --------------------------------------------------------------------------
    # managing restart file
    # --------------------------------------------------------------------------
    # If we have to restart, we need to manage things to make sure the files 
    # there are appropriate. 
    # If we're restarting, copy the last previous successfull file over
    # to be the start of this run.
    if restart and is_root:
        print("in second loop")
        max_scale = "a0.0000"  # keep as string to avoid parsing. comparison still work
        # use a regular expression to get what we need to compare to
        regex = re.compile("[a][0-1]\.[\d]{4}")
        for file in halos_dir.iterdir():
            match = regex.search(file.name)
            if match:
                scale = match.group(0)
                if scale > max_scale:
                    max_scale = scale
        if max_scale == "a0.0000":
            raise RuntimeError("Rockstar restart failed")
        print(max_scale)

        # we can then find all the files that have this and copy them over
        for file in halos_dir.iterdir():
            if max_scale in file.name:
                new_file = rockstar_dir / file.name.replace(max_scale, "0")
                shutil.copy2(file, new_file)


        # Then modify the config file
        restart_old = rockstar_dir / "restart.cfg"
        restart_new = rockstar_dir / "restart.cfg.temp"

        # we need the actual scale factor
        ds_last = yt.load(str(sim_dir / f"continuous_{max_scale}.art"))
        scale = ds_last.scale_factor
        del ds_last

        with open(restart_old, "r") as old:
            with open(restart_new, "w") as new:
                for line in old:
                    split = line.split()
                    key = split[0]
                    value = split[-1]
                    # we need to change the restart number
                    if key == "RESTART_SNAP":
                        line = line.replace(value, "1")
                    # # and the scale factor
                    # if key == "SCALE_NOW":
                    #     line = line.replace(value, str(scale))
                    new.write(line)

        # then copy the file over
        shutil.move(restart_new, restart_old)

    # --------------------------------------------------------------------------
    # then we can run the halo finder on these files
    # --------------------------------------------------------------------------
    ts = yt.DatasetSeries(chunk)

    # check what kind of particles are present
    if ('N-BODY_0', 'MASS') in ts[0].derived_field_list:
        particle_type = "N-BODY_0"
    else:
        particle_type = "N-BODY"

    rh = RockstarHaloFinder(ts, num_readers=readers, num_writers=writers, 
                            outbase=str(rockstar_dir), 
                            particle_type=particle_type)
    rh.run(restart=restart)

    # to clear memory, explicity delete and garbage collect this. Not sure if
    # this helps, but it can't hurt.
    for sim in ts:
        del sim
    del ts
    del rh
    gc.collect()

    # --------------------------------------------------------------------------
    # then move the output files to the correct directory
    # --------------------------------------------------------------------------
    if is_root:
        for sim_idx, sim in enumerate(chunk):
            old_list = rockstar_dir / f"out_{sim_idx}.list"
            old_list.rename(sim_to_out(sim))
            for write_idx in range(writers):
                for suffix in ["ascii", "bin"]:
                    old_halo = Path(rockstar_dir / f"halos_{sim_idx}.{write_idx}.{suffix}")
                    old_halo.rename(sim_to_halo(sim, write_idx, suffix))

        # then delete all the other stuff manually
        (rockstar_dir / "auto-rockstar.cfg").unlink()
        (rockstar_dir / "datasets.txt").unlink()
        (rockstar_dir / "rockstar.cfg").unlink()
        shutil.rmtree(rockstar_dir / "profiling")

    # now that we've done one loop, we can restart from the next set
    restart = True

# assert that every simulation has a matching halo file
for sim in all_sims:
    assert sim_to_halo(sim, 0).exists()

# update the sentinel file
# if is_root:
#     pathlib.Path(out_dir + "sentinel.txt").touch()
