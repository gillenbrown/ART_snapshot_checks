import sys
from pathlib import Path
from collections import defaultdict

from tqdm import tqdm

sentinel = Path(sys.argv[1]).resolve()
tidal_dir = sentinel.parent
log_dir = tidal_dir.parent / "log"

# First go through and find all the tidal output files.
all_tidal = []
for d in log_dir.iterdir():
    for t in (d / "tidal").iterdir():
        all_tidal.append(t)

# then determine which files are associated with the same star particle.
star_files = defaultdict(list)
for tidal_file in all_tidal:
    # the star ID is the name of the file
    star_files[tidal_file.stem].append(tidal_file)

# then combine these particles.
for star, file_list in tqdm(star_files.items()):
    out_file = tidal_dir / f"{star}.txt"
    # first go through and determine which output files have already been
    # written to this file
    pre_existing_runs = []
    if out_file.exists():  # can't read if it doesn't exist
        with open(out_file, "r") as old_file:
            for line in old_file:
                if line.startswith("# run name"):
                    pre_existing_runs.append(line.split()[-1])

    # then go through and append to the file, if the current source has not
    # already been written.
    with open(out_file, "a") as full_tidal_file:
        for individual_tidal_file in sorted(file_list):
            # only write this if we haven't already added this log file
            log_name = individual_tidal_file.parent.parent.name
            if log_name not in pre_existing_runs:
                # note that we're writing from this file
                full_tidal_file.write(f"# run name {log_name}\n")
                # then write all the tidal info from this file
                with open(individual_tidal_file, "r") as file:
                    for l in file:
                        full_tidal_file.write(l)

sentinel.touch()
