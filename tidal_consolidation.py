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
    with open(tidal_dir / f"{star}.txt", "w") as out:
        for filename in sorted(file_list):
            with open(filename, "r") as file:
                for l in file:
                    out.write(l)

sentinel.touch()
