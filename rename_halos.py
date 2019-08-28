import sys
import os
import shutil

halo_output_file = os.path.abspath(sys.argv[1])
halo_output_dir = os.path.dirname(halo_output_file)
rockstar_dir = halo_output_dir.replace("/halos", "/rockstar_halos")
datasets_txt = rockstar_dir + os.sep + "datasets.txt"

scale_factor_str = halo_output_file[-12:-6]

with open(datasets_txt) as datasets:
    found = False
    for row in datasets:
        if row.startswith("#"):
            continue
        ds_name = row.split()[0]
        idx = row.split()[1]
        
        # parse the filename to get the scale factor
        this_scale_factor = ds_name[-10:-4]
        
        if scale_factor_str == this_scale_factor:
            found = True
            break
    
    if not found:
        raise ValueError("Dataset {} not found in {}".format(scale_factor_str, datasets_txt))
# now we have the right index
old_prefix = "halos_{}.".format(idx)
new_prefix = "halos_a{}.".format(scale_factor_str)

all_rockstar_outputs = os.listdir(rockstar_dir)
for rockstar_output in all_rockstar_outputs:
    if old_prefix in rockstar_output:
        old_path = rockstar_dir + os.sep + rockstar_output
        new_name = rockstar_output.replace(old_prefix, new_prefix)
        new_path = halo_output_dir + os.sep + new_name
        shutil.copy(old_path, new_path)
