#!/bin/zsh

# Script to copy plots to my personal machine from shangrila and Stampede2

local_dir=/Users/gillenb/code/ART_snapshot_checks/comparison_plots
rm -r $local_dir

# shangrila first
scp -r gillenb@shangrila.astro.lsa.umich.edu:code/mine/ART_snapshot_checks/comparison_plots $local_dir

# then Stampede2. This uses the the globus command line interface
globus transfer -r --label plots --notify failed ceea5ca0-89a9-11e7-a97f-22000a92523b:/work2/06912/tg862118/stampede2/ART_snapshot_checks_analysis/comparison_plots d7dd2e72-e628-11e9-893b-0ea49c25678c:$local_dir