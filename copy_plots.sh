#!/bin/zsh

# Script to copy plots to my personal machine from shangrila and Stampede2

local_dir=/Users/gillenb/code/ART_snapshot_checks/comparison_plots
rm -r $local_dir

# shangrila first
scp -r gillenb@shangrila.astro.lsa.umich.edu:code/mine/ART_snapshot_checks/comparison_plots $local_dir

# then Stampede2
