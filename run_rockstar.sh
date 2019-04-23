# Run this command in the directory where you want to deposit the halo catalogs.
# the first parameter is the location of the directory holding the output files.
# an optional second parameter can only be "remove". If passed, this will remove old
#     halo catalogs. If not provided, none will be deleted. This is just a check 
#     to make sure you really want to do this.
mpirun -n 6 --mca btl ^openib python /u/home/gillenb/code/mine/halo_scripts/halo_finding_rockstar.py $1 $2 $3 $4
