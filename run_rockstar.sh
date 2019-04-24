# Run this command in the directory where you want to deposit the halo catalogs.
# the first parameter is the location of the directory holding the output files.
# an optional second parameter can only be "remove". If passed, this will remove old
#     halo catalogs. If not provided, none will be deleted. This is just a check 
#     to make sure you really want to do this.
# Lou analysis machines have Skylake nodes with 20 cores
mpirun -n 20 --mca btl ^openib python ./halo_finding_rockstar.py $1 $2 $3 $4
