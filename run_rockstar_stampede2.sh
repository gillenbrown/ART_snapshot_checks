# This script takes two parameters
# the first parameter is the location of the directory holding the output files
# the second is the directory to put the rockstar outputs
mpiexec -np 68 python3 ./halo_finding_rockstar.py $1 $2 68
