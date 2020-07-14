# This script takes two parameters
# the first parameter is the location of the directory holding the output files
# the second is the directory to put the rockstar outputs
mpiexec -np 4 python ./halo_finding_rockstar_single_output.py $1 $2
