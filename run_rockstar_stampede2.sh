# This script takes two parameters
# the first parameter is the location of the directory holding the output files
# the second is the directory to put the rockstar outputs.
# Here we use 47 cores (with an offset of 1) to use all the processes not used
# by the master make process
ibrun -n 47 -o 1 python3 ./halo_finding_rockstar.py $1 $2 47
