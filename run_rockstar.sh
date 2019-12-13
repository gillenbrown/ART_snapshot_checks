# This script takes three parameters
# the first parameter is the location of the directory holding the output files
# the second is the directory to put the rockstar outputa
# an optional third parameter can only be "silent". If passed, this will 
#     suppress printed output statements 
# Lou analysis machines have Skylake nodes with 20 cores
mpiexec -np 4 python ./halo_finding_rockstar.py $1 $2 $3
