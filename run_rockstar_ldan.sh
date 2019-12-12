# see ./halo_finding_rockstar.py for the parameters this takes
# Set some environment variables so MPI is sure to be launched properly
export SEVERAL_TRIES_NTRIES=10
export MPI_LAUNCH_TIMEOUT=300 
# Lou analysis machines have Skylake nodes with 20 cores
/pleiades/u/scicon/tools/bin/several_tries mpiexec -np 3 python ./halo_finding_rockstar.py $1 $2 $3
