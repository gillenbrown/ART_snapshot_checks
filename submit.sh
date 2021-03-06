#!/bin/bash
#SBATCH --job-name=snapshot_checks
#SBATCH --account=TG-AST200017
#SBATCH --mail-user=gillenb@umich.edu
#SBATCH --mail-type=ALL
#SBATCH --no-requeue
#SBATCH --partition=skx-normal
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --cpus-per-task=1
#SBATCH --output=stdout_%x_%j

# This took a lot of messing with to get right. Do NOT use `module reset`. This
# will stop conda from recognizing my distribution somehow.
# We do need to source bashrc, as it holds the code to activate conda and
# activate the correct env
module load remora
source ~/.bashrc
# Make sure the system python is not accessible
unset PYTHONPATH

# double check these variables
echo $(which python)
echo $PYTHONPATH

# go to the directory I submitted the job from
cd $SLURM_SUBMIT_DIR


# halos uses all the cores, so it needs to be done one at a time. The number
# used thereafter depends on the memory, it's just empirical. For halos we have
# one master process for make. The halo finding script will use ibrun to spawn
# processes on the rest of the cores, although we have to restrict based on
# memory considerations
make dirs
remora ibrun -n 1 -o 0 make halos
remora make -j10
# &>> $SLURM_JOB_NAME.stdout.$SLURM_JOB_ID