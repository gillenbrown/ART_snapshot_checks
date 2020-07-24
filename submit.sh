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

source ~/.bashrc

module reset
module load remora python3/3.7.0

# go to the directory I submitted the job from
cd $SLURM_SUBMIT_DIR

# halos uses all the cores, so it needs to be done one at a time. The number
# used thereafter depends on the memory, it's just empirical. For halos we have
# one master process for make. The halo finding script will use ibrun to spawn
# processes on the rest of the cores
make dirs
remora ibrun -n 1 -o 0 make halos
remora make -j48
# &>> $SLURM_JOB_NAME.stdout.$SLURM_JOB_ID