#!/bin/bash
#SBATCH --job-name=snapshot_checks
#SBATCH --account=TG-AST200017
#SBATCH --mail-user=gillenb@umich.edu
#SBATCH --mail-type=ALL
#SBATCH --no-requeue
#SBATCH --partition=normal
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=68
#SBATCH --output=stdout_%x_%j

source ~/.bashrc

module reset
module load remora python3/3.7.0

# go to the directory I submitted the job from
cd $SLURM_SUBMIT_DIR

# so far I have one production run. I'll need 5 cores for that halo finding.
# that leaves the rest for other things, so I can use 64 parallel ranks
make dirs
remora make -j64
# &>> $SLURM_JOB_NAME.stdout.$SLURM_JOB_ID