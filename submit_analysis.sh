#!/bin/bash
#SBATCH --job-name=analysis
#SBATCH --account=TG-AST200017
#SBATCH --mail-user=gillenb@umich.edu
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH --partition=skx-dev
#SBATCH --time=0:10:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --cpus-per-task=1
#SBATCH --output=stdout_%x_%j

# This took a lot of messing with to get right. Do NOT use `module reset`. This
# will stop conda from recognizing my distribution somehow.
# We do need to source bashrc, as it holds the code to activate conda and
# activate the correct env
source ~/.bashrc

# double check these variables
# echo $(which python)
# echo $PYTHONPATH

# I do run the halo finding, but just as a check. Nothing should ever run in reality.
make dirs
ibrun -n 1 -o 0 make halos
make -j5
# &>> $SLURM_JOB_NAME.stdout.$SLURM_JOB_ID
