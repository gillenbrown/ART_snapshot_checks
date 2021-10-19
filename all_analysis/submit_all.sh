#!/bin/bash
#SBATCH --job-name=all_analysis
#SBATCH --account=TG-AST200017
#SBATCH --mail-user=gillenb@umich.edu
#SBATCH --mail-type=ALL
#SBATCH --no-requeue
#SBATCH --partition=skx-normal
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --cpus-per-task=1
#SBATCH --output=stdout_%x_%j

# This took a lot of messing with to get right. Do NOT use `module reset`. This
# will stop conda from recognizing my distribution somehow.
# We do need to source bashrc, as it holds the code to activate conda and
# activate the correct env
# module load remora
source ~/.bashrc
# Make sure the system python is not accessible
unset PYTHONPATH
# double check these variables
# echo $(which python)
# echo $PYTHONPATH

# slurm setup, copied from example at $LAUNCHER_DIR/extras/batch-scripts
export LAUNCHER_RMI=SLURM
export LAUNCHER_PLUGIN_DIR=$LAUNCHER_DIR/plugins
# setup of my tasks
export LAUNCHER_WORKDIR="/work/06912/tg862118/stampede2/ART_snapshot_checks/all_analysis"
export LAUNCHER_JOB_FILE="job_list.txt"
$LAUNCHER_DIR/paramrun
