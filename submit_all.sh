#!/bin/bash
#SBATCH --job-name=analysis_all
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
# -u stops output buffering
python -u analysis_all.py

