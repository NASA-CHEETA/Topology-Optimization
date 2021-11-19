#!/bin/bash
##############################SLURM SCRIPT ######################################



#SBATCH --time=4:00:00
#SBATCH --job-name="hex"
#SBATCH --mail-user=gkdas2@illinois.edu
#SBATCH --mail-type=BEGIN,END
#SBATCH --output=outputprint.dat
#SBATCH --error=FAILURE.dat
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --partition=secondary-Eth
#####################################


# Change to the directory from which the batch job was submitted
export SLURM_SUBMIT_DIR=/home/gkdas2/scratch/workspace/TOPOPT_BENCHMARK_STUDIES/ccxVsmatlab/ccx
# Change to the job directory
cd $SLURM_SUBMIT_DIR

export OMP_NUM_THREADS=20
time python main.py
