#!/bin/bash
##############################SLURM SCRIPT ######################################



#SBATCH --time=4:00:00
#SBATCH --job-name="raeghrf0"
#SBATCH --mail-user=gkdas2@illinois.edu
#SBATCH --mail-type=BEGIN,END
#SBATCH --output=outputprint.dat
#SBATCH --error=FAILURE.dat
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --partition=secondary-eth
#####################################


# Change to the directory from which the batch job was submitted
export SLURM_SUBMIT_DIR=/home/gkdas2/scratch/workspace/RAE_AEROELASTIC/RAE_ghrvc_fo
# Change to the job directory
cd $SLURM_SUBMIT_DIR

export OMP_NUM_THREADS=20

#time ccx_2.15_mod -i cant3d.inp -p 3 -r 1.5 -v 0.12
time python TO_config.py
