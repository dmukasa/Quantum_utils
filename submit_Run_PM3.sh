#!/bin/bash

# Submit this script with: sbatch <this-filename>

#SBATCH --time=99:00:00   # walltime
#SBATCH --ntasks=16   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=10G   # memory per CPU core
#SBATCH -J "Dimers"   # job name
#SBATCH --mail-user=dmukasa@caltech.edu   # email address

# Notify at the beginning, end of job and on failure.
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module load openmpi/2.1.6

# RUN FILE
sh /groups/GaoGroup/dmukasa/ML_DFT/Dimer_editor/Run_PM3.sh
