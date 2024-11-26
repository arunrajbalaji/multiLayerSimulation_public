#!/bin/bash

#SBATCH -J N117Long
#SBATCH -D /home/mani/abalaji/multiLayerEPNP/src/
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 72:00:00
#SBATCH -o %x_%A.out
#SBATCH -w compute-3-33
#SBATCH --mail-type=NONE
#SBATCH --mail-user=abalaji@stanford.edu

module load apps/matlab/R2021a

##########################
### Beginning of Execution
##########################

echo
echo The master node of this job is `hostname`
echo The working directory is `pwd`
echo
echo -------------------------
echo Execution Starting Time:
echo `date`
echo -------------------------
echo

matlab -batch "singleNodeParallel"

echo
echo -------------------------
echo Execution Ending Time:
echo `date`
echo -------------------------
echo

