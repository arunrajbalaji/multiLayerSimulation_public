#!/bin/bash

#SBATCH -J lowV
#SBATCH -D /home/mani/abalaji/multiLayerEPNP/src/
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 2:00:00
#SBATCH -o %x_%A.out
#SBATCH --mail-type=NONE
#SBATCH --mail-user=abalaji@stanford.edu

module load apps/matlab/R2023a

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

matlab -batch "main('/fastscratch/abalaji/shuttlingPaper/h2so4_fbm_koh_refined/output/appliedVoltage_0_/','inputFile_141875.json')"

echo
echo -------------------------
echo Execution Ending Time:
echo `date`
echo -------------------------
echo

