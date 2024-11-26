#!/bin/bash -l

#SBATCH -J fbm
#SBATCH -D /home/mani/abalaji/multiLayerEPNP/src/
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 3:00:00
#SBATCH -o O_%x_%A_%a.out
#SBATCH --mail-type=NONE
#SBATCH --mail-user=abalaji@stanford.edu
#SBATCH --array 1-101

module load apps/matlab/R2023a

##########################
### Beginning of Execution
##########################

echo
echo The master node of this job is `hostname`
echo The working directory is `pwd`
echo The slurm task ID is $SLURM_ARRAY_TASK_ID
echo
echo -------------------------
echo Execution Starting Time:
echo `date`
echo -------------------------
echo

FOLDER_1="/fastscratch/abalaji/shuttlingPaper/hcl_fbm_koh_withGroups/input.json"
LOOPSTRING="{'appliedVoltage', 0.00, 0.01, 1.00}"

matlab -batch "runParallelJob('${FOLDER_1}', $SLURM_ARRAY_TASK_ID, $SLURM_ARRAY_JOB_ID, ${LOOPSTRING})"

echo
echo -------------------------
echo Execution Ending Time:
echo `date`
echo -------------------------
echo

