#!/bin/bash -l

#SBATCH -J FBMC_koh0
#SBATCH -D /home/mani/abalaji/multiLayerEPNP/src/
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 2:00:00
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

#FOLDER_1="/fastscratch/abalaji/june2023_correctPorosity/alpha_1.8/koh0/input.json"
#FOLDER_2="/fastscratch/abalaji/june2023_correctPorosity/alpha_1.8/koh0.25/input.json"
FOLDER_3="/fastscratch/abalaji/FBPM_curveFitting/koh0/ref_koh1/input.json"
FOLDER_4="/fastscratch/abalaji/FBPM_curveFitting/koh0/L1e-3/input.json"
FOLDER_5="/fastscratch/abalaji/FBPM_curveFitting/koh0/L1.5e-3/input.json"
LOOPSTRING="{'appliedVoltage', -0.50, 0.01, 1.00}"

#matlab -batch "runParallelJob('${FOLDER_1}', $SLURM_ARRAY_TASK_ID, $SLURM_ARRAY_JOB_ID, ${LOOPSTRING})"
#matlab -batch "runParallelJob('${FOLDER_2}', $SLURM_ARRAY_TASK_ID, $SLURM_ARRAY_JOB_ID, ${LOOPSTRING})"
matlab -batch "runParallelJob('${FOLDER_3}', $SLURM_ARRAY_TASK_ID, $SLURM_ARRAY_JOB_ID, ${LOOPSTRING})"
matlab -batch "runParallelJob('${FOLDER_4}', $SLURM_ARRAY_TASK_ID, $SLURM_ARRAY_JOB_ID, ${LOOPSTRING})"
matlab -batch "runParallelJob('${FOLDER_5}', $SLURM_ARRAY_TASK_ID, $SLURM_ARRAY_JOB_ID, ${LOOPSTRING})"

echo
echo -------------------------
echo Execution Ending Time:
echo `date`
echo -------------------------
echo

