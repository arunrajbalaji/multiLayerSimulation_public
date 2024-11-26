#!/bin/bash -l

#SBATCH -J eps_R
#SBATCH -D /home/mani/abalaji/multiLayerEPNP/src/
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 00:10:00
#SBATCH -o MLEPNP_%A_%a.out
#SBATCH --mail-type=NONE
#SBATCH --mail-user=abalaji@stanford.edu
#SBATCH --array 1-11

module load apps/matlab/R2021a

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

INPUT="/fastscratch/abalaji/reservoirPorosityVariation_fumasep/poro_1_noRxn/input.json"
LOOPARGSTRING="{'appliedVoltage', -1.0, 0.1, 0.0}"
matlab -batch "runParallelJob('${INPUT}', $SLURM_ARRAY_TASK_ID, ${LOOPARGSTRING})"

echo
echo -------------------------
echo Execution Ending Time:
echo `date`
echo -------------------------
echo
