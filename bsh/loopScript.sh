echo "Beginning job submission script."

jobsList=($(seq 1 1 170))

for jobNum in "${jobsList[@]}"; do
  sbatch ~/multiLayerEPNP/src/run.sh
done

