echo "Beginning job submission script."

voltages=($(seq -1.0 0.025 0.0))

for appliedVoltage in "${voltages[@]}"; do
    
    echo "Voltage submitted is $appliedVoltage"

    sbatch ~/multiLayerEPNP/src/run.sh 2.0e-10 1.0e-5 $appliedVoltage

done
