echo "Beginning job submission script."

voltages=($(seq -2.0 0.025 0.0))

for appliedVoltage in "${voltages[@]}"; do
    
    echo "Voltage submitted is $appliedVoltage"

    sbatch ~/multiLayerEPNP/src/run.sh 1.0e-11 1.0e-5 $appliedVoltage

done
