#!/bin/bash

#SBATCH --job-name SEPIIDA_argo_final_electron_71320.3keV_100deg_100000particles
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 40
#SBATCH --time 1-00:00:00
#SBATCH --output /projects/jucl6426/SEPIIDA/results/SEPIIDA_argo_final_electron_71320.3keV_100deg_100000particles.log
#SBATCH --qos=preemptable
#SBATCH --exclude=bhpc-c5-u7-20,bhpc-c5-u7-21,bhpc-c5-u7-22,bmem-rico1
#SBATCH --requeue
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=julia.claxton@colorado.edu

# Terminate on any non-zero exit status
set -e


# Print job ID
echo "Job ID: $SLURM_JOB_ID"
echo "Node ID: $SLURMD_NODENAME"

# Load modules
module purge
module load gcc/14.2.0

# Run simulation
set -x
cd /projects/jucl6426/SEPIIDA/build/
./SEPIIDA 100000 e- 71320.3 100 -magnetic_model jrm33 -atmosphere_filename ARGO_Jupiter.csv -injection_altitude 500.0 -backscatter_altitude 501.0 -brem_splitting 10 -min_energy_eV 1000 -lat 85 -cache_radius_km 1.0 -prefix argo_final

# Copy results to safe folder
cp /projects/jucl6426/SEPIIDA/build/results/argo_final*electron_input*71320.3keV_100deg_100000particles* /projects/jucl6426/SEPIIDA/results
set +x
