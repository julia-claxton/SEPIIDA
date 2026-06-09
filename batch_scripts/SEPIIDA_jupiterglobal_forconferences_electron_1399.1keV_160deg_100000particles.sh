#!/bin/bash

#SBATCH --job-name SEPIIDA_jupiterglobal_forconferences_electron_1399.1keV_160deg_100000particles
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 40
#SBATCH --time 1-00:00:00
#SBATCH --output /projects/jucl6426/SEPIIDA/results/SEPIIDA_jupiterglobal_forconferences_electron_1399.1keV_160deg_100000particles.log
#SBATCH --qos=preemptable
#SBATCH --exclude=bhpc-c5-u7-19,bhpc-c5-u7-22,bmem-rico1
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
./SEPIIDA 100000 e- 1399.1 160 -magnetic_model jrm33 -atmosphere_filename jupiter_gram.csv -injection_altitude 990.0 -backscatter_altitude 991.0 -brem_splitting 500 -min_energy_eV 10 -lat 85 -cache_radius_km 1.0 -prefix jupiterglobal_forconferences

# Copy results to safe folder
cp /projects/jucl6426/SEPIIDA/build/results/jupiterglobal_forconferences*electron_input*1399.1keV_160deg_100000particles* /projects/jucl6426/SEPIIDA/results
set +x
