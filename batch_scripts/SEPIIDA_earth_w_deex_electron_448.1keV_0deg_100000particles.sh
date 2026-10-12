#!/bin/bash

#SBATCH --job-name SEPIIDA_earth_w_deex_electron_448.1keV_0deg_100000particles
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 40
#SBATCH --time 1-00:00:00
#SBATCH --output /projects/jucl6426/SEPIIDA/results/SEPIIDA_earth_w_deex_electron_448.1keV_0deg_100000particles.log
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
./SEPIIDA 100000 e- 448.1 0 -magnetic_model igrf2025 -atmosphere_filename msis_earth.csv -injection_altitude 350.0 -backscatter_altitude 351.0 -brem_splitting 5000 -min_energy_eV 1000 -cache_radius_km 1.0 -prefix earth_w_deex

# Copy results to safe folder
cp /projects/jucl6426/SEPIIDA/build/results/earth_w_deex*electron_input*448.1keV_0deg_100000particles* /projects/jucl6426/SEPIIDA/results
set +x
