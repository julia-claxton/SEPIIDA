#!/bin/bash

#SBATCH --job-name SEPIIDA_segfault_test_electron_10000.0keV_0.0deg_500particles
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 40
#SBATCH --time 1-00:00:00
#SBATCH --output /projects/jucl6426/SEPIIDA/results/SEPIIDA_segfault_test_electron_10000.0keV_0.0deg_500particles.log
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
cd /projects/jucl6426/SEPIIDA/build/
./SEPIIDA 500 e- 10000.0 0.0 -magnetic_model igrf2025 -atmosphere_filename msis_earth_atmosphere_profile.csv -backscatter_altitude 451.0 -brem_splitting 1 -min_energy_eV 10 -prefix segfault_test

# Copy results to safe folder
cp /projects/jucl6426/SEPIIDA/build/results/segfault_test*electron_input*10000.0keV_0.0deg_500particles* /projects/jucl6426/SEPIIDA/results
