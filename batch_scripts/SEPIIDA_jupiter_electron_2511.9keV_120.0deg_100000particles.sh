#!/bin/bash

#SBATCH --job-name SEPIIDA_jupiter_electron_2511.9keV_120.0deg_100000particles
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 40
#SBATCH --time 1-00:00:00
#SBATCH --output /projects/jucl6426/SEPIIDA/results/SEPIIDA_jupiter_electron_2511.9keV_120.0deg_100000particles.log
#SBATCH --qos=preemptable
#SBATCH --exclude=bhpc-c5-u7-19,bhpc-c5-u7-22
#SBATCH --requeue
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=julia.claxton@colorado.edu

# Terminate on any non-zero exit status
set -e

# Print job ID
echo "Job ID: $SLURM_JOB_ID"

# Load modules
module purge
module load gcc/14.2.0

# Run simulation
cd /projects/jucl6426/SEPIIDA/build/
./SEPIIDA 100000 e- 2511.9 120.0 -brem_splitting 100 -magnetic_model jrm33 -atmosphere_filename jupiter_atmosphere_profile.csv -prefix jupiter

# Copy results to safe folder
cp /projects/jucl6426/SEPIIDA/build/results/jupiter*electron_input*2511.9keV_120.0deg_100000particles* /projects/jucl6426/SEPIIDA/results
