#!/bin/bash

#SBATCH --job-name SEPIIDA_patest_electron_10000.0keV_5deg_1000particles
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 40
#SBATCH --time 1-00:00:00
#SBATCH --output /projects/jucl6426/SEPIIDA/results/SEPIIDA_patest_electron_10000.0keV_5deg_1000particles.log
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
./SEPIIDA 1000 e- 10000.0 5 -magnetic_model jrm33 -atmosphere_filename jupiter_atmosphere_profile.csv -prefix patest

# Copy results to safe folder
cp /projects/jucl6426/SEPIIDA/build/results/patest*electron_input*10000.0keV_5deg_1000particles* /projects/jucl6426/SEPIIDA/results
