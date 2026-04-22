#!/bin/bash

#SBATCH --job-name SEPIIDA_wtf_electron_10000.0keV_0.1deg_10000particles
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 40
#SBATCH --time 1-00:00:00
#SBATCH --output /projects/jucl6426/SEPIIDA/results/SEPIIDA_wtf_electron_10000.0keV_0.1deg_10000particles.log
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
./SEPIIDA 10000 e- 10000.0 0.1 -magnetic_model jrm33 -atmosphere_filename jupiter_atmosphere_profile.csv -backscatter_altitude 451.0 -brem_splitting 100 -min_energy_eV 10 -prefix wtf

# Copy results to safe folder
cp /projects/jucl6426/SEPIIDA/build/results/wtf*electron_input*10000.0keV_0.1deg_10000particles* /projects/jucl6426/SEPIIDA/results
