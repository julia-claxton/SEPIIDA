#!/bin/bash

#SBATCH --job-name SEPIIDA_jlc_electron_1000.0keV_102deg_1000particles
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 40
#SBATCH --time 1-00:00:00
#SBATCH --output /projects/jucl6426/SEPIIDA/results/SEPIIDA_jlc_electron_1000.0keV_102deg_1000particles.log
#SBATCH --qos=preemptable
#SBATCH --exclude=bhpc-c5-u7-19,bhpc-c5-u7-22,bmem-rico1
#SBATCH --requeue
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=julia.claxton@colorado.edu

# Terminate on any non-zero exit status
set -e

# Echo commands to log
set -x

# Print job ID
echo "Job ID: $SLURM_JOB_ID"
echo "Node ID: $SLURMD_NODENAME"

# Load modules
module purge
module load gcc/14.2.0

# Run simulation
cd /projects/jucl6426/SEPIIDA/build/
./SEPIIDA 1000 e- 1000.0 102 -magnetic_model jrm33 -atmosphere_filename jupiter_gram.csv -injection_altitude 990.0 -backscatter_altitude 991.0 -brem_splitting 1 -min_energy_eV 10 -lat 85 -prefix jlc

# Copy results to safe folder
cp /projects/jucl6426/SEPIIDA/build/results/jlc*electron_input*1000.0keV_102deg_1000particles* /projects/jucl6426/SEPIIDA/results
