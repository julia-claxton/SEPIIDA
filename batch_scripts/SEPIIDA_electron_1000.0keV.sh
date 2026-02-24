#!/bin/bash

#SBATCH --job-name SEPIIDA_electron_1000.0keV
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 40
#SBATCH --time 1-00:00:00
#SBATCH --output /projects/jucl6426/SEPIIDA/results/log_SEPIIDA_electron_1000.0keV.out
#SBATCH --qos=preemptable
#SBATCH --exclude=bhpc-c5-u7-19,bhpc-c5-u7-22,bhpc-c5-u7-23
#SBATCH --requeue
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jucl6426@colorado.edu

# Terminate on any non-zero exit status
set -e

# Run simulation
cd /projects/jucl6426/SEPIIDA/build/
./SEPIIDA 100000 e- 1000.0

# Copy results to safe folder
cp /projects/jucl6426/SEPIIDA/build/results/lat_45deg_input_450km/electron_input_1000.0keV_100000particles_*_spectra.csv /projects/jucl6426/SEPIIDA/results

