#!/bin/bash

#SBATCH --job-name AG4_proton_83915824.8keV
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 40
#SBATCH --time 3-00:00:00
#SBATCH --output /projects/jucl6426/SEPIIDA/results/log_AG4_proton_83915824.8keV.out
#SBATCH --qos=blanca-lair
#SBATCH --exclude=bhpc-c5-u7-19,bhpc-c5-u7-22,bhpc-c5-u7-23
#SBATCH --requeue
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jucl6426@colorado.edu

# Terminate on any non-zero exit status
set -e

# Run simulation
cd /projects/jucl6426/SEPIIDA/build/
./SEPIIDA 100000 proton 83915824.8

# Copy results to safe folder
cp /projects/jucl6426/SEPIIDA/build/results/mlat_45deg_input_450km/proton_input_83915824.8keV_100000particles_*_spectra.csv /projects/jucl6426/SEPIIDA/results

