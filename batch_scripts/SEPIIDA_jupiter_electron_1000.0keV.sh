#!/bin/bash

#SBATCH --job-name SEPIIDA_jupiter_electron_1000.0keV
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 40
#SBATCH --time 1-00:00:00
#SBATCH --output /projects/jucl6426/SEPIIDA/results/log_electron_input_1000.0keV_10000particles.log
#SBATCH --qos=preemptable
#SBATCH --exclude=bhpc-c5-u7-19,bhpc-c5-u7-22
#SBATCH --requeue
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jucl6426@colorado.edu

# Terminate on any non-zero exit status
set -e

# Run simulation
cd /projects/jucl6426/SEPIIDA/build/
./SEPIIDA 10000 e- 1000.0 0.0 -magnetic_model jrm33 -atmosphere_filename crude_jupiter_atmosphere_profile.csv -lat 70 -brem_splitting 20 -altitude_offset 0 -injection_altitude 500 -backscatter_altitude 505

# Copy results to safe folder
find /projects/jucl6426/SEPIIDA/build/results/lat_70deg_input_500km/ -name 'electron_input_1000.0keV_0.0deg_10000particles_*.csv' -exec sh -c 'mv {} jupiter_$(basename {})' \;
cp /projects/jucl6426/SEPIIDA/build/results/lat_70deg_input_500km/jupiter_electron_input_1000.0keV_0.0deg_10000particles_* /projects/jucl6426/SEPIIDA/results

