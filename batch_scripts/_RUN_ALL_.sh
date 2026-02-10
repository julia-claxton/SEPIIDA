#!/bin/bash

for i in /projects/$USER/SEPIIDA/source/batch_scripts/*keV.sh; do
  sbatch --quiet "$i"
done
watch -n1 squeue --format=\"%.18i %.15P %.30j %.10u %.10T %.13M %.13l %.8D %R\" --me