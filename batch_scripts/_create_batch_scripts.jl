using Statistics, LinearAlgebra
using Glob
using Printf

function main()
  number_of_particles = 100_000  # Number of particles to input
  particle = "e-" # "e-" = electrons, "proton" = protons, "gamma" = photons
  energies_to_simulate = 1e3

  # Create shell scripts
  rm.(glob("*keV.sh", @__DIR__))
  written = 0
  skipped = 0
  lair = 0

  for E in energies_to_simulate
    input_particle_longname = particle == "e-" ? "electron" : particle
    energy_string = @sprintf "%.1f" E
    job_name = "SEPIIDA_$(input_particle_longname)_$(energy_string)keV"
    qos = "preemptable"
    time_limit = "1-00:00:00"

    # Send long-runtime beams to blanca-lair
    if E ≥ 2e7
      qos = "blanca-lair"
      time_limit = "3-00:00:00"
      lair += 1
    end

    file = open("$(@__DIR__)/$(job_name).sh", "w")
    println(file,
    """
    #!/bin/bash

    #SBATCH --job-name $(job_name)
    #SBATCH --nodes 1
    #SBATCH --ntasks-per-node 40
    #SBATCH --time $(time_limit)
    #SBATCH --output /projects/jucl6426/SEPIIDA/results/log_$(job_name).out
    #SBATCH --qos=$(qos)
    #SBATCH --exclude=bhpc-c5-u7-19,bhpc-c5-u7-22,bhpc-c5-u7-23
    #SBATCH --requeue
    #SBATCH --mail-type=BEGIN,END,FAIL
    #SBATCH --mail-user=jucl6426@colorado.edu

    # Terminate on any non-zero exit status
    set -e

    # Run simulation
    cd /projects/jucl6426/SEPIIDA/build/
    ./SEPIIDA $(number_of_particles) $(particle) $(energy_string)

    # Copy results to safe folder
    cp /projects/jucl6426/SEPIIDA/build/results/mlat_45deg_input_450km/$(input_particle_longname)_input_$(energy_string)keV_$(number_of_particles)particles_* /projects/jucl6426/SEPIIDA/results
    """
    )
    close(file)

    written += 1
  end

  println("$(written) files written.")
  println("$(skipped) files skipped.")
  println("$(lair) files sent to blanca-lair.")
end

main()