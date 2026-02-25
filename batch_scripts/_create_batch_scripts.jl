using Statistics, LinearAlgebra
using Glob
using Printf

function main()
  number_of_particles = 10_000 #100_000  # Number of particles to input
  particle = "e-" # "e-" = electrons, "proton" = protons, "gamma" = photons
  energies_to_simulate = 1e3
  pitch_angles_to_simulate = 0

  magnetic_models = ["earth_tilted_dipole", "jrm33"]
  atmospheres = ["msis_earth_atmosphere_profile.csv", "crude_jupiter_atmosphere_profile.csv"]

  # Create shell scripts
  rm.(glob("*keV.sh", @__DIR__))
  written = 0
  skipped = 0
  lair = 0

  for i in eachindex(magnetic_models)
    
    i == 1 ? planet = "earth" : planet = "jupiter"

    input_particle_longname = particle == "e-" ? "electron" : particle
    energy_string = "1000.0" #@sprintf "%.1f" E
    job_name = "SEPIIDA_$(planet)_$(input_particle_longname)_$(energy_string)keV"
    qos = "preemptable"
    time_limit = "1-00:00:00"

    # Send long-runtime beams to blanca-lair
    #=
    if E ≥ 2e7
      qos = "blanca-lair"
      time_limit = "3-00:00:00"
      lair += 1
    end
    =#

    println("TODO ADD TIMESTAMPS TO RESULTS AND LOGS")


    file = open("$(@__DIR__)/$(job_name).sh", "w")
    println(file,
    """
    #!/bin/bash

    #SBATCH --job-name $(job_name)
    #SBATCH --nodes 1
    #SBATCH --ntasks-per-node 40
    #SBATCH --time $(time_limit)
    #SBATCH --output /projects/jucl6426/SEPIIDA/results/log_$(input_particle_longname)_input_$(energy_string)keV_$(number_of_particles)particles.log
    #SBATCH --qos=$(qos)
    #SBATCH --exclude=bhpc-c5-u7-19,bhpc-c5-u7-22
    #SBATCH --requeue
    #SBATCH --mail-type=BEGIN,END,FAIL
    #SBATCH --mail-user=jucl6426@colorado.edu

    # Terminate on any non-zero exit status
    set -e

    # Run simulation
    cd /projects/jucl6426/SEPIIDA/build/
    ./SEPIIDA 10000 e- 1000.0 0.0 -magnetic_model $(magnetic_models[i]) -atmosphere_filename $(atmospheres[i]) -lat 70 -brem_splitting 20 -altitude_offset 0 -injection_altitude 500 -backscatter_altitude 505

    # Copy results to safe folder
    find /projects/jucl6426/SEPIIDA/build/results/lat_70deg_input_500km/ -name '$(input_particle_longname)_input_1000.0keV_0.0deg_$(number_of_particles)particles_*.csv' -exec sh -c 'mv {} $(planet)_\$(basename {})' \\;
    cp /projects/jucl6426/SEPIIDA/build/results/lat_70deg_input_500km/$(planet)_$(input_particle_longname)_input_1000.0keV_0.0deg_$(number_of_particles)particles_* /projects/jucl6426/SEPIIDA/results
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