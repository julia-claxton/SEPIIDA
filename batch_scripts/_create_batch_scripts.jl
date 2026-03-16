using Statistics, LinearAlgebra
using Glob
using Printf

function main()
    # Remove existing jobscripts
    rm.(glob("*keV*.sh", @__DIR__))
    
    # Write new scripts
    write_job_script("preemptable", 1e5, "e-", 1000, 180, 
        prefix = "jupiter",
        flags = "
            -brem_splitting 100
            -magnetic_model jrm33
            -atmosphere_filename jupiter_atmosphere_profile.csv
        "
    )
    write_job_script("preemptable", 1e5, "e-", 1000, 0, 
        prefix = "earth",
        flags = "
            -brem_splitting 100
            -magnetic_model igrf2025
            -atmosphere_filename msis_earth_atmosphere_profile.csv
        "
    )
end

function write_job_script(qos, n_particles, input_particle, energy, pa; prefix = "", flags = "")
    # Prepare input
    input_particle_longname = input_particle == "e-" ? "electron" : input_particle
    energy_string = @sprintf "%.1f" energy
    if (prefix != "") && (!contains(flags, "-result_prefix"))
        flags = "$(flags) -prefix $(prefix)"
    end
    flags = replace(flags, "\n" => " ")
    qos == "blanca-lair" ? time_limit = "4-00:00:00" : time_limit = "1-00:00:00"

    # Remove double+ and leading whitespaces from flags
    while contains(flags, "  ")
        flags = replace(flags, "  " => " ")
    end
    if startswith(flags, " ")
        flags = flags[2:end]
    end

    # Write file
    job_name = "SEPIIDA_$(prefix)_$(input_particle_longname)_$(energy_string)keV_$(pa)deg_$(n_particles)particles"
    file = open("$(@__DIR__)/$(job_name).sh", "w")
    print(file,
    """
    #!/bin/bash

    #SBATCH --job-name $(job_name)
    #SBATCH --nodes 1
    #SBATCH --ntasks-per-node 40
    #SBATCH --time $(time_limit)
    #SBATCH --output /projects/jucl6426/SEPIIDA/results/$(job_name).log
    #SBATCH --qos=$(qos)
    #SBATCH --exclude=bhpc-c5-u7-19,bhpc-c5-u7-22
    #SBATCH --requeue
    #SBATCH --mail-type=BEGIN,END,FAIL
    #SBATCH --mail-user=julia.claxton@colorado.edu

    # Terminate on any non-zero exit status
    set -e

    # Print job ID
    echo "Job ID: \$SLURM_JOB_ID"

    # Load modules
    module purge
    module load gcc/14.2.0

    # Run simulation
    cd /projects/jucl6426/SEPIIDA/build/
    ./SEPIIDA $(n_particles) $(input_particle) $(energy_string) $(pa) $(flags)

    # Copy results to safe folder
    cp /projects/jucl6426/SEPIIDA/build/results/$(prefix)*$(input_particle_longname)_input*$(energy_string)keV_$(pa)deg_$(n_particles)particles* /projects/jucl6426/SEPIIDA/results
    """
    )
    close(file)
    println("Wrote batch_scripts/$(job_name).sh.")
end

main()