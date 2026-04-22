using Statistics, LinearAlgebra
using Glob
using Printf

function write_job_script(qos, n_particles, input_particle, energy, pa; prefix = "", flags = "")
    # Prepare input
    input_particle_longname = input_particle == "e-" ? "electron" : input_particle
    energy_string = @sprintf "%.1f" energy
    n_particles_string = @sprintf "%i" n_particles
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
    job_name = "SEPIIDA_$(prefix)_$(input_particle_longname)_$(energy_string)keV_$(pa)deg_$(n_particles_string)particles"
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
    #SBATCH --exclude=bhpc-c5-u7-19,bhpc-c5-u7-22,bmem-rico1
    #SBATCH --requeue
    #SBATCH --mail-type=BEGIN,END,FAIL
    #SBATCH --mail-user=julia.claxton@colorado.edu

    # Terminate on any non-zero exit status
    set -e

    # Print job ID
    echo "Job ID: \$SLURM_JOB_ID"
    echo "Node ID: \$SLURMD_NODENAME"

    # Load modules
    module purge
    module load gcc/14.2.0

    # Run simulation
    cd /projects/jucl6426/SEPIIDA/build/
    ./SEPIIDA $(n_particles_string) $(input_particle) $(energy_string) $(pa) $(flags)

    # Copy results to safe folder
    cp /projects/jucl6426/SEPIIDA/build/results/$(prefix)*$(input_particle_longname)_input*$(energy_string)keV_$(pa)deg_$(n_particles_string)particles* /projects/jucl6426/SEPIIDA/results
    """
    )
    close(file)
    println("Wrote batch_scripts/$(job_name).sh.")
end

# Remove existing jobscripts
rm.(glob("*.sh", @__DIR__))

write_job_script("preemptable", 1e3, "e-", 10_000, 0.0, 
    prefix = "segfault_test",
    flags = "
        -backscatter_altitude 451.0
        -brem_splitting 100
        -min_energy_eV 10
    "
)



#=
for lat in -90:30:90
    for energy in [100.0, 1_000.0, 10_000.0]
        for pitch_angle in [0.1, 45, 70, 90]
            write_job_script("preemptable", 1e5, "e-", energy, pitch_angle, 
                prefix = "gamma_variation_lat$(lat)",
                flags = "
                    -magnetic_model jrm33
                    -atmosphere_filename jupiter_atmosphere_profile.csv
                    -backscatter_altitude 451.0
                    -brem_splitting 100
                    -min_energy_eV 10
                    -lat $(lat)
                "
            )
        end
    end
end
=#