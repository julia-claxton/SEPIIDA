using Statistics, LinearAlgebra
using Glob
using Printf

include("$(dirname(@__DIR__))/postprocessing_scripts/SEPIIDA_Reader.jl")

function write_job_script(qos, n_particles, input_particle, energy, pa; prefix = "", flags = "")
    # Prepare input
    input_particle_longname = input_particle == "e-" ? "electron" : input_particle
    energy_string = @sprintf "%.1f" energy
    n_particles_string = @sprintf "%i" n_particles
    if (prefix != "")
        flags = "$(flags) -prefix $(prefix)"
    end
    flags = replace(flags, "\n" => " ")
    qos == "blanca-lair" ? time_limit = "5-00:00:00" : time_limit = "1-00:00:00"

    # Remove double+ and leading whitespaces from flags
    while contains(flags, "  ")
        flags = replace(flags, "  " => " ")
    end
    if startswith(flags, " ")
        flags = flags[2:end]
    end

    # Write file
    job_prefix = prefix == "" ? "" : "$(prefix)_"
    job_name = "SEPIIDA_$(job_prefix)$(input_particle_longname)_$(energy_string)keV_$(pa)deg_$(n_particles_string)particles"
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
    set -x
    cd /projects/jucl6426/SEPIIDA/build/
    ./SEPIIDA $(n_particles_string) $(input_particle) $(energy_string) $(pa) $(flags)

    # Copy results to safe folder
    cp /projects/jucl6426/SEPIIDA/build/results/$(prefix)*$(input_particle_longname)_input*$(energy_string)keV_$(pa)deg_$(n_particles_string)particles* /projects/jucl6426/SEPIIDA/results
    set +x
    """
    )
    close(file)
    println("Wrote batch_scripts/$(job_name).sh")
end

function status()
    jobs = glob("*.sh", @__DIR__)
    println("\nThere are $(length(jobs)) jobs.")
end

# Remove existing jobscripts
rm.(glob("*.sh", @__DIR__))

beamlist = get_beamlist("/Users/luna/Research/geant4/SEPIIDA/results/2026-07-20--15.14_jupiterfinal")
existing_e = round.(energy_list(beamlist), digits = 1)
existing_pa = pitch_angle_list(beamlist)
existing_beams = collect(zip(existing_e, existing_pa))

# Write new jobs
for E in logrange(30, 1e5, 20)
    for pa in [105, 110:10:140..., 180]
        if (round(E, digits = 1), pa) ∈ existing_beams; continue; end

        N = 1e5
        split_factor = E > 1e3 ? 1000 : 10
        write_job_script("blanca-lair", N, "e-", E, pa, 
            prefix = "jupiterglobal_forconferences",
            flags = "
                -magnetic_model jrm33
                -atmosphere_filename jupiter_gram.csv
                -injection_altitude 990.0
                -backscatter_altitude 991.0
                -brem_splitting $(split_factor)
                -min_energy_eV 10
                -lat 85
                -cache_radius_km 1.0
            "
        )
    end
end


status()