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
rm.(glob("*keV*.sh", @__DIR__))

write_job_script("preemptable", 100000, "e-", 10000.0, 0.1, 
    prefix = "earth_xray_for_bob",
    flags = "
        -magnetic_model igrf2025
        -atmosphere_filename msis_earth_atmosphere_profile.csv
        -backscatter_altitude 451.0
        -brem_splitting 100
        -min_energy_eV 10
    "
)



#=
# Write new scripts
alpha_energies = 1e6 .* [1.467029285016976, 1.7817875218749928, 2.1434580789902244, 2.5531026851640726, 3.0195499666152683, 3.5436236292778935, 4.125627948346392, 4.774463797977579, 5.490273042738451, 6.282166640382038, 7.159519110086783, 8.122386626142525, 9.170597549523674, 10.313477444726429, 11.560505948882511, 12.911525353931168, 14.376071249205166, 15.963767866189537, 17.674457592576562, 20.53486106783365, 26.080040217819565, 35.821267113820014, 52.667749722531134, 78.32854513709988, 114.50220914157535, 165.98409421382254]
proton_energies = 1e6 .* [0.4919968561663142, 0.6202564919316679, 0.7632187244759833, 0.9246328837346689, 1.1043104921012188, 1.3019751267372308, 1.521935705828299, 1.7640765698686003, 2.032968417343461, 2.3285487124854978, 2.650701655257517, 3.0041418809864373, 3.388794404104349, 3.809461740784525, 4.271004102390247, 4.773387063190187, 5.316554218799155, 5.9053912772825905, 6.544811688653951, 7.234780943620419, 7.98022843236318, 8.78609801629206, 9.652363544291253, 11.097336568805867, 13.890406704078622, 18.78300886683073, 27.226294269353872, 40.07138487336281, 58.168090056375284, 83.9158247545856]

for i in eachindex(alpha_energies)
    write_job_script("preemptable", 1e5, "alpha", alpha_energies[i], 0.1, 
        prefix = "GCR_muons",
        flags = "
            -magnetic_model igrf2025
            -atmosphere_filename msis_earth_atmosphere_profile.csv
            -backscatter_altitude 451.0
            -brem_splitting 1
        "
    )
end

for i in eachindex(proton_energies)
    write_job_script("preemptable", 1e5, "proton", proton_energies[i], 0.1, 
        prefix = "GCR_muons",
        flags = "
            -magnetic_model igrf2025
            -atmosphere_filename msis_earth_atmosphere_profile.csv
            -backscatter_altitude 451.0
            -brem_splitting 1
        "
    )
end
=#