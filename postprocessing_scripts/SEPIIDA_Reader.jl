using Statistics, LinearAlgebra          # Core math
using BenchmarkTools, Profile, TickTock  # Debugging
using NPZ, DelimitedFiles                # File interactions
using Glob
using Printf
using BasicInterpolators
using Base.Threads
using Plots, Plots.PlotMeasures
    default(
        dpi = 300,
        framestyle = :box,
        tickdirection = :out,
        label = false
    )

const TOP_LEVEL = dirname(@__DIR__)
include("$(@__DIR__)/General_Functions.jl")

"""
docstring

Written by Julia Claxton. Contact: julia.claxton@colorado.edu
Released under MIT License.
"""

struct BeamInfo
    prefix::String
    particle::String
    injection_alt::Float64
    backscatter_alt::Float64
    latitude::Float64
    energy::Float64
    pitch_angle::Float64
    n_particles::Int64
    dir::String
    base_filename::String
end

struct Beam
    info::BeamInfo

    energy_deposition::Vector{Float64}
    electron_production::Vector{Base.Float64}
    backscatter::Dict{String, Any}
    spectra::Dict{String, Array{Float64}}

    altitude_bin_edges::Vector{Float64}
    altitude_bin_means::Vector{Float64}

    energy_bin_edges::Vector{Float64}
    energy_bin_means::Vector{Float64}

    pitch_angle_bin_edges::Vector{Float64}
    pitch_angle_bin_means::Vector{Float64}
end

# =====================================
# Helpful functions (frontend)
# =====================================
import Base.print
function print(beaminfo::BeamInfo)
    beaminfo.prefix == "" ? prefix_to_print = "No run prefix" : prefix_to_print = "Run prefix: `$(beaminfo.prefix)`"
    print("""
    $(prefix_to_print)
        $(beaminfo.n_particles) $(beaminfo.particle)s:
        Energy: $(beaminfo.energy) keV
        Pitch Angle: $(beaminfo.pitch_angle)º
        Latitude: $(beaminfo.latitude)º
        Injection: $(beaminfo.injection_alt) km
        Backscatter: $(beaminfo.backscatter_alt) km

    Beam files located at $(beaminfo.dir)/$(beaminfo.base_filename)*.csv
    """
    )
end

function print(beam::Beam)
    print(beam.info)
end

# =====================================
# Main getter functions (frontend)
# =====================================
function get_available_beams(dir_to_search)
    regex_any = "(.*)"
    regex_float = "([0-9]+[\\.]?[0-9]*)" # With optional decimal point
    regex_int = "([0-9]+)"
    search_string = 
        "$(regex_any)_input_" *
        "inject$(regex_float)km_" *
        "lat$(regex_float)deg_" *
        "$(regex_float)keV_" *
        "$(regex_float)deg_" * 
        "$(regex_int)particles"
    regex_search = Regex(search_string)

    filenames = basename.(glob("*.csv", dir_to_search))
    matches = match.(regex_search, filenames)
    base_filenames = unique([el.match for el in matches])

    beams = Vector{BeamInfo}()
    for base_filename in base_filenames
        captures = match(regex_search, base_filename).captures
        length(captures) == 6 || error("Beam filename `$(filename)` in unexpected format. 6 matches expected, $(length(this_files_matches)) matches found")
        
        # Parse prefix, if there is one
        prefix_exists = count(Regex("$(regex_any)_"), captures[1]) > 0
        if prefix_exists
            prefix = match(Regex("$(regex_any)_"), captures[1]).captures[1]
            particle = match(Regex("$(prefix)_$(regex_any)"), captures[1]).captures[1]
        else
            prefix = ""
            particle = captures[1]
        end

        # Parse other captures
        injection_altitude = parse(Float64, captures[2])
        latitude = parse(Float64, captures[3])
        energy = parse(Float64, captures[4])
        pitch_angle = parse(Float64, captures[5])
        n_particles = parse(Int64, captures[6])

        # Find backscatter altitude
        backscatter_search = [match(Regex("$(base_filename)_backscatter_$(regex_float)km.csv"), filename) for filename in filenames]
        backscatter_altitude = backscatter_search[backscatter_search .≠ nothing][1].captures[1]
        backscatter_altitude = parse(Float64, backscatter_altitude)

        # Construct BeamInfo object
        push!(beams, BeamInfo(
            prefix, 
            particle, 
            injection_altitude, 
            backscatter_altitude, 
            latitude, 
            energy,
            pitch_angle, 
            n_particles, 
            dir_to_search,
            base_filename
        ))
    end
    return beams
end

function load_beam(beaminfo::BeamInfo)
    # ==============================
    # Get backscatter
    # ==============================
    backscatter = get_backscatter(beaminfo)

    # ==============================
    # Get energy deposition
    # ==============================
    energy_deposition_data = get_energy_deposition(beaminfo)

    # ==============================
    # Get spectra
    # ==============================
    prebaked_path = "$(beaminfo.dir)/preprocessed_spectra/$(beaminfo.base_filename).npz"
    if isfile(prebaked_path)
        spectra = npzread(prebaked_path)
    else
        prebake_beam(beaminfo)
        return load_beam(beaminfo)
    end

    # Assign bin edge variables
    energy_deposition_data["altitude_bin_edges"] == spectra["altitude_bin_edges"] || error("Altitude bin edge mismatch between spectra and energy deposition!")
    altitude_bin_edges = spectra["altitude_bin_edges"]
    energy_bin_edges = spectra["energy_bin_edges"]
    pitch_angle_bin_edges = spectra["pitch_angle_bin_edges"]

    # ==============================
    # Construct and return
    # ==============================
    return Beam(
        beaminfo,

        energy_deposition_data["energy_deposition"],
        energy_deposition_data["electron_production"],
        backscatter,
        spectra,

        altitude_bin_edges,
        edges_to_means(altitude_bin_edges),

        energy_bin_edges,
        edges_to_means(energy_bin_edges),

        pitch_angle_bin_edges,
        edges_to_means(pitch_angle_bin_edges)
    )
end

# =====================================
# Prebake functions
# =====================================
function prebake_directory(dir_to_prebake; overwrite = true)
    if overwrite && isdir("$(dir_to_prebake)/preprocessed_spectra")
        rm("$(dir_to_prebake)/preprocessed_spectra", recursive = true)
    end

    beams = get_available_beams(dir_to_prebake)
    println("Prebaking $(length(beams)) beams in $(dir_to_prebake)...")
    print_progress_bar(0)

    stdout_lock = ReentrantLock()
    completed_beams = 0
    Threads.@threads for i in eachindex(beams)
        overwrite ? prebake_beam(beams[i]) : load_beam(beams[i])

        # Print status message
        lock(stdout_lock) do 
            print_progress_bar(completed_beams/length(beams))
            completed_beams += 1
        end
    end
    return
end

function prebake_beam(beaminfo::BeamInfo)
    # ==============================
    # Get spectra
    # ==============================
    altitude_bin_edges = nothing
    energy_bin_edges = nothing
    pitch_angle_bin_edges = nothing

    species_to_get = ["e-", "gamma", "proton", "alpha", "neutron"]
    spectra = Dict{String, Array{Float64}}()
    for species in species_to_get
        species_data = read_spectrum(beaminfo, species)

        if isnothing(altitude_bin_edges)
            spectra["altitude_bin_edges"] = species_data["altitude_bin_edges"]
            spectra["energy_bin_edges"] = species_data["energy_bin_edges"]
            spectra["pitch_angle_bin_edges"] = species_data["pitch_angle_bin_edges"]
        else
            spectra["altitude_bin_edges"] == species_data["altitude_bin_edges"]       || error("Mismatch in axis labels across spectra files!")
            spectra["energy_bin_edges"] == species_data["energy_bin_edges"]           || error("Mismatch in axis labels across spectra files!")
            spectra["pitch_angle_bin_edges"] == species_data["pitch_angle_bin_edges"] || error("Mismatch in axis labels across spectra files!")
        end

        spectra[species] = species_data["counts"] ./ beaminfo.n_particles
    end

    # ==============================
    # Write prebaked spectra
    # ==============================
    prebaked_dir = "$(beaminfo.dir)/preprocessed_spectra"
    if isdir(prebaked_dir) == false
        mkdir(prebaked_dir)
    end
    npzwrite("$(prebaked_dir)/$(beaminfo.base_filename).npz", spectra)
    return
end

# =====================================
# Spectrum reader functions
# =====================================
function read_spectrum(beaminfo::BeamInfo, species::String)
    path = "$(beaminfo.dir)/$(beaminfo.base_filename)_spectra_$(species).csv"
    file = open(path, "r")
    contents = read(file, String)
    lines = split(contents, "\n", keepempty = false)
    lines = lines[endswith.(lines, ";")] # Discard headers

    altitude_bin_edges    = parse_csv(lines[1])
    energy_bin_edges      = parse_csv(lines[2])
    pitch_angle_bin_edges = parse_csv(lines[3])
    counts                = parse_csv(lines[4])

    return Dict(
        "altitude_bin_edges" => altitude_bin_edges, 
        "energy_bin_edges" => energy_bin_edges, 
        "pitch_angle_bin_edges" => pitch_angle_bin_edges, 
        "counts" => counts
    )
end

function parse_csv(contents; delimiter = ";")
    # Detect number of dimensions in csv
    n_dimensions = 1
    while occursin(repeat(delimiter, n_dimensions+1), contents)
        n_dimensions += 1
    end

    # Extract data
    nested_result = recursive_split(contents, n_dimensions)
    dims = get_dims_from_nested_vectors(nested_result)
    result = fill("0.000000", Tuple(dims))
    fill_value(result, nested_result, dims, n_dimensions)
    return parse.(Float64, result)

end

function recursive_split(contents, n_dimensions)
    result = split(contents, repeat(";", n_dimensions), keepempty = false)
    if n_dimensions == 1; return result; end

    [recursive_split(el, n_dimensions - 1) for el in result]
end

function get_dims_from_nested_vectors(v; current_dims = [])
    if typeof(v) == SubString{String}; return current_dims; end # Used to be length(v) == 1
    get_dims_from_nested_vectors(v[1], current_dims = [current_dims..., length(v)])
end

function fill_value(result, nested_result, dims, n_dimensions; current_position = (0 .* dims) .+ 1)
    if length(dims) == 0
        if nested_result == "0.000000"; return; end
        result[current_position...] = nested_result
        return
    end

    for idx in 1:dims[1]
        current_position[n_dimensions - length(dims) + 1] = idx
        fill_value(result, nested_result[idx], dims[2:end], n_dimensions, current_position = current_position)
    end
end

# =====================================
# Ionization reader functions
# =====================================
# TODO get ionization function

function get_energy_deposition(beaminfo::BeamInfo)
    path = "$(beaminfo.dir)/$(beaminfo.base_filename)_energydeposition_ioncount.csv"
    data = readdlm(path, ',', skipstart = 1)

    # Parse altitude labels
    altitude = data[:, 1]
    altitude_bin_edges = zeros(length(altitude)+1)

    for i in eachindex(altitude)
        bin_edges = split(altitude[i], "-")
        altitude_bin_min = parse(Float64, split(bin_edges[1], "km", keepempty = false)[1])
        altitude_bin_max = parse(Float64, split(bin_edges[2], "km", keepempty = false)[1])

        if i == 1
            altitude_bin_edges[begin] = altitude_bin_min
            altitude_bin_edges[2] = altitude_bin_max
            continue
        end

        altitude_bin_edges[i] == altitude_bin_min  || error("Mismatch in altitude bin edges in energy deposition file!")
        
        altitude_bin_edges[i+1] = altitude_bin_max
    end

    # Parse data
    energy_deposition = Float64.(data[:,2])
    electron_production = Float64.(data[:,4])

    # Return
    return Dict(
        "altitude_bin_edges" => altitude_bin_edges,
        "energy_deposition" => energy_deposition ./ beaminfo.n_particles,
        "electron_production" => electron_production ./ beaminfo.n_particles
    )
end

function get_W_from_atmosphere_profile(atmosphere_path)
    # Energy per ion pair data by species. Units: eV/pair
    W = Dict{String, Float64}(
        "Ar" => 27,
        "He" => 32.5,
        "H2" => 38.0,
        "N2" => 35.8,
        "O2" => 32.2,
        "O"  => 36 #?????? unsure on value
    )

    # Atomic masses. Units: kg
    m = Dict{String, Float64}(
        "O"  => 2.6567e-26,
        "N"  => 2.3259e-26,
        "He" => 6.646477e-27,
        "Ar" => 6.6335e-26,
        "H"  => 1.674e-27,
    )
    # Diatomics
    m["O2"] = 2 * m["O"]
    m["N2"] = 2 * m["N"]
    m["H2"] = 2 * m["H"]

    data = readdlm(atmosphere_path, ',')
    header = data[1, :]
    composition_data = data[2:end, :]
    altitude = composition_data[:,1]

    constituents = ["O", "N2", "O2", "He", "Ar", "H", "N", "H2"]
    density_data = Dict{String, Vector{Float64}}(
        "total" => zeros(length(altitude))
    )
    for species in constituents
        species_column = findfirst(header .== "$species (kg/m3)")
        density_data[species] = Float64.(composition_data[:, species_column]) ./ m[species] # Divide by atomic mass to get number density
            # Units: number m⁻³
        density_data["total"] .+= density_data[species]
    end

    W_net = zeros(length(altitude))
    for species in constituents
        if species ∉ keys(W); continue; end
        W_net .+= W[species] .* density_data[species] ./ density_data["total"]
    end

    return LinearInterpolator(altitude, W_net)
        # Units: eV pair⁻¹
end

# =====================================
# Backscatter reader functions
# =====================================
function get_backscatter(beaminfo::BeamInfo)
    backscatter_alt = @sprintf "%.2f" beaminfo.backscatter_alt
    path = "$(beaminfo.dir)/$(beaminfo.base_filename)_backscatter_$(backscatter_alt)km.csv"
    data = readdlm(path, ',', skipstart = 1)

    result = Dict{String, Any}()
    column_names = ["particle_names", "weights", "energies_keV", "pitch_angles_deg", "x_unit_momenta", "y_unit_momenta", "z_unit_momenta", "x_positions_m", "y_positions_m", "z_positions_m"]
    for i in eachindex(column_names)
        i == 1 ? to_write = String.(data[:,i]) : to_write = Float64.(data[:,i])
        result[column_names[i]] = to_write
    end

    backscattered_species = unique(result["particle_names"])
    for species in backscattered_species
        result["$(species)_mask"] = result["particle_names"] .== species
    end

    return result
end


results_dir = "$(dirname(TOP_LEVEL))/results/2026-03-18--11.46"

#prebake_directory(results_dir)

beaminfos = get_available_beams(results_dir)
beams = load_beam.(beaminfos)



for beam in beams
    omnidirectional = dropdims(sum(beam.spectra["e-"], dims = 3), dims = 3)
    heatmap(log10.(beam.energy_bin_means), beam.altitude_bin_edges, log10.(omnidirectional),
        title = "Omnidirectional",
        xlabel = "Log10 Energy, keV",
        ylabel = "Altitude, km",
        bg=:black,
        colorbar_title = "Log10 counts",
        clims = (-log10.(beam.info.n_particles)-2, 0),
        xticks = (-2:8),
        ylims = (0, 500)
    )
    display(plot!())
end
