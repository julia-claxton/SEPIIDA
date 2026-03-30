using Statistics, LinearAlgebra          # Core math
using BenchmarkTools, Profile, TickTock  # Debugging
using NPZ, DelimitedFiles                # File interactions
using Glob
using Printf
using BasicInterpolators
using Base.Threads
using Plots
    default(
        dpi = 300,
        framestyle = :box,
        tickdirection = :out,
        label = false,
        colormap = :bone_1,
        linecolor = :black,
        linewidth = 1.2
    )        

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
    log_path::String
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
function load_dir(dir;
    prefix = nothing,
    particle = nothing,
    energy = nothing,
    pitch_angle = nothing,
    sort_by = "energy"
    )

    return load_beam(get_beamlist(dir,
        prefix = prefix,
        particle = particle,
        energy = energy,
        pitch_angle = pitch_angle,
        sort_by = sort_by
    ))
end

function get_beamlist(dir_to_search;
    prefix = nothing,
    particle = nothing,
    energy = nothing,
    pitch_angle = nothing,
    sort_by = "energy"
    )
    regex_any = "(.*)"
    regex_float = "([\\-]?[0-9]+[\\.]?[0-9]*)" # With optional sign and decimal point
    regex_int = "([\\-]?[0-9]+)" # With optional sign
    
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
            parsed_prefix = match(Regex("$(regex_any)_"), captures[1]).captures[1]
            parsed_particle = match(Regex("$(parsed_prefix)_$(regex_any)"), captures[1]).captures[1]
        else
            parsed_prefix = ""
            parsed_particle = captures[1]
        end

        # Parse other captures
        injection_altitude = parse(Float64, captures[2])
        latitude = parse(Float64, captures[3])
        parsed_energy = parse(Float64, captures[4])
        parsed_pitch_angle = parse(Float64, captures[5])
        n_particles = parse(Int64, captures[6])

        # Find backscatter altitude
        backscatter_filename = glob("$(base_filename)_backscatter_*", dir_to_search)[1]
        backscatter_match = match(Regex("$(base_filename)_backscatter_$(regex_float)km.csv"), backscatter_filename)
        backscatter_altitude = parse(Float64, backscatter_match.captures[1])

        # Get log path
        parsed_prefix == "" ? log_prefix = "" : log_prefix = "$(parsed_prefix)_"
        log_path = "$(dir_to_search)/SEPIIDA_$(captures[1])_$(captures[4])keV_$(captures[5])deg_$(captures[6])particles.log"
        if isfile(log_path) == false; log_path = ""; end

        # Pass over beams that don't match user-sepcified particle/prefix, if specified
        if (prefix ≠ nothing) && (parsed_prefix ≠ prefix); continue; end
        if (particle ≠ nothing) && (parsed_particle ≠ particle); continue; end
        if (energy ≠ nothing) && (parsed_energy ≠ energy); continue; end
        if (pitch_angle ≠ nothing) && (parsed_pitch_angle ≠ pitch_angle); continue; end

        # Construct BeamInfo object
        push!(beams, BeamInfo(
            parsed_prefix, 
            parsed_particle, 
            injection_altitude, 
            backscatter_altitude, 
            latitude, 
            parsed_energy,
            parsed_pitch_angle, 
            n_particles, 
            dir_to_search,
            base_filename,
            log_path
        ))
    end

    # Sort result
    if sort_by == "energy"
        vec_to_sort_by = energy_list(beams)
    elseif sort_by == "pitch angle"
        vec_to_sort_by = pitch_angle_list(beams)
    else
        error("Sorting option \"$(sort_by)\" not recognized")
    end
    sortvec = sortperm(vec_to_sort_by)
    return beams[sortvec]
end

function load_beam(beaminfo::BeamInfo;
    load_backscatter = true
    )
    # ==============================
    # Get backscatter
    # ==============================
    backscatter = get_backscatter(beaminfo, return_empty = !load_backscatter)

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

function load_beam(beaminfo::Vector{BeamInfo},
    load_backscatter = true
)
    return [load_beam(el, load_backscatter = load_backscatter) for el in beaminfo]
end

function energy_list(beamlist::Vector{BeamInfo})
    return energy_list.(beamlist)
end

function energy_list(beams::Vector{Beam})
    return energy_list.(beams)
end

function energy_list(beam::Beam)
    return energy_list(beam.info)
end

function energy_list(beaminfo::BeamInfo)
    return beaminfo.energy
end

function pitch_angle_list(beamlist::Vector{BeamInfo})
    return pitch_angle_list.(beamlist)
end

function pitch_angle_list(beams::Vector{Beam})
    return pitch_angle_list.(beams)
end

function pitch_angle_list(beam::Beam)
    return pitch_angle_list(beam.info)
end

function pitch_angle_list(beaminfo::BeamInfo)
    return beaminfo.pitch_angle
end



# =====================================
# Visualization (frontend)
# =====================================
function quicklook(beam::Beam)
    primary_species = replace(beam.info.particle, "electron" => "e-")
    primary_species == "gamma" ? secondary_species = "e-" : secondary_species = "gamma"

    plot_energy_spectrum(beam, 
        species = primary_species, 
        show_title = false, 
        show_plot = false
    )
    p1 = plot!(title = primary_species)

    plot_energy_spectrum(beam, 
        species = secondary_species, 
        show_title = false, 
        show_plot = false
    )
    p2 = plot!(title = secondary_species)

    p3 = plot_pitch_angle_spectrum(beam, 
        species = primary_species, 
        show_title = false, 
        show_plot = false
    )
    p4 = plot_pitch_angle_spectrum(beam, 
        species = secondary_species, 
        show_title = false, 
        show_plot = false
    )
    p5 = plot_energy_deposition(beam,
        show_title = false, 
        show_plot = false
    )

    layout = @layout [
    [grid(2,2)] b{0.25w}
    ]

    plot(p1, p2, p3, p4, p5,
        suptitle = "$(beam.info.prefix) $(beam.info.energy) keV $(beam.info.pitch_angle)º",
        layout = layout,
        size = (2, 1.1) .* 600
    )
    display(plot!())
    return
end

function plot_energy_spectrum(beam::Beam; species = "e-", show_title = true, show_plot = true)
    to_plot = dropdims(sum(beam.spectra[species], dims = 3), dims = 3)
    
    show_title ? title = "$(beam.info.prefix) $(beam.info.energy) keV $(beam.info.pitch_angle)º" : title = ""
    
    xmin = beam.energy_bin_edges[begin]
    xmax = beam.energy_bin_edges[end]

    ymin = beam.altitude_bin_edges[begin]
    ymax = beam.altitude_bin_edges[end]

    cmin = log10(minimum(to_plot))
    cmax = max(log10(maximum(to_plot)), 0)
    
    heatmap(beam.energy_bin_means, beam.altitude_bin_edges, log10.(to_plot),
        title = title,
        
        xlabel = "Energy (keV)",
        xlims = (xmin, xmax),
        xticks = 10.0 .^ (log10(xmin):1:log10(xmax)),
        xscale = :log10,
        
        ylabel = "Altitude (km)",
        ylims = (ymin, ymax),
        yticks = ymin:50:ymax,

        colorbar_title = "Log10 $(species)/Input Particle",
        clims = (cmin, cmax),
        
        bg_color_inside = :black
    )
    box_aspect!(1)
    if show_plot; display(plot!()); end
    return plot!()
end

function plot_pitch_angle_spectrum(beam::Beam; species = "e-", show_title = true, show_plot = true)
    to_plot = dropdims(sum(beam.spectra[species], dims = 2), dims = 2)
    
    show_title ? title = "$(beam.info.prefix) $(beam.info.energy) keV $(beam.info.pitch_angle)º" : title = ""

    xmin = beam.pitch_angle_bin_edges[begin]
    xmax = beam.pitch_angle_bin_edges[end]

    ymin = beam.altitude_bin_edges[begin]
    ymax = beam.altitude_bin_edges[end]

    cmin = log10(minimum(to_plot))
    cmax = max(log10(maximum(to_plot)), 0)
    
    heatmap(beam.pitch_angle_bin_means, beam.altitude_bin_edges, log10.(to_plot),
        title = title,
        
        xlabel = "Pitch Angle (deg)",
        xlims = (xmin, xmax),
        xticks = xmin:30:xmax,
        
        ylabel = "Altitude (km)",
        ylims = (ymin, ymax),
        yticks = ymin:50:ymax,

        colorbar_title = "Log10 $(species)/Input Particle",
        clims = (cmin, cmax),
        
        bg_color_inside = :black
    )
    box_aspect!(1)
    if show_plot; display(plot!()); end
    return plot!()
end

function quicklook(beaminfo::BeamInfo)
    return quicklook(load_beam(beaminfo, load_backscatter = false))
end

function plot_energy_spectrum(beaminfo::BeamInfo)
    return plot_energy_spectrum(load_beam(beaminfo, load_backscatter = false), species = species, show_title = show_title, show_plot = show_plot)
end

function plot_pitch_angle_spectrum(beaminfo::BeamInfo; species = "e-", show_title = true, show_plot = true)
    return plot_pitch_angle_spectrum(load_beam(beaminfo, load_backscatter = false), species = species, show_title = show_title, show_plot = show_plot)
end


function plot_energy_deposition(beam::Beam; show_title = true, show_plot = true)
    show_title ? title = "$(beam.info.prefix) $(beam.info.energy) keV $(beam.info.pitch_angle)º" : title = ""
    
    plot_mask = (beam.energy_deposition .≠ 0) .&& isfinite.(beam.energy_deposition)
    xmin = minimum(beam.energy_deposition[plot_mask])
    xmax = maximum(beam.energy_deposition[plot_mask])

    ymin = beam.altitude_bin_edges[begin]
    ymax = beam.altitude_bin_edges[end]

    plot(beam.energy_deposition[plot_mask], beam.altitude_bin_means[plot_mask],
        title = title,
        
        xlabel = "Energy Deposition (keV/input particle)",
        xlims = (xmin, 2.5*xmax),
        xticks = 10.0 .^ (ceil(log10(xmin)):2:floor(log10(xmax))),
        xminorticks = 1,
        xminorgrid = true,
        xscale = :log10,
        
        ylabel = "Altitude (km)",
        ylims = (ymin, ymax),
        yticks = ymin:50:ymax,
    )
    box_aspect!(2)
    if show_plot; display(plot!()); end
    return plot!()
end

# =====================================
# Prebake functions
# =====================================
function prebake_directory(dir_to_prebake; overwrite = true)
    if overwrite && isdir("$(dir_to_prebake)/preprocessed_spectra")
        rm("$(dir_to_prebake)/preprocessed_spectra", recursive = true)
    end

    beams = get_beamlist(dir_to_prebake)
    println("Prebaking $(length(beams)) beams in $(dir_to_prebake)...")
    print_progress_bar(0)

    stdout_lock = ReentrantLock()
    completed_beams = 0
    Threads.@threads for i in eachindex(beams)
        overwrite ? prebake_beam(beams[i]) : load_beam(beams[i])

        # Print status message
        lock(stdout_lock) do 
            completed_beams += 1
            print_progress_bar(completed_beams/length(beams))
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

    # Detect what species were recorded
    spectra_files = glob("$(beaminfo.base_filename)_spectra_*.csv", beaminfo.dir)
    species_to_get = String.([match(r"spectra_(.*).csv", el).captures[1] for el in spectra_files])

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
function get_backscatter(beaminfo::BeamInfo; return_empty = false)
    backscatter_alt = @sprintf "%.2f" beaminfo.backscatter_alt
    path = "$(beaminfo.dir)/$(beaminfo.base_filename)_backscatter_$(backscatter_alt)km.csv"
    column_names = ["particle_names", "weights", "energies_keV", "pitch_angles_deg", "x_unit_momenta", "y_unit_momenta", "z_unit_momenta", "x_positions_m", "y_positions_m", "z_positions_m"]
    
    # Exit without reading if needed
    if return_empty == true
        result = Dict{String, Any}()
        [result[colname] = [] for colname in column_names]
        return result
    end

    # Read data
    file_contents = readdlm(path, ',')

    # If no backscatter data, return empty dict
    if size(file_contents)[1] == 1
        result = Dict{String, Any}()
        [result[colname] = [] for colname in column_names]
        return result
    end

    data = file_contents[2:end, :]
    result = Dict{String, Any}()
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