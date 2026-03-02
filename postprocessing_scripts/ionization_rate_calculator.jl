using Statistics, LinearAlgebra          # Core math
using BenchmarkTools, Profile, TickTock  # Debugging
using NPZ, DelimitedFiles                # File interactions
using Glob
using Plots, Plots.PlotMeasures
    default(
        dpi = 300,
        framestyle = :box,
        tickdirection = :out,
        label = false
    )

const TOP_LEVEL = dirname(@__DIR__)


function get_W_from_atmosphere_profile(atmosphere_path)
    # Energy per ion pair data by species
    W = Dict{String, Float64}(
        "Ar" => 27,
        "He" => 32.5,
        "H2" => 38.0,
        "N2" => 35.8,
        "O2" => 32.2,
        "O"  => 36 #?????? unsure on value
    )

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
        density_data[species] = Float64.(composition_data[:, species_column])
        density_data["total"] .+= density_data[species]
    end

    W_net = zeros(length(altitude))
    for species in constituents
        if species ∉ keys(W); continue; end
        W_net .+= W[species] .* density_data[species] ./ density_data["total"]
    end

    return W_net
        # Units: eV pair⁻¹
end

function get_energy_deposition()
    dir = glob("*", "$(dirname(TOP_LEVEL))/build/results/")[1]
    path = glob("*energydeposition*", dir)[1]

    data = readdlm(path, ',', skipstart = 1)

    # Parse altitude labels
    altitude = data[:, 1]
    altitude_bins_min = zeros(length(altitude))
    altitude_bins_max = zeros(length(altitude))

    for i in eachindex(altitude)
        bin_edges = split(altitude[i], "-")
        altitude_bins_min[i] = parse(Float64, split(bin_edges[1], "km", keepempty = false)[1])
        altitude_bins_max[i] = parse(Float64, split(bin_edges[2], "km", keepempty = false)[1])
    end

    energy_deposition = Float64.(data[:,2])

    return altitude_bins_min, altitude_bins_max, energy_deposition
end


W = get_W_from_atmosphere_profile("$(TOP_LEVEL)/atmosphere_profile.csv")
z_min, z_max, energy_deposition = get_energy_deposition()
mean_altitude = (z_min .+ z_max) ./ 2

ionization = (energy_deposition .* 1e3) ./ W

plot(W, mean_altitude,
    xlabel = "W (eV pair⁻¹)",
    xlims = (37, 47),

    ylabel = "Altitude (km)",
    ylims = (0, 1000)
)
display(plot!())

plot(log10.(ionization), mean_altitude,
    xlabel = "Ionization (TODO units)",

    ylabel = "Altitude (km)",
    ylims = (0, 1000)
)
display(plot!())