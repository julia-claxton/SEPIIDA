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
include("$(@__DIR__)/General_Functions.jl")

function parse_csv(contents; delimiter = ";")
    # Detect number of dimensions in csv
    n_dimensions = 1
    while occursin(repeat(delimiter, n_dimensions+1), contents)
        n_dimensions += 1
    end

    # Extract data
    nested_result = recursive_split(contents, n_dimensions)
    dims = get_dims_from_nested_vectors(nested_result)
    result = fill("", Tuple(dims))
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
        result[CartesianIndex(Tuple(current_position))] = nested_result
        return
    end

    for idx in 1:dims[1]
        current_position[n_dimensions - length(dims) + 1] = idx
        fill_value(result, nested_result[idx], dims[2:end], n_dimensions, current_position = current_position)
    end
end


prefix = ""

path = glob("$(prefix)*spectra_e-*", "$(dirname(TOP_LEVEL))/build/results")[1]
file = open(path, "r")
contents = read(file, String)
lines = split(contents, "\n", keepempty = false)
lines = lines[endswith.(lines, ";")] # Discard headers

altitude = parse_csv(lines[1])
energy   = parse_csv(lines[2])
pa       = parse_csv(lines[3])
counts   = parse_csv(lines[4])

omnidirectional = dropdims(sum(counts, dims = 3), dims = 3)
heatmap(log10.(energy[begin:end-1]), altitude, log10.(omnidirectional),
    title = "Omnidirectional",
    xlabel = "Log10 Energy, keV",
    ylabel = "Altitude, km",
    bg=:black,
    colorbar_title = "Log10 counts",
    clims = (-2, 3),
    ylims = (0, 1000)
)
display(plot!())

function my_histogram(data, edges)
    weights = exact_1d_histogram(data, edges)
    plot(edges, [weights..., weights[end]],
        linetype = :steppost
    )
    return plot!()
end


path = glob("$(prefix)*backscatter*", "$(dirname(TOP_LEVEL))/build/results")[1]
data = readdlm(path, ',', skipstart = 1)

particle_name = data[:,1]
mask = particle_name .== "gamma"
e = data[mask,3]
pa = data[mask,4]
weight = data[mask,2]

e_bin_edges = 10.0 .^ (1:.1:3)
my_histogram(e, e_bin_edges)
display(plot!(
    xscale = :log10,
    xlabel = "Energy (keV)",
    ylabel = "Gamma counts"
))


pa_bin_edges = 0:5:180
Ω = 2π .* [cosd(pa_bin_edges[i]) - cosd(pa_bin_edges[i+1]) for i in 1:length(pa_bin_edges)-1]
pa_weights = exact_1d_histogram(pa, pa_bin_edges) ./ Ω
plot(pa_bin_edges, [pa_weights..., pa_weights[end]],
    xlabel = "Pitch Angle (deg)",
    ylabel = "Gamma couns (#/str)",
    xticks = 0:30:180,
    linetype = :steppost
)
display(plot!())


pa_bin_edges_rad = pa_bin_edges .* (π/180)

plot(pa_bin_edges_rad, [pa_weights..., pa_weights[end]],
    linetype = :steppost,
    projection = :polar
)