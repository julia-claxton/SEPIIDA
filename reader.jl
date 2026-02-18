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


path = glob("*spectra_e-*", "/Users/luna/Research/geant4/SEPIIDA/build/results/mlat_45deg_input_500km")[1]
file = open(path, "r")
contents = read(file, String)
lines = split(contents, "\n", keepempty = false)
lines = lines[endswith.(lines, ";")] # Discard headers

altitude = parse_csv(lines[1])
energy   = parse_csv(lines[2])
pa       = parse_csv(lines[3])
counts   = parse_csv(lines[4])

#=
path = glob("*energydeposition*", "/Users/luna/Research/geant4/SEPIIDA/results/mlat_45deg_input_450km/")[1]
data = readdlm(path, ',', skipstart = 1)
edep = data[:,2]
ion = data[:,3]


plot(log10.(ion./1e5), eachindex(ion./1e5),
ylims = (0, 200)
)

plot(log10.(edep./1e5), eachindex(edep./1e5),
ylims = (0, 200)
)

error()


for i in 1:size(counts)[3]
    heatmap(log10.(energy[begin:end-1]), altitude, log10.(counts[:, :, i]),
        title = "$(pa[i])º - $(pa[i+1])º",
        xlabel = "Log10 Energy, keV",
        ylabel = "Altitude, km",
        bg=:black,
        colorbar_title = "Log10 counts",
        clims = (-2, 3),
        ylims = (0, 500)
    )
    display(plot!())
end
=#

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
