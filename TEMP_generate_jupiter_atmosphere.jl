using Statistics, LinearAlgebra          # Core math
using BenchmarkTools, Profile, TickTock  # Debugging
using NPZ, DelimitedFiles                # File interactions
using BasicInterpolators
using Glob
using Plots, Plots.PlotMeasures
    default(
        dpi = 300,
        framestyle = :box,
        tickdirection = :out,
        label = false
    )


""" 
Extremely crude model of the top 1000 km of Jupiter's atmosphere. This is for placeholder
purposes until we can get access to a more detailed model. This just uses an exponential
pressure profile, isothermal assumption, ideal gas law, and a constant mixing ratio of 90/10 H2/He.
"""

molar_mass_h2 = 2.016e-3   # kg mol⁻¹
molar_mass_he = 4.00260e-3 # kg mol⁻¹
R = 8.31446261815324       # J K⁻¹ mol⁻¹
R_h2 = R / molar_mass_h2   # J K⁻¹ kg⁻¹
R_he = R / molar_mass_he   # J K⁻¹ kg⁻¹

igl_density(P, R, T) = P / R * T # Ideal gas law

# Get T-P profile
data = readdlm("$(@__DIR__)/TEMP_jupiter_TP_profile_Knížek_2025.csv", ',', skipstart = 0)

profile_T = data[:,1]
profile_P = data[:,2] .* 100000 # Convert bar to Pa

sortvec = sortperm(profile_P)
profile_T = profile_T[sortvec]
profile_P = profile_P[sortvec]

temperature = LinearInterpolator(profile_P, profile_T, NoBoundaries())


# Get pressure profile
H = 27 # Scale height [km]
P0 = 1e5 # Surface pressure, Pa
z = 0:1:999 # Altitudes to sample [km]

P = P0 .* exp.(-z./H)
T = temperature.(P)

ρ_h2 = 0.9 .* igl_density.(P, R_h2, T)
ρ_he = 0.1 .* igl_density.(P, R_he, T)

file = open("$(@__DIR__)/atmosphere_profile.csv", write = true)
columns = ["Altitude (km)", "O (kg/m3)", "N2 (kg/m3)", "O2 (kg/m3)", "Total (kg/m3)", "Neutral Temp. (K)", "He (kg/m3)", "Ar (kg/m3)", "H (kg/m3)", "N (kg/m3)", "H2 (kg/m3)"]

# Write header
for i in eachindex(columns)
    write(file, columns[i])
    i == length(columns) ? write(file, "\n") : write(file, ",")
end

for altitude_idx = eachindex(z)
    write(file, "$(z[altitude_idx]),$(repeat("0.0,",4))$(T[altitude_idx]),$(ρ_he[altitude_idx]),$(repeat("0.0,",3))$(ρ_h2[altitude_idx])\n")
end

close(file)
println("Done")