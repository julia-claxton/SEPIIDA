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


""" 
Extremely crude model of the top 1000 km of Jupiter's atmosphere. This is for placeholder
purposes until we can get access to a more detailed model. This just uses an exponential
P, T profile, ideal gas law, and a constant mixing ratio of 90/10 H2/He.
"""

molar_mass_h2 = 2.016e-3   # kg mol⁻¹
molar_mass_he = 4.00260e-3 # kg mol⁻¹
R = 8.31446261815324       # J K⁻¹ mol⁻¹
R_h2 = R / molar_mass_h2   # J K⁻¹ kg⁻¹
R_he = R / molar_mass_he   # J K⁻¹ kg⁻¹

igl_density(P, R, T) = R * T / P # Ideal gas law

H = 27 # Scale height [km]
P0 = 1e5 # Surface pressure, Pa
T0 = 165 # Surface temperature, K

z = 0:1:999
P = P0 .* exp.(-z./H)
T = T0 .* exp.(-z./H)

ρ = igl_density.(P, R_h2, T)

ρ_h2 = 0.9ρ
ρ_he = 0.1ρ

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