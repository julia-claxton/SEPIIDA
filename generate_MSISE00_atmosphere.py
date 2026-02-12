
'''
# This script calls MSISe-00 at a given time, latitude, and longitude 
# over a given altitude range. Although this library can generate 
# n-dimensional arrays over all these abscissa, it will not write to
# the file correctly unless a single time, lat., and long. are chosen.
# The script overwrites the profile every time it's called. 
#
# Units of atmospheric species are meters^-3, temperature in Kelvin, total
# mass density in kilogram/meters^3, and altitude in kilometers. MASS 
# DENSITIES are written to csv file in kilograms/meter^3.
#
# To call:
#
#   python generate_atmosphere_profile.py arg1
# 
# if arg1 is anything, script will plot data in addition to writing file
# if arg1 is nothing, script will only write file
#
# GB
'''

import msise00
from datetime import datetime
import numpy as np
import sys

################################
##Input parameters to MSISe-00##
################################

altitudeStepSize = 1   # km
lowerAlt         = 0    # km
upperAlt         = 999  # km

# Datetime for MSIS model run
year  = 2018
month = 1
day   = 1
hour  = 3
time = datetime(year, month, day, hour);

# Geographic latitude , longitude
g_lat_lon = [65. , -148.]

################################
################################
################################

# Atomic masses in kg
mOxygen   = 2.6567e-26;
mNitrogen = 2.3259e-26;
mHelium   = 6.646477e-27;
mArgon    = 6.6335e-26;
mHydrogen = 1.674e-27;
mOxygen2  = mOxygen * 2.;
mNitrogen2= mNitrogen * 2.;


# Construct altitude array
altitudeArr = np.linspace(lowerAlt, 
    upperAlt, 
    int((upperAlt+altitudeStepSize)/altitudeStepSize)
);

# Call MSIS model
atmos = msise00.run(
    time = time, 
    altkm = altitudeArr, 
    glat = g_lat_lon[0], 
    glon = g_lat_lon[1]
)

# If program gets any arguments it'll plot the profiles
if len(sys.argv) > 1: 
    import matplotlib.pyplot as plt
    fig = plt.figure();
    ax1 = fig.add_subplot(131)
    for var in atmos.species:
        if var != 'Total':
            atmos[var].plot(y='alt_km', label=var, ax=ax1)
    plt.xscale('log'); plt.legend(fontsize=8);
    plt.xlabel('Species number density [m$^{-3}$]')
    plt.ylabel('Altitude [km]');
    plt.title('Number Density Profile')
    plt.grid();

    ax2 = fig.add_subplot(132)
    atmos['Tn'].plot(y='alt_km', ax=ax2);
    plt.xlabel('Neutral Temperature [K]')
    plt.ylabel('');
    plt.title('Neutral Temperature Profile')
    plt.grid(); 
    
    ax3 = fig.add_subplot(133)
    atmos['Total'].plot(y='alt_km', ax=ax3);
    plt.xscale('log');
    plt.xlabel(r'Total Mass Density [kg/m$^3$]')
    plt.ylabel('');
    plt.title('Density Profile')
    plt.grid(); 
   
    plt.suptitle('MSISe-00 Results: %s , %.1f$^\circ$N %.1f$^\circ$W' % (time, g_lat_lon[0], g_lat_lon[1]))
    plt.show();     

# Write order defined in src/DetectorConstruction.cc
write_order = ['alt_km','O','N2','O2','Total','Tn','He','Ar','H','N']
header = ["Altitude (km)", "O (kg/m3)", "N2 (kg/m3)","O2 (kg/m3)","Total (kg/m3)","Neutral Temp. (K)","He (kg/m3)","Ar (kg/m3)","H (kg/m3)","N (kg/m3)","H2 (kg/m3)"]
mass_mult   = [1, mOxygen, mNitrogen2, mOxygen2, 1, 1, mHelium, mArgon, mHydrogen, mNitrogen];

# Write to atmosphere file
with open('atmosphere_profile.csv', 'w') as f:
    for i in range(len(header)):
        f.write(header[i])
        f.write(',') if i != len(header)-1 else f.write('\n')
    
    for alt in altitudeArr:
        # Query all variables at an altitude slice
        # Outputs are in number per cubic meter
        line = atmos.sel(
            time = time,
            alt_km = alt,
            lat = g_lat_lon[0],
            lon = g_lat_lon[1]
        )

        for index, item in enumerate(write_order):
            f.write(str(line[item].data * mass_mult[index]))
            f.write(',') if item != write_order[-1] else f.write(',0.0\n') # The 0 is for the H2 column, which this MSIS library doesn't output
    # file closes when scope is left

