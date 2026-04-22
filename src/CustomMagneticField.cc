

#include "CustomMagneticField.hh"
#include "CustomMagneticFieldMessenger.hh"
#include "ANSIColors.h"
#include <numeric>
#include <functional>
#include <iostream>
#include <unistd.h>
#include "GlobalFunctions.h"

/* The following class calculates the magnetic field strength 
 * and direction according to a model of user's choice. This class inherits
 * from G4MagneticField. GetFieldValue() is a Geant virtual method that 
 * is called to obtain magnetic and electric field values for particle 
 * propagation purposes.
 */

CustomMagneticField::CustomMagneticField()
: G4MagneticField(),
  fMagneticFieldMessenger(0),
  fLAT_degrees(-999.0),  // Units: deg
  fFieldModel("throw an error"), 
  fRe(6371e3),           // Units: m
  fu0(1.257e-6),         // Units: N * A^-2
  fpi(3.14159265358979)
{
  fMagneticFieldMessenger = new CustomMagneticFieldMessenger(this);
}

CustomMagneticField::~CustomMagneticField(){}

void CustomMagneticField::GetFieldValue(const G4double Point[4],G4double *Bfield) const
{
  // Point is a spacetime 4-vector: Point[0..3] = (x, y, z, t)
  
  // Guard
  if(fLAT_degrees == -999.0){
    G4cout << ANSI_RED <<
      __FILE__ << ": " << __FUNCTION__ << "\n" <<
      "ERROR: Magnetic latitude not set. You should not see this." <<
    ANSI_NOCOLOR << G4endl;
    throw;
  }

  // Guard
  std::vector<G4String> availableModels = {
    "igrf2025",
    "jrm33"
  };
  if(std::find(availableModels.begin(), availableModels.end(), fFieldModel) == availableModels.end()){
    G4cout << "\n" <<
      ANSI_RED <<
      __FILE__ << ": " << __FUNCTION__ << "\n" <<
      "ERROR: Magnetic field mode \"" << fFieldModel << "\" not recognized" <<
      ANSI_NOCOLOR <<
    G4endl;
    throw;
  }

  // Set planetary radius based on B-field model
  G4double Rplanet;
  G4double Re = 6371200.0; // Earth radius, m
  G4double Rj = 71492000.0;	// Jupiter equatorial radius. Units: m

  if(fFieldModel == "igrf2025"){
    Rplanet = Re;
  }
  else if(fFieldModel == "jrm33"){
    Rplanet = Rj;
  }

  // Get position vectors for B field calculation

  // Position vector from planet center to world origin. Add 500 due to origin of simulation being 500 km above the bottom of the simulation.
  // Units: m
  // Frame: Planet-centered (Z parallel to spin axis, X/Y in equatorial plane
  G4double LAT_radians = fLAT_degrees * fpi / 180.0;
  G4double r_planetCenter_to_origin[3] = {
    (Rplanet) * std::cos(LAT_radians),
    0,
    (Rplanet + 500e3) * std::sin(LAT_radians),
  }; 

  // Position vector between world origin and particle in world coordinates.
  // Units: m
  // Frame: Planet-centered
  std::vector<G4double> r_origin_to_particle = 
    G4world_to_planetCentric(
      Point[0]/m,
      Point[1]/m,
      Point[2]/m
  );

  // Position vector between planet center and particle.
  // Units: m
  // Frame: Planet-centered
  G4double r[3] = {
    r_planetCenter_to_origin[0] + r_origin_to_particle[0], 
    r_planetCenter_to_origin[1] + r_origin_to_particle[1], 
    r_planetCenter_to_origin[2] + r_origin_to_particle[2]
  };

  // Get B-field values
  // For all models, position r is expected in meters in a cartesian planet-centered frame
  // where Z is parallel to the axis from which latitude is measured (usually the spin axis).
  // Output is in nanotesla.
  G4double Bx_nT_planetCentered, By_nT_planetCentered, Bz_nT_planetCentered;
  if(fFieldModel == "igrf2025"){
    igrf2025Field(
      r[0]/Re, r[1]/Re, r[2]/Re, 
      &Bx_nT_planetCentered, &By_nT_planetCentered, &Bz_nT_planetCentered
    );
  }
  else if(fFieldModel == "jrm33"){
    jrm33Field(
      r[0]/Rj, r[1]/Rj, r[2]/Rj, 
      &Bx_nT_planetCentered, &By_nT_planetCentered, &Bz_nT_planetCentered
    );
  }
  // Don't need a guard here for unrecognized models since we already guarded earlier in the function

  // Rotate back into world coordinates
  // Units: nT
  // Frame: G4 world
  std::vector<G4double> B_worldFrame = 
    planetCentric_to_G4world(
      Bx_nT_planetCentered, 
      By_nT_planetCentered, 
      Bz_nT_planetCentered
  );

  // Assign values
  Bfield[0] = B_worldFrame[0] * 1e-9 * tesla; // Bx
  Bfield[1] = B_worldFrame[1] * 1e-9 * tesla; // By
  Bfield[2] = B_worldFrame[2] * 1e-9 * tesla; // Bz
  Bfield[3] = 0; // Ex
  Bfield[4] = 0; // Ey
  Bfield[5] = 0; // Ez
}

std::vector<G4double> CustomMagneticField::planetCentric_to_G4world(G4double x_siii, G4double y_siii, G4double z_siii) const {
  // This works for any planet-centric frame with Z on the spin axis and X-Y in the equatorial plane

  // We arbitrarily define the world origin as being in the X-Z plane of planet frame for ease of calculation.
  // So at 0º LAT, Z of the G4 world is parallel to X of planetframe, and X of G4 world is anti-parallel to Z of planetframe.
  // Y is constant between the two systems.
  // And at 90º LAT, the coordinate systems are aligned on all three axes.
  // So to rotate between planetframe to G4 world, we need a positive rotation (defined right-handedly) of planetframe about its Y-axis by 90º - LAT.
  // And rotating the frame by θ means rotating vectors by -θ, so we rotate the B-field output by -(90º - LAT) about the Y-axis
  G4double rotation_angle_deg = -(90 - fLAT_degrees); // Units: deg
  G4double rotation_angle_rad = rotation_angle_deg * fpi / 180.0; // Units: rad

  G4double x_g4world = (x_siii * std::cos(rotation_angle_rad)) + (z_siii * std::sin(rotation_angle_rad));
  G4double y_g4world = y_siii;
  G4double z_g4world = (-x_siii * std::sin(rotation_angle_rad)) + (z_siii * std::cos(rotation_angle_rad));

  return std::vector<G4double> {x_g4world, y_g4world, z_g4world};
}

std::vector<G4double> CustomMagneticField::G4world_to_planetCentric(G4double x_g4world, G4double y_g4world, G4double z_g4world) const {
  // This is the reverse process of planetframe -> G4world, so we simply need to rotate vectors by the same amount defined there,
  // just in the opposite direction.
  G4double rotation_angle_deg = 90 - fLAT_degrees; // Units: deg
  G4double rotation_angle_rad = rotation_angle_deg * fpi / 180.0; // Units: rad

  G4double x_siii = (x_g4world * std::cos(rotation_angle_rad)) + (z_g4world * std::sin(rotation_angle_rad));
  G4double y_siii = y_g4world;
  G4double z_siii = (-x_g4world * std::sin(rotation_angle_rad)) + (z_g4world * std::cos(rotation_angle_rad));

  return std::vector<G4double> {x_siii, y_siii, z_siii};
}

G4double CustomMagneticField::vectorMagnitude(std::vector<G4double> v) const {
  return sqrt( (v[0]*v[0]) + (v[1]*v[1]) + (v[2]*v[2]) );
}