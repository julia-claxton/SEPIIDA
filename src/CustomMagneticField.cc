

#include "CustomMagneticField.hh"
#include "CustomMagneticFieldMessenger.hh"
#include <numeric>
#include <functional>
#include <iostream>
#include <unistd.h>

/* The following class calculates the magnetic field strength 
 * and direction according to a model of user's choice. This class inherits
 * from G4MagneticField. GetFieldValue() is a Geant virtual method that 
 * is called to obtain magnetic and electric field values for particle 
 * propagation purposes.
 */

CustomMagneticField::CustomMagneticField()
: G4MagneticField(),
  fMagneticFieldMessenger(0),
  fDipoleMoment(6.4e22), // Earth magnetic moment, A * m^2. Value source: https://sciencedemonstrations.fas.harvard.edu/presentations/earths-magnetic-field
  fLAT_degrees(65.77),  // Units: deg
  fFieldModel("none"), 
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
  // Bfield is a pointer to a 6x1 array of E- and B-field components
  if(fFieldModel == "earth_tilted_dipole"){
    earthFieldDipole(Point, Bfield); // For centered dipole at Earth
  }
  else if(fFieldModel == "jrm33"){
    getJrm33Field(Point, Bfield);
  }
  else {
    G4cout << "\n" <<
      "\033[0;31m" <<
      __FILE__ << ": " << __FUNCTION__ << "\n" <<
      "ERROR: Magnetic field mode \"" << fFieldModel << "\" not recognized" <<
      "\033[0m" <<
    G4endl;
    throw;
  }
}

void CustomMagneticField::earthFieldDipole(const G4double Point[4], G4double *Bfield) const {
  // Calculate field components using centered dipole model
  G4double LAT_radians = fLAT_degrees * fpi / 180.0;
  G4double magMoment[3] = {0, -1 * fDipoleMoment * std::cos(LAT_radians), -1 * fDipoleMoment * std::sin(LAT_radians)}; // Magnetic moment of Earth in world coordinates
  G4double r_earthCenter_to_origin[3] = {0, 0, fRe + 500e3}; // Units: m. Add 500 due to origin of simulation being 500 km above sea level.
  G4double r_origin_to_particle[3] = {Point[0]/m, Point[1]/m, Point[2]/m}; // Position vector between world origin and particle in world coordinates. Units: m
  G4double r[3] = {
    r_earthCenter_to_origin[0] + r_origin_to_particle[0],
    r_earthCenter_to_origin[1] + r_origin_to_particle[1],
    r_earthCenter_to_origin[2] + r_origin_to_particle[2]
  }; // Units: m

  // Get dot product of m and r
  G4double dotProd = 0;
  for(int i = 0; i < 3; i++){
    dotProd += magMoment[i] * r[i];
  }

  // Get magnitude of r
  G4double rMag = std::sqrt((r[0]*r[0]) + (r[1]*r[1]) + (r[2]*r[2]));

  // Get each component of field strength
  G4double B[3];
  for(int i = 0; i < 3; i++){
    B[i] = (fu0/(4*fpi)) * ( ((3*dotProd*r[i])/pow(rMag, 5)) - (magMoment[i]/pow(rMag, 3)) );
  }

  // Assign values
  // x = East direction
  // y = North direction
  // z = Up direction, radially out from Earth 
  Bfield[0] = B[0] * tesla; // Bx
  Bfield[1] = B[1] * tesla; // By
  Bfield[2] = B[2] * tesla; // Bz
  Bfield[3] = 0; // Ex
  Bfield[4] = 0; // Ey
  Bfield[5] = 0; // Ez
}

void CustomMagneticField::getJrm33Field(const G4double Point[4],G4double *Bfield) const {
  G4double Rj = 71492000.0;	// Jupiter equatorial radius. Units: m
  G4double LAT_radians = fLAT_degrees * fpi / 180.0;

  // Wilson (2022), Space Sci. Rev. "The coordinate system in use is the right handed System III (1965), which deﬁned Jupiter’s rotation rate as 870.536º per day"
  // https://lasp.colorado.edu/mop/files/2015/02/CoOrd_systems12.pdf
  // SIII: X = 0º latitude, corotating with Jupiter. Z = Jupiter spin axis.

  // Position vector from Jupiter center to world origin. Add 500 due to origin of simulation being 500 km above the bottom of the simulation.
  // Units: m
  // Frame: SIII
  G4double r_planetCenter_to_origin[3] = {
    (Rj + 500.0) * std::cos(LAT_radians),
    0,
    (Rj + 500.0) * std::sin(LAT_radians),
  }; 

  // Position vector between world origin and particle in world coordinates.
  // Units: m
  // Frame: G4 world
  G4double r_origin_to_particle_g4WorldFrame[3] = {
    Point[0]/m,
    Point[1]/m,
    Point[2]/m
  };

  // Position vector between world origin and particle in world coordinates.
  // Units: m
  // Frame: SIII
  std::vector<G4double> r_origin_to_particle = G4world_to_SIII(
    r_origin_to_particle_g4WorldFrame[0],
    r_origin_to_particle_g4WorldFrame[1],
    r_origin_to_particle_g4WorldFrame[2]
  );
  
  // Position vector between Jupiter center and particle.
  // Units: m
  // Frame: SIII
  G4double r[3] = {
    r_planetCenter_to_origin[0] + r_origin_to_particle[0], 
    r_planetCenter_to_origin[1] + r_origin_to_particle[1], 
    r_planetCenter_to_origin[2] + r_origin_to_particle[2]
  };

  // Get field value
  // Units: nT
  // Frame: SIII
  G4double Bx_nT_siii, By_nT_siii, Bz_nT_siii;
  jrm33Field(r[0]/Rj, r[1]/Rj, r[2]/Rj, &Bx_nT_siii, &By_nT_siii, &Bz_nT_siii);

  // Rotate back into world coordinates
  std::vector<G4double> B_worldFrame = SIII_to_G4world(Bx_nT_siii, By_nT_siii, Bz_nT_siii);

  // Assign values
  Bfield[0] = B_worldFrame[0] * 1e-9 * tesla; // Bx
  Bfield[1] = B_worldFrame[1] * 1e-9 * tesla; // By
  Bfield[2] = B_worldFrame[2] * 1e-9 * tesla; // Bz
  Bfield[3] = 0; // Ex
  Bfield[4] = 0; // Ey
  Bfield[5] = 0; // Ez
}

std::vector<G4double> CustomMagneticField::SIII_to_G4world(G4double x_siii, G4double y_siii, G4double z_siii) const {
  // We arbitrarily define the world origin as being in the X-Z plane of SIII for ease of calculation.
  // So at 0º LAT, Z of the G4 world is parallel to X of SIII, and X of G4 world is anti-parallel to Z of SIII.
  // And at 90º LAT, the coordinate systems are aligned on all three axes.
  // So to rotate between SIII to G4 world, we need a positive rotation (defined right-handedly) of SIII about its Y-axis by 90º - LAT.
  // And rotating the frame by θ means rotating vectors by -θ, so we rotate the B-field output by -(90º - LAT) about the Y-axis
  G4double rotation_angle_deg = -(90 - fLAT_degrees); // Units: deg
  G4double rotation_angle_rad = rotation_angle_deg * fpi / 180.0; // Units: rad

  G4double x_g4world = (x_siii * std::cos(rotation_angle_rad)) + (z_siii * std::sin(rotation_angle_rad));
  G4double y_g4world = y_siii;
  G4double z_g4world = (-x_siii * std::sin(rotation_angle_rad)) + (z_siii * std::cos(rotation_angle_rad));

  return std::vector<G4double> {x_g4world, y_g4world, z_g4world};
}

std::vector<G4double> CustomMagneticField::G4world_to_SIII(G4double x_g4world, G4double y_g4world, G4double z_g4world) const {
  // This is the reverse process of SIII -> G4world, so we simply need to rotate vectors by the same amount defined there,
  // just in the opposite direction.
  G4double rotation_angle_deg = 90 - fLAT_degrees; // Units: deg
  G4double rotation_angle_rad = rotation_angle_deg * fpi / 180.0; // Units: rad

  G4double x_siii = (x_g4world * std::cos(rotation_angle_rad)) + (z_g4world * std::sin(rotation_angle_rad));
  G4double y_siii = y_g4world;
  G4double z_siii = (-x_g4world * std::sin(rotation_angle_rad)) + (z_g4world * std::cos(rotation_angle_rad));

  return std::vector<G4double> {x_siii, y_siii, z_siii};
}