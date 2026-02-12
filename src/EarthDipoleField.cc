

#include "EarthDipoleField.hh"
#include "EarthDipoleFieldMessenger.hh"
#include <numeric>
#include <functional>
#include <iostream>
#include <unistd.h>

#include "/Users/luna/Downloads/test/libjupitermag/include/jupitermag.h"



/* The following class calculates the Earth's magnetic field strength 
 * and direction according to a tilted dipole model. This class inherits
 * from G4MagneticField. GetFieldValue() is a Geant virtual method that 
 * is called to obtain magnetic and electric field values for particle 
 * propagation purposes.
 *
 */

EarthDipoleField::EarthDipoleField()
: G4MagneticField(),
  fDipoleFieldMessenger(0),
  fDipoleMoment(6.4e22), // Earth magnetic moment, A * m^2. Value source: https://sciencedemonstrations.fas.harvard.edu/presentations/earths-magnetic-field
  fMLAT_degrees(65.77),  // Units: deg
  fRe(6371e3),           // Units: m
  fu0(1.257e-6),         // Units: N * A^-2
  fpi(3.14159265358979)
{
  fDipoleFieldMessenger = new EarthDipoleFieldMessenger(this);



  InternalModel jupiterMagModel;

  jupiterMagModel = InternalModel();
  jupiterMagModel.SetModel("jrm33");





}


EarthDipoleField::~EarthDipoleField(){}

void EarthDipoleField::GetFieldValue(const G4double Point[4],G4double *Bfield) const
{
  // Point is a spacetime 4-vector: Point[0..3] = (x, y, z, t)
  // Bfield is a pointer to a 6x1 array of E- and B-field components
  // Calculate field components using centered dipole model
  G4double MLAT_radians = fMLAT_degrees * fpi / 180.0;
  G4double magMoment[3] = {0, -1 * fDipoleMoment * std::cos(MLAT_radians), -1 * fDipoleMoment * std::sin(MLAT_radians)}; // Magnetic moment of Earth in world coordinates
  G4double r_earthCenter_to_origin[3] = {0, 0, fRe + 500e3}; // Units: m. Add 500 due to origin of simulation being 500 km above sea level.
  G4double r_origin_to_particle[3] = {Point[0]/m, Point[1]/m, Point[2]/m}; // Position vector between world origin and particle in world coordinates. Units: m
  G4double r[3] = {r_earthCenter_to_origin[0]+r_origin_to_particle[0], r_earthCenter_to_origin[1]+r_origin_to_particle[1], r_earthCenter_to_origin[2]+r_origin_to_particle[2]}; // Units: m

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
    B[i] = B[i] * tesla;
  }





  /*
  double x = 10.0, y = 0.0, z = 0.0;
  double Bx, By, Bz;
  jupiterMagModel.Field(x,y,z,&Bx,&By,&Bz);
  G4cout << Bx << ", " << By << ", " << Bz << G4endl;
  */







  // Assign values
  // x = East direction
  // y = North direction
  // z = Up direction, or radially out from Earth 
  Bfield[0] = B[0]; // Bx
  Bfield[1] = B[1]; // By
  Bfield[2] = B[2]; // Bz
  Bfield[3] = 0;    // Ex
  Bfield[4] = 0;    // Ey
  Bfield[5] = 0;    // Ez
}