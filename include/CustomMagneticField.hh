#ifndef CUSTOMMAGNETICFIELD_h
#define CUSTOMMAGNETICFIELD_h 1

#include "G4MagneticField.hh"
#include "G4SystemOfUnits.hh"
#include "GlobalFunctions.h"
#include "globals.hh"

#include "jupitermag.h" // Jupiter magnetic field model. Source: https://github.com/mattkjames7/libjupitermag
  // Note: To compile on Mac M1, I needed to inline FluxCan and FluxDip, and comment out definition of M_PI in the header for this library

class CustomMagneticFieldMessenger;

class CustomMagneticField : public G4MagneticField
{
public:
  CustomMagneticField();
  virtual ~CustomMagneticField() override;

  virtual void GetFieldValue(const G4double Point[4], G4double *Bfield) const override;

  std::vector<G4double> planetCentric_to_G4world(G4double x_siii, G4double y_siii, G4double z_siii) const;
  std::vector<G4double> G4world_to_planetCentric(G4double x_g4world, G4double y_g4world, G4double z_g4world) const;

  // Messenger methods
  void SetLAT(G4double newLat_degrees);
  void SetFieldModel(G4String newFieldModell);

  G4double GetLAT();
  G4String GetFieldModel();

  // Misc
  G4double vectorMagnitude(std::vector<G4double> v) const;
  
  mutable G4double lat_degrees = -999.0;
  mutable G4String fieldModel = "throw an error";

  private:
  CustomMagneticFieldMessenger* fMagneticFieldMessenger;
  G4double fRe;
  G4double fu0;
  G4double fpi;
  void Guard() const;
};

#endif
