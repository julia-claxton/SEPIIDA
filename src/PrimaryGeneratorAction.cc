//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: PrimaryGeneratorAction.cc 94307 2015-11-11 13:42:46Z gcosmo $
//
/// \file PrimaryGeneratorAction.cc
/// \brief Creates particles that will be propagated in the simulation

#include "PrimaryGeneratorAction.hh"
#include "PrimaryGeneratorMessenger.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
//#include "G4GeneralParticleSource.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"
#include "G4TransportationManager.hh"
#include "G4FieldManager.hh"
#include "G4MagneticField.hh"
#include <math.h>


PrimaryGeneratorAction::PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
  fParticleGun(0),
  fPrimaryMessenger(0),
  fBeamEnergy(100.),
  fBeamPitchAngle_deg(0.0),
  fInitialParticleAlt(450.0),
  fPI(3.14159265359),
  fRad2Deg(180.0 / 3.14159265359),
  fSourceType("e-")
{
  fPrimaryMessenger = new PrimaryGeneratorMessenger(this);
}


PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
  delete fPrimaryMessenger;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // Select input particle type
  fParticleGun  = new G4ParticleGun();
  G4ParticleDefinition* inputParticle = G4ParticleTable::GetParticleTable()->FindParticle(fSourceType); // Electron = "e-", proton = "proton", photon = "gamma"
  fParticleGun->SetParticleDefinition(inputParticle);
  
  // Set up container for particle properties
  ParticleSample* r = new ParticleSample();

  // Set particle energy
  r->energy = fBeamEnergy * keV;

  // Initial position vector
  G4ThreeVector x0(0.0, 0.0, (fInitialParticleAlt - 500.0)*km); // Subtract 500 from z due to coordinate axis location in middle of world volume

  // Get initial velocity vector. We do this by getting the B field vector, finding an
  // orthogonal vector to it, then rotating by the desired pitch angle about the orthogonal
  // rotation vector.

  // Get B field vector
  G4double spacetimePoint[4] = {x0[0], x0[1], x0[2], 0};

  G4double emComponents[6];
  G4FieldManager* fieldManager = G4TransportationManager::GetTransportationManager()->GetFieldManager();
  fieldManager->GetDetectorField()->GetFieldValue(spacetimePoint, emComponents);
  G4double B[3] = {emComponents[0], emComponents[1], emComponents[2]};
  G4double normB = std::sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);

  // Rotate B vector by the desired input pitch angle to get input velocity vector.
  // We first need to find a vector that is orthogonal to B that we can rotate about.
  // We will do this by crossing B with the X-axis.

  // It sucks that we have to do this for every particle rather than pre-calculating the velocity
  // direction, but GetFieldValue() isn't alive when this script starts, so whatever. Oh well.

  // In the extremely rare case B is perfectly aligned with the x-axis, we'll throw an error. This should never happen
  // in this simulation as it is now, so I'll deal with it in the moment if this happens. Somehow.
  G4ThreeVector unitB(B[0]/normB, B[1]/normB, B[2]/normB);
  G4ThreeVector x(1.0, 0.0, 0.0);

  if(std::abs(1 - unitB.dot(x)) < 1e-10){
    // If this error happens, switch to a different vector than x. y or z would be fine.
    G4cout << "\n" <<
      "\033[0;31m" <<
      __FILE__ << ": " << __FUNCTION__ << "\n" <<
      "ERROR: B is parallel to X-axis at primary generation point. Cannot generate orthogonal vector." 
      "\033[0m" <<
    G4endl;
  }

  // TODO this ^^ is actually very easy to fix by switching to y unit vector if the condition is true

  // Get rotation axis
  G4ThreeVector rotationAxis = unitB.cross(x);
  rotationAxis = rotationAxis / rotationAxis.mag(); // Convert to unit vector

  // Euler's finite rotation formula
  G4double fBeamPitchAngle_rad = fBeamPitchAngle_deg * fPI / 180.0;
  G4ThreeVector v0 = 
    (unitB * std::cos(fBeamPitchAngle_rad))
    + (rotationAxis.cross(unitB) * std::sin(fBeamPitchAngle_rad))
    + (rotationAxis * rotationAxis.dot(unitB) * (1 - std::cos(fBeamPitchAngle_rad)))
  ;
  v0 = v0 / v0.mag();

  // Safety check: Verify that pitch angle generation is correct
  G4double generatedPitchAngle_deg;
  G4double for_acos = unitB.dot(v0) / (unitB.mag() * v0.mag());
  if(std::abs(for_acos - 1.0) < 1e-10){generatedPitchAngle_deg = 0;} // Float precision near 0º and 90º can cause out-of-domain errors resulting in NaNs
  else if(std::abs(for_acos) < 1e-10){generatedPitchAngle_deg = 90;}
  else {generatedPitchAngle_deg = std::acos(for_acos) * 180/fPI;}

  if(std::abs(generatedPitchAngle_deg - fBeamPitchAngle_deg) > 1e-5){
    G4cout << "\n" <<
      "\033[0;31m" <<
      __FILE__ << ": " << __FUNCTION__ << "\n" <<
      "ERROR: Primary generated with incorrect pitch angle.\n" <<
      "You should never see this. Please email julia.claxton@colorado.edu with this error and the conditions that produced it.\n" <<
      "\tDesired pitch angle: " << fBeamPitchAngle_deg << "º\n" <<
      "\tGenerated pitch angle: " << generatedPitchAngle_deg << "º\n" <<
      "\tv0 = (" << v0[0] << ", " << v0[1] << ", " << v0[2] << ")\n" << 
      "\tB = (" << unitB[0] << ", " << unitB[1] << ", " << unitB[2] << ")" << 
      "\033[0m" <<
    G4endl;
    throw;
  }
 
  // Assign position & velocity
  r->xPos = x0[0];
  r->yPos = x0[1];
  r->zPos = x0[2];
  r->xDir = v0[0];
  r->yDir = v0[1];
  r->zDir = v0[2];
  
  // Communicate parameters to particle gun
  fParticleGun->SetParticlePosition(G4ThreeVector(r->xPos, r->yPos, r->zPos)); 
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(r->xDir, r->yDir, r->zDir));
  fParticleGun->SetParticleEnergy(r->energy);
  
  // Geant method to create initial particle with the above properties 
  fParticleGun->GeneratePrimaryVertex(anEvent);

  // Free memory from ParticleSample struct
  delete r;
}

std::vector<G4double> PrimaryGeneratorAction::randDowngoingDirection(){
  // Use rejection sampling to generate point on the downgoing unit hemisphere
  while(true){
    G4double x = (2 * G4UniformRand()) - 1; // [-1, 1]
    G4double y = (2 * G4UniformRand()) - 1; // [-1, 1]
    G4double z = (1 * G4UniformRand()) - 1; // [-1, 0] Only look in the downgoing hemisphere

    G4double radius = std::sqrt( (x*x) + (y*y) + (z*z) );

    if( (radius <= 1) && (radius > 0) ){
      x /= radius;
      y /= radius;
      z /= radius;

      // Safety check
      if( std::abs(1 - std::sqrt((x*x) + (y*y) + (z*z))) > 1e-12 ){
        G4cout << "You should never ever ever ever see this." << G4endl;
        G4cout << "(x, y, z) = (" << x << ", " << y << ", " << z << ")" << G4endl;
        throw;
      }

      // Return breaks the while(true) loop
      return std::vector<G4double> {x, y, z};
    }
  }
}
