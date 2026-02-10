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
  fBeamPitchAngle(40.),
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
 
  // Set starting location
  r->xPos = 0; 
  r->yPos = 0;
  r->zPos = (fInitialParticleAlt - 500.0)*km; // Subtraction due to coordinate axis location in middle of world volume

  /* Keeping this block because it's useful. TODO Consider adding an isotropic switch for the user.
  // Isotropic downgoing incidence - choose point on downgoing unit hemisphere for momentum
  std::vector<G4double> momentum = randDowngoingDirection();
  G4double vx0 = momentum.at(0);
  G4double vy0 = momentum.at(1);
  G4double vz0 = momentum.at(2);
  */

  G4double vx0 = 0;
  G4double vy0 = 0;
  G4double vz0 = -1; // Field-aligned

  // Get B vector to rotate velocity to be isotropic downgoing in pitch angle rather than in worldspace
  G4double spacetimePoint[4] = {r->xPos, r->xPos, r->zPos, 0};
  G4double emComponents[6];

  G4FieldManager* fieldManager = G4TransportationManager::GetTransportationManager()->GetFieldManager();
  fieldManager->GetDetectorField()->GetFieldValue(spacetimePoint, emComponents);
  G4double B[3] = {emComponents[0], emComponents[1], emComponents[2]};
  G4double normB = std::sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);
  
  // Get necessary tilt angle
  G4double tilt_angle_rad = std::acos(std::abs(B[2]/normB));

  // Rotate coordinates into world frame
  r->xDir = vx0;
  r->yDir = (std::cos(tilt_angle_rad) * vy0) - (std::sin(tilt_angle_rad) * vz0);
  r->zDir = (std::sin(tilt_angle_rad) * vy0) + (std::cos(tilt_angle_rad) * vz0);

  // Verify that pitch angle generation is correct
  double normMomentum = std::sqrt(pow(r->xDir, 2) + pow(r->yDir, 2) + pow(r->zDir, 2));
  double dotProd = (r->xDir * B[0]) + (r->yDir * B[1]) + (r->zDir * B[2]);
  double generatedPitchAngle_deg = std::acos(dotProd / (normMomentum * normB)) * 180/3.14159265358979;

  if( generatedPitchAngle_deg > 90){
    G4cout << "** ERROR: Primary generated with incorrect pitch angle. You should never see this. Please email julia.claxton@colorado.edu with this error and the conditions that produced it." << G4endl;
    throw;
  }
  
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
