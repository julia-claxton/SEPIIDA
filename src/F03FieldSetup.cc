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
/// \file field/field03/src/F03FieldSetup.cc
/// \brief Implementation of the F03FieldSetup class
//
//
// $Id: F03FieldSetup.cc 104351 2017-05-26 07:23:04Z gcosmo $
//
//
//   Field Setup class implementation.
//
//

#include "F03FieldSetup.hh"
#include "F03FieldMessenger.hh"

#include "G4MagneticField.hh"
#include "G4UniformMagField.hh"
#include "CustomMagneticField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4Transportation.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4ChordFinder.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

#include "G4ExplicitEuler.hh"
#include "G4ImplicitEuler.hh"
#include "G4SimpleRunge.hh"
#include "G4SimpleHeum.hh"
#include "G4ClassicalRK4.hh"
#include "G4HelixExplicitEuler.hh"
#include "G4HelixImplicitEuler.hh"
#include "G4HelixSimpleRunge.hh"
#include "G4CashKarpRKF45.hh"
#include "G4RKG3_Stepper.hh"
#include "G4DormandPrince745.hh"
#include "G4ExactHelixStepper.hh"
#include "G4HelixHeum.hh"
#include "G4HelixMixedStepper.hh"
#include "G4NystromRK4.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "GlobalFunctions.h"
#include "ANSIColors.h"

F03FieldSetup::F03FieldSetup()
 : fFieldManager(0),
   fChordFinder(0),
   fEquation(0),
   fStepper(0),
   fStepperType("unset"),
   fFieldMessenger(0),
   cacheRadius(HACKY_CACHE_RADIUS * km)
{
  Guard();

  // Constants
  fStepperType = "G4DormandPrince745"; // Set the stepper here
  fMinStep = 0.01 * km;

  // Field options
  nonCachedMagneticField = new CustomMagneticField(); 
  G4MagneticField* fieldToUse = nonCachedMagneticField;
  CACHED_MAGNETIC_FIELD = false;
  
  if(cacheRadius > 0.0){
    cachedMagneticField = new G4CachedMagneticField(nonCachedMagneticField,  cacheRadius * km);
    fieldToUse = cachedMagneticField;
    CACHED_MAGNETIC_FIELD = true;
  }

  // Field setup stuff
  fFieldManager = G4TransportationManager::GetTransportationManager()->GetFieldManager();
  fFieldMessenger = new F03FieldMessenger(this);
  fEquation = new G4Mag_UsualEqRhs(fieldToUse);

  // 1. Clean up previous state
  delete fChordFinder;
  fChordFinder = nullptr;

  // 2. Create the steppers (Note: this also deletes the previous ones.)
  SetStepper();

  // 3. Create the chord finder(s)
  fChordFinder = new G4ChordFinder(fieldToUse, fMinStep, fStepper);
  fFieldManager->SetChordFinder(fChordFinder);

  // 4. Ensure that the field is updated (in Field manager & equation)
  fFieldManager->SetDetectorField(fieldToUse);
}

F03FieldSetup::~F03FieldSetup()
{
  delete fChordFinder;
  delete fStepper;
  delete fFieldMessenger;
}


void F03FieldSetup::SetStepper()
{
  delete fStepper;
  fStepper = nullptr;
  bool reportStepper = true; // Whether to print the stepper beign used to terminal at start of simulation

  if(fStepperType == "G4ExplicitEuler"){
    fStepper = new G4ExplicitEuler( fEquation );
  }
  else if(fStepperType == "G4ImplicitEuler"){
    fStepper = new G4ImplicitEuler( fEquation );
  }
  else if(fStepperType == "G4SimpleRunge"){
    fStepper = new G4SimpleRunge( fEquation );
  }
  else if(fStepperType == "G4SimpleHeum"){
    fStepper = new G4SimpleHeum( fEquation );
  }
  else if(fStepperType == "G4ClassicalRK4"){
    fStepper = new G4ClassicalRK4( fEquation );
  }
  else if(fStepperType == "G4HelixExplicitEuler"){
    fStepper = new G4HelixExplicitEuler( fEquation );
  }
  else if(fStepperType == "G4HelixImplicitEuler"){
    fStepper = new G4HelixImplicitEuler( fEquation );
  }
  else if(fStepperType == "G4HelixSimpleRunge"){
    fStepper = new G4HelixSimpleRunge( fEquation );
  }
  else if(fStepperType == "G4CashKarpRKF45"){
    fStepper = new G4CashKarpRKF45( fEquation );
  }
  else if(fStepperType == "G4RKG3_Stepper"){
    fStepper = new G4RKG3_Stepper( fEquation );
  }
  else if(fStepperType == "G4DormandPrince745"){
    fStepper = new G4DormandPrince745( fEquation );
  }
  else if(fStepperType == "G4ExactHelixStepper"){
    fStepper = new G4ExactHelixStepper( fEquation );
  }
  else if(fStepperType == "G4HelixHeum"){
    fStepper = new G4HelixHeum( fEquation );
  }
  else if(fStepperType == "G4HelixMixedStepper"){
    fStepper = new G4HelixMixedStepper( fEquation );
  }
  else if (fStepperType == "G4NystromRK4"){
    fStepper = new G4NystromRK4( fEquation );
  }
  else {
    G4cout
      << __FILE__ << ": " << __FUNCTION__ << "\n"
      << "Stepper " << fStepperType << " not recognized."
    << G4endl;
    throw;
  }
  if(reportStepper){G4cout << "Using stepper " << fStepperType << G4endl;}
}

G4FieldManager* F03FieldSetup::GetGlobalFieldManager(){
  return G4TransportationManager::GetTransportationManager()->GetFieldManager();
}

void F03FieldSetup::SetCacheRadius(G4double newCacheRadius){
  println("Setting radius");
  cacheRadius = newCacheRadius;
}

void F03FieldSetup::Guard(){
  if(cacheRadius == -999.0){
    G4cout 
      << ANSI_RED << "\n"
      << __FILE__ << ": " << __FUNCTION__ << "\n"
      << "Cache radius unset. You should never see this."
    << G4endl;
    throw;
  }
}