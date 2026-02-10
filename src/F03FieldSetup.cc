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
#include "EarthDipoleField.hh"
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

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"


F03FieldSetup::F03FieldSetup()
 : fFieldManager(0),
   fChordFinder(0),
   fEquation(0),
   fMagneticField(0),
   fStepper(0),
   fFieldMessenger(0)
{
  fMagneticField = new EarthDipoleField(); 
  fFieldMessenger = new F03FieldMessenger(this);
  fEquation = new G4Mag_UsualEqRhs(fMagneticField);

  // Default values
  fMinStep     = 0.01*km ; 
  fStepperType = 10;
  fFieldManager = GetGlobalFieldManager();

  UpdateField();
}


F03FieldSetup::~F03FieldSetup()
{
  delete fMagneticField;
  delete fChordFinder;
  delete fStepper;
  delete fFieldMessenger;
}


void F03FieldSetup::UpdateField()
{
  // It must be possible to call 'again' - e.g. to choose an alternative stepper
  // has been chosen, or in case other changes have been made.

  // 1. Clean up previous state
  delete fChordFinder;
  fChordFinder= nullptr;

  // 2. Create the steppers ( Note: this also deletes the previous ones. )
  SetStepper();

  // 3. Create the chord finder(s)
  fChordFinder = new G4ChordFinder(fMagneticField, fMinStep, fStepper);
  fFieldManager->SetChordFinder(fChordFinder);

  // 4. Ensure that the field is updated (in Field manager & equation)
  fFieldManager->SetDetectorField(fMagneticField);
}

void F03FieldSetup::SetStepper()
{
  delete fStepper;
  fStepper= nullptr;
  bool reportStepper = false; // Whether to print the stepper beign used to terminal at start of simulation

  switch ( fStepperType )
  {
    case 0:
      fStepper = new G4ExplicitEuler( fEquation );
      if(reportStepper){G4cout<<"G4ExplicitEuler is called"<<G4endl;}
      break;
    case 1:
      fStepper = new G4ImplicitEuler( fEquation );
      if(reportStepper){G4cout<<"G4ImplicitEuler is called"<<G4endl;}
      break;
    case 2:
      fStepper = new G4SimpleRunge( fEquation );
      if(reportStepper){G4cout<<"G4SimpleRunge is called"<<G4endl;}
      break;
    case 3:
      fStepper = new G4SimpleHeum( fEquation );
      if(reportStepper){G4cout<<"G4SimpleHeum is called"<<G4endl;}
      break;
    case 4:
      fStepper = new G4ClassicalRK4( fEquation );
      if(reportStepper){G4cout<<"G4ClassicalRK4 (default) is called"<<G4endl;}
      break;
    case 5:
      fStepper = new G4HelixExplicitEuler( fEquation );
      if(reportStepper){G4cout<<"G4HelixExplicitEuler is called"<<G4endl;}
      break;
    case 6:
      fStepper = new G4HelixImplicitEuler( fEquation );
      if(reportStepper){G4cout<<"G4HelixImplicitEuler is called"<<G4endl;}
      break;
    case 7:
      fStepper = new G4HelixSimpleRunge( fEquation );
      if(reportStepper){G4cout<<"G4HelixSimpleRunge is called"<<G4endl;}
      break;
    case 8:
      fStepper = new G4CashKarpRKF45( fEquation );
      if(reportStepper){G4cout<<"G4CashKarpRKF45 is called"<<G4endl;}
      break;
    case 9:
      fStepper = new G4RKG3_Stepper( fEquation );
      if(reportStepper){G4cout<<"G4RKG3_Stepper is called"<<G4endl;}
      break;
    case 10:
      fStepper = new G4DormandPrince745( fEquation );
      if(reportStepper){G4cout<<"ode45 is called"<<G4endl;}
      break;
    default: fStepper = 0;
  }
}

void F03FieldSetup::SetFieldValue(G4double fieldStrength)
{
  G4ThreeVector fieldSetVec(0.0, 0.0, fieldStrength);
  SetFieldValue( fieldSetVec );
}

void F03FieldSetup::SetFieldValue(G4ThreeVector fieldVector)
{
  if(fMagneticField) delete fMagneticField;

  if(fieldVector != G4ThreeVector(0.,0.,0.))
  {
    fMagneticField = new EarthDipoleField();
  }
  else
  {
    // If the new field's value is Zero, then
    // setting the pointer to zero ensures
    // that it is not used for propagation.
    fMagneticField = 0;
  }
  

  // Either
  //   - UpdateField() to reset all (ChordFinder, Equation);
  // UpdateField();
  //     or simply update the field manager & equation of motion
  //     with pointer to new field
  GetGlobalFieldManager()->SetDetectorField(fMagneticField);
  fEquation->SetFieldObj( fMagneticField );
}

G4FieldManager* F03FieldSetup::GetGlobalFieldManager(){
  return G4TransportationManager::GetTransportationManager()->GetFieldManager();
}

