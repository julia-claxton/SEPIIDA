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

/*
# To build the SEPIIDA executable: (in zsh)
cd path/to/SEPIIDA-build # You'll need to create this directory if it doesn't exist already
rm -r /path/to/SEPIIDA-build/ *  # Remove space between / and * at the end before running (C++ gets angry if there's a slash-star in a block comment). Removes any existing files or old builds
cmake -DCMAKE_INSTALL_PREFIX="/path/to/geant4-install" -DGeant4_DIR="/path/to/geant4-install/lib" /path/to/SEPIIDA-source
make # Build executable
chmod +x ./RUN_ALL.sh # Make the runall script executable by anyone

# To run the simulation
./SEPIIDA <PARTICLE NAME> <BEAM ENERGY IN KEV> <BEAM PITCH ANGLE IN DEG>
# or
./RUN_ALL.sh
*/

/// \file atmosphericEPP_main.cc
/// \brief Main function to run Geant4 EPP simulation

// Base simulation building classes
#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"
#include "RunAction.hh"

// Multithreaded run manager
#include "G4MTRunManager.hh"
#include "G4Threading.hh"

// Physics lists
#include "G4UImanager.hh"
#include "FTFP_BERT.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4PhysListFactory.hh"
#include "G4StepLimiterPhysics.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4EmParameters.hh"
#include "G4HadronicProcessStore.hh"

#include <chrono>
#include <filesystem>

// For Printing statistic from Transporation process(es)
#include "G4Electron.hh"
#include "G4Transportation.hh"
#include "G4CoupledTransportation.hh"


int main(int argc,char** argv)
{
  // We need 4 arguments provided: a number of particles, particle type, energy, and pitch angle to run.
  // Error out if we don't get those
  if(argc != 4)
  {
    std::cout << "Incorrect number of command line arguments provided. " << argc-1 << " given, 3 required. Format: ./SEPIIDA <number of particles> <particle name> <particle energy>" << std::endl;
    throw;
  }
  
  // Start simulation timer
  auto t_start = std::chrono::high_resolution_clock::now();

  // Choose the Random engine
  G4Random::setTheEngine(new CLHEP::RanecuEngine);
  
  // Set the seeds
  long seeds[2];
  time_t systime = time(NULL);
  
  // Seed built in c-rand engine
  srand(systime);

  // Geant rand engine
  seeds[0] = (long) systime;
  seeds[1] = (long) (systime*G4UniformRand());
  G4Random::setTheSeeds(seeds);

  // Construct the run manager in multithreaded mode
  G4MTRunManager* runManager = new G4MTRunManager;
  runManager->SetNumberOfThreads(G4Threading::G4GetNumberOfCores()); // Use maximum number of possible cores

  // Physics list
  G4PhysListFactory factory;
  G4VModularPhysicsList* physicsList = factory.GetReferencePhysList("QBBC"); // QBBC uses EM v1. Do we need more updated EM model? TODO
  
  //G4VModularPhysicsList* physicsList = new FTFP_BERT;
  physicsList->RegisterPhysics(new G4StepLimiterPhysics());
  physicsList->SetVerboseLevel(0);
  runManager->SetUserInitialization(new DetectorConstruction());
  runManager->SetUserInitialization(physicsList);
  runManager->SetUserInitialization(new ActionInitialization());

  G4double lowLimit = 250.0 * eV;
  G4double highLimit = 100.0 * GeV;
  G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(lowLimit, highLimit);

  // Suppress large verbosity from EM & hadronic processes
  G4EmParameters::Instance()->SetVerbose(-1);
  G4HadronicProcessStore::Instance()->SetVerbose(0);

  // Get the pointer to the user interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  // Set the simulation parameters
  UImanager->ApplyCommand("/control/execute EDIT_THIS_FILE.mac");
  
  // Set input particle number
  G4String nParticles = argv[1];
  UImanager->ApplyCommand("/control/alias NUMBER_OF_PARTICLES " + nParticles);

  // Set particle definition variable, uses Geant4's particle names: https://fismed.ciemat.es/GAMOS/GAMOS_doc/GAMOS.5.1.0/x11519.html
  G4String particle = argv[2];
  UImanager->ApplyCommand("/control/alias BEAM_PARTICLE " + particle); 
  
  // Set particle longname - what the result file will call the input particle. This is just for clarity to 
  // the end user on what each result file represents, as I think "photon" and "electron" are clearer than
  // G4's internal names "gamma" and "e-".
  G4String longname;
  if(particle == "e-")         {longname = "electron";}
  else if(particle == "gamma") {longname = "photon";}
  else                         {longname = particle;}
  UImanager->ApplyCommand("/control/alias BEAM_PARTICLE_LONGNAME " + longname);

  // Set beam energy
  G4String energy = argv[3];
  UImanager->ApplyCommand("/control/alias BEAM_ENERGY_KEV " + energy);

  // Print status block
  std::cout << "=====================================================================" << std::endl;
  std::time_t t = std::time(nullptr);
  std::tm tm = *std::localtime(&t);
  std::cout << "Starting Simulation: " << std::put_time(&tm, "%F %T") << std::endl;
  std::cout << std::endl;

  std::cout << "Multithreading Active" << std::endl;
  std::cout << "    " << G4Threading::G4GetNumberOfCores() << " cores available" << std::endl;
  std::cout << "    " << runManager->GetNumberOfThreads() << " threads active" << std::endl;
  std::cout << std::endl;

  std::cout << "Beam Parameters" << std::endl;
  std::cout << 
    "    Input:       " << nParticles  << " " << longname << "s " << std::endl <<
    "    Energy:      " << energy      << " keV" << std::endl <<
  std::endl;
  std::cout << "=====================================================================" << std::endl;
  std::cout << std::endl;

  // Execute run
  UImanager->ApplyCommand("/control/execute run_beam.mac");

  // End run
  delete runManager;
  auto t_end = std::chrono::high_resolution_clock::now();

  // Report run statistics
  double elapsed_time_ms = std::chrono::duration<double, std::milli>(t_end-t_start).count();
  t = std::time(nullptr);
  tm = *std::localtime(&t);

  std::cout << "=====================================================================" << std::endl;
  std::cout << "Simulation completed in " << elapsed_time_ms/1000.0 << " seconds (" << nParticles << " " << longname << "s @ " << energy << " keV)" << std::endl;
  std::cout << "Simulation Finish: " << std::put_time(&tm, "%F %T") << std::endl;
  std::cout << "=====================================================================" << std::endl << std::endl;

  return 0;
}