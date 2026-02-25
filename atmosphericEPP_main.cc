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

# To see syntax & flags, ./SEPIIDA -help
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
#include <map>

#include "G4Electron.hh"
#include "G4Transportation.hh"
#include "G4CoupledTransportation.hh"

// These function declarations feel like bad practice... Oh well.
extern void printHelpScreen();
extern void println(G4String line);

int main(int argc,char** argv)
{
  // Check for help flag
  if((argc == 2) && (std::strcmp(argv[1], "-help") == 0)){
    printHelpScreen();
    return 0;
  }
  
  // We need at least 4 arguments provided: a number of particles, particle type, energy, and pitch angle to run.
  // Error out if we don't get those.
  if(argc < 5){
    G4cout << 
      "Incorrect number of command line arguments provided. " << argc-1 << " given, at least 4 required. Format: ./SEPIIDA <number of particles> <particle name> <particle energy> <particle pitch angle> <optional arguments>\n" << 
      "Call `/path/to/SEPIIDA -help` for help." <<
    G4endl;
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
  UImanager->ApplyCommand("/control/execute EDIT_THIS_FILE.mac");

  // Verbosity off
  UImanager->ApplyCommand("/control/verbose 0");
  UImanager->ApplyCommand("/run/verbose 0");
  UImanager->ApplyCommand("/event/verbose 0");
  UImanager->ApplyCommand("/tracking/verbose 0");
  UImanager->ApplyCommand("/run/particle/verbose 0");
  UImanager->ApplyCommand("/geometry/navigator/verbose 0");
  UImanager->ApplyCommand("/process/had/verbose 0");

  // Let's go, lesbians!
  UImanager->ApplyCommand("/run/initialize");

  // ==========================================
  // Parse command line arguments
  // ==========================================
  
  // Required args:
  // Set input particle number
  G4String nParticles = argv[1];
  UImanager->ApplyCommand("/control/alias NUMBER_OF_PARTICLES " + nParticles);

  // Set particle definition variable, uses Geant4's particle names: https://fismed.ciemat.es/GAMOS/GAMOS_doc/GAMOS.5.1.0/x11519.html
  G4String particle = argv[2];
  UImanager->ApplyCommand("/control/alias BEAM_PARTICLE " + particle); 
  UImanager->ApplyCommand("/beamParameters/setBeamParticle {BEAM_PARTICLE}");
  
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
  UImanager->ApplyCommand("/beamParameters/setBeamEnergy {BEAM_ENERGY_KEV}");
  
  // Set beam pitch angle
  G4String pitchAngle = argv[4];
  UImanager->ApplyCommand("/control/alias BEAM_PITCH_ANGLE_DEG " + pitchAngle);
  UImanager->ApplyCommand("/beamParameters/setBeamPitchAngle {BEAM_PITCH_ANGLE_DEG}");

  // Set default values for optional arguments
  std::map<G4String, G4String> optionalFlags = {
    {"-magnetic_model",  "earth_tilted_dipole"},   // What magnetic field model to use. Current options: "earth-tilted-dipole", "jrm33"
    {"-lat",                  "67.0"},                  // Magnetic latitude [deg]
    {"-atmosphere_filename", "atmosphere_profile.csv"}, // Filename for atmospheric profile
    {"-brem_splitting",      "100"},                    // Number of times to split bremsstrahlung photons
    {"-altitude_offset",     "0.0"},                    // Amount by which to offset altitude axis labels [km] TODO not implemented
    {"-injection_altitude",  "450.0"},                  // Altitude to inject particles at [km]
  };
  // Add backscatter argument after map is made so we can reference the injection altitude for its default value
  optionalFlags.insert(
    {"-backscatter_altitude", std::to_string(std::stod(optionalFlags["-injection_altitude"]) + 1.0)} // Altitude to track backscattered particles at [km]
  );

  // Get length of longest flag for printing purposes later on
  G4int flagsMaxLength = 0;
  for (auto const& el : optionalFlags){
    if(el.first.length() > flagsMaxLength){flagsMaxLength = el.first.length();}
  }

  // Assign optional arguments
  std::vector<G4String> flagsUserChanged;
  for(int argIdx = 5; argIdx < argc; argIdx += 2){
    G4String flagName = argv[argIdx];
    G4String flagValue = argv[argIdx+1];

    // Error for repeated flags
    if(std::find(flagsUserChanged.begin(), flagsUserChanged.end(), flagName) != flagsUserChanged.end()){
      G4cout <<
        "\033[0;31m" <<
        "ERROR: Flag \"" << flagName << "\" repeated." <<
        "\033[0m" <<
      G4endl;
      throw;
    }

    // Error if an unrecognized flag is provided
    if(optionalFlags.find(flagName) == optionalFlags.end()){
      G4cout <<
        "\033[0;31m" <<
        "ERROR: Flag \"" << flagName << "\" not found. Check spelling or run `path/to/SEPIIDA -help` for a list of available flags." <<
        "\033[0m" <<
      G4endl;
      throw;
    }

    // Do things
    flagsUserChanged.push_back(flagName);
    optionalFlags[flagName] = flagValue;
  }

  // Pass optional arguments to simulation
  // B field mode
  UImanager->ApplyCommand("/fieldParameters/setFieldModel " + optionalFlags["-magnetic_model"]);

  // Latitude
  UImanager->ApplyCommand("/fieldParameters/setLAT " + optionalFlags["-lat"]);
  UImanager->ApplyCommand("/control/alias LAT_DEGREES " + optionalFlags["-lat"]);

  // Atmosphere filename
  UImanager->ApplyCommand("/dataCollection/setAtmosFileName " + optionalFlags["-atmosphere_filename"]);

  // Brem splitting
  UImanager->ApplyCommand("/process/em/setSecBiasing eBrem world " + optionalFlags["-brem_splitting"] + " 100 MeV"); // Syntax: /process/em/setSecBiasing processName Region factor energyUpperLimit energyUnit

  // Altitude offset
  UImanager->ApplyCommand("/dataCollection/setAltitudeOffset " + optionalFlags["-altitude_offset"]);

  // Injection altitude
  UImanager->ApplyCommand("/beamParameters/setParticleStartingAltitude " + optionalFlags["-injection_altitude"]);
  UImanager->ApplyCommand("/control/alias INJECTION_ALTITUDE_KM " + optionalFlags["-injection_altitude"]);

  // Backscatter altitude
  UImanager->ApplyCommand("/dataCollection/setCollectionAltitude " + optionalFlags["-backscatter_altitude"]);


  // Print status block
  G4cout << "=====================================================================" << G4endl;
  std::time_t t = std::time(nullptr);
  std::tm tm = *std::localtime(&t);
  G4cout << "Starting Simulation: " << std::put_time(&tm, "%F %T") << G4endl;
  G4cout << G4endl;

  G4cout << "Multithreading Active" << G4endl;
  G4cout << "    " << G4Threading::G4GetNumberOfCores() << " cores available" << G4endl;
  G4cout << "    " << runManager->GetNumberOfThreads() << " threads active" << G4endl;
  G4cout << G4endl;

  G4cout << "Beam Parameters" << G4endl;
  G4cout << 
    "    Input:       " << nParticles  << " " << longname << "s " << G4endl <<
    "    Energy:      " << energy      << " keV" << G4endl <<
    "    Pitch Angle: " << pitchAngle  << " deg" << G4endl <<
  G4endl;

  G4cout << "Simulation Configuration:" << G4endl;
  for (auto const& el : optionalFlags){
    G4cout << "    " << el.first << std::string(flagsMaxLength - el.first.length() + 1, ' ') << "| " << el.second << G4endl;
  }
  G4cout << "=====================================================================" << G4endl;
  G4cout << G4endl;

  // Execute run
  UImanager->ApplyCommand("/control/execute run_beam.mac");

  // End run
  delete runManager;
  auto t_end = std::chrono::high_resolution_clock::now();

  // Report run statistics
  double elapsed_time_ms = std::chrono::duration<double, std::milli>(t_end-t_start).count();
  t = std::time(nullptr);
  tm = *std::localtime(&t);

  G4cout << "=====================================================================" << G4endl;
  G4cout << "Simulation completed in " << elapsed_time_ms/1000.0 << " seconds (" << nParticles << " " << longname << "s @ " << energy << " keV, " << pitchAngle << "º)" << G4endl;
  G4cout << "Simulation Finish: " << std::put_time(&tm, "%F %T") << G4endl;
  G4cout << "=====================================================================" << G4endl << G4endl;

  return 0;
}

void printHelpScreen(){
  println("");
  println("-----------------------------------------------------------------------------------------------------");
  println("Help Page:");
  println("  Simulation of Energetic Particle Incidence, Ionization, and Dynamics in an Atmosphere (SEPIIDA)");
  println("  By Julia Claxton (she/they), based on work by Grant Berland");
  println("  Questions: julia.claxton@colorado.edu");
  println("");
  println("To run executable:");
  println("path/to/SEPIIDA <number of particles> <input particle species> <input particle energy (keV)> <input particle pitch angle (deg)>");
  println("");
  println("Optional flags:");
  println("  -magnetic_model");
  println("      Magnetic field configuration");
  println("      Default: earth_tilted_dipole");
  println("      Options: earth_tilted_dipole, jrm33");
  println("");
  println("  -lat");
  println("      Latitude [deg]");
  println("      Default: 67.0");
  println("");
  println("  -atmosphere_filename");
  println("      Atmosphere profile filename");
  println("      Default: atmosphere_profile.csv");
  println("");
  println("  -brem_splitting");
  println("      Bremsstrahlung photon splitting");
  println("      Default: 100");
  println("");
  println("  -altitude_offset");
  println("      Amount to offset world altitude labels by [km]");
  println("      Default: 0.0");
  println("");
  println("  -injection_altitude");
  println("      Particle injection altitude [km]");
  println("      Default: 450.0");
  println("");
  println("  -backscatter_altitude");
  println("      Backscatter recording altitude [km]");
  println("      Default: injection_altitude + 1.0");
  println("");
  println("  -help");
  println("      Print help screen if value is 1");
  println("      Default: 0");
  println("");
  println("Press enter to continue...");
  println("-----------------------------------------------------------------------------------------------------");
  getchar();
}

void println(G4String line){
  G4cout << line << G4endl;
}