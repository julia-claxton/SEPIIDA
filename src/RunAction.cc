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
// $Id: RunAction.cc 99560 2016-09-27 07:03:29Z gcosmo $
//
/// \file RunAction.cc
/// \brief Implementation of the RunAction class

#include "RunAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4AccumulableManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "G4ParticleDefinition.hh"
#include "G4Transportation.hh"
#include "G4CoupledTransportation.hh"
#include "G4Electron.hh"

#include "RunActionMessenger.hh"

#include <fstream>
#include <regex>
#include <filesystem>

RunAction::RunAction():
  G4UserRunAction(),
  fRunActionMessenger(),
  fBaseResultPath() 
{
  // Set killing energies
  fWarningEnergy = 0.01 * keV; // Particles below this energy are killed after 1 step. Value arbitrary 
  fImportantEnergy = 0.1 * keV; // Particles above this energy are killed after fNumberOfTrials if they are looping. Value arbitrary 
  fNumberOfTrials = 1000; // Number of trials before a looping 'important' particle is killed. Value arbitrary

  fRunActionMessenger = new RunActionMessenger(this); 

  // Create energy spectra sample planes
  sampleAltitudes_km = linspace(fMinSampleAltitude_km, fMaxSampleAltitude_km, fNumberOfSamplePlanes);

  // Create energy spectra energy bins
  energyBinEdges_keV = linspace(std::log10(fEnergyMinkeV), std::log10(fEnergyMaxkeV), fNumberOfEnergyBins + 1); 
  for(int i = 0; i < fNumberOfEnergyBins+1; i++){
    energyBinEdges_keV.at(i) = pow(10, energyBinEdges_keV.at(i));
  }

  // Create bin edges
  energyDepositionBinEdges = linspace(fMinSampleAltitude_km, fMaxSampleAltitude_km, fNumberOfSamplePlanes); // This is the same value as sampleAltitudes_km. Copied to a new name for clarity with different types of data recording.
  pitchAngleBinEdges_deg = linspace(fPitchAngleMin, fPitchAngleMax, fNumberOfPitchAngleBins+1);

  // Initialize result histograms
  mainSpectrum.resize(particlesToRecord.size(), std::vector<std::vector<std::vector<G4double>>>(fNumberOfSamplePlanes, std::vector<std::vector<G4double>>(fNumberOfEnergyBins, std::vector<G4double>(fNumberOfPitchAngleBins, 0))));
  totalEnergyDeposition.resize(fNumberOfSamplePlanes-1, 0);
  ionizingEnergyDeposition.resize(fNumberOfSamplePlanes-1, 0);
  ionCounts.resize(fNumberOfSamplePlanes-1, 0);

  std::vector<std::string> fBackscatteredParticleNames;
  std::vector<double> fBackscatteredTrackWeights;
  std::vector<double> fBackscatteredEnergieskeV;
  std::vector<double> fBackscatteredPitchAnglesDeg;
  std::vector<std::array<double,3>> fBackscatterDirections;
  std::vector<std::array<double,3>> fBackscatterPositions;

  // Precalculate things
  energySpectrumHistogramFactor = pow(10, (std::log10(fEnergyMaxkeV) - std::log10(fEnergyMinkeV)) / fNumberOfEnergyBins);
  altitudeSpacing_km = std::abs(sampleAltitudes_km[1] - sampleAltitudes_km[0]);
  pitchAngleBinSize_deg = (fPitchAngleMax - fPitchAngleMin) / fNumberOfPitchAngleBins;
}

RunAction::~RunAction()
{
  delete fRunActionMessenger;
}

void RunAction::BeginOfRunAction(const G4Run*)
{
  int threadID = G4Threading::G4GetThreadId();

  // If we are a worker thread (not the main thread)
  if(threadID != -1){
    printTimestamp();
    G4cout <<"STARTING: Thread " << threadID << G4endl;
    ChangeLooperParameters( G4Electron::Definition() );
    return;
  }

  // If we are here, then we are the main thread
  // First, make sure that the build directory set by the user is correct
  std::filesystem::path baseResultsPath = fBaseResultPath.c_str();
  std::filesystem::path buildDirectory = baseResultsPath.parent_path().parent_path();
  if(std::filesystem::is_directory(buildDirectory) == false)
  {
    G4cout << "\n" <<
      "\033[0;31m" <<
      __FILE__ << ": " << __FUNCTION__ << "\n" <<
      "ERROR: User-specified build directory " << buildDirectory << " does not exist. This path is user-specified in set_simulation_parameters.mac. Check that SEPIIDA_BUILD_DIR in EDIT_THIS_FILE.mac matches your build directory and does not have a slash at the end."
      "\033[0m" <<
    G4endl;
    throw;
  }

  // Change parameters for looping particles
  ChangeLooperParameters( G4Electron::Definition() );
}

void RunAction::ChangeLooperParameters(const G4ParticleDefinition* particleDef)
{
  if(particleDef == nullptr)
    particleDef = G4Electron::Definition();
  auto transportPair= findTransportation(particleDef);
  auto transport = transportPair.first;
  auto coupledTransport = transportPair.second;

  if(transport != nullptr)
  { 
    // Change the values of the looping particle parameters of Transportation 
    if(fWarningEnergy >= 0.0)
      transport->SetThresholdWarningEnergy(  fWarningEnergy ); 
    if(fImportantEnergy >= 0.0)
      transport->SetThresholdImportantEnergy(  fImportantEnergy ); 
    if(fNumberOfTrials > 0)
      transport->SetThresholdTrials( fNumberOfTrials );
  }
  else if(coupledTransport != nullptr)
  { 
    // Change the values for Coupled Transport
    if(fWarningEnergy >= 0.0)
      coupledTransport->SetThresholdWarningEnergy(fWarningEnergy); 
    if(fImportantEnergy >= 0.0)
      coupledTransport->SetThresholdImportantEnergy(fImportantEnergy); 
    if(fNumberOfTrials > 0)
      coupledTransport->SetThresholdTrials(fNumberOfTrials);
  }
}

void RunAction::EndOfRunAction(const G4Run*){
  // Get thread ID to see if we are main thread or not
  int threadID = G4Threading::G4GetThreadId();

  // If we are not the main thread, write results to file and exit
  if(threadID != -1)
  {
    printTimestamp();
    G4cout << "\033[0;36m WRITING: Thread " << threadID << "\033[0m" << G4endl;
    threadWriteSpectra(threadID);
    threadWriteEnergyDeposition(threadID);
    threadWriteIonProduction(threadID);
    threadWriteBackscatter(threadID);
    printTimestamp();
    G4cout << "\033[0;32mFINISHED: Thread " << threadID << "\033[0m" << G4endl;
    return;
  }
  // If we are the main thread, merge datafiles from each thread. Main thread ends after workers are done, so this is the end of the simulation
  // TODO metadata headers
  G4cout << G4endl;
  mergeEnergySpectra();
  mergeEnergyDeposition();
  mergeBackscatter();
}

void RunAction::threadWriteSpectra(int threadID){
  for(int particleIndex = 0; particleIndex < particlesToRecord.size(); particleIndex++){
    std::string spectraPath = 
      fBaseResultPath 
      + "_spectra_" + particlesToRecord[particleIndex]
      + "_thread" + std::to_string(threadID)
      + ".csv"
    ;

    // No header for thread-specific files, go right to writing
    std::string data = "";
    for(int altitudeIndex = 0; altitudeIndex < fNumberOfSamplePlanes; altitudeIndex++){
      for(int energyIndex = 0; energyIndex < fNumberOfEnergyBins; energyIndex++){
        for(int paIndex = 0; paIndex < fNumberOfPitchAngleBins; paIndex++){
          data.append(std::to_string(mainSpectrum[particleIndex][altitudeIndex][energyIndex][paIndex]));
          data.append(";");
        }
        data.append(";");
      }
      data.append(";");
    }

    std::ofstream dataFile;
    dataFile.open(spectraPath, std::ios_base::out); // Open file in write mode to overwrite any previous results
    dataFile << data;
    dataFile.close();
  }
}

void RunAction::threadWriteEnergyDeposition(int threadID){
  std::string depositionPath = 
    fBaseResultPath 
    + "_energy_deposition_thread" + std::to_string(threadID)
    + ".csv"
  ;
  std::ofstream dataFile;
  dataFile.open(depositionPath, std::ios_base::out); // Open file in write mode to overwrite any previous results

  // Write results
  for(int altitudeIndex = 0; altitudeIndex < fNumberOfSamplePlanes-1; altitudeIndex++){
    dataFile << totalEnergyDeposition[altitudeIndex] << "\n";
  }
  dataFile.close();
}

void RunAction::threadWriteIonProduction(int threadID){
  std::string path = 
    fBaseResultPath 
    + "_ion_counts_thread" + std::to_string(threadID)
    + ".csv"
  ;
  std::ofstream dataFile;
  dataFile.open(path, std::ios_base::out); // Open file in write mode to overwrite any previous results

  // Write results
  for(int altitudeIndex = 0; altitudeIndex < fNumberOfSamplePlanes-1; altitudeIndex++){
    dataFile << ionCounts[altitudeIndex] << "\n";
  }
  dataFile.close();
}

void RunAction::threadWriteBackscatter(int threadID){
  std::string backscatterPath = 
    fBaseResultPath 
    + "_backscatter_"
    + std::to_string(fCollectionAltitude) + "km"
    + "_thread" + std::to_string(threadID)
    + ".csv"
  ;
  std::ofstream dataFile;
  dataFile.open(backscatterPath, std::ios_base::out); // Open file in write mode to overwrite any previous results

  // Make sure we don't have missing data for any backscattered particle
  int n = fBackscatteredParticleNames.size();
  if((fBackscatteredEnergieskeV.size() != n) || ((fBackscatteredTrackWeights.size() != n)) || (fBackscatterDirections.size() != n) ||(fBackscatterPositions.size() != n) || (fBackscatteredPitchAnglesDeg.size() != n))
  {
    G4cout << "\n" <<
      "\033[0;31m" <<
      __FILE__ << ": " << __FUNCTION__ << "\n" <<
      "ERROR: Incomplete backscatter data! You should never see this." <<
      "\033[0m" <<
    G4endl;
    throw;
  }

  for(int i = 0; i < n; i++)
  {
    dataFile << 
      fBackscatteredParticleNames.at(i)  << "," <<
      fBackscatteredTrackWeights.at(i)   << "," <<
      fBackscatteredEnergieskeV.at(i)    << "," <<
      fBackscatteredPitchAnglesDeg.at(i) << "," <<

      fBackscatterDirections.at(i)[0] << "," <<
      fBackscatterDirections.at(i)[1] << "," <<
      fBackscatterDirections.at(i)[2] << "," <<

      fBackscatterPositions.at(i)[0] << "," <<
      fBackscatterPositions.at(i)[1] << "," <<
      fBackscatterPositions.at(i)[2] << "\n"
    ;
  }
  dataFile.close();
}

void RunAction::printTimestamp(){
  // Pad with spaces to have consistent print location
  int threadID = G4Threading::G4GetThreadId();
  int nThreads = G4Threading::GetNumberOfRunningWorkerThreads();
  int paddingLength = std::to_string(nThreads-1).length() - std::to_string(threadID).length();
  
  // Get current time to print
  std::time_t t = std::time(nullptr);
  std::tm tm = *std::localtime(&t);

  // Print message
  G4cout << std::string(paddingLength, ' ') << "(" << std::put_time(&tm, "%F %T") <<") ";

}

void RunAction::mergeEnergySpectra(){
  // Result already allocated in beginning of run action (mainSpectrum)
  // TODO make the others use already-allocated vars as well

  G4cout << "Merging energy-pitch angle spectra..." << G4endl;

  // Progress bar variables
  int mergedCount = 0;
  int barLength = 30;

  // Loop over threads and add up histograms
  int nThreads = G4Threading::GetNumberOfRunningWorkerThreads();
  for(int thread = 0; thread < nThreads; thread++){
    for(int particleIndex = 0; particleIndex < particlesToRecord.size(); particleIndex++){
      std::string threadFilepath =
        fBaseResultPath
        + "_spectra_" + particlesToRecord[particleIndex]
        + "_thread" + std::to_string(thread)
        + ".csv"
      ;

      // Get data from this thread and add it to the result
      addThreadSpectraToMainHistogram(threadFilepath, particleIndex);

      // Progress bar
      mergedCount++;
      double fraction = mergedCount / static_cast<double>(nThreads*particlesToRecord.size());
      printProgressBar(fraction, barLength);

      // Remove this thread's results file
      std::remove(threadFilepath.c_str());
    }
  }
  G4cout << G4endl;

  for(int particleIndex = 0; particleIndex < particlesToRecord.size(); particleIndex++){
    // Write summed result to string
    std::string mainFilename =
      fBaseResultPath
      + "_spectra_"
      + particlesToRecord[particleIndex]
      + ".csv"
    ;
    std::string data = "Counts:\n";
    for(int altitudeIndex = 0; altitudeIndex < fNumberOfSamplePlanes; altitudeIndex++){
      for(int energyIndex = 0; energyIndex < fNumberOfEnergyBins; energyIndex++){
        for(int paIndex = 0; paIndex < fNumberOfPitchAngleBins; paIndex++){
          data.append(std::to_string(mainSpectrum[particleIndex][altitudeIndex][energyIndex][paIndex]));
          data.append(";");
        }
        data.append(";");
      }
      data.append(";");
    }

    // Write to file
    std::ofstream dataFile;
    dataFile.open(mainFilename, std::ios_base::out); // Open file in write mode to overwrite any previous results
    writeAxisLabelHeader(dataFile);
    dataFile << data;
    dataFile.close();
  }
}

void RunAction::writeAxisLabelHeader(std::ofstream& file){
  file << "Recording Altitudes (km):" << "\n";
  for(int altitudeIndex = 0; altitudeIndex < fNumberOfSamplePlanes; altitudeIndex++){
    file << sampleAltitudes_km[altitudeIndex] + fAltitudeOffset << ";";
  }
  file << "\n";

  file << "Energy Bin Edges (keV):" << "\n";
      for(int energyIndex = 0; energyIndex < fNumberOfEnergyBins+1; energyIndex++){
    file << energyBinEdges_keV[energyIndex] << ";";
  }
  file << "\n";

  file << "Pitch Angle Bin Edges (deg):" << "\n";
      for(int paIndex = 0; paIndex < fNumberOfPitchAngleBins+1; paIndex++){
    file << pitchAngleBinEdges_deg[paIndex] << ";";
  }
  file << "\n";
}

void RunAction::mergeEnergyDeposition(){
  G4cout << "Merging energy deposition..." << G4endl;

  // Progress bar variables
  int mergedCount = 0;
  int barLength = 30;

  // Loop over threads and add up histograms
  int nThreads = G4Threading::GetNumberOfRunningWorkerThreads();
  for(int thread = 0; thread < nThreads; thread++){
    std::string threadEnergyDepFilepath =
      fBaseResultPath 
      + "_energy_deposition_thread" + std::to_string(thread)
      + ".csv"
    ;
    std::string threadIonCountFilepath =
      fBaseResultPath 
      + "_ion_counts_thread" + std::to_string(thread)
      + ".csv"
    ;

    // Get data from this thread and add it to the result
    std::vector<G4double> eDepThreadData;
    std::vector<G4double> ionCountThreadData;

    eDepThreadData.resize(fNumberOfSamplePlanes-1, 0);
    ionCountThreadData.resize(fNumberOfSamplePlanes-1, 0);

    eDepThreadData = read1Dcsv(threadEnergyDepFilepath, eDepThreadData);
    ionCountThreadData = read1Dcsv(threadIonCountFilepath, ionCountThreadData);

    totalEnergyDeposition = add1DAltitudeVectors(totalEnergyDeposition, eDepThreadData);
    ionCounts = add1DAltitudeVectors(ionCounts, ionCountThreadData);

    // Progress bar
    mergedCount++;
    double fraction = mergedCount / static_cast<double>(nThreads);
    printProgressBar(fraction, barLength);

    // Delete this thread-specific file
    std::remove(threadEnergyDepFilepath.c_str());
    std::remove(threadIonCountFilepath.c_str());
  }
  G4cout << G4endl;

  // Write summed result to file
  std::string mainFilename =
    fBaseResultPath
    + "_energydeposition_ionization.csv"
  ;
  std::ofstream dataFile;
  dataFile.open(mainFilename, std::ios_base::out); // Open file in write mode to overwrite any previous results

  // Write header
  dataFile << "Altitude (km),Energy Deposition (keV),Electron Production (counts)" << "\n";

  // Write rows
  for(int altitudeIndex = 0; altitudeIndex < fNumberOfSamplePlanes-1; altitudeIndex++){
    dataFile 
      << energyDepositionBinEdges[altitudeIndex] + fAltitudeOffset << "km-" 
      << energyDepositionBinEdges[altitudeIndex+1] + fAltitudeOffset << "km,"
      << totalEnergyDeposition[altitudeIndex] << ","
      << ionCounts[altitudeIndex]
      << "\n"
    ;
  }
  dataFile.close();
}

void RunAction::mergeBackscatter(){
  G4cout << "Merging backscatter..." << G4endl;

  // Progress bar variables
  int mergedCount = 0;
  int barLength = 30;

  // Allocate result
  std::vector<std::vector<std::string>> result;
  result.resize(10, std::vector<std::string>(0, ""));

  // Loop over threads
  int nThreads = G4Threading::GetNumberOfRunningWorkerThreads();
  for(int thread = 0; thread < nThreads; thread++){
    std::string threadFilepath =
      fBaseResultPath 
      + "_backscatter_"
      + std::to_string(fCollectionAltitude) + "km"
      + "_thread" + std::to_string(thread)
      + ".csv"
    ;

    // Get data from this thread and append it to the result
    result = read2DcsvBackscatter(threadFilepath, result);

    // Progress bar
    mergedCount++;
    double fraction = mergedCount/static_cast<double>(nThreads);
    printProgressBar(fraction, barLength);

    // Delete this thread-specific file
    std::remove(threadFilepath.c_str());
  }
  G4cout << G4endl;

  // Write summed result to file
  // First, format backscatter collection altitude first for brevity in filename
  char buffer[100];
  snprintf(buffer, 100, "%.2f", fCollectionAltitude);
  std::string collectionAltitudeString = buffer;

  std::string mainFilename =
    fBaseResultPath 
    + "_backscatter_"
    + collectionAltitudeString + "km"
    + ".csv"
  ;
  std::ofstream dataFile;
  dataFile.open(mainFilename, std::ios_base::out); // Open file in write mode to overwrite any previous results

  // Write header
  dataFile
    << "Particle Name,"
    << "Weight,"
    << "Energy (keV),"
    << "Pitch Angle (deg),"
    << "Unit Momentum x,"
    << "Unit Momentum y,"
    << "Unit Momentum z,"
    << "Position x (m),"
    << "Position y (m),"
    << "Position z (m)"
    << "\n"
  ;

  // Write rows
  for(int i = 0; i < result.at(0).size(); i++){
    for(int j = 0; j < result.size(); j++){
      dataFile << result.at(j).at(i);
      if(j != result.size()-1){dataFile << ",";}
    }
    dataFile << "\n";
  }
  dataFile.close();
}

void RunAction::writeRunMetadata(G4String filepath)
{
  // run date
  // lat /& dip angle?
  // injection altitude
  // msis run info
  // bs only, collection altitude
  // brem splitting

  //getenv(

}

std::pair<G4Transportation*, G4CoupledTransportation*> RunAction::findTransportation(const G4ParticleDefinition* particleDef, bool reportError)
{
  const auto *partPM = particleDef->GetProcessManager();
    
  G4VProcess* partTransport = partPM->GetProcess("Transportation");
  auto transport= dynamic_cast<G4Transportation*>(partTransport);

  partTransport = partPM->GetProcess("CoupledTransportation");
  auto coupledTransport = dynamic_cast<G4CoupledTransportation*>(partTransport);

  if(reportError && !transport && !coupledTransport)
  {
    G4cerr << "Unable to find Transportation process for particle type "
           << particleDef->GetParticleName()
           << "  ( PDG code = " << particleDef->GetPDGEncoding() << " ) "
    << G4endl;
  }
  
  return std::make_pair( transport, coupledTransport );
}

void RunAction::addThreadSpectraToMainHistogram(std::string path, int particleIndex){
  // Open file
  std::ifstream file;
  file.open(path, std::ifstream::in);

  // Read data
  int altitudeIndex = 0;
  int energyIndex = 0;
  int paIndex = 0;
  for(std::string word; std::getline(file, word, ';');){
    if(word == ""){continue;}

    // Add to histogram
    mainSpectrum[particleIndex][altitudeIndex][energyIndex][paIndex] += std::stod(word);

    // Increment iterators
    paIndex++;
    if(paIndex == (fNumberOfPitchAngleBins)){
      paIndex = 0;
      energyIndex++;
    }
    if(energyIndex == (fNumberOfEnergyBins)){
      energyIndex = 0;
      altitudeIndex++;
    }
    if(altitudeIndex == (fNumberOfSamplePlanes)){
      altitudeIndex = 0;
    }
  }
}

std::vector<std::vector<std::string>> RunAction::read2DcsvBackscatter(std::string path, std::vector<std::vector<std::string>> result){
  // Open file
  std::ifstream file;
  file.open(path, std::ifstream::in);

  // Parse lines
  int dim1Index = 0;
  int dim2Index = 0;

  std::string line;
  std::string token;
  while( std::getline(file, line) ){
    std::istringstream word(line);
    while ( std::getline(word, token, ',') ){
      result.at(dim2Index).push_back(token);
      dim2Index++;
    }
    dim1Index++;
    dim2Index = 0;
  }
  file.close();
  return result;
}

std::vector<G4double> RunAction::read1Dcsv(std::string path, std::vector<G4double> result){
  // Open file
  std::ifstream file;
  file.open(path, std::ifstream::in);

  // Parse lines
  int dim1Index = 0;

  std::string line;
  std::string token;
  while( std::getline(file, line) ){
    std::istringstream word(line);
    while ( std::getline(word, token) ){
      result[dim1Index] = std::stod(token);
    }
    dim1Index++;
  }
  file.close();
  return result;
}

std::vector<std::vector<G4double>> RunAction::addEnergySpectraVectors(std::vector<std::vector<G4double>> v1, std::vector<std::vector<G4double>> v2){
  std::vector<std::vector<G4double>> result;
  result.resize(fNumberOfSamplePlanes, std::vector<G4double>(fNumberOfEnergyBins, 0));

  for(int altitudeIndex = 0; altitudeIndex < fNumberOfSamplePlanes; altitudeIndex++){
    for(int energyIndex = 0; energyIndex < fNumberOfEnergyBins; energyIndex++){
      result[altitudeIndex][energyIndex] = v1[altitudeIndex][energyIndex] + v2[altitudeIndex][energyIndex];
    }
  }
  return result;
}

std::vector<G4double> RunAction::add1DAltitudeVectors(std::vector<G4double> v1, std::vector<G4double> v2){
  std::vector<G4double> result;
  result.resize(fNumberOfSamplePlanes-1, 0);

  for(int altitudeIndex = 0; altitudeIndex < fNumberOfSamplePlanes-1; altitudeIndex++){
    result[altitudeIndex] = v1[altitudeIndex] + v2[altitudeIndex];
  }
  return result;
}

std::vector<G4double> RunAction::linspace(G4double start, G4double stop, int n){
  std::vector<G4double> result;
  G4double stepSize = (stop - start) / (n-1);

  for(int i = 0; i < n; i++){
    result.push_back(start + (i * stepSize));
  }
  return result;
}

void RunAction::printProgressBar(double fraction, int barLength){
  int filled = std::floor(barLength * fraction);
  G4cout
    << "["
    << std::string(filled, '=')
    << std::string(barLength - filled, ' ')
    << "] "
    << 100.0 * std::round(fraction * 1000.0) / 1000 << "%"
    << "\r"
    << std::flush
  ;
  if(fraction == 1.0){
    G4cout << G4endl;
  }
}