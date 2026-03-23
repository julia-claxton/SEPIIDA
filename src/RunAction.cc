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

#include "ANSIColors.h"

#include <fstream>
#include <regex>
#include <filesystem>

RunAction::RunAction():
  G4UserRunAction(),
  fRunActionMessenger(),
  fBaseResultPath(),
  fCollectionAltitude(-999.0)
{
  // Set killing energies
  fWarningEnergy = 10 * keV; // Looping particles above this energy generate a warning when killed.
  fImportantEnergy = 1 * keV; // Looping particles above this energy are given fNumberOfTrials to become unlooping before being killed. Below this, no trials are given. 
  fNumberOfTrials = 1000; // Number of trials before a looping particles above fImportantEnergy is killed.

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
  // Guard
  if(fCollectionAltitude == -999.0){
    G4cout 
      << ANSI_RED << "\n"
      << __FILE__ << ": " << __FUNCTION__ << "\n"
      << "Backscatter altitude unset. You should never see this."
    << G4endl;
    throw;
  }

  // If we are a worker thread (not the main thread)
  int threadID = G4Threading::G4GetThreadId();
  if(threadID != -1){
    printTimestamp();
    G4cout <<"STARTING: Thread " << threadID << G4endl;
    simStartTimestamp = std::time(nullptr);
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
    // Get simulation wallclock time
    simEndTimestamp = std::time(nullptr);
    std::time_t simTimeSeconds = simEndTimestamp - simStartTimestamp;

    // Write thread data
    printTimestamp();
    G4cout << ANSI_CYAN << " WRITING: Thread " << threadID << ANSI_NOCOLOR << " (Simulation in " << simTimeSeconds << " seconds)" << G4endl;
    threadWriteSpectra(threadID);
    threadWriteEnergyDepositionAndIonCount(threadID);
    threadWriteBackscatter(threadID);
    
    // Get writing wallclock time
    std::time_t writeEndTimestamp = std::time(nullptr);
    std::time_t writeTimeSeconds = writeEndTimestamp - simEndTimestamp;

    // Print finish message
    printTimestamp();
    G4cout << ANSI_GREEN << "FINISHED: Thread " << threadID << ANSI_NOCOLOR << " (Write in " << writeTimeSeconds << " seconds)" << G4endl;
    return;
  }
  // If we are the main thread, merge datafiles from each thread. Main thread ends after workers are done, so this is the end of the simulation
  // TODO metadata headers
  G4cout << G4endl;
  mergeEnergySpectra();
  mergeEnergyDepositionAndIonCount();
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
          G4String countsAsString = std::to_string(mainSpectrum[particleIndex][altitudeIndex][energyIndex][paIndex]);

          // Remove trailing zeroes after decimal point for filesize considerations
          std::string::size_type trimAfter = countsAsString.find_last_not_of('0');
          G4String toWrite = countsAsString.substr(0, trimAfter+1);
          if(toWrite == "0."){toWrite = "0";}
          
          data.append(toWrite);
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

void RunAction::threadWriteEnergyDepositionAndIonCount(int threadID){
  std::string depositionPath = 
    fBaseResultPath 
    + "_energydeposition_ioncount_thread" + std::to_string(threadID)
    + ".csv"
  ;
  std::ofstream dataFile;
  dataFile.open(depositionPath, std::ios_base::out); // Open file in write mode to overwrite any previous results

  // Write results
  for(int altitudeIndex = 0; altitudeIndex < fNumberOfSamplePlanes-1; altitudeIndex++){
    dataFile << totalEnergyDeposition[altitudeIndex] << "," << ionizingEnergyDeposition[altitudeIndex] << "," << ionCounts[altitudeIndex] << "\n";
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
      ANSI_RED <<
      __FILE__ << ": " << __FUNCTION__ << "\n" <<
      "ERROR: Incomplete backscatter data! You should never see this." <<
      ANSI_NOCOLOR <<
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
          G4String countsAsString = std::to_string(mainSpectrum[particleIndex][altitudeIndex][energyIndex][paIndex]);

          // Remove trailing zeroes after decimal point for filesize considerations
          std::string::size_type trimAfter = countsAsString.find_last_not_of('0');
          G4String toWrite = countsAsString.substr(0, trimAfter+1);
          if(toWrite == "0."){toWrite = "0";}
          
          data.append(toWrite);
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

void RunAction::mergeEnergyDepositionAndIonCount(){
  G4cout << "Merging energy deposition and ion count..." << G4endl;

  // Progress bar variables
  int mergedCount = 0;
  int barLength = 30;

  // Loop over threads and add up histograms
  int nThreads = G4Threading::GetNumberOfRunningWorkerThreads();
  for(int thread = 0; thread < nThreads; thread++){
    std::string resultsFilepath =
      fBaseResultPath 
      + "_energydeposition_ioncount_thread" + std::to_string(thread)
      + ".csv"
    ;

    // Get data from this thread and add it to the result
    std::vector<std::vector<G4String>> eDepThreadDataString;

    // Get size of energy deposition file
    std::tuple<int, int> eDepVectorSize = get2DcsvSize(resultsFilepath);
    int dim1Size = std::get<0>(eDepVectorSize);
    int dim2Size = std::get<1>(eDepVectorSize);

    eDepThreadDataString.resize(dim2Size, std::vector<G4String>(0, ""));
    eDepThreadDataString = append2DcsvToString(resultsFilepath, eDepThreadDataString);
    std::vector<std::vector<G4double>> eDepThreadData = convert2DVectorStringToDouble(eDepThreadDataString);
    
    // Allocate vectors for individual columns
    std::vector<G4double> threadTotalEnergyDeposition = eDepThreadData.at(0);
    std::vector<G4double> threadIonizingEnergyDeposition = eDepThreadData.at(1);
    std::vector<G4double> threadIonCounts = eDepThreadData.at(2);

    // Add to main histograms
    totalEnergyDeposition = add1DAltitudeVectors(totalEnergyDeposition, threadTotalEnergyDeposition);
    ionizingEnergyDeposition = add1DAltitudeVectors(ionizingEnergyDeposition, threadIonizingEnergyDeposition);
    ionCounts = add1DAltitudeVectors(ionCounts, threadIonCounts);

    // Progress bar
    mergedCount++;
    double fraction = mergedCount / static_cast<double>(nThreads);
    printProgressBar(fraction, barLength);

    // Delete this thread-specific file
    std::remove(resultsFilepath.c_str());
  }
  G4cout << G4endl;

  // Write summed result to file
  std::string mainFilename =
    fBaseResultPath
    + "_energydeposition_ioncount.csv"
  ;
  std::ofstream dataFile;
  dataFile.open(mainFilename, std::ios_base::out); // Open file in write mode to overwrite any previous results

  // Write header
  dataFile << "Altitude (km),Total Energy Deposition (keV),Ionizing Energy Deposition (keV),Electron Production (counts)" << "\n";

  // Write rows
  for(int altitudeIndex = 0; altitudeIndex < fNumberOfSamplePlanes-1; altitudeIndex++){
    dataFile 
      << energyDepositionBinEdges[altitudeIndex] + fAltitudeOffset << "km-" 
      << energyDepositionBinEdges[altitudeIndex+1] + fAltitudeOffset << "km,"
      << totalEnergyDeposition[altitudeIndex] << ","
      << ionizingEnergyDeposition[altitudeIndex] << ","
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
  std::vector<std::vector<G4String>> result;
  result.resize(10, std::vector<G4String>(0, ""));

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
    result = append2DcsvToString(threadFilepath, result);

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

std::vector<std::vector<G4String>> RunAction::append2DcsvToString(std::string path, std::vector<std::vector<G4String>> result){
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

std::vector<std::vector<G4double>> RunAction::convert2DVectorStringToDouble(std::vector<std::vector<G4String>> input){
  int dim1Size = input[0].size();
  int dim2Size = input.size();

  std::vector<std::vector<G4double>> result = std::vector<std::vector<G4double>>(dim2Size, std::vector<G4double>(dim1Size, 0));
  for(int i = 0; i < dim1Size; i++){
    for(int j = 0; j < dim2Size; j++){
      result[j][i] = std::stod(input[j][i]);
    }
  }
  return result;
}

std::tuple<int, int> RunAction::get2DcsvSize(std::string path){
  // Open file
  std::ifstream file;
  file.open(path, std::ifstream::in);

  // Parse lines
  int dim1Index = 0;
  int dim2Index = 0;

  // Result variables
  int dim1Size = 0;
  int dim2Size = 0;

  std::string line;
  std::string token;
  while( std::getline(file, line) ){
    std::istringstream word(line);
    while ( std::getline(word, token, ',') ){
      dim2Index++;
      if(dim2Index > dim2Size){dim2Size = dim2Index;}
    }
    dim1Index++;
    dim1Size = dim1Index;
    dim2Index = 0;
  }
  return {dim1Size, dim2Size};
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
  G4double stepSize = (stop - start) / (static_cast<G4double>(n)-1.0);

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
    << static_cast<int>(100.0 * std::round(fraction * 1000.0) / 1000) << "%"
    << "\r"
    << std::flush
  ;
  if(fraction == 1.0){
    G4cout << G4endl;
  }
}