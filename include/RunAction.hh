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
// $Id: RunAction.hh 99560 2016-09-27 07:03:29Z gcosmo $
//
/// \file RunAction.hh
/// \brief Definition of the RunAction class

#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "G4Accumulable.hh"
#include "globals.hh"
#include <vector>
#include "myHistogram.hh"

class G4ParticleDefinition;
class G4Transportation;
class G4CoupledTransportation;
class G4Run;
class SteppingAction;
class RunActionMessenger;

class RunAction : public G4UserRunAction
{
  public:
    RunAction();
    virtual ~RunAction();
    virtual void BeginOfRunAction(const G4Run*);
    virtual void   EndOfRunAction(const G4Run*);
    void ChangeLooperParameters(const G4ParticleDefinition* particleDef ); // Helper method to change the Transportation's 'looper' parameters 
    std::pair<G4Transportation*, G4CoupledTransportation*> findTransportation(const G4ParticleDefinition * particleDef, bool reportError= true); // Helper method to find the Transportation process for a particle type 
    
    // Messengers
    void SetCollectionAltitude(G4double collectionAltitude){fCollectionAltitude = collectionAltitude;};
    void SetBaseResultPath(G4String name){fBaseResultPath=name;};

    // Data writers
    void writeAxisLabelHeader(std::ofstream& file);

    void threadWriteSpectra(int threadID);
    void threadWriteEnergyDeposition(int threadID);
    void threadWriteIonProduction(int threadID);
    void threadWriteBackscatter(int threadID);

    void mergeEnergySpectra();
    void mergeEnergyDeposition();
    void mergeBackscatter();

    void writeRunMetadata(G4String filepath);

    // Data readers
    void addThreadSpectraToMainHistogram(std::string path, int particleIndex);
    std::vector<std::vector<std::string>> read2DcsvBackscatter(std::string path, std::vector<std::vector<std::string>> result);
    std::vector<G4double> read1Dcsv(std::string path, std::vector<G4double> result);
    
    // Misc
    std::vector<std::vector<G4double>> addEnergySpectraVectors(std::vector<std::vector<G4double>> v1, std::vector<std::vector<G4double>> v2);
    std::vector<G4double> add1DAltitudeVectors(std::vector<G4double> v1, std::vector<G4double> v2);
    std::vector<G4double> linspace(G4double start, G4double stop, int n);
    void printProgressBar(double fraction, int barLength);
    void printTimestamp();

  public:
    void     SetNumberOfTrials( G4int val ){fNumberOfTrials  = val;}
    void     SetWarningEnergy( double val ){fWarningEnergy   = val;}
    void     SetImportantEnergy( double val ){fImportantEnergy = val;}   
    G4int    GetNumberOfTrials(){ return fNumberOfTrials; }
    G4double GetWarningEnergy(){ return fWarningEnergy; }
    G4double GetImportantEnergy(){ return fImportantEnergy; } 

  public:
    G4double fCollectionAltitude;

    // Histogram limits
    // Sample planes are also bin edges for energy deposition histogram
    static constexpr G4double fMinSampleAltitude_km = 0.0;
    static constexpr G4double fMaxSampleAltitude_km = 1000.0;
    static constexpr G4int fNumberOfSamplePlanes = 1001; //1001; // 1 plane per km

    static constexpr G4double fEnergyMinkeV = 1e-3; // 1 eV
    static constexpr G4double fEnergyMaxkeV = 100e6; // 100 GeV
    static constexpr G4int fNumberOfEnergyBins = 180; // 20 bins per decade

    static constexpr G4double fPitchAngleMin = 0.0;
    static constexpr G4double fPitchAngleMax = 180.0;
    static constexpr G4int fNumberOfPitchAngleBins = 36; // 5º per bin

    // Axis labels
    std::vector<G4double> sampleAltitudes_km;
    std::vector<G4double> energyDepositionBinEdges; 
    std::vector<G4double> energyBinEdges_keV;
    std::vector<G4double> pitchAngleBinEdges_deg;

    // Histograms
    std::vector<std::string> particlesToRecord = {"proton", "e-", "alpha", "gamma", "neutron"};
    std::vector<std::vector<std::vector<std::vector<G4double>>>> mainSpectrum; // Dimensions: particle species, altitude, energy, pitch angle
    std::vector<double> totalEnergyDeposition; 
    std::vector<double> ionizingEnergyDeposition; 
    std::vector<double> ionCounts;

    std::vector<std::string> fBackscatteredParticleNames;
    std::vector<double> fBackscatteredTrackWeights;
    std::vector<double> fBackscatteredEnergieskeV;
    std::vector<double> fBackscatteredPitchAnglesDeg;
    std::vector<std::array<double,3>> fBackscatterDirections;
    std::vector<std::array<double,3>> fBackscatterPositions;

    // Precalculated for speed
    G4double energySpectrumHistogramFactor; // Factor we multiply by for each successive energy histogram edge
    G4double altitudeSpacing_km; // Space between sample altitudes in km
    G4double pitchAngleBinSize_deg; // Size of pitch angle bins

  private:
    RunActionMessenger* fRunActionMessenger;
    G4String fBaseResultPath;
    G4String fEnergyDepositionFileName;

    // Values for initialising 'loopers' parameters of Transport process
    G4int    fNumberOfTrials  =  0;    // Default will not overwrite
    G4double fWarningEnergy   = -1.0;  // Default values - non operational 
    G4double fImportantEnergy = -1.0;  // Default - will not overwrite

    int theVerboseLevel = 0;
};

#endif