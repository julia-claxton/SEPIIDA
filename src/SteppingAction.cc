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
// $Id: SteppingAction.cc 74483 2013-10-09 13:37:06Z gcosmo $
//
/// \file SteppingAction.cc
/// \brief Implementation of the SteppingAction class

#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "RunAction.hh"
#include "G4LogicalVolume.hh"
#include "G4SystemOfUnits.hh"
#include "SteppingActionMessenger.hh"
#include "G4AutoLock.hh"
#include "G4TransportationManager.hh"
#include "G4FieldManager.hh"
#include "G4MagneticField.hh"
#include <chrono>

SteppingAction::SteppingAction(EventAction* eventAction, RunAction* RuAct)
: G4UserSteppingAction(),
  fEventAction(eventAction),
  fRunAction(RuAct),
  fBackscatterFilename(),
  fSteppingMessenger()
{
  fSteppingMessenger = new SteppingActionMessenger(this);
}

SteppingAction::~SteppingAction(){delete fSteppingMessenger;}

void SteppingAction::UserSteppingAction(const G4Step* step)
{
  // Dividing by a unit outputs data in that unit, so divisions by keV result in outputs in keV
  // https://geant4-internal.web.cern.ch/sites/default/files/geant4/collaboration/working_groups/electromagnetic/gallery/units/SystemOfUnits.html
  G4Track* track = step->GetTrack();
  const G4int trackID = track->GetTrackID();
  const G4VProcess* parentProcess = step->GetTrack()->GetCreatorProcess();
  const G4ThreeVector position = track->GetPosition();
  const G4ThreeVector momentumDirection = track->GetMomentumDirection();
  const G4String particleName = track->GetDynamicParticle()->GetDefinition()->GetParticleName();
  const G4double preStepKineticEnergy = step->GetPreStepPoint()->GetKineticEnergy();
  const G4double postStepKineticEnergy = step->GetPostStepPoint()->GetKineticEnergy();
  const G4double trackWeight = track->GetWeight();
  
  // Altitude information
  const G4double preStepAlt_km  = (step->GetPreStepPoint()->GetPosition().z()/km) + 500.0; // Add 500 because the coordinate origin is in center of the 1000 km atmospheric column
  const G4double postStepAlt_km = (step->GetPostStepPoint()->GetPosition().z()/km) + 500.0;

  // Get altitude indices (float) of start and stop point
  const G4double preStepAltitudeIndex = (preStepAlt_km -  fRunAction->fMinSampleAltitude_km) / fRunAction->altitudeSpacing_km;
  const G4double postStepAltitudeIndex = (postStepAlt_km -  fRunAction->fMinSampleAltitude_km) / fRunAction->altitudeSpacing_km;

  // ===========================
  // Guard Block
  // ===========================

  // Check for NaN energy
  if(std::isnan(postStepKineticEnergy))
  {  
    G4cout << "WARNING: Killed " << particleName << "at " << postStepKineticEnergy/keV << " keV. Reason: NaN energy. Process: " << step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() << G4endl;
    track->SetTrackStatus(fStopAndKill);
  }
  // Check for exceeding 1 second of simulation time
  if(track->GetProperTime()/second > 1)
  {
    G4cout << "WARNING: Killed " << particleName << "at " << postStepKineticEnergy/keV << " keV. Reason: Exceeded 1s simulation time. Process: " << step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() << G4endl;
    track->SetTrackStatus(fStopAndKill);
  }
  // Check for stuck photons. Occassionally they seem to get 'wedged' between atmospheric layers and stop propagating without being automatically killed, hanging the program forever
  if((step->GetStepLength()/m < 1e-12) && particleName == "gamma"){
    G4cout << "WARNING: Killed " << particleName << " at " << postStepKineticEnergy/keV << " keV. Reason: Stuck gamma. Current step length: " << step->GetStepLength()/m << " m" << G4endl;
    track->SetTrackStatus(fStopAndKill);
  }

  // ===========================
  // Data recording
  // ===========================
  logEnergySpectra(step, fRunAction,
    particleName,
    trackWeight,
    preStepKineticEnergy,
    postStepKineticEnergy,
    preStepAltitudeIndex,
    postStepAltitudeIndex,
    preStepAlt_km,
    postStepAlt_km
  );
  logEnergyDeposition(step, fRunAction,
    trackWeight,
    preStepAltitudeIndex,
    postStepAltitudeIndex
  );
  logBackscatter(step, fRunAction,
    particleName,
    trackWeight,
    preStepAlt_km,
    postStepAlt_km,
    position,
    momentumDirection,
    postStepKineticEnergy
  );
  if((parentProcess != nullptr) && (particleName == "e-") && (ionizationTracks[trackID] == false)){
    ionizationTracks[trackID] = true;
    logIonProduction(fRunAction, preStepAltitudeIndex, trackWeight);
  }
}

void SteppingAction::logEnergySpectra(const G4Step* step, RunAction* fRunAction,
  const G4String particleName,
  const G4double trackWeight,
  const G4double preStepKineticEnergy,
  const G4double postStepKineticEnergy,
  const G4double preStepAltitudeIndex,
  const G4double postStepAltitudeIndex,
  const G4double preStepAlt_km,
  const G4double postStepAlt_km
){
  // Kick out if particle isn't one we care about
  auto particleIdx = std::find(fRunAction->particlesToRecord.begin(), fRunAction->particlesToRecord.end(), particleName) - fRunAction->particlesToRecord.begin();
  if(particleIdx == fRunAction->particlesToRecord.size()){return;} // Equalling the size means the index is out of bounds because of zero-indexing

  // Kick out if step is entirely outside the altitudes we care about
  bool preStepInRange = (0 <= preStepAltitudeIndex) && (preStepAltitudeIndex < (fRunAction->fNumberOfSamplePlanes-1));
  bool postStepInRange = (0 <= postStepAltitudeIndex) && (postStepAltitudeIndex < (fRunAction->fNumberOfSamplePlanes-1));
  if( (preStepInRange == false) && (postStepInRange == false) ){ return; }
  
  // Get bounding indices of planes that have been crossed
  int startIdx = std::ceil(std::min(preStepAltitudeIndex, postStepAltitudeIndex));
  int stopIdx = std::floor(std::max(preStepAltitudeIndex, postStepAltitudeIndex));

  // Kick out particles that didn't cross any planes
  if(startIdx > stopIdx){return;}

  // Precalculations for efficiency in the next loop
  const G4ThreeVector prePosition = step->GetPreStepPoint()->GetPosition();
  const G4ThreeVector postPosition = step->GetPostStepPoint()->GetPosition();
  const G4ThreeVector positionDiff = postPosition - prePosition;

  const G4ThreeVector preMomentum = step->GetPreStepPoint()->GetMomentumDirection();
  const G4ThreeVector postMomentum = step->GetPostStepPoint()->GetMomentumDirection();
  const G4ThreeVector momentumDiff = postMomentum - preMomentum;

  // Loop over crossed planes and add to energy spectra
  for(int altitudeIndex = startIdx; altitudeIndex <= stopIdx; altitudeIndex++){
    // Kick out invalid indices
    if( (altitudeIndex < 0) || (altitudeIndex > (fRunAction->fNumberOfSamplePlanes-1)) ){continue;}

    // Do an interpolation to get approximate energy at the plane crossing
    G4double t = (fRunAction->sampleAltitudes_km[altitudeIndex] - preStepAlt_km) / (postStepAlt_km - preStepAlt_km);
    G4double crossingEnergy = preStepKineticEnergy + (t * (postStepKineticEnergy - preStepKineticEnergy)); // Linear interpolation

    // Interpolate to get approximate pitch angle at time of crossing
    const G4ThreeVector interpolatedPosition = prePosition + (positionDiff * t);
    const G4ThreeVector interpolatedMomentum = preMomentum + (momentumDiff * t);
    G4double interpolatedPitchAngleDeg = getPitchAngle(interpolatedPosition, interpolatedMomentum);

    // Find the energy bin this particle resides in utilizing regular spacing to directly calculate the index.
    int energyIndex = std::floor(logbase(fRunAction->energySpectrumHistogramFactor, (crossingEnergy/keV)/(fRunAction->fEnergyMinkeV)));
    if( (energyIndex < 0) || (energyIndex > (fRunAction->fNumberOfEnergyBins-1)) ){continue;} // Don't record energy if particle energy is out of range of histogram

    // Find pitch angle bin
    int paIndex = std::floor((interpolatedPitchAngleDeg - fRunAction->fPitchAngleMin) / fRunAction->pitchAngleBinSize_deg);
    
    // Add to histogram
    fRunAction->mainSpectrum[particleIdx][altitudeIndex][energyIndex][paIndex] += 1 * trackWeight;
  }
}

void SteppingAction::logIonProduction(RunAction* fRunAction,
  const G4double preStepAltitudeIndex, 
  const G4double trackWeight
){
  fRunAction->ionCounts.at(std::floor(preStepAltitudeIndex)) += 1 * trackWeight;
}

void SteppingAction::logEnergyDeposition(const G4Step* step, RunAction* fRunAction,
  const G4double trackWeight,
  const G4double preStepAltitudeIndex,
  const G4double postStepAltitudeIndex
){
  // Kick out if step is entirely outside the altitudes we care about
  bool preStepInRange = (0 <= preStepAltitudeIndex) && (preStepAltitudeIndex < (fRunAction->fNumberOfSamplePlanes-1));
  bool postStepInRange = (0 <= postStepAltitudeIndex) && (postStepAltitudeIndex < (fRunAction->fNumberOfSamplePlanes-1));
  if( (preStepInRange == false) && (postStepInRange == false) ){ return; }

  // Track energy deposition, dividing the energy deposition between the crossed bins
  G4double weightedEnergyDeposition_keV = (step->GetTotalEnergyDeposit() * trackWeight) / keV;

  G4double startIdx = std::min(preStepAltitudeIndex, postStepAltitudeIndex);
  G4double stopIdx = std::max(preStepAltitudeIndex, postStepAltitudeIndex);
  G4double traversedDistanceIdx = stopIdx - startIdx;

  // Catch steps with no movement to avoid NaNs
  if(startIdx == stopIdx){
    fRunAction->energyDeposition.at(std::floor(startIdx)) += weightedEnergyDeposition_keV;
    return;
  }

  for(int idx = std::floor(startIdx); idx <= std::floor(stopIdx); idx++){
    G4double fractionOfStepInThisIndex = overlap(startIdx, stopIdx, idx, idx+1) / traversedDistanceIdx;
    fRunAction->energyDeposition.at(idx) += weightedEnergyDeposition_keV * fractionOfStepInThisIndex;
  }
}

void SteppingAction::logBackscatter(const G4Step* step, RunAction* fRunAction,
  const G4String particleName,
  const G4double trackWeight,
  const G4double preStepAlt_km,
  const G4double postStepAlt_km,
  const G4ThreeVector position,
  const G4ThreeVector momentumDirection,
  const G4double postStepKineticEnergy
){
  // Kick out non-backscattering particles
  bool backscattering = (preStepAlt_km < fRunAction->fCollectionAltitude) && (postStepAlt_km > fRunAction->fCollectionAltitude);
  if(backscattering == false){return;}

  G4double pitchAngleDeg = getPitchAngle(position, momentumDirection);

  // Write particle parameters to memory
  fRunAction->fBackscatteredParticleNames.push_back(particleName);
  fRunAction->fBackscatteredTrackWeights.push_back(trackWeight);
  fRunAction->fBackscatteredEnergieskeV.push_back(postStepKineticEnergy/keV);
  fRunAction->fBackscatteredPitchAnglesDeg.push_back(pitchAngleDeg);
  fRunAction->fBackscatterDirections.push_back({momentumDirection.x(), momentumDirection.y(), momentumDirection.z()});
  fRunAction->fBackscatterPositions.push_back({position.x()/m, position.y()/m, (position.z()/m) + 500000.0}); // Shift z-axis so we are writing altitude above sea level to file rather than the world coordinates. Remember: World origin is at 500 km above sea level
}

G4double SteppingAction::getPitchAngle(G4ThreeVector position, G4ThreeVector momentumDirection){
  G4double spacetimePoint[4] = {position.x(), position.y(), position.z(), 0};
  G4double emComponents[6];

  G4FieldManager* fieldManager = G4TransportationManager::GetTransportationManager()->GetFieldManager();
  fieldManager->GetDetectorField()->GetFieldValue(spacetimePoint, emComponents);
  G4double B[3] = {emComponents[0], emComponents[1], emComponents[2]};
  G4double normB = std::sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);

  G4double normMomentum = std::sqrt(pow(momentumDirection.x(), 2) + pow(momentumDirection.y(), 2) + pow(momentumDirection.z(), 2));
  G4double dotProd = (momentumDirection.x() * B[0]) + (momentumDirection.y() * B[1]) + (momentumDirection.z() * B[2]);

  const G4double pitchAngleDeg = std::acos(dotProd / (normMomentum * normB)) * 180/3.14159265358979;
  return pitchAngleDeg;
}

G4double SteppingAction::logbase(G4double base, G4double x){
  // Log base of x, log_base(x)
  return log2(x)/log2(base);
}

G4double SteppingAction::overlap(G4double a1, G4double a2, G4double b1, G4double b2){
  // min of maxes - max of mins
  if(a1 > a2){G4cout << "SteppingAction::overlapFraction: Ranges must be sorted!" << G4endl; throw;}
  if(b1 > b2){G4cout << "SteppingAction::overlapFraction: Ranges must be sorted!" << G4endl; throw;}

  return std::min(a2, b2) - std::max(a1, b1);
}

void SteppingAction::printVector(G4ThreeVector v, G4double unit){
  G4cout << "(" << v.x()/unit << ", " << v.y()/unit << ", " << v.z()/unit << ")";
}