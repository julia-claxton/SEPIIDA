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
#include "F03FieldSetup.hh"
#include "G4LogicalVolume.hh"
#include "G4SystemOfUnits.hh"
#include "SteppingActionMessenger.hh"
#include "G4AutoLock.hh"
#include "G4TransportationManager.hh"
#include "G4FieldManager.hh"
#include "G4MagneticField.hh"
#include <chrono>
#include "ANSIColors.h"
#include "GlobalFunctions.h"
#include "TrackInformation.hh"

SteppingAction::SteppingAction(EventAction* eventAction, RunAction* RuAct)
: G4UserSteppingAction(),
  fEventAction(eventAction),
  fRunAction(RuAct),
  fBackscatterFilename(),
  fSteppingMessenger()
{
  fSteppingMessenger = new SteppingActionMessenger(this);
  uncachedField = new CustomMagneticField();
}

SteppingAction::~SteppingAction(){delete fSteppingMessenger;}

void SteppingAction::UserSteppingAction(const G4Step* step){
  // Dividing by a unit outputs data in that unit, so divisions by keV result in outputs in keV
  // https://geant4-internal.web.cern.ch/sites/default/files/geant4/collaboration/working_groups/electromagnetic/gallery/units/SystemOfUnits.html
  G4Track* track = step->GetTrack();
  TrackInformation* trackInfo = (TrackInformation*) track->GetUserInformation();
  const G4String particleName = track->GetDynamicParticle()->GetDefinition()->GetParticleName();
  const G4int trackID = track->GetTrackID();
  const G4int eventID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
  const G4VProcess* parentProcess = track->GetCreatorProcess();
  
  const G4ThreeVector position = track->GetPosition();
  const G4ThreeVector momentumDirection = track->GetMomentumDirection();

  const G4ThreeVector prePosition = step->GetPreStepPoint()->GetPosition();
  const G4ThreeVector preMomentumDirection = step->GetPreStepPoint()->GetMomentumDirection();
  const G4double preStepKineticEnergy = step->GetPreStepPoint()->GetKineticEnergy();

  const G4ThreeVector postPosition = step->GetPostStepPoint()->GetPosition();
  const G4ThreeVector postMomentumDirection = step->GetPostStepPoint()->GetMomentumDirection();
  const G4double postStepKineticEnergy = step->GetPostStepPoint()->GetKineticEnergy();
  
  const G4double trackWeight = track->GetWeight();
  const G4String particleIdentifier = std::to_string(eventID) + ":" + std::to_string(trackID);
  
  // Altitude information
  const G4double preStepAlt_km  = (prePosition.z()/km) + 500.0; // Add 500 because the coordinate origin is in center of the 1000 km atmospheric column
  const G4double postStepAlt_km = (postPosition.z()/km) + 500.0;

  // Get altitude indices (float) of start and stop point
  const G4double preStepAltitudeIndex = (preStepAlt_km -  fRunAction->fMinSampleAltitude_km) / fRunAction->altitudeSpacing_km;
  const G4double postStepAltitudeIndex = (postStepAlt_km -  fRunAction->fMinSampleAltitude_km) / fRunAction->altitudeSpacing_km;

  // Attach my user information to the track if it doesn't have it yet, then retrieve it
  if(!trackInfo){
    G4VUserTrackInformation* newTrackInfo = new TrackInformation;
    track->SetUserInformation(newTrackInfo);
    trackInfo = (TrackInformation*) track->GetUserInformation();
  }

  // ===========================
  // Guard Block
  // ===========================
  // Check for NaN momentum
  if(std::isnan(momentumDirection[0]) || std::isnan(momentumDirection[1]) || std::isnan(momentumDirection[2])){
    G4cout << ANSI_YELLOW <<
      __FILE__ << ": " << __FUNCTION__ << "\n" <<
      "WARNING: Killed " << particleName << " found with NaN momentum." <<
    ANSI_NOCOLOR << G4endl;
    track->SetTrackStatus(fStopAndKill);
  }

  // Check for NaN energy
  if(std::isnan(postStepKineticEnergy)){  
    G4cout << ANSI_YELLOW <<
      __FILE__ << ": " << __FUNCTION__ << "\n" <<
      "WARNING: Killed " << particleName << " found with NaN energy." <<
    ANSI_NOCOLOR << G4endl;
    track->SetTrackStatus(fStopAndKill);
  }

  // Check for exceeding 1 second of simulation time
  if(track->GetProperTime()/second > 1){
    G4cout << ANSI_YELLOW <<
      __FILE__ << ": " << __FUNCTION__ << "\n" <<
      "WARNING: Killed " << particleName << " at " << postStepKineticEnergy/keV << " keV exceeding 1s simulation time." <<
    ANSI_NOCOLOR << G4endl;
    track->SetTrackStatus(fStopAndKill);
  }

  // Check for stuck photons. Occassionally they seem to get 'wedged' between atmospheric layers and stop propagating without being automatically killed, hanging the program forever
  if((step->GetStepLength()/m < 1e-12) && particleName == "gamma"){
    G4cout << ANSI_YELLOW <<
      __FILE__ << ": " << __FUNCTION__ << "\n" <<
      "WARNING: Killed " << particleName << " at " << postStepKineticEnergy/keV << " keV. Reason: Stuck gamma. Current step length: " << step->GetStepLength()/m << " m" <<
    ANSI_NOCOLOR << G4endl;
    track->SetTrackStatus(fStopAndKill);
  }

  // ===========================
  // Account for cached B-fields
  // ===========================
  // If the magnetic field is cached, we need to manually adjust the pitch angle to produce mirroring
  // We use adiabatic pitch angle change (sin2(alpha1)/sin2(alpha2) = B1/B2) for this, which is nearly
  // exact (with minor inaccuracy from CSDA/continuous processes, but I don't think that will cause any
  // problems).
  if(CACHED_MAGNETIC_FIELD){
    G4bool firstStep = step->GetPreStepPoint()->GetProcessDefinedStep() ? false : true;
    G4String preProcessName = firstStep ? "None" : step->GetPreStepPoint()->GetProcessDefinedStep()->GetProcessName();
    G4String postProcessName = step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
    G4bool pureTransportationStep = (preProcessName == "Transportation") && (postProcessName == "Transportation");
    
    G4bool applyNudge = (track->GetDynamicParticle()->GetCharge() != 0) && (firstStep || pureTransportationStep); // Particle is charged and we are either on the first step of track's life or there were no processes during this step
    if(!applyNudge){needToUpdateBPre = true;}
    if(applyNudge){
      applyAdiabaticPitchAngleChange(track, step, prePosition, postPosition, preMomentumDirection, postMomentumDirection);
    }
  }

  // ===========================
  // Data recording
  // ===========================
  logSpectra(step, fRunAction,
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
    trackInfo,
    particleName,
    trackWeight,
    preStepAlt_km,
    postStepAlt_km,
    position,
    momentumDirection,
    postStepKineticEnergy
  );
  if((parentProcess != nullptr) && (particleName == "e-")){
    logIonProduction(fRunAction, trackInfo, preStepAltitudeIndex, trackWeight);
  }
}

void SteppingAction::logSpectra(const G4Step* step, RunAction* fRunAction,
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

  const G4ThreeVector preMomentum = step->GetPreStepPoint()->GetMomentumDirection();
  const G4ThreeVector postMomentum = step->GetPostStepPoint()->GetMomentumDirection();
  
  const G4double prePitchAngle = getPitchAngle(prePosition, preMomentum, uncachedField);
  const G4double postPitchAngle = getPitchAngle(postPosition, postMomentum, uncachedField);
  const G4double pitchAngleDiff = postPitchAngle - prePitchAngle;

  // Loop over crossed planes and add to energy spectra
  for(int altitudeIndex = startIdx; altitudeIndex <= stopIdx; altitudeIndex++){
    // Kick out invalid indices
    if( (altitudeIndex < 0) || (altitudeIndex > (fRunAction->fNumberOfSamplePlanes-1)) ){continue;}

    // Do an interpolation to get approximate energy at the plane crossing
    G4double t = (fRunAction->sampleAltitudes_km[altitudeIndex] - preStepAlt_km) / (postStepAlt_km - preStepAlt_km); // Interpolant, not time
    if(std::isnan(t)){t = 0;} // If we didn't travel in altitude, avoid killing the whole program with a NaN value of t

    // Perform linear interpolation for energy
    G4double crossingEnergy = preStepKineticEnergy + ((postStepKineticEnergy - preStepKineticEnergy) * t); 
    if( ((crossingEnergy/keV) < fRunAction->fEnergyMinkeV) || ((crossingEnergy/keV) > (fRunAction->fEnergyMaxkeV)) ){continue;} // Don't record energy if particle energy is out of range of histogram

    // Perform linear interpolation for pitch angle
    const G4double interpolatedPitchAngleDeg = prePitchAngle + (pitchAngleDiff * t);
    if( (interpolatedPitchAngleDeg < fRunAction->fPitchAngleMin) || (interpolatedPitchAngleDeg > (fRunAction->fPitchAngleMax)) ){continue;} // Don't record pitch angle if it's out of range

    // Find the energy bin this particle resides in utilizing regular spacing to directly calculate the index.
    int energyIndex = std::floor(logbase(fRunAction->energySpectrumHistogramFactor, (crossingEnergy/keV)/(fRunAction->fEnergyMinkeV)));

    // Find pitch angle bin
    int paIndex = std::floor((interpolatedPitchAngleDeg - fRunAction->fPitchAngleMin) / fRunAction->pitchAngleBinSize_deg);

    // Clip top of pitch angle range
    if(paIndex == fRunAction->fNumberOfPitchAngleBins){
      paIndex -= 1; // Technically should be the start of the first out-of-range bin, but that obviously causes problems, so we make the last bin inclusive on both sides of its range
    }
    
    // Add to histogram
    fRunAction->mainSpectrum[particleIdx][altitudeIndex][energyIndex][paIndex] += 1 * trackWeight;
  }
}

void SteppingAction::logIonProduction(RunAction* fRunAction,
  TrackInformation* trackInfo,
  const G4double preStepAltitudeIndex, 
  const G4double trackWeight
){
  if((std::floor(preStepAltitudeIndex) >= fRunAction->ionCounts.size()) || trackInfo->GetIonProductionLogFlag()){return;}
  trackInfo->SetIonProductionLogFlag(true);
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
  G4double weightedNonIonizingEnergyDeposition_keV = (step->GetNonIonizingEnergyDeposit() * trackWeight) / keV;
  G4double weightedIonizingEnergyDeposition_keV = weightedEnergyDeposition_keV - weightedNonIonizingEnergyDeposition_keV;

  G4double startIdx = std::min(preStepAltitudeIndex, postStepAltitudeIndex);
  G4double stopIdx = std::max(preStepAltitudeIndex, postStepAltitudeIndex);
  G4double traversedDistanceIdx = stopIdx - startIdx;

  // Catch steps with no movement to avoid NaNs
  if(startIdx == stopIdx){
    fRunAction->totalEnergyDeposition.at(std::floor(startIdx)) += weightedEnergyDeposition_keV;
    fRunAction->ionizingEnergyDeposition.at(std::floor(startIdx)) += weightedIonizingEnergyDeposition_keV;
    return;
  }

  for(int idx = std::floor(startIdx); idx <= std::floor(stopIdx); idx++){
    G4double fractionOfStepInThisIndex = overlap(startIdx, stopIdx, idx, idx+1) / traversedDistanceIdx;
    fRunAction->totalEnergyDeposition[idx] += weightedEnergyDeposition_keV * fractionOfStepInThisIndex;
    fRunAction->ionizingEnergyDeposition[idx] += weightedIonizingEnergyDeposition_keV * fractionOfStepInThisIndex;
  }
}

void SteppingAction::logBackscatter(const G4Step* step, RunAction* fRunAction,
  TrackInformation* trackInfo,
  const G4String particleName,
  const G4double trackWeight,
  const G4double preStepAlt_km,
  const G4double postStepAlt_km,
  const G4ThreeVector position,
  const G4ThreeVector momentumDirection,
  const G4double postStepKineticEnergy
){
  // Kick out non-backscattering particles and particles that have already been recorded
  bool backscattering = 
    (preStepAlt_km < fRunAction->fCollectionAltitude) 
    && (postStepAlt_km > fRunAction->fCollectionAltitude) 
    && (trackInfo->GetBackscatterLogFlag() == false)
  ;
  if(backscattering == false){return;}

  // If we've reached here, the particle is backscattering and should be logged
  trackInfo->SetBackscatterLogFlag(true);
  
  // Write particle parameters to memory
  G4double pitchAngleDeg = getPitchAngle(position, momentumDirection, uncachedField);
  fRunAction->fBackscatteredParticleNames.push_back(particleName);
  fRunAction->fBackscatteredTrackWeights.push_back(trackWeight);
  fRunAction->fBackscatteredEnergieskeV.push_back(postStepKineticEnergy/keV);
  fRunAction->fBackscatteredPitchAnglesDeg.push_back(pitchAngleDeg);
  fRunAction->fBackscatterDirections.push_back({momentumDirection.x(), momentumDirection.y(), momentumDirection.z()});
  fRunAction->fBackscatterPositions.push_back({position.x()/m, position.y()/m, (position.z()/m) + 500000.0}); // Shift z-axis so we are writing altitude above sea level to file rather than the world coordinates. Remember: World origin is at 500 km above sea level
}

void SteppingAction::applyAdiabaticPitchAngleChange(
  G4Track*      track,
  const G4Step* step,
  G4ThreeVector prePosition,
  G4ThreeVector postPosition,
  G4ThreeVector preMomentumDirection,
  G4ThreeVector postMomentumDirection
){
  // Get pre/post B vectors
  if(needToUpdateBPre){Bpre = getB(uncachedField, prePosition);}
  G4ThreeVector BvecPre = Bpre;
  G4ThreeVector BvecPost = getB(uncachedField, postPosition);
  
  // Calculate pitch angle at the end of this step under adiabatic theory
  G4double Bpre = BvecPre.mag();
  G4double Bpost = BvecPost.mag();
  G4double paPre_rad = getPitchAngle(preMomentumDirection, BvecPre) * (M_PI/180);    // Units: Radian
  G4double paPost_rad = getPitchAngle(postMomentumDirection, BvecPost) * (M_PI/180); // Units: Radian
  
  // Calculate desired pitch angle
  G4double for_asin = std::sqrt(Bpost/Bpre) * std::sin(paPre_rad);    
  for_asin = std::clamp(for_asin, -1.0, 1.0); // Account for float precision
  G4double desiredPaPost_rad = std::asin(for_asin);
  
  // Resolve quadrant ambiguity
  if(paPre_rad > (M_PI/2)){desiredPaPost_rad = M_PI - desiredPaPost_rad;} 
  G4double dPa_deg = (desiredPaPost_rad - paPost_rad) * (180/M_PI);

  // Get new momentum vector
  G4ThreeVector rotateAbout = BvecPost.cross(postMomentumDirection);

  // Apply the nudge and return
  G4ThreeVector newMomentumDirection = rotateVector(postMomentumDirection, rotateAbout, dPa_deg);
  step->GetPostStepPoint()->SetMomentumDirection(newMomentumDirection);
  track->SetMomentumDirection(newMomentumDirection);
  return;
}

G4ThreeVector SteppingAction::getB(G4MagneticField* field, G4ThreeVector position){
  G4double _position[4] = {position[0], position[1], position[2], 0};
  G4double EMField[6];
  field->GetFieldValue(_position, EMField);
  G4ThreeVector B(EMField[0], EMField[1], EMField[2]);

  return B;
}

G4ThreeVector SteppingAction::getB(G4ThreeVector position){
  // Uses global detector field if no field is specified
  G4double _position[4] = {position[0], position[1], position[2], 0};
  G4double EMField[6];
  
  G4FieldManager* fieldManager = G4TransportationManager::GetTransportationManager()->GetFieldManager();
  fieldManager->GetDetectorField()->GetFieldValue(_position, EMField);
  G4ThreeVector B(EMField[0], EMField[1], EMField[2]);

  return B;
}

G4double SteppingAction::getBMagnitude(G4MagneticField* field, G4ThreeVector position){
  G4ThreeVector B = getB(field, position);
  return B.mag();
}

G4double SteppingAction::getBMagnitude(G4ThreeVector position){
  G4ThreeVector B = getB(position);
  return B.mag();
}

G4double SteppingAction::getPitchAngle(G4ThreeVector position, G4ThreeVector momentumDirection, G4MagneticField* field){  
  G4ThreeVector B = getB(field, position);
  return getPitchAngle(momentumDirection, B);
}

G4double SteppingAction::getPitchAngle(G4ThreeVector momentumDirection, G4ThreeVector B){ 
  G4double pitchAngleDeg;
  G4double for_acos = momentumDirection.dot(B) / (momentumDirection.mag() * B.mag());
  
  // Throw if for_acos is out-of-range by significantly more than float precision
  if(std::abs(for_acos) > (1.0 + 1e-5)){
    G4cout << ANSI_RED <<
      __FILE__ << ": " << __FUNCTION__ << "\n" <<
      "ERROR: Value provided to std::acos() out of range. std::abs(for_acos) - 1.0 = " << std::abs(for_acos) - 1.0 <<
    ANSI_NOCOLOR << G4endl;
    throw;
  }

  // Float precision near 0º and 180º can cause out-of-domain errors resulting in NaNs, clip values near the edges of the domain
  if (for_acos > 1.0)      {pitchAngleDeg = 0.0;}
  else if(for_acos < -1.0) {pitchAngleDeg = 180.0;}
  else                     {pitchAngleDeg = std::acos(for_acos) * 180/M_PI;}

  // Guard
  if(std::isnan(pitchAngleDeg)){
    G4cout << ANSI_RED <<
      __FILE__ << ": " << __FUNCTION__ << "\n" <<
      "ERROR: Pitch angle is nan. This should never happen." <<
    ANSI_NOCOLOR << G4endl;
    G4cout <<
      "Debug Info:" << G4endl <<
      "    v = " << momentumDirection[0] << ", " << momentumDirection[1] << ", " << momentumDirection[2] << G4endl <<
      "    B = " << B[0] << ", " << B[1] << ", " << B[2] << G4endl <<
      "    pitchAngleDeg = " << pitchAngleDeg << G4endl <<
      "    for_acos = " << for_acos <<
    G4endl;
    throw;
  }

  // Return
  return pitchAngleDeg;
}

G4double SteppingAction::logbase(G4double base, G4double x){
  // Log base of x, log_base(x)
  return log2(x)/log2(base);
}

G4double SteppingAction::overlap(G4double a1, G4double a2, G4double b1, G4double b2){
  // min of maxes - max of mins
  if(a1 > a2){
    G4cout << ANSI_RED <<
      __FILE__ << ": " << __FUNCTION__ << "\n" <<
      "Ranges must be sorted!" <<
    ANSI_NOCOLOR << G4endl;
    throw;
  }
  if(b1 > b2){
    G4cout << ANSI_RED <<
      __FILE__ << ": " << __FUNCTION__ << "\n" <<
      "Ranges must be sorted!" <<
    ANSI_NOCOLOR << G4endl;
    throw;
  }

  return std::min(a2, b2) - std::max(a1, b1);
}