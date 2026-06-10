#include "TrackInformation.hh"

TrackInformation::TrackInformation() : 
  G4VUserTrackInformation(),
  hasBeenBackscatterLogged(false),
  hasBeenIonProductionLogged(false)
{}

TrackInformation::~TrackInformation(){}