#ifndef INCLUDE_TRACKINFORMATION_H 
#define INCLUDE_TRACKINFORMATION_H 1

#include "G4Types.hh"
#include "G4VUserTrackInformation.hh"

class TrackInformation : public G4VUserTrackInformation {
public:
  TrackInformation(); 
  virtual ~TrackInformation();

  void Print() const {};

  void SetBackscatterLogFlag(G4bool value) {hasBeenBackscatterLogged = value;}
  G4bool GetBackscatterLogFlag() const {return hasBeenBackscatterLogged;}

  void SetIonProductionLogFlag(G4bool value) {hasBeenIonProductionLogged = value;}
  G4bool GetIonProductionLogFlag() const {return hasBeenIonProductionLogged;}

private:
  G4bool hasBeenBackscatterLogged;
  G4bool hasBeenIonProductionLogged;
};

#endif
