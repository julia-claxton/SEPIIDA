
#ifndef DetectorMessenger_h
#define DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class DetectorConstruction;

class G4UIdirectory;
class G4UIcmdWithAString;

class DetectorMessenger : public G4UImessenger
{
public:
  DetectorMessenger(DetectorConstruction* );
  virtual ~DetectorMessenger();

  virtual void SetNewValue(G4UIcommand*, G4String);


private:
  DetectorConstruction*      fDetectorMessenger;
  G4UIdirectory*             fPrimDir;
  G4UIcmdWithAString*        fcmd1;

};

#endif
