#ifndef RunActionMessenger_h
#define RunActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"


class RunAction;
class G4UIdirectory;
class G4UIcmdWithAString;
class RunActionMessenger : public G4UImessenger
{
public:
  RunActionMessenger(RunAction* );
  virtual ~RunActionMessenger();
  virtual void SetNewValue(G4UIcommand*, G4String);

private:
  RunAction*                 fRunAction;
  G4UIdirectory*             fPrimDir;
  G4UIcmdWithAString*        fcmd1;
  G4UIcmdWithAString*        fcmd2;
};

#endif