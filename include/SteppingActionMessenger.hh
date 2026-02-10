#ifndef SteppingActionMessenger_h
#define SteppingActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"


class SteppingAction;
class G4UIdirectory;
class G4UIcmdWithADouble;
class G4UIcmdWithAnInteger;
class G4UIcmdWithAString;
class SteppingActionMessenger : public G4UImessenger
{
public:
  SteppingActionMessenger(SteppingAction* );
  virtual ~SteppingActionMessenger();

  virtual void SetNewValue(G4UIcommand*, G4String);


private:
  SteppingAction* 	       fSteppingAction;
  G4UIdirectory*           fSteppingActionDir;
  G4UIcmdWithAString*      fcmd1;
  G4UIcmdWithADouble*      fcmd2;

};

#endif
