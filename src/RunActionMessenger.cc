#include "RunActionMessenger.hh"
#include "RunAction.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIdirectory.hh"


RunActionMessenger::RunActionMessenger(RunAction* runAct)
  : G4UImessenger(),
    fRunAction(runAct) 
{
  fPrimDir = new G4UIdirectory("/dataCollection/");

  fcmd1 = new G4UIcmdWithAString("/dataCollection/setBaseResultPath", this);
  fcmd1->AvailableForStates(G4State_PreInit, G4State_Idle);

  fcmd2 = new G4UIcmdWithAString("/dataCollection/setCollectionAltitude", this);
  fcmd2->SetDefaultValue("451.0");

  fcmd3 = new G4UIcmdWithAString("/dataCollection/setAltitudeOffset", this);
  fcmd3->SetDefaultValue("0.0");
}

RunActionMessenger::~RunActionMessenger()
{
  delete fPrimDir;
  delete fcmd1;
  delete fcmd2;
  delete fcmd3;
}

void RunActionMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if(command == fcmd1){fRunAction->SetBaseResultPath(newValue);}
  if(command == fcmd2){fRunAction->SetCollectionAltitude(std::stod(newValue));}
  if(command == fcmd3){fRunAction->SetAltitudeOffset(std::stod(newValue));}
}