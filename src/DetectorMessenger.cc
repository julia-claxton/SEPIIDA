

#include "DetectorMessenger.hh"

#include "DetectorConstruction.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIdirectory.hh"

DetectorMessenger::DetectorMessenger(DetectorConstruction* detCon)
  : G4UImessenger(),
    fDetectorMessenger(detCon) 
{
  fPrimDir = new G4UIdirectory("/atmosphere/");

  fcmd1 = new G4UIcmdWithAString("/atmosphere/setFilename", this);
  fcmd1->AvailableForStates(G4State_PreInit, G4State_Init, G4State_Idle);
}

DetectorMessenger::~DetectorMessenger()
{
  delete fPrimDir;
  delete fcmd1;
}

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if(command == fcmd1){fDetectorMessenger->SetAtmosphereFilename(newValue);}    	  

}
