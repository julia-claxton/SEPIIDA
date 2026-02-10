

#include "DetectorMessenger.hh"

#include "DetectorConstruction.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIdirectory.hh"

DetectorMessenger::DetectorMessenger(DetectorConstruction* detCon)
  : G4UImessenger(),
    fDetectorMessenger(detCon) 
{
  fPrimDir = new G4UIdirectory("/dataCollection/");
  fPrimDir->SetGuidance("Choose file to generate atmosphere.");

  fcmd1 = new G4UIcmdWithAString("/dataCollection/setAtmosFileName",this);
  
  fcmd1->SetParameterName("Enter atmosphere file name.",true);
  fcmd1->SetDefaultValue("MSISE00_atmosphere.csv");
  fcmd1->AvailableForStates(G4State_PreInit);
}

DetectorMessenger::~DetectorMessenger()
{
  delete fPrimDir;
  delete fcmd1;
}

void DetectorMessenger::SetNewValue(G4UIcommand* command, 
					    G4String newValue)
{
  if(command == fcmd1){
    fDetectorMessenger->SetAtmosFilename(newValue);
  }    	  

}
