

#include "PrimaryGeneratorMessenger.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIdirectory.hh"


PrimaryGeneratorMessenger::PrimaryGeneratorMessenger(PrimaryGeneratorAction* prim)
  :G4UImessenger(),
  fPrimaryGenerator(prim) 
{
  fPrimDir = new G4UIdirectory("/beamParameters/");
  fPrimDir->SetGuidance("Select particle distributions.");
  
  fcmd3 = new G4UIcmdWithAString("/beamParameters/setBeamParticle",this);
  fcmd3->SetParameterName("Particle to input. e- = electron, proton = proton, gamma = photon", true);
  fcmd3->SetDefaultValue("e-");
  fcmd3->AvailableForStates(G4State_PreInit, G4State_Idle);
  
  fDcmd = new G4UIcmdWithADouble("/beamParameters/setBeamEnergy",this);
  fDcmd->SetParameterName("Beam energy [keV]",true);
  fDcmd->SetDefaultValue(100.0);
  fDcmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  
  fDcmd2 = new G4UIcmdWithADouble("/beamParameters/setBeamPitchAngle",this);
  fDcmd2->SetParameterName("Beam pitch angle [deg]",true);
  fDcmd2->SetDefaultValue(0.0);
  fDcmd2->AvailableForStates(G4State_PreInit, G4State_Idle);

  fDcmd3 = new G4UIcmdWithADouble("/beamParameters/setParticleStartingAltitude",this);
  fDcmd3->SetParameterName("Particle initial altitude [km]",true);
  fDcmd3->SetDefaultValue(450.0);
  fDcmd3->AvailableForStates(G4State_PreInit, G4State_Idle);
}

PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger()
{
  delete fPrimDir;
  delete fcmd3;
  delete fDcmd;
  delete fDcmd2;
  delete fDcmd3;
}
void PrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if(command == fDcmd) {fPrimaryGenerator->SetBeamEnergy(std::stod(newValue));}
  if(command == fDcmd2){fPrimaryGenerator->SetBeamPitchAngle(std::stod(newValue));}
  if(command == fDcmd3){fPrimaryGenerator->SetParticleInitialAlt(std::stod(newValue));}
  if(command == fcmd3) {fPrimaryGenerator->SetInputParticleType(newValue);}
}
