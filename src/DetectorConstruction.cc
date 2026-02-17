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
// $Id: DetectorConstruction.cc 94307 2015-11-11 13:42:46Z gcosmo $
//
/// \file DetectorConstruction.cc
/// \brief Creates the atmosphere for EPP simulation from an MSIS-generated file

#include "DetectorConstruction.hh"

#include "F03FieldSetup.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4TransportationManager.hh"
#include "G4GeometryManager.hh"

#include "G4AutoDelete.hh"
#include "G4SDManager.hh"

#include "DetectorMessenger.hh"


DetectorConstruction::DetectorConstruction():
  G4VUserDetectorConstruction(),
  fAtmosphereFilename("atmosphere_profile.csv"),
  fDetectorMessenger(),
  fTableSize(0),
  fLogicWorld(0)
{
  fDetectorMessenger = new DetectorMessenger(this);
}


DetectorConstruction::~DetectorConstruction()
{
  delete fDetectorMessenger;
}


G4VPhysicalVolume* DetectorConstruction::Construct()
{
  G4GeometryManager::GetInstance()->SetWorldMaximumExtent(1000*km);

  // Material: Vacuum
  G4Material* vacuum_material = new G4Material(
    "Vacuum",
    1.0,
    1.01*g/mole,
    1.0E-25*g/cm3,
    kStateGas,
    2.73*kelvin,
    3.0E-18*pascal
  );

  G4Material* low_density_material = new G4Material(
    "Low_dens",
    1.0,
    1.01*g/mole,
    1.0E-10*g/cm3,
    kStateGas, 2.73*kelvin,
    3.0E-12*pascal
  );

  // Option to switch on/off checking of volumes overlaps
  G4bool checkOverlaps = false; // Set to true if you need to debug. Set false as there are no issues right now and it's very verbose.

  // World
  G4double world_sizeXY = 1000.0*km;
  G4double world_sizeZ  = 1000.0*km;

  G4Tubs* solidWorld = new G4Tubs(
    "World",          // its name
    0.,  			        // inner radius
    world_sizeXY,  	  // outer radius
    0.5*world_sizeZ,  // z half length
    0.,			          // starting phi
    360.*deg  	    	// segment angle
  );

  fLogicWorld = new G4LogicalVolume(
    solidWorld,       // its solid
    vacuum_material,  // its material
    "World"           // its name
  );            

  G4VPhysicalVolume* physWorld = new G4PVPlacement(
    0,                // no rotation
    G4ThreeVector(),  // at (0,0,0)
    fLogicWorld,      // its logical volume
    "World",          // its name
    0,                // its mother  volume
    false,            // no boolean operation
    0,                // copy number
    checkOverlaps     // overlaps checking
  );

  /*
  msisAtmosTable columns:
  
  [0]  alt [km] 
  [1]  O   * 
  [2]  N2  *
  [3]  O2  *
  [4]  Total mass density [kg/m^3] 
  [5]  Temp [K] 
  [6]  He  *
  [7]  Ar  *
  [8]  H   *
  [9]  N   *
  [10] H2  *

  all species' mass density in [kg/m^3]
  */

  fTableSize = GetMSIStableSize(fAtmosphereFilename);

  // Cast to const for table instantiation
  unsigned const int tableSize = fTableSize;
  G4double msisAtmosTable[tableSize][11];
 
  // Populate array with MSIS atmosphere table
  GetMSIStable(msisAtmosTable, fAtmosphereFilename, tableSize);

  // Atmospheric material definitions
  G4Element* O  = new G4Element("Oxygen",   "O",  8.0,  16.0   * g/mole);
  G4Element* N  = new G4Element("Nitrogen", "N",  7.0,  14.0   * g/mole);
  G4Element* Ar = new G4Element("Argon",    "Ar", 18.0, 39.948 * g/mole);
  G4Element* He = new G4Element("Helium",   "He", 2.0,  4.0    * g/mole);
  G4Element* H  = new G4Element("Hydrogen", "H",  1.0,  1.0078 * g/mole);

  // Layers are the size needed to fill the 1000 km column
  G4double layerThickness = (1000.0 / tableSize) * km;
  G4double layerLocation;
  
  G4Tubs* atmosphereLayer = new G4Tubs(
    "AtmosphereLayer",                 // its name
    0.0,                                // inner radius
    world_sizeXY-1.*mm,                // outer radius
    (0.5*layerThickness) - (0.5*um),   // z half length
		0.0,                                // starting phi
		360.0*deg                           // segment angle
  );

  G4Material*      layerMaterial;
  G4LogicalVolume* logicLayer;
  G4double         pressure;
  G4double         R_gas_constant_air = 287.0;  // J/kg-K air
  G4double         totalLayerMassDensity;
  G4double         zeroThreshold = 1e-21; // approximate minimum number density [cm^-3] that Geant will tolerate
  G4int            nComponents;

  for(int i = 0; i < fTableSize; i++)
  {
    // Ideal gas law for atmospheric pressure
    // P [Pa] = R [J/kg-K air] * rho [kg/m^3] * T [K]
    pressure = R_gas_constant_air * msisAtmosTable[i][4] * msisAtmosTable[i][5];
    // TODO: replace R gas constant with layer mean mass! 

    // Material definitions for non-elements
    G4Material* N2;
    G4Material* O2;
    G4Material* H2;

    if(msisAtmosTable[i][2] > zeroThreshold){
      N2 = new G4Material(
        "N2-Layer"+std::to_string(i),
        msisAtmosTable[i][2]*kg/m3,
        1,
        kStateGas,
        msisAtmosTable[i][5]*kelvin,
        pressure*pascal
      );
      N2->AddElement(N, 2);
    }
    if(msisAtmosTable[i][3] > zeroThreshold){
      O2 = new G4Material(
        "O2-Layer"+std::to_string(i),
        msisAtmosTable[i][3]*kg/m3,
        1,
        kStateGas,
        msisAtmosTable[i][4]*kelvin,
        pressure*pascal
      );
      O2->AddElement(O, 2);
    }
    if(msisAtmosTable[i][10] > zeroThreshold){
      H2 = new G4Material(
        "H2-Layer"+std::to_string(i),
        msisAtmosTable[i][10]*kg/m3,
        1,
        kStateGas,
        msisAtmosTable[i][10]*kelvin,
        pressure*pascal
      );
      H2->AddElement(H, 2);
    }

    // Get number of components & total mass density in layer
    std::vector<int> columnsWithDensityValues = {1, 2, 3, 6, 7, 8, 9, 10};
    nComponents = 0;
    totalLayerMassDensity = 0;
    for(int j = 0; j < columnsWithDensityValues.size(); j++){
      int columnIdx = columnsWithDensityValues[j];
      if(msisAtmosTable[i][columnIdx] < zeroThreshold){continue;}

      nComponents++;
      totalLayerMassDensity += msisAtmosTable[i][columnIdx];
    }
    
    // Create material for layer
    layerMaterial = new G4Material(
      "AirLayer"+std::to_string(i), // name
      totalLayerMassDensity*kg/m3,  // density
      nComponents,                  // number of components
      kStateGas,                    // state
      msisAtmosTable[i][4]*kelvin,  // temperature
      pressure*pascal   	          // pressure
    );
        
    // I'd love to replace this with a for loop but it doesn't seem worth trying to figure out how
    // to put both G4Element and G4Material in the same vector. Too bad!
    if(msisAtmosTable[i][1] > zeroThreshold)  // O
      layerMaterial->AddElement(O, msisAtmosTable[i][1]/totalLayerMassDensity);
   
    if(msisAtmosTable[i][2] > zeroThreshold) // N2
      layerMaterial->AddMaterial(N2, msisAtmosTable[i][2]/totalLayerMassDensity);
   
    if(msisAtmosTable[i][3] > zeroThreshold) // O2
      layerMaterial->AddMaterial(O2, msisAtmosTable[i][3]/totalLayerMassDensity);
    
    if(msisAtmosTable[i][6] > zeroThreshold) // He 
      layerMaterial->AddElement(He, msisAtmosTable[i][6]/totalLayerMassDensity);
    
    if(msisAtmosTable[i][7] > zeroThreshold) // Ar 
      layerMaterial->AddElement(Ar, msisAtmosTable[i][7]/totalLayerMassDensity);
   
    if(msisAtmosTable[i][8] > zeroThreshold) // H 
      layerMaterial->AddElement(H, msisAtmosTable[i][8]/totalLayerMassDensity);
  
    if(msisAtmosTable[i][9] > zeroThreshold) // N
      layerMaterial->AddElement(N, msisAtmosTable[i][9]/totalLayerMassDensity);

    if(msisAtmosTable[i][10] > zeroThreshold) // H2
      layerMaterial->AddMaterial(H2, msisAtmosTable[i][10]/totalLayerMassDensity);

    logicLayer = new G4LogicalVolume(atmosphereLayer, layerMaterial, "AtmosphereLayer"+std::to_string(i));
    layerLocation = (i-fTableSize/2)*layerThickness + layerThickness/2.;
    
    // Create layer
    new G4PVPlacement(
      0,                                     // rotation
		  G4ThreeVector(0., 0., layerLocation),  // location
		  logicLayer,                            // its logical volume
		  "AtmosphereLayer"+std::to_string(i),   // its name
		  fLogicWorld,                           // its mother volume
		  false,                                 // no boolean operation
		  i,                                     // copy number
		  checkOverlaps                          // overlaps checking
    );
  }
  return physWorld;
}

void DetectorConstruction::ConstructSDandField()
{
  if(!fEmFieldSetup.Get())
  {    
    F03FieldSetup* emFieldSetup = new F03FieldSetup();    
    fEmFieldSetup.Put(emFieldSetup);
      
    G4AutoDelete::Register(emFieldSetup);
  }
  fLogicWorld->SetFieldManager(fEmFieldSetup.Get()->GetGlobalFieldManager(), true); 
}

void DetectorConstruction::GetMSIStable(G4double tableEntry[][11], G4String filename, unsigned int tableSize)
{
  // Fill data array with data
  std::ifstream file;
  file.open(filename, std::ifstream::in);

  // Parse lines
  int numberOfHeaderLines = 1;
  int dim1Index = -1 * numberOfHeaderLines; // Negative to offset header line
  int dim2Index = 0;

  std::string line;
  std::string token;
  while( std::getline(file, line) ){
    if(dim1Index < 0){dim1Index++; continue;} // Skip header lines

    std::istringstream word(line); // Read line
    while ( std::getline(word, token, ',') ){
      tableEntry[dim1Index][dim2Index] = std::stod(token);
      dim2Index++;
    }
    dim1Index++;
    dim2Index = 0;
  }
  file.close();
}

G4int DetectorConstruction::GetMSIStableSize(G4String filename)
{
  std::ifstream filePtr(filename, std::ios::in);
  G4String line;
  G4int counter = 0;

  while(getline(filePtr, line)){counter++;}
  filePtr.close();

  int numberOfHeaderLines = 1;
  return counter - numberOfHeaderLines;
}