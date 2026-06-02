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
/// \brief Creates the atmosphere for EPP simulation from an atmosphere file

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
#include "ANSIColors.h"
#include "GlobalFunctions.h"

#include <fstream>
#include <regex>

DetectorConstruction::DetectorConstruction():
  G4VUserDetectorConstruction(),
  atmosphereRelPath("none"),
  fDetectorMessenger(),
  fLogicWorld(0)
{
  fDetectorMessenger = new DetectorMessenger(this);
}

DetectorConstruction::~DetectorConstruction()
{
  delete fDetectorMessenger;
}

G4VPhysicalVolume* DetectorConstruction::Construct(){
  // Option to switch on/off checking of volumes overlaps
  G4bool checkOverlaps = false; // Set to true if you need to debug. Set false as there are no issues right now and it's very verbose.
  G4double layerGap = 1.0 * um; // Gap between air layers

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

  // World
  G4double world_sizeXY = 1000.0*km;
  G4double world_sizeZ  = 1000.0*km;

  G4Tubs* solidWorld = new G4Tubs(
    "World",           // its name
    0.0,  			       // inner radius
    world_sizeXY,  	   // outer radius
    0.5 * world_sizeZ, // z half length
    0.0,	             // starting phi
    360.0 *deg  	     // segment angle
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
    0,                // its mother volume
    false,            // no boolean operation
    0,                // copy number
    checkOverlaps     // overlaps checking
  );

  // Get atmosphere data array size
  const int nLayers = getNumberOfAtmosphereLayers(atmosphereRelPath);
  std::vector<G4String> tableColNames = readAtmosphereHeader(atmosphereRelPath);
  G4int nCols = tableColNames.size();

  // Read atmosphere data into vector
  std::vector<std::vector<G4double>> atmosphereData(nLayers, std::vector<G4double>(nCols, -999.0));
  readAtmosphereData(atmosphereData, atmosphereRelPath, nLayers); // Populate array with atmosphere table data

  // Get NIST material database for most accurate material definitions
  G4NistManager* nistManager = G4NistManager::Instance();
  
  // Create all the materials we will need. This also strips the units from the density values in the header vector for easier lookup later
  std::vector<G4bool> columnHasDensity(nCols, false);
  for(int i = 0; i < nCols; i++){
    if(tableColNames[i] == "Total (kg/m3)"){continue;}
    
    // Strip units from density fields
    G4String toErase = " (kg/m3)";
    G4int eraseFrom = tableColNames[i].find(toErase);
    if(eraseFrom == -1){continue;} // Skip over non-density fields
    tableColNames[i].erase(eraseFrom, toErase.length()); // Strip units
    columnHasDensity[i] = true;

    // Create material
    createMaterialFromChemicalSymbol(nistManager, tableColNames[i]);
  }

  // Layers are the size needed to fill the 1000 km column
  // TODO dynamic simulation space size
  G4double layerThickness = (1000.0 / static_cast<G4double>(nLayers)) * km;
  G4double layerLocation;
  
  G4Tubs* atmosphereLayer = new G4Tubs(
    "AtmosphereLayer",                   // its name
    0.0,                                 // inner radius
    world_sizeXY-1.*mm,                  // outer radius
    (0.5*layerThickness) - (layerGap/2), // z half length
		0.0,                                 // starting phi
		360.0*deg                            // segment angle
  );

  G4Material*      layerMaterial;
  G4LogicalVolume* logicLayer;
  G4double R = 8.31446261815324; // Universal gas constant. J /(kg mol)
  G4double R_gas_constant_air = 287.0;  // J/kg-K air <- TODO update this for / mol
  G4double zeroThreshold = 1e-21; // approximate minimum number density [cm^-3] that Geant will tolerate

  // Get the column with temperature data
  auto temperatureColIdx = std::find(tableColNames.begin(), tableColNames.end(), "Neutral Temp. (K)");
  if(temperatureColIdx == tableColNames.end()){
    G4cout << ANSI_RED <<
      __FILE__ << ": " << __FUNCTION__ << "\n" <<
      "ERROR: Atmosphere file does not contain column with label `Neutral Temp. (K)`" <<
    ANSI_NOCOLOR << G4endl;
    throw;
  }
  G4int temperatureIdxAsInt = std::distance(tableColNames.begin(), temperatureColIdx);

  // Iterate over every layer in the atmosphere table
  for(int i = 0; i < nLayers; i++){
    // Get number of components & total mass density in layer
    G4double totalLayerMassDensity = 0;
    G4int nComponents = 0;
    for(int j = 0; j < nCols; j++){
      if(columnHasDensity[j] == false){continue;} // Skip non-density fields

      G4double constituentDensity = atmosphereData[i][j];
      if(constituentDensity < zeroThreshold){continue;}

      totalLayerMassDensity += constituentDensity;
      nComponents++;
    }

    // Ideal gas law for layer pressure
    // P [Pa] = R [J/kg-K air] * rho [kg/m^3] * T [K]
    G4double temperature = atmosphereData[i][temperatureIdxAsInt];
    G4double pressure = R_gas_constant_air * totalLayerMassDensity * temperature; // TODO specific gas constant https://en.wikipedia.org/wiki/Gas_constant#Specific_gas_constant
    
    // Create material for layer
    layerMaterial = new G4Material(
      "LayerMaterial_"+std::to_string(i), // name
      totalLayerMassDensity*kg/m3,  // density
      nComponents,                  // number of components
      kStateGas,                    // state
      temperature * kelvin,         // temperature
      pressure * pascal   	        // pressure
    );

    // Loop over each molecule in the atmosphere composition table and add it to the layer material (if needed)
    for(int j = 0; j < nCols; j++){
      if(columnHasDensity[j] == false){continue;} // Skip non-density fields

      G4double constituentDensity = atmosphereData[i][j];
      if(constituentDensity < zeroThreshold){continue;}
      G4String chemSymbol = tableColNames[j];

      G4Material* materialToAdd = nistManager->BuildMaterialWithNewDensity(
        chemSymbol + "_layer" + std::to_string(i), // Name attached to the material
        chemSymbol,                                // Name of base material to modify. See all materials with manager->ListMaterials("all");
        constituentDensity * kg/m3,                // Density
        temperature * kelvin,                      // Temperature
        pressure * pascal                          // Pressure
      );
      layerMaterial->AddMaterial(materialToAdd, constituentDensity/totalLayerMassDensity);
    }

    // Create layer
    layerLocation = (i-nLayers/2)*layerThickness + layerThickness/2.0;
    logicLayer = new G4LogicalVolume(
      atmosphereLayer, 
      layerMaterial, 
      "AtmosphereLayer"+std::to_string(i)
    );
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
  if(!fEmFieldSetup.Get()){    
    F03FieldSetup* emFieldSetup = new F03FieldSetup();
    fEmFieldSetup.Put(emFieldSetup);

    G4AutoDelete::Register(emFieldSetup);
  }
  __DEBUG_PING__;
  fLogicWorld->SetFieldManager(fEmFieldSetup.Get()->GetGlobalFieldManager(), true); 
  __DEBUG_PING__;
}

std::vector<G4String> DetectorConstruction::readAtmosphereHeader(G4String path){
  std::ifstream file;
  file.open(path, std::ifstream::in);

  std::string line;
  std::string token;
  std::getline(file, line);

  std::vector<G4String> result;
  std::istringstream word(line); 
  while ( std::getline(word, token, ',') ){    
    result.push_back(token);
  }
  return result;
}

void DetectorConstruction::readAtmosphereData(std::vector<std::vector<G4double>> &atmosphereData, G4String filename, unsigned int nLayers){
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
      atmosphereData.at(dim1Index).at(dim2Index) = std::stod(token);
      dim2Index++;
    }
    dim1Index++;
    dim2Index = 0;
  }
  file.close();
}

G4Material* DetectorConstruction::createMaterialFromChemicalSymbol(G4NistManager* nistManager, G4String chemSymbol){
  // Returns material at STP

  // Get atom groups with their numbers
  std::regex regexSplitAtoms("[A-Z][a-z]*[0-9]*");
  std::vector<G4String> atoms = regexParse(chemSymbol, regexSplitAtoms);

  // Loop through atom groups
  std::vector<G4Element*> elements;
  std::vector<G4int> elementNumbers;
  for(int i = 0; i < atoms.size(); i++){
    // Parse atom
    std::regex captureAtom("[A-Z][a-z]*");
    G4String atomName = regexParse(atoms[i], captureAtom)[0];

    // Parse number
    std::regex captureAtomNumber("[0-9]+");
    std::vector<G4String> nAtomsResult = regexParse(atoms[i], captureAtomNumber);
    G4int nAtoms = nAtomsResult.size() == 0 ? 1 : std::stoi(nAtomsResult[0]);

    // Add to vectors
    elements.push_back(nistManager->FindOrBuildElement(atomName));
    elementNumbers.push_back(nAtoms);
  }

  // Create composite material
  G4int nComponents = elements.size();
  G4Material* result = new G4Material(
    chemSymbol,      // Material name (must be unique)
    1e-1 * kg/m3,    // Density
    nComponents  // Number of components in material
  );
  for(int i = 0; i < nComponents; i++){
    result->AddElement(elements[i], elementNumbers[i]);
  }

  return result;
}

std::vector<G4String> DetectorConstruction::regexParse(G4String s, std::regex r){
  std::vector<G4String> result;
  for(
    std::sregex_iterator i = std::sregex_iterator(s.begin(), s.end(), r);
    i != std::sregex_iterator();
    ++i
  ){
    std::smatch match = *i;
    result.push_back(match.str());
  }
  return result;
}

G4int DetectorConstruction::getNumberOfAtmosphereLayers(G4String filename)
{
  std::ifstream filePtr(filename, std::ios::in);

  if(filePtr.is_open() == false){
    G4cout << ANSI_RED <<
      __FILE__ << ": " << __FUNCTION__ << "\n" <<
      "ERROR: Atmosphere file \"" <<
      filename <<
      "\" could not be opened. Are you running the executable from the build directory? This function uses relative paths for platform agnosticity."
    ANSI_NOCOLOR << G4endl;
    throw;
  }

  G4String line;
  G4int counter = 0;

  while(getline(filePtr, line)){counter++;}
  filePtr.close();

  G4int numberOfHeaderLines = 1;
  return counter - numberOfHeaderLines;
}