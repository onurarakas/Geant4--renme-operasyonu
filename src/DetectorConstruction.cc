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
//
/// \file B1/src/DetectorConstruction.cc
/// \brief Implementation of the B1::DetectorConstruction class

#include "DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4SDManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4SystemOfUnits.hh"



//ok
#include "G4RotationMatrix.hh"
#include "G4Tubs.hh"
#include "G4Types.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"


#include <fstream>
#include <iostream>
#include <vector>

namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
   

  // Envelope parameters
  //
  G4double env_sizeXY = 1*m, env_sizeZ = 1*m;
  G4Material* env_mat = nist->FindOrBuildMaterial("G4_Galactic");

  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;
  
  //
  // World
  //
  G4double world_sizeXY = env_sizeXY;
  G4double world_sizeZ  = env_sizeZ;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_Galactic");
  

  auto solidWorld = new G4Box("World",                           // its name
     world_sizeXY,  world_sizeXY,  world_sizeZ);  // its size

  auto logicWorld = new G4LogicalVolume(solidWorld,  // its solid
    world_mat,                                       // its material
    "World");   
    
    G4VisAttributes* world_renk = new G4VisAttributes(renk);
  
    scint1_renk->SetVisibility(true);                                     // its name

  auto physWorld = new G4PVPlacement(nullptr,  // no rotation
    G4ThreeVector(),                           // at (0,0,0)
    logicWorld,                                // its logical volume
    "World",                                   // its name
    nullptr,                                   // its mother  volume
    false,                                     // no boolean operation
    0,                                         // copy number
    checkOverlaps);                            // overlaps checking


  // Envelope
  //
  auto solidEnv = new G4Box("Envelope",                    // its name
    env_sizeXY, env_sizeXY, env_sizeZ);  // its size

  auto logicEnv = new G4LogicalVolume(solidEnv,  // its solid
    env_mat,                                     // its material
    "Envelope");                                 // its name

  new G4PVPlacement(nullptr,  // no rotation
    G4ThreeVector(),          // at (0,0,0)
    logicEnv,                 // its logical volume
    "Envelope",               // its name
    logicWorld,               // its mother  volume
    false,                    // no boolean operation
    0,                        // copy number
    checkOverlaps);           // overlaps checking
 
  //
  // Shape 1
  //
  /*
  G4Material* shape1_mat = nist->FindOrBuildMaterial("G4_A-150_TISSUE");
  G4ThreeVector pos11 = G4ThreeVector(0, 2*cm, -7*cm);

  // Conical section shape
  G4double shape1_rmina =  0.*cm, shape1_rmaxa = 2.*cm;
  G4double shape1_rminb =  0.*cm, shape1_rmaxb = 4.*cm;
  G4double shape1_hz = 3.*cm;
  G4double shape1_phimin = 0.*deg, shape1_phimax = 360.*deg;
  auto solidShape1 = new G4Cons("Shape1", shape1_rmina, shape1_rmaxa, shape1_rminb, shape1_rmaxb,
    shape1_hz, shape1_phimin, shape1_phimax);

  auto logicShape1 = new G4LogicalVolume(solidShape1,  // its solid
    shape1_mat,                                        // its material
    "Shape1");                                         // its name

  new G4PVPlacement(nullptr,  // no rotation
    pos11,                     // at position
    logicShape1,              // its logical volume
    "Shape1",                 // its name
    logicEnv,                 // its mother  volume
    false,                    // no boolean operation
    0,                        // copy number
    checkOverlaps);           // overlaps checking

  //
  // Shape 2
  //
  G4Material* shape2_mat = nist->FindOrBuildMaterial("G4_BONE_COMPACT_ICRU");
  G4ThreeVector pos22 = G4ThreeVector(0, -1*cm, 7*cm);

  // Trapezoid shape
  G4double shape2_dxa = 12*cm, shape2_dxb = 12*cm;
  G4double shape2_dya = 10*cm, shape2_dyb = 16*cm;
  G4double shape2_dz  = 6*cm;
  auto solidShape2 = new G4Trd("Shape2",  // its name
    0.5 * shape2_dxa, 0.5 * shape2_dxb, 0.5 * shape2_dya, 0.5 * shape2_dyb,
    0.5 * shape2_dz);  // its size

  auto logicShape2 = new G4LogicalVolume(solidShape2,  // its solid
    shape2_mat,                                        // its material
    "Shape2");                                         // its name

  new G4PVPlacement(nullptr,  // no rotation
    pos22,                     // at position
    logicShape2,              // its logical volume
    "Shape2",                 // its name
    logicEnv,                 // its mother  volume
    false,                    // no boolean operation
    0,                        // copy number
    checkOverlaps);           // overlaps checking

  // Set Shape2 as scoring volume
  //
  fScoringVolume = logicShape2;

  //
  //always return the physical World
  //
  
  */
  
  G4Material* scintillator = nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
  
  G4Material* mylar = nist->FindOrBuildMaterial("G4_MYLAR");
  
  G4Material* Shield_mat = nist->FindOrBuildMaterial("G4_Pb");
  
  G4ThreeVector pos1 = G4ThreeVector(0*m, 0*m, 40*cm);
  G4ThreeVector pos2 = G4ThreeVector(0*m, 0*m, 70*cm);
  G4ThreeVector pos3 = G4ThreeVector(0*m, 0*m, 100*cm);
  G4ThreeVector pos_pb = G4ThreeVector(0*m, 0*m, 90 *cm);
  
  
//rotasyon
  G4RotationMatrix* rotationMatrix = new G4RotationMatrix();
    
  rotationMatrix->rotateY(CLHEP::pi/2. );
  
  
  G4int y1 = 20 * mm;
  
  G4int y2 = 100 * mm;
  
  G4int y3 = 20 * mm; 
  
  G4int width = 300 * mm;
  
  G4int z = 500 * mm;
  
  G4cout << "7" <<G4endl;
  
  //1.sint başı
  
  G4VSolid* scint1 = new G4Box("det1",  //name
    y1,
    width,
    z); // segment angles
    
  //scint 1 için kaplama solid
    
  G4VSolid* scint1_kaplama_1 = new G4Box("scint1_kaplama_kalın",  //name
    y1 + 0.5*mm,
    width + 0.5*mm,
    z + 0.5*mm); 
    
  G4VSolid* scint1_kaplama_2 = new G4Box("scint1_kaplama_ince",  //name
    y1,
    width,
    z); // segment angles
    
  G4VSolid * scint1_kaplama = new G4SubtractionSolid("scint1_kaplama_kalın - scint1_kaplama_ince", scint1_kaplama_1, scint1_kaplama_2);
  
  G4cout << "8" <<G4endl;
  
  //scint1 kaplamasının logic
  
   G4LogicalVolume* logic_scint1_kaplama = new G4LogicalVolume(scint1_kaplama, mylar, "logic_scint1_kaplama");
   
   
  G4Colour renk1(0, 255, 0,125);

  G4VisAttributes* kaplama1_renk = new G4VisAttributes(renk1);
  
  kaplama1_renk->SetVisibility(true);
  
  logic_scint1_kaplama->SetVisAttributes(kaplama1_renk);
   
    new G4PVPlacement(rotationMatrix,  // no rotation
    pos1,                     // at position
    logic_scint1_kaplama,              // its logical volume
    "PVscint1_kaplama",                 // its name
    logicEnv,                 // its mother  volume
    false,                    // no boolean operation
    0,                        // copy number
    checkOverlaps);
  
  //scint1 kaplamasının logic sonu

  
    
   G4LogicalVolume* logic_scint1 = new G4LogicalVolume(scint1, scintillator, "logic_scint1");
   
    G4Colour renk(255, 0, 0, 125);

    G4VisAttributes* scint1_renk = new G4VisAttributes(renk);
  
    scint1_renk->SetVisibility(true);
  
    logic_scint1->SetVisAttributes(scint1_renk);
   
    new G4PVPlacement(rotationMatrix,  // no rotation
    pos1,                     // at position
    logic_scint1,              // its logical volume
    "PVscint1",                 // its name
    logicEnv,                 // its mother  volume
    false,                    // no boolean operation
    0,                        // copy number
    checkOverlaps);           // overlaps checking
    

   
    
 //1.sint sonu 
 
 
 //2.sint başı 
  G4VSolid* scint2 = new G4Box("det2",  //name
    y2,
    width,
    z); // segment angles
    
   G4LogicalVolume* logic_scint2 = new G4LogicalVolume(scint2, scintillator, "logic_scint3");
   
     
    
    new G4PVPlacement(rotationMatrix,  // no rotation
    pos2,                     // at position
    logic_scint2,              // its logical volume
    "PVscint2",                 // its name
    logicEnv,                 // its mother  volume
    false,                    // no boolean operation
    0,                        // copy number
    checkOverlaps);           // overlaps checking
    
    //2.sint sonu
    
    
    
    
    //3.sint başı
    
    G4VSolid* scint3 = new G4Box("det3",  //name
    y3,
    width,
    z); // segment angles
    
   G4LogicalVolume* logic_scint3 = new G4LogicalVolume(scint3, scintillator, "logic_scint3");
    
    new G4PVPlacement(rotationMatrix,  // no rotation
    pos3,                     // at position
    logic_scint3,              // its logical volume
    "PVscint3",                 // its name
    logicEnv,                 // its mother  volume
    false,                    // no boolean operation
    0,                        // copy number
    checkOverlaps);           // overlaps checking
    
 //3.sint sonu
 
 //kurşun plaka başı 
 
    G4VSolid* kursun = new G4Box("kursun",  //name
    y3,
    width,
    z); // segment angles
    
   G4LogicalVolume* logic_kursun = new G4LogicalVolume(kursun, Shield_mat, "logic_kursun");
   
    G4Colour sa(0.7, 40, 0.1);
   
    G4VisAttributes* visAttributes = new G4VisAttributes(sa);
  
    visAttributes->SetVisibility(true);
  
    logic_kursun->SetVisAttributes(visAttributes);
     
     
    
    new G4PVPlacement(rotationMatrix,  // no rotation
    pos_pb,                     // at position
    logic_kursun,              // its logical volume
    "PVscint3",                 // its name
    logicEnv,                 // its mother  volume
    false,                    // no boolean operation
    0,                        // copy number
    checkOverlaps);           // overlaps checking
 
 //son
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
