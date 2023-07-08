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


#include "TrackerSD.hh"
#include "TrackerHit.hh"

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



//ok
#include "G4RotationMatrix.hh"
#include <G4Tubs.hh>
#include <G4Types.hh>

#include <fstream>
#include <iostream>
#include <vector>



namespace B1
{


//varsayılan birim metredir


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


G4VPhysicalVolume* DetectorConstruction::Construct()
{	
  //varsayılan birim metredir
  G4double h1 = 1.31*m;
  G4double h2 = 1.50*m;
  G4double R = 0.5*m;
  G4double d =5.0*m;
  //G4cout<< h1 <<G4endl;
  
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();

  // Envelope parameters
  //
  G4double env_sizeZ = 2*(h1 + d + h2), env_sizeY = (2*R), env_sizeX = (6.45*R);
  G4Material* env_mat = nist->FindOrBuildMaterial("G4_Galactic");

  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  //
  // World
  //
  //G4double world_sizeXY = 1.2*env_sizeXY;
  //G4double world_sizeZ  = 1.2*env_sizeZ;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_Galactic");

  auto solidWorld = new G4Box("World",                           // its name
    2.5*env_sizeX, 2.5*env_sizeY, 1*env_sizeZ);  // its size

  auto logicWorld = new G4LogicalVolume(solidWorld,  // its solid
    world_mat,                                       // its material
    "World");                                        // its name

  auto physWorld = new G4PVPlacement(nullptr,  // no rotation
    G4ThreeVector(),                           // at (0,0,0)
    logicWorld,                                // its logical volume
    "World",                                   // its name
    nullptr,                                   // its mother  volume
    false,                                     // no boolean operation
    0,                                         // copy number
    checkOverlaps);                            // overlaps checking

  //
  // Envelope
  //
  auto solidEnv = new G4Box("Envelope",                    // its name
    2 * env_sizeX, 2 * env_sizeY,  1*env_sizeZ);  // its size

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
/*
  //
  // Shape 1
  //
  G4Material* shape1_mat = nist->FindOrBuildMaterial("G4_A-150_TISSUE");
  G4ThreeVector pos1 = G4ThreeVector(0, 2*cm, -7*cm);

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
    pos1,                     // at position
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
  G4ThreeVector pos2 = G4ThreeVector(0, -1*cm, 7*cm);

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
    pos2,                     // at position
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
  

// ok
  
  G4Material* scintillator = nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
  G4ThreeVector pos3 = G4ThreeVector(0*m, 0*m, d);
  //G4ThreeVector pos4 = G4ThreeVector(0*m, 0*m, -d/2);
  
  
// det1
  G4VSolid* tubeSolid1 =
new G4Tubs("det1",  //name
0* cm, // inner radius
R, // outer radius
h1, // height
0.0 * deg, 360.0 * deg); // segment angles
  
  
  auto logicShape3 = new G4LogicalVolume(tubeSolid1,  // its solid
    scintillator,                                        // its material
    "LVDet1");                                         // its name


  
  new G4PVPlacement(nullptr,  // no rotation
    pos3,                     // at position
    logicShape3,              // its logical volume
    "PVDet1",                 // its name
    logicEnv,                 // its mother  volume
    false,                    // no boolean operation
    0,                        // copy number
    checkOverlaps);           // overlaps checking

  //fScoringVolume = logicShape3;
  
  
//det2



  G4VSolid* tubeSolid2 =
new G4Tubs("det2",  //name
0* m, // inner radius
R, // outer radius
h2, // height
0.0 * deg, 360.0 * deg); // segment angles

auto logicShape4 = new G4LogicalVolume(tubeSolid2,  // its solid
    scintillator,                                        // its material
    "LVDet2");     
    
  G4int detnumber = 1;
  for(G4int icount = 0; icount < detnumber; icount++)
  {
    
    G4ThreeVector pos4 = G4ThreeVector(d*std::cos((2 * 3.1415/detnumber) * icount) , 0*m, d*std::sin((2 * 3.1415/detnumber) * icount) );
    
    G4RotationMatrix* rotationMatrix = new G4RotationMatrix();
    
    rotationMatrix->rotateY(CLHEP::pi/2. + icount * 2*3.1415/detnumber);
    
    //G4Transform3D* transform = new G4Transform3D(rotationMatrix,pos4);
  
  
    new G4PVPlacement(rotationMatrix,
      pos4,
      logicShape4,              // its logical volume
      "PVDet2",                 // its name
      logicEnv,                 // its mother  volume
      false,                    // no boolean operation
      icount,                        // copy number
      checkOverlaps);           // overlaps checking
      
      
      
      G4String trackerChamberSD2name = "/TrackerSD2_copy_no_" + std::to_string(icount);
      auto aTrackerSD2 = new TrackerSD(trackerChamberSD2name, "DetCollection2");
      G4SDManager::GetSDMpointer()->AddNewDetector(aTrackerSD2);
      SetSensitiveDetector(logicShape4->GetName(), aTrackerSD2, true);
     
      
      }

  fScoringVolume = logicShape4;
  
  //GUN VÜCUDU. 
  //Bu vücudu parçacıklar cidden çembersel bir şekilde çıkıyor mu anlamak için koyuyorum sınırlardan taşıyorlarsa demek ki bir terslik var.
  G4double innerRadius = 0.0*m; // Inner radius of the circle

  G4double outerRadius = 0.2*m; // Outer radius of the circle
  G4double height = 1.0*mm; // Height of the cylinder
  G4double startAngle = 0.0 * deg; // Start angle of the circle
  G4double spanningAngle = 360.0*deg; // Spanning angle of the circle

  G4Sphere* gunVolume = new G4Sphere("GunVolume", innerRadius, outerRadius, 0, 2 * M_PI, 0, M_PI);

  // Create logical volume
  G4LogicalVolume* gunLogical = new G4LogicalVolume(gunVolume, world_mat,"GunLogical");
  
  G4VisAttributes* visAttributes = new G4VisAttributes(G4Colour(0, 14.0, 27.0));
  
  visAttributes->SetVisibility(true);
  
  gunLogical->SetVisAttributes(visAttributes);
  


  new G4PVPlacement(0, G4ThreeVector(0,0,0),   gunLogical, "GunPhysical", logicEnv,false,checkOverlaps);

 //fScoringVolume = gunLogical;
  
  return physWorld;
}

void DetectorConstruction::ConstructSDandField()
{

  // Sensitive detectors

  /*
  G4String trackerChamberSD2name = "/TrackerSD2";
  auto aTrackerSD2 = new TrackerSD(trackerChamberSD2name, "DetCollection2");
  G4SDManager::GetSDMpointer()->AddNewDetector(aTrackerSD2);
  SetSensitiveDetector("LVDet2", aTrackerSD2, true);
  */
  
  G4String trackerChamberSD1name = "/TrackerSD1";
  auto aTrackerSD1 = new TrackerSD(trackerChamberSD1name, "DetCollection1");
  G4SDManager::GetSDMpointer()->AddNewDetector(aTrackerSD1);
  SetSensitiveDetector("LVDet1", aTrackerSD1, true);


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// yapmak istediğim şeylerin verisi olan bir makale bul ve aynı sonucu üretmeye çalış
//https://www.iap.kit.edu/corsika/70.php
}

