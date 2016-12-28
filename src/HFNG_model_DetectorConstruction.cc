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
// $Id: HFNG_model_DetectorConstruction.cc 94307 2015-11-11 13:42:46Z gcosmo $
//
/// \file HFNG_model_DetectorConstruction.cc
/// \brief Implementation of the HFNG_model_DetectorConstruction class

#include "HFNG_model_DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4GenericTrap.hh"

#include "G4SubtractionSolid.hh"

#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HFNG_model_DetectorConstruction::HFNG_model_DetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HFNG_model_DetectorConstruction::~HFNG_model_DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* HFNG_model_DetectorConstruction::Construct()
{  
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
  
  // Envelope parameters
  //
  G4double env_sizeXY = 10*m, env_sizeZ = 10*m;
  G4Material* env_mat = nist->FindOrBuildMaterial("G4_AIR");
   
  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  //     
  // World
  //
  G4double world_sizeXY = 1.2*env_sizeXY;
  G4double world_sizeZ  = 1.2*env_sizeZ;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");
  
  G4Box* solidWorld =    
    new G4Box("World",                       //its name
       0.8*world_sizeXY, 0.8*world_sizeXY, 0.8*world_sizeZ);     //its size
      
  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "World");            //its name
                                   
  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking
                     
  //     
  // Envelope
  //  
  G4Box* solidEnv =    
    new G4Box("Envelope",                    //its name
        0.5*env_sizeXY, 0.5*env_sizeXY, 0.5*env_sizeZ); //its size
      
  G4LogicalVolume* logicEnv =                         
    new G4LogicalVolume(solidEnv,            //its solid
                        env_mat,             //its material
                        "Envelope");         //its name
               
  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),         //at (0,0,0)
                    logicEnv,                //its logical volume
                    "Envelope",              //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
 
  //
  // Shrould
  //

  // The shrould is built by forming a cylinder, subtracting out two inner
  // cylinders and inner block, and subtracting off blocks on each side.
  // Subtraction is done as a few boolean operations in series. Note: The
  // solids subtracted away have a larger z dimension to them than the
  // shrould they are subtracting away from.

  // Assume an aluminum material for the shrould.
  G4Material* shrould_mat = nist->FindOrBuildMaterial("G4_Al");
  G4ThreeVector pos_shrould = G4ThreeVector();

  // Define the main cylinder solid volume.
  G4Tubs* main_cylinder =
    new G4Tubs("Main Cylinder",
               0, 1*m, 1*m, 0, twopi);       // r:    0 m -> 1 m
                                             // z:   -1 m -> 1 m
                                             // phi:    0 -> 2 pi

  // Define the inner two cylinders using their placement for subtraction.
  // The solid volume only needs to be defined for one and the placement
  // has to be defined for two. Assume each inner cylinder has one-third
  // of the radius of the outer cylinder and the same height.
  G4ThreeVector pos_innercyl_1(-50*cm, 0, 0);
  G4ThreeVector pos_innercyl_2(50*cm, 0, 0);

  G4Tubs* inner_cylinder =
    new G4Tubs("Inner Cylinder",
               0, 33*cm, 2*m, 0, twopi);     // r:    0 m -> 0.66 m
                                             // z:   -1 m -> 1 m
                                             // phi:    0 -> 2 pi

  // Define the inner block and its placement for subtraction.
  G4ThreeVector pos_inner_block(0, 0, 0);

  G4Box* inner_block =
    new G4Box("Inner Block",         // Its name
              50*cm, 33*cm, 2*m);  // Its size

  // Define the two outer blocks and their placement for subtraction.
  G4ThreeVector pos_outer_block_1(0, 1*m, 0);
  G4ThreeVector pos_outer_block_2(0, -1*m, 0);

  G4Box* outer_block =
    new G4Box("Outer Block",        // Its name
              2*m, 33*cm, 2*m);    // Its size


  // Define the shrould subtraction solids via iterative subtraction operations.
  G4SubtractionSolid* shrould_1 =
    new G4SubtractionSolid("Shrould 1",           // Name of the 1st shrould
                           main_cylinder,         // Main cylinder
                           inner_cylinder,        // Inner cylinder 1
                           0,                     // No rotation of inner cylinder 1
                           pos_innercyl_1);       // Position of inner cylinder 1

  G4SubtractionSolid* shrould_2 =
    new G4SubtractionSolid("Shrould 2",           // Name of the 2nd shrould
                           shrould_1,             // 1st shrould
                           inner_cylinder,        // Inner cylinder 2
                           0,                     // No rotation of inner cylinder 2
                           pos_innercyl_2);       // Position of inner cylinder 2

  G4SubtractionSolid* shrould_3 =
    new G4SubtractionSolid("Shrould 3",           // Name of the 3rd shrould
                           shrould_2,             // 2nd shrould
                           inner_block,           // Inner block
                           0,                     // No rotation of the inner block
                           pos_inner_block);      // Position of inner block

  G4SubtractionSolid* shrould_4 =
    new G4SubtractionSolid("Shrould 4",           // Name of the 4th shrould
                           shrould_3,             // 3rd shrould
                           outer_block,           // Outer block 1
                           0,                     // No rotation of outer block 1
                           pos_outer_block_1);    // Position of outer block 1


  G4SubtractionSolid* shrould_5 =
    new G4SubtractionSolid("Shrould 5",           // Name of the 5th shrould
                           shrould_4,             // 4rd shrould
                           outer_block,           // Outer block 2
                           0,                     // No rotation of outer block 2
                           pos_outer_block_2);    // Position of outer block 2

  G4LogicalVolume* logical_shrould =
    new G4LogicalVolume(shrould_5,            // Its solid volume
                        shrould_mat,          // Its material
                        "Logical Shrould");   // Its name

  new G4PVPlacement(0,                    // No rotation
                    pos_shrould,          // At position
                    logical_shrould,      // Its logical volume
                    "Physical Shrould",   // Its name
                    logicEnv,             // Its mother volume
                    true,                 // Boolean operation
                    0,                    // Copy number
                    checkOverlaps);       // Overlaps checking

  // Set the shrould as the scoring volume.
  //
  fScoringVolume = logical_shrould;

  //
  //always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
