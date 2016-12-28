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
  // Geometry of interest
  //

  // Place SOLID VOLUME code for item of interest here.
  //

  // Various shapes for the solid volume are placed here as example references.
  //
  // NOTE: Each shape requires pointers to the header file for that shape. A
  // call to the header file MUST be made at the beginning of this source file.
  // The calls to the shapes in the examples below have been included.

  // Box shape
  G4double x_half_dimension = 0.5*cm;
  G4double y_half_dimension = 0.5*cm;
  G4double z_half_dimension = 0.5*cm;

  G4Box* solid_volume_of_interest =
    new G4Box("Geometry",             // Its name
              x_half_dimension,       // Half the intended x-dimension
              y_half_dimension,       // Half the intended y-dimension
              z_half_dimension);      // Half the intended z-dimension

  // Cylinderical (tubular) shape
  G4double r_min = 1*cm;
  G4double r_max = 1*cm;

  G4double half_length = 0.5*cm;

  G4double start_angle = 0;
  G4double end_angle = twopi;

  G4Tubs* solid_volume_of_interest =
    new G4Tubs("Geometry",            // Its name
               r_min,                 // Inner radius
               r_max,                 // Outer radius
               half_length,           // Half the intended length
               start_angle,           // Start angle of the cylindrical segment
               end_angle);            // End angle of the cylindrical segment

  // Conical shape
  G4double r_min_at_-z_dim = 0;
  G4double r_max_at_-z_dim = twopi;

  G4double r_min_at_+z_dim = 0;
  G4double r_max_at_+z_dim = twopi;

  G4double z_dim = 0.5*cm;

  G4double start_angle = 0;
  G4double end_angle = twopi;

  G4Cons* solid_volume_of_interest =
    new G4Cons("Geometry",            // Its name
               r_min_at_-z_dim,       // First inner radius
               r_max_at_-z_dim,       // First outer radius
               r_min_at_+z_dim,       // Second inner radius
               r_max_at_+z_dim,       // Second outer radius
               z_dim,                 // Half the intended length
               start_angle,           // Start angle of the conical segment
               end_angle);            // End angle of the conical segment


  // Spherical shape
  G4double r_min = 0;
  G4double r_max = 1*cm;
  G4double phi_min = 0;
  G4double phi_max = twopi;
  G4double theta_min = 0;
  G4double theta_max = twopi;

  G4Sphere* solid_volume_of_interest =
    new G4Sphere("Geometry",          // Its name
                 r_min,               // Inner radius
                 r_max,               // Outer radius
                 phi_min,             // Start phi angle
                 phi_max,             // End phi angle
                 theta_min,           // Start theta angle
                 theta_max);          // End theta angle

  // Triangular prism shape (uses general trapezoid constructor)



  // Place LOGICAL VOLUME code for item of interest here.
  //

  // Place PHYSICAL VOLUME code for item of interest here.
  //

  // Set the shrould as the scoring volume.
  //
  fScoringVolume = logical_shrould;

  //
  // Always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
