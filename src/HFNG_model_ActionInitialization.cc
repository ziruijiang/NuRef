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
// $Id: HFNG_model_ActionInitialization.cc 68058 2013-03-13 14:47:43Z gcosmo $
//
/// \file HFNG_model_ActionInitialization.cc
/// \brief Implementation of the HFNG_model_ActionInitialization class

#include "HFNG_model_ActionInitialization.hh"
#include "HFNG_model_PrimaryGeneratorAction.hh"
#include "HFNG_model_RunAction.hh"
#include "HFNG_model_EventAction.hh"
#include "HFNG_model_SteppingAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HFNG_model_ActionInitialization::HFNG_model_ActionInitialization()
 : G4VUserActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HFNG_model_ActionInitialization::~HFNG_model_ActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HFNG_model_ActionInitialization::BuildForMaster() const
{
  HFNG_model_RunAction* runAction = new HFNG_model_RunAction;
  SetUserAction(runAction);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HFNG_model_ActionInitialization::Build() const
{
  SetUserAction(new HFNG_model_PrimaryGeneratorAction);

  HFNG_model_RunAction* runAction = new HFNG_model_RunAction;
  SetUserAction(runAction);
  
  HFNG_model_EventAction* eventAction = new HFNG_model_EventAction(runAction);
  SetUserAction(eventAction);
  
  SetUserAction(new HFNG_model_SteppingAction(eventAction));
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
