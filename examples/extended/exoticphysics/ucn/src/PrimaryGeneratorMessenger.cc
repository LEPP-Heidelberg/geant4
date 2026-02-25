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
// Please cite the following paper if you use this software
// Nucl.Instrum.Meth.B260:20-27, 2007

#include "PrimaryGeneratorMessenger.hh"
#include "ExUCNPrimaryGeneratorAction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

PrimaryGeneratorMessenger::PrimaryGeneratorMessenger(ExUCNPrimaryGeneratorAction* Gun)
:fAction(Gun)
//PrimaryGeneratorMessenger::PrimaryGeneratorMessenger()
{ 
  //fGunDir = new G4UIdirectory("/beam/");i
  //fGunDir->SetGuidance("beam control");

  /*
   * fEmissionCmd = new G4UIcmdWithAnInteger("/beam/emission",this);
  fEmissionCmd->SetGuidance("Select emission model.");
  fEmissionCmd->SetParameterName("emission",false);
  fEmissionCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
*/

  gunDirectory = new G4UIdirectory("/gun/");
  gunDirectory->SetGuidance("Particle Gun control commands.");

  setGunPositionCmd = new G4UIcmdWith3VectorAndUnit("/gun/gunPosition",this);
  setGunPositionCmd->SetGuidance(" Set coord. of the gun position.");
  setGunPositionCmd->SetParameterName("X","Y","Z",true,true);
  setGunPositionCmd->SetDefaultUnit("mm");
  //setGunPositionCmd->SetDefaultValue(0.0*mm, 0.0*mm, 0.0*mm) ; 

  setxGunDirectionCmd = new G4UIcmdWithADouble("/gun/xGunDirection",this);
  setxGunDirectionCmd->SetGuidance(" Set x coord. of the gun direction.");
  setxGunDirectionCmd->SetParameterName("xv",true);
  setxGunDirectionCmd->SetDefaultValue(0.0) ;

  setyGunDirectionCmd = new G4UIcmdWithADouble("/gun/yGunDirection",this);
  setyGunDirectionCmd->SetGuidance(" Set y coord. of the gun direction.");
  setyGunDirectionCmd->SetParameterName("yv",true);
  setyGunDirectionCmd->SetDefaultValue(0.0) ;

  setzGunDirectionCmd = new G4UIcmdWithADouble("/gun/zGunDirection",this);
  setzGunDirectionCmd->SetGuidance(" Set z coord. of the gun direction.");
  setzGunDirectionCmd->SetParameterName("zv",true);
  setzGunDirectionCmd->SetDefaultValue(0.0) ;

  setGunEnergy_neVCmd = new G4UIcmdWithADouble("/gun/gunEnergyneV",this);
  setGunEnergy_neVCmd->SetGuidance(" Set energy of the gun in neV");
  setGunEnergy_neVCmd->SetParameterName("e",true);
  setGunEnergy_neVCmd->SetDefaultValue(0.0) ;

  setGunDirectionRangeCmd = new G4UIcmdWithADouble("/gun/gunDirectionRange",this);
  setGunDirectionRangeCmd->SetGuidance(" Set gun direction range (centered around the gun direction)");
  setGunDirectionRangeCmd->SetParameterName("xv",true);
  setGunDirectionRangeCmd->SetDefaultValue(0.0) ;

  setGunMaxRadCmd = new G4UIcmdWithADoubleAndUnit("/gun/gunMaxRad",this);
  setGunMaxRadCmd->SetGuidance(" Set radius of the gun.");
  setGunMaxRadCmd->SetParameterName("zv",true);
  //setGunMaxRadCmd->SetDefaultValue(0.0*mm) ;

  setEspecCmd = new G4UIcmdWithAString("/gun/spectrumfile",this);
  setEspecCmd->SetGuidance(" Set name of the energy spectrum file");
  setEspecCmd->SetParameterName("spectrumfile",true);
  setEspecCmd->SetDefaultValue("spec.dat") ;

    setAngspecCmd = new G4UIcmdWithAString("/gun/angspectrumfile",this);
  setAngspecCmd->SetGuidance(" Set name of the angular spectrum file");
  setAngspecCmd->SetParameterName("angspectrumfile",true);

    setAngspeclogCmd = new G4UIcmdWithAString("/gun/angspectrumlog",this);
  setAngspeclogCmd->SetGuidance(" Set name of the angular spectrum log file");
  setAngspeclogCmd->SetParameterName("angspectrumlog",true);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger()
{

  delete gunDirectory;
  delete setGunPositionCmd;
  delete setxGunDirectionCmd;
  delete setyGunDirectionCmd;
  delete setzGunDirectionCmd;
  delete setGunEnergy_neVCmd;
  delete setGunDirectionRangeCmd;
  delete setGunMaxRadCmd;
  delete setEspecCmd;
  delete setAngspecCmd;
  delete setAngspeclogCmd;


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void PrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{ 
//  if (command == fEmissionCmd)
//  {fAction->SetEmission(fEmissionCmd->GetNewIntValue(newValue));}   

	/*
  if( command == setGunPositionCmd)
   { fAction->SetGunPosition(setGunPositionCmd->GetNew3VectorValue(newValue));}
  if( command == setxGunDirectionCmd)
   { fAction->SetxGunDirection(setxGunDirectionCmd->GetNewDoubleValue(newValue));}
  if( command == setyGunDirectionCmd)
   { fAction->SetyGunDirection(setyGunDirectionCmd->GetNewDoubleValue(newValue));}
  if( command == setzGunDirectionCmd)
   {fAction->SetzGunDirection(setzGunDirectionCmd->GetNewDoubleValue(newValue));}
 */
   if( command == setGunEnergy_neVCmd)
   { fAction->SetGunEnergy_neV(setGunEnergy_neVCmd->GetNewDoubleValue(newValue));
   }
/*
   if( command == setGunDirectionRangeCmd)
   { fAction->SetGunDirectionRange(setGunDirectionRangeCmd->GetNewDoubleValue(newValue));}
  if( command == setGunMaxRadCmd)
   { fAction->SetGunMaxRad(setGunMaxRadCmd->GetNewDoubleValue(newValue));}
*/
   if( command == setEspecCmd)
         //G4cout << " new val " << newValue << G4endl;
    {
	   // SetGunEnergy_neV(G4double e)
	     fAction->SetEspec(G4String(newValue));
    }
  if( command == setAngspecCmd)
         //G4cout << " new val " << newValue << G4endl;
   {// fAction->SetAngspec(G4String(newValue));
   }
  if( command == setAngspeclogCmd)
         //G4cout << " new val " << newValue << G4endl;
   {// fAction->SetAngspeclog(G4String(newValue));
   }
}
