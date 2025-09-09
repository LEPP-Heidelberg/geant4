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

#include "ExUCNDetectorConstructionMessenger.hh"
#include "ExUCNDetectorConstruction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

ExUCNDetectorConstructionMessenger::ExUCNDetectorConstructionMessenger(ExUCNDetectorConstruction* Gun)
:fAction(Gun)
//PrimaryGeneratorMessenger::PrimaryGeneratorMessenger()
{
G4cout << "create messenger for the UCN DetectorConstruction class " << G4endl;
  gunDirectory = new G4UIdirectory("/det/");
  gunDirectory->SetGuidance("Geometry control commands.");

  setxGunDirectionCmd1 = new G4UIcmdWithADouble("/det/optpot1",this);
  setxGunDirectionCmd1->SetGuidance("Set optical potential 1");
  setxGunDirectionCmd1->SetParameterName("xv1",true);
  setxGunDirectionCmd1->SetDefaultValue(0.0) ;

  setxGunDirectionCmd2 = new G4UIcmdWithADouble("/det/optpot2",this);
  setxGunDirectionCmd2->SetGuidance("Set optical potential 2");
  setxGunDirectionCmd2->SetParameterName("xv2",true);
  setxGunDirectionCmd2->SetDefaultValue(0.0) ;

  setxGunDirectionCmd3 = new G4UIcmdWithADouble("/det/optpot3",this);
  setxGunDirectionCmd3->SetGuidance("Set optical potential 3");
  setxGunDirectionCmd3->SetParameterName("xv3",true);
  setxGunDirectionCmd3->SetDefaultValue(0.0) ;

    setxGunDirectionCmd4 = new G4UIcmdWithADouble("/det/optpot4",this);
  setxGunDirectionCmd4->SetGuidance("Set optical potential 4");
  setxGunDirectionCmd4->SetParameterName("xv4",true);
  setxGunDirectionCmd4->SetDefaultValue(0.0) ;

    setxGunDirectionCmd5 = new G4UIcmdWithADouble("/det/eta1",this);
  setxGunDirectionCmd5->SetGuidance("Set eta 1");
  setxGunDirectionCmd5->SetParameterName("xv5",true);
  setxGunDirectionCmd5->SetDefaultValue(0.0) ;

    setxGunDirectionCmd6 = new G4UIcmdWithADouble("/det/eta2",this);
  setxGunDirectionCmd6->SetGuidance("Set eta 2");
  setxGunDirectionCmd6->SetParameterName("xv6",true);
  setxGunDirectionCmd6->SetDefaultValue(0.0) ;

    setxGunDirectionCmd7 = new G4UIcmdWithADouble("/det/eta3",this);
  setxGunDirectionCmd7->SetGuidance("Set eta 3");
  setxGunDirectionCmd7->SetParameterName("xv7",true);
  setxGunDirectionCmd7->SetDefaultValue(0.0) ;

    setxGunDirectionCmd8 = new G4UIcmdWithADouble("/det/eta4",this);
  setxGunDirectionCmd8->SetGuidance("Set eta 4");
  setxGunDirectionCmd8->SetParameterName("xv8",true);
  setxGunDirectionCmd8->SetDefaultValue(0.0) ;

    setxGunDirectionCmd9 = new G4UIcmdWithADouble("/det/diff1",this);
  setxGunDirectionCmd9->SetGuidance("Set diff 1");
  setxGunDirectionCmd9->SetParameterName("xv9",true);
  setxGunDirectionCmd9->SetDefaultValue(0.0) ;

     setxGunDirectionCmd10 = new G4UIcmdWithADouble("/det/diff2",this);
  setxGunDirectionCmd10->SetGuidance("Set diff 2");
  setxGunDirectionCmd10->SetParameterName("xv10",true);
  setxGunDirectionCmd10->SetDefaultValue(0.0) ;

     setxGunDirectionCmd11 = new G4UIcmdWithADouble("/det/diff3",this);
  setxGunDirectionCmd11->SetGuidance("Set diff 3");
  setxGunDirectionCmd11->SetParameterName("xv11",true);
  setxGunDirectionCmd11->SetDefaultValue(0.0) ;

     setxGunDirectionCmd12 = new G4UIcmdWithADouble("/det/diff4",this);
  setxGunDirectionCmd12->SetGuidance("Set diff 4");
  setxGunDirectionCmd12->SetParameterName("xv12",true);
  setxGunDirectionCmd12->SetDefaultValue(0.0) ;

     setxGunDirectionCmd13 = new G4UIcmdWithADouble("/det/valveheight",this);
  setxGunDirectionCmd13->SetGuidance("Set height of the valve in mm");
  setxGunDirectionCmd13->SetParameterName("xv13",true);
  setxGunDirectionCmd13->SetDefaultValue(0.0) ;

  setxGunDirectionCmd14 = new G4UIcmdWithADouble("/det/stepsize",this);
  setxGunDirectionCmd14->SetGuidance("Set step size limit in mm in mm");
  setxGunDirectionCmd14->SetParameterName("xv14",true);
  setxGunDirectionCmd14->SetDefaultValue(0.0) ;



  //optpot1  
  /*
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
*/
G4cout << "generated detector messenger class " << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExUCNDetectorConstructionMessenger::~ExUCNDetectorConstructionMessenger()
{

  delete gunDirectory;
  delete setxGunDirectionCmd1;  
  delete setxGunDirectionCmd2;
  delete setxGunDirectionCmd3;
  delete setxGunDirectionCmd4;
  delete setxGunDirectionCmd5;
  delete setxGunDirectionCmd6;
  delete setxGunDirectionCmd7;
  delete setxGunDirectionCmd8;
  delete setxGunDirectionCmd9;
  delete setxGunDirectionCmd10;
  delete setxGunDirectionCmd11;
  delete setxGunDirectionCmd12;
  delete setxGunDirectionCmd13;
  delete setxGunDirectionCmd14;

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void ExUCNDetectorConstructionMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{

  if( command == setxGunDirectionCmd1)
   { fAction->SetProperty1(setxGunDirectionCmd1->GetNewDoubleValue(newValue));}
  if( command == setxGunDirectionCmd2)
   { fAction->SetProperty2(setxGunDirectionCmd2->GetNewDoubleValue(newValue));}
  if( command == setxGunDirectionCmd3)
   { fAction->SetProperty3(setxGunDirectionCmd3->GetNewDoubleValue(newValue));}
  if( command == setxGunDirectionCmd4)
   { fAction->SetProperty4(setxGunDirectionCmd4->GetNewDoubleValue(newValue));}
  if( command == setxGunDirectionCmd5)
   { fAction->SetProperty5(setxGunDirectionCmd5->GetNewDoubleValue(newValue));}
  if( command == setxGunDirectionCmd6)
   { fAction->SetProperty6(setxGunDirectionCmd6->GetNewDoubleValue(newValue));}
  if( command == setxGunDirectionCmd7)
   { fAction->SetProperty7(setxGunDirectionCmd7->GetNewDoubleValue(newValue));}
  if( command == setxGunDirectionCmd8)
   { fAction->SetProperty8(setxGunDirectionCmd8->GetNewDoubleValue(newValue));}
  if( command == setxGunDirectionCmd9)
   { fAction->SetProperty9(setxGunDirectionCmd9->GetNewDoubleValue(newValue));}
  if( command == setxGunDirectionCmd10)
   { fAction->SetProperty10(setxGunDirectionCmd10->GetNewDoubleValue(newValue));}
  if( command == setxGunDirectionCmd11)
   { fAction->SetProperty11(setxGunDirectionCmd11->GetNewDoubleValue(newValue));}
  if( command == setxGunDirectionCmd12)
   { fAction->SetProperty12(setxGunDirectionCmd12->GetNewDoubleValue(newValue));}
  if( command == setxGunDirectionCmd13)
   { fAction->SetProperty13(setxGunDirectionCmd13->GetNewDoubleValue(newValue));}
  if( command == setxGunDirectionCmd14)
   { fAction->SetProperty14(setxGunDirectionCmd14->GetNewDoubleValue(newValue));}


}

