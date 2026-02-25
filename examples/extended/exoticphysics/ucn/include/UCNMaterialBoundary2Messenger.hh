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

#ifndef UCNMaterialBoundary2Messenger_h
#define UCNMaterialBoundary2Messenger_h 1

#include "G4UImessenger.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UImessenger.hh"
#include "UCNMaterialBoundary2.hh"

class UCNMaterialBoundary2;

///////
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADouble;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWith3Vector;
class G4UIcmdWith3VectorAndUnit;



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class UCNMaterialBoundary2Messenger: public G4UImessenger
{
  public:

   UCNMaterialBoundary2Messenger(UCNMaterialBoundary2*);
//     PrimaryGeneratorMessenger();
   ~UCNMaterialBoundary2Messenger();

    void SetNewValue(G4UIcommand*, G4String);

  private:

    UCNMaterialBoundary2* fAction;

   //G4string especfile;
    G4UIcmdWith3VectorAndUnit* setGunPositionCmd;
    G4UIcmdWithADouble* setxGunDirectionCmd;
    G4UIcmdWithADouble* setyGunDirectionCmd;
    G4UIcmdWithADouble* setzGunDirectionCmd;
    G4UIcmdWithADouble* setGunEnergy_neVCmd;
    G4UIcmdWithADouble* setGunDirectionRangeCmd;
    G4UIcmdWithADoubleAndUnit* setGunMaxRadCmd;
    G4UIdirectory*     gunDirectory;
    G4UIcmdWithAString* setEspecCmd;
    G4UIcmdWithAString* setAngspecCmd;
    
    G4UIcmdWithAString* setAngspeclogCmd;
    
        G4UIcmdWithAString* setFileName1Cmd;    
	G4UIcmdWithAString* setFileName2Cmd;
	G4UIcmdWithAString* setFileName3Cmd;
         G4UIcmdWithADouble* setstepsCmd;
         G4UIcmdWithADouble* setlifeCmd;
	 G4UIcmdWithADouble* sethelifeCmd;
         G4UIcmdWithADouble* setfoilCmd; 
         G4UIcmdWithADouble* settimeCmd;
	 G4UIdirectory*     matDirectory;
    //G4UIcmdWithAString*
         G4UIcmdWithADouble* setDiff1Cmd;
         G4UIcmdWithADouble* setDiff2Cmd;
         G4UIcmdWithADouble* setDiff3Cmd;
         G4UIcmdWithADouble* setDiff4Cmd;
         G4UIcmdWithADouble* setDiff5Cmd;
          G4UIcmdWithADouble* setOpt1Cmd;
         G4UIcmdWithADouble* setOpt2Cmd;
         G4UIcmdWithADouble* setOpt3Cmd;
         G4UIcmdWithADouble* setOpt4Cmd;
         G4UIcmdWithADouble* setOpt5Cmd;
         G4UIcmdWithADouble* setEta1Cmd;
         G4UIcmdWithADouble* setEta2Cmd;
         G4UIcmdWithADouble* setEta3Cmd;
         G4UIcmdWithADouble* setEta4Cmd;
         G4UIcmdWithADouble* setEta5Cmd;
          G4UIcmdWithADouble* setCCmd;
         G4UIcmdWithADouble* setVCmd;
         G4UIcmdWithADouble* setf_aCmd;
         G4UIcmdWithADouble* setf_hCmd;
  	 G4UIcmdWithADouble* setlambdaCmd;
         G4UIcmdWithADouble* setGamma1aCmd;
         G4UIcmdWithADouble* setGamma1hCmd;
         G4UIcmdWithADouble* setGamma2aCmd;
         G4UIcmdWithADouble* setGamma2hCmd;
	 G4UIcmdWithADouble* setGamma3aCmd;
         G4UIcmdWithADouble* setGamma3hCmd;
         G4UIcmdWithADouble* setu_minCmd;
         G4UIcmdWithADouble* setu_maxCmd;
         // Excited States
	 G4UIcmdWithADouble* setexcitedFracCmd;
	 G4UIcmdWithADouble* setmeanDownShiftLifetimeCmd;
	 G4UIcmdWithADouble* setE_recoilCmd;
	 G4UIcmdWithADouble* setmaxGlobalTimeCmd; 
	 // End Excited States
         G4UIcmdWithAString* setChopperFileCmd;
    
};

#endif

