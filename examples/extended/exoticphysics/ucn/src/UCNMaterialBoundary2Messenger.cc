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

#include "UCNMaterialBoundary2Messenger.hh"
#include "UCNMaterialBoundary2.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

UCNMaterialBoundary2Messenger::UCNMaterialBoundary2Messenger(UCNMaterialBoundary2* Gun)
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

//  G4cout << "create messenger for the UCN Material Boundary 2  class " << G4endl;
  matDirectory = new G4UIdirectory("/mat/");
  matDirectory->SetGuidance("UCN Material Boundary 2  control commands.");

  setFileName1Cmd = new G4UIcmdWithAString("/mat/outputfile1",this);
  setFileName1Cmd->SetGuidance(" Set name of the output file 1");
  setFileName1Cmd->SetParameterName("outputfile1",true);
  setFileName1Cmd->SetDefaultValue("output1.dat") ;

  setFileName2Cmd = new G4UIcmdWithAString("/mat/outputfile2",this);
  setFileName2Cmd->SetGuidance(" Set name of the output file 2");
  setFileName2Cmd->SetParameterName("outputfile2",true);
  setFileName2Cmd->SetDefaultValue("output2.dat") ;



   setstepsCmd = new G4UIcmdWithADouble("/mat/numofsteps",this);
  setstepsCmd->SetGuidance(" Set the number of steps before we kill the particle");
  setstepsCmd->SetParameterName("steps",true);
  setstepsCmd->SetDefaultValue(10) ;

    setlifeCmd = new G4UIcmdWithADouble("/mat/nlifetime",this);
  setlifeCmd->SetGuidance(" Set the neutron lifetime");
  setlifeCmd->SetParameterName("life",true);
  setlifeCmd->SetDefaultValue(10) ;

     sethelifeCmd = new G4UIcmdWithADouble("/mat/helifetime",this);
  sethelifeCmd->SetGuidance(" Set the neutron life time inside helium");
  sethelifeCmd->SetParameterName("helife",true);
  sethelifeCmd->SetDefaultValue(10) ;


     setfoilCmd = new G4UIcmdWithADouble("/mat/foil",this);
  setfoilCmd->SetGuidance(" Set the vacuum separation foil absorption");
  setfoilCmd->SetParameterName("foil",true);
  setfoilCmd->SetDefaultValue(1) ;

settimeCmd = new G4UIcmdWithADouble("/mat/timesteps",this);
  settimeCmd->SetGuidance(" Use of time steps in ms");
  settimeCmd->SetParameterName("1",true);
  settimeCmd->SetDefaultValue(1) ;

  // optical potentials optpot1, optpot2, optpot3
setOpt1Cmd = new G4UIcmdWithADouble("/mat/optpot1",this);
  setOpt1Cmd->SetParameterName("optpot1",true);
  setOpt1Cmd->SetDefaultValue(1) ;
setOpt2Cmd = new G4UIcmdWithADouble("/mat/optpot2",this);
  setOpt2Cmd->SetParameterName("optpot2",true);
  setOpt2Cmd->SetDefaultValue(1) ;
setOpt3Cmd = new G4UIcmdWithADouble("/mat/optpot3",this);
  setOpt3Cmd->SetParameterName("optpot3",true);
  setOpt3Cmd->SetDefaultValue(1) ;
setOpt4Cmd = new G4UIcmdWithADouble("/mat/optpot4",this);
  setOpt4Cmd->SetParameterName("optpot4",true);
  setOpt4Cmd->SetDefaultValue(1) ;
setOpt5Cmd = new G4UIcmdWithADouble("/mat/optpot5",this);
  setOpt5Cmd->SetParameterName("optpot5",true);
  setOpt5Cmd->SetDefaultValue(1) ;


  // eta1, eta2, eta3
  setEta1Cmd = new G4UIcmdWithADouble("/mat/eta1",this);
  setEta1Cmd->SetParameterName("eta1",true);
  setEta1Cmd->SetDefaultValue(1) ;
setEta2Cmd = new G4UIcmdWithADouble("/mat/eta2",this);
  setEta2Cmd->SetParameterName("eta2",true);
  setEta2Cmd->SetDefaultValue(1) ;
setEta3Cmd = new G4UIcmdWithADouble("/mat/eta3",this);
  setEta3Cmd->SetParameterName("eta3",true);
  setEta3Cmd->SetDefaultValue(1) ;
setEta4Cmd = new G4UIcmdWithADouble("/mat/eta4",this);
  setEta4Cmd->SetParameterName("eta4",true);
  setEta4Cmd->SetDefaultValue(1) ;
setEta5Cmd = new G4UIcmdWithADouble("/mat/eta5",this);
  setEta5Cmd->SetParameterName("eta5",true);
  setEta5Cmd->SetDefaultValue(1) ;

  // diffuse scattering diff1, diff2 ...
setDiff1Cmd = new G4UIcmdWithADouble("/mat/diff1",this);
  setDiff1Cmd->SetParameterName("diff1",true);
  setDiff1Cmd->SetDefaultValue(1) ;
setDiff2Cmd = new G4UIcmdWithADouble("/mat/diff2",this);
  setDiff2Cmd->SetParameterName("diff2",true);
  setDiff2Cmd->SetDefaultValue(1) ;
setDiff3Cmd = new G4UIcmdWithADouble("/mat/diff3",this);
  setDiff3Cmd->SetParameterName("diff3",true);
  setDiff3Cmd->SetDefaultValue(1) ;
setDiff4Cmd = new G4UIcmdWithADouble("/mat/diff4",this);
  setDiff4Cmd->SetParameterName("diff4",true);
  setDiff4Cmd->SetDefaultValue(1) ;
setDiff5Cmd = new G4UIcmdWithADouble("/mat/diff5",this);
  setDiff5Cmd->SetParameterName("diff5",true);
  setDiff5Cmd->SetDefaultValue(1) ;
  
  
  // theory parameters C,V, u_0, u_min, u_max, lambda, f_a, Gamma_1a, f_h, Gamma1h, Gamma2a, Gamma3a, Gamma4a, Gamma2h,Gamma3h, Gamma4h

setCCmd = new G4UIcmdWithADouble("/mat/C",this);
  setCCmd->SetParameterName("C",true);
  setCCmd->SetDefaultValue(1) ;
setVCmd = new G4UIcmdWithADouble("/mat/V",this);
  setVCmd->SetParameterName("V",true);
  setVCmd->SetDefaultValue(1) ;
setu_minCmd = new G4UIcmdWithADouble("/mat/u_min",this);
  setu_minCmd->SetParameterName("u_min",true);
  setu_minCmd->SetDefaultValue(1) ;
setu_maxCmd = new G4UIcmdWithADouble("/mat/u_max",this);
  setu_maxCmd->SetParameterName("u_max",true);
  setu_maxCmd->SetDefaultValue(1) ;
setlambdaCmd = new G4UIcmdWithADouble("/mat/lambda",this);
  setlambdaCmd->SetParameterName("lambda",true);
  setlambdaCmd->SetDefaultValue(1) ;
setf_aCmd = new G4UIcmdWithADouble("/mat/f_a",this);
  setf_aCmd->SetParameterName("f_a",true);
  setf_aCmd->SetDefaultValue(1) ;
setf_hCmd = new G4UIcmdWithADouble("/mat/f_h",this);
  setf_hCmd->SetParameterName("f_h",true);
  setf_hCmd->SetDefaultValue(1) ;

setGamma1hCmd = new G4UIcmdWithADouble("/mat/Gamma1_h",this);
  setGamma1hCmd->SetParameterName("Gamma1_h",true);
  setGamma1hCmd->SetDefaultValue(1) ;
setGamma1aCmd = new G4UIcmdWithADouble("/mat/Gamma1_a",this);
  setGamma1aCmd->SetParameterName("Gamma1_a",true);
  setGamma1aCmd->SetDefaultValue(1) ;
setGamma2hCmd = new G4UIcmdWithADouble("/mat/Gamma2_h",this);
  setGamma2hCmd->SetParameterName("Gamma2_h",true);
  setGamma2hCmd->SetDefaultValue(1) ;
setGamma2aCmd = new G4UIcmdWithADouble("/mat/Gamma2_a",this);
  setGamma2aCmd->SetParameterName("Gamma2_a",true);
  setGamma2aCmd->SetDefaultValue(1) ;
setGamma3hCmd = new G4UIcmdWithADouble("/mat/Gamma3_h",this);
  setGamma3hCmd->SetParameterName("Gamma3_h",true);
  setGamma3hCmd->SetDefaultValue(1) ;
setGamma3aCmd = new G4UIcmdWithADouble("/mat/Gamma3_a",this);
  setGamma3aCmd->SetParameterName("Gamma3_a",true);
  setGamma3aCmd->SetDefaultValue(1) ;

// Excited States
setexcitedFracCmd = new G4UIcmdWithADouble("/mat/excitedFrac",this);
  setexcitedFracCmd->SetParameterName("excitedFrac",true);
  setexcitedFracCmd->SetDefaultValue(1);

setmeanDownShiftLifetimeCmd = new G4UIcmdWithADouble("/mat/meanDownShiftLifetime",this);
  setmeanDownShiftLifetimeCmd->SetParameterName("meanDownShiftLifetime",true);
  setmeanDownShiftLifetimeCmd->SetDefaultValue(1);

setE_recoilCmd = new G4UIcmdWithADouble("/mat/E_recoil",this);
  setE_recoilCmd->SetParameterName("E_recoil",true);
  setE_recoilCmd->SetDefaultValue(1);

setmaxGlobalTimeCmd = new G4UIcmdWithADouble("/mat/maxGlobalTime",this);
setmaxGlobalTimeCmd->SetParameterName("maxGlobalTime",true);
  setmaxGlobalTimeCmd ->SetDefaultValue(1);
// End Excited States


setChopperFileCmd = new G4UIcmdWithAString("/mat/chopperfile",this);
  setChopperFileCmd->SetGuidance(" Set name of the chopper function file ");
  setChopperFileCmd->SetParameterName("chopperfile",true);
  setChopperFileCmd->SetDefaultValue("chopperfunction.dat") ;

  /*

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

G4cout << "generated primary messenger class " << G4endl;  
*/
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

UCNMaterialBoundary2Messenger::~UCNMaterialBoundary2Messenger()
{

  delete matDirectory;
  delete setFileName1Cmd;
  delete setFileName2Cmd;
  delete setFileName3Cmd;
  delete setstepsCmd;
    delete setfoilCmd;
    delete settimeCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



void UCNMaterialBoundary2Messenger::SetNewValue(G4UIcommand* command, G4String newValue)
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
   if( command == setFileName1Cmd)
   {
	   //G4cout << " set file name 1 " << G4endl;
	   fAction->SetFileName1(G4String(newValue));
   }
 if( command == setFileName2Cmd)
   {
           fAction->SetFileName2(G4String(newValue));
   }

 if( command == setFileName3Cmd)
   {
           fAction->SetFileName3(G4String(newValue));
   }

if( command == setstepsCmd)
   {
           fAction->SetSteps(setstepsCmd->GetNewDoubleValue(newValue));
   }
if( command == setlifeCmd)
   {
           fAction->SetLifeTime(setlifeCmd->GetNewDoubleValue(newValue));
   }
if( command == sethelifeCmd)
   {
           fAction->SetHeLifeTime(sethelifeCmd->GetNewDoubleValue(newValue));
   }
   if( command == setfoilCmd)
   {
           fAction->SetFoil(setfoilCmd->GetNewDoubleValue(newValue));
   }
  if( command == settimeCmd)
   {
           fAction->SetTime(settimeCmd->GetNewDoubleValue(newValue));
   }


// 
if( command == setDiff1Cmd)  {      fAction->SetDiff1(setDiff1Cmd->GetNewDoubleValue(newValue));}
if( command == setDiff2Cmd)  {      fAction->SetDiff2(setDiff2Cmd->GetNewDoubleValue(newValue));}
if( command == setDiff3Cmd)  {      fAction->SetDiff3(setDiff3Cmd->GetNewDoubleValue(newValue));}
if( command == setDiff4Cmd)  {      fAction->SetDiff4(setDiff4Cmd->GetNewDoubleValue(newValue));}
if( command == setDiff5Cmd)  {      fAction->SetDiff5(setDiff5Cmd->GetNewDoubleValue(newValue));}

if( command == setOpt1Cmd)  {      fAction->SetOpt1(setOpt1Cmd->GetNewDoubleValue(newValue));}
if( command == setOpt2Cmd)  {      fAction->SetOpt2(setOpt2Cmd->GetNewDoubleValue(newValue));}
if( command == setOpt3Cmd)  {      fAction->SetOpt3(setOpt3Cmd->GetNewDoubleValue(newValue));}
if( command == setOpt4Cmd)  {      fAction->SetOpt4(setOpt4Cmd->GetNewDoubleValue(newValue));}
if( command == setOpt5Cmd)  {      fAction->SetOpt5(setOpt5Cmd->GetNewDoubleValue(newValue));}
if( command == setEta1Cmd)  {      fAction->SetEta1(setEta1Cmd->GetNewDoubleValue(newValue));}
if( command == setEta2Cmd)  {      fAction->SetEta2(setEta2Cmd->GetNewDoubleValue(newValue));}
if( command == setEta3Cmd)  {      fAction->SetEta3(setEta3Cmd->GetNewDoubleValue(newValue));}
if( command == setEta4Cmd)  {      fAction->SetEta4(setEta4Cmd->GetNewDoubleValue(newValue));}
if( command == setEta5Cmd)  {      fAction->SetEta5(setEta5Cmd->GetNewDoubleValue(newValue));}
if( command == setCCmd)  {      fAction->SetC(setCCmd->GetNewDoubleValue(newValue));}
if( command == setVCmd)  {      fAction->SetV(setVCmd->GetNewDoubleValue(newValue));}
if( command == setf_aCmd)  {      fAction->Setf_a(setf_aCmd->GetNewDoubleValue(newValue));}
if( command == setf_hCmd)  {      fAction->Setf_h(setf_hCmd->GetNewDoubleValue(newValue));}
if( command == setlambdaCmd)  {      fAction->Setlambda(setlambdaCmd->GetNewDoubleValue(newValue));}
if( command == setGamma1aCmd)  {      fAction->SetGamma1a(setGamma1aCmd->GetNewDoubleValue(newValue));}
if( command == setGamma1hCmd)  {      fAction->SetGamma1h(setGamma1hCmd->GetNewDoubleValue(newValue));}
if( command == setGamma2aCmd)  {      fAction->SetGamma2a(setGamma2aCmd->GetNewDoubleValue(newValue));}
if( command == setGamma2hCmd)  {      fAction->SetGamma2h(setGamma2hCmd->GetNewDoubleValue(newValue));}
if( command == setGamma3aCmd)  {      fAction->SetGamma3a(setGamma3aCmd->GetNewDoubleValue(newValue));}
if( command == setGamma3hCmd)  {      fAction->SetGamma3h(setGamma3hCmd->GetNewDoubleValue(newValue));}
if( command == setu_minCmd)  {      fAction->Setu_min(setu_minCmd->GetNewDoubleValue(newValue));}
if( command == setu_maxCmd)  {      fAction->Setu_max(setu_maxCmd->GetNewDoubleValue(newValue));}
// Excited States
if( command == setexcitedFracCmd)  {      fAction->SetexcitedFrac(setexcitedFracCmd->GetNewDoubleValue(newValue));}
if( command == setmeanDownShiftLifetimeCmd)  {      fAction->SetmeanDownShiftLifetime(setmeanDownShiftLifetimeCmd->GetNewDoubleValue(newValue));}
if( command == setE_recoilCmd)  {      fAction->SetE_recoil(setE_recoilCmd->GetNewDoubleValue(newValue));}
if( command == setmaxGlobalTimeCmd)  {      fAction->SetmaxGlobalTime(setmaxGlobalTimeCmd->GetNewDoubleValue(newValue));}
// End Excited States


//if( command == setChopperFileCmd)  {      fAction->SetChopperFunction(G4String(newValue));}
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
         //    fAction->SetEspec(G4String(newValue));
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
           
