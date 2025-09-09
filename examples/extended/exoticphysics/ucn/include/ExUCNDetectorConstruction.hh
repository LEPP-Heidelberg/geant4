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
/// \file exoticphysics/ucn/include/ExUCNDetectorConstruction.hh
/// \brief Definition of the ExUCNDetectorConstruction class
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef ExUCNDetectorConstruction_h
#define ExUCNDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
//#include "ExUCNDetectorConstructionMessenger.hh"

class G4Material;
class ExUCNDetectorConstructionMessenger;
class G4UniformGravityField;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class ExUCNDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    ExUCNDetectorConstruction();
    virtual ~ExUCNDetectorConstruction();
    virtual void SetProperty1(G4double e);
    virtual void SetProperty2(G4double e);
    virtual void SetProperty3(G4double e);
    virtual void SetProperty4(G4double e);
    virtual void SetProperty5(G4double e);
    virtual void SetProperty6(G4double e);
    virtual void SetProperty7(G4double e);
    virtual void SetProperty8(G4double e);
    virtual void SetProperty9(G4double e);
    virtual void SetProperty10(G4double e);
    virtual void SetProperty11(G4double e);
    virtual void SetProperty12(G4double e);
   virtual void SetProperty13(G4double e);
 virtual void SetProperty14(G4double e);

  public:
    virtual G4VPhysicalVolume* Construct();
    virtual void ConstructSDandField();

//       ExUCNDetectorConstructionMessenger* fMessenger;

  private:
    G4Material* fVacuum;
    G4Material* fGuideMaterial;
    G4Material* DetMat;
    G4Material* WorldMaterial;
    G4Material* GuideMaterial1;
    G4Material* GuideMaterial2;
    G4Material* GuideMaterial3;
    G4Material* GuideMaterial4;
    G4Material* GuideMaterial5;


    ExUCNDetectorConstructionMessenger* fMessenger;

    static G4ThreadLocal G4UniformGravityField* fField;

  private:
    void DefineMaterials();
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
