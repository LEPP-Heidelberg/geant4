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
/// \file exoticphysics/ucn/src/ExUCNDetectorConstruction.cc
/// \brief Implementation of the ExUCNDetectorConstruction class
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "ExUCNDetectorConstruction.hh"
#include "ExUCNDetectorConstructionMessenger.hh"
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4FieldManager.hh"
#include "G4GeometryManager.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4PhysicalConstants.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4RepleteEofM.hh"
#include "G4SolidStore.hh"
#include "G4SystemOfUnits.hh"
#include "G4TransportationManager.hh"
#include "G4Tubs.hh"
#include "G4Torus.hh"
#include "G4Cons.hh"
#include "G4UCNMaterialPropertiesTable.hh"
#include "G4UniformGravityField.hh"
#include "G4UserLimits.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VisAttributes.hh"
// #include "G4EqGravityField.hh"
#include "G4CashKarpRKF45.hh"
#include "G4ChordFinder.hh"
#include "G4ClassicalRK4.hh"
#include "G4DormandPrince745.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4PropagatorInField.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"

G4double optpot1 = 0.;
G4double optpot2 = 0.;
G4double optpot3 = 0.;
G4double optpot4 = 0.;
G4double optpot5 = 0.;

G4double eta1 = 0.;
G4double eta2 = 0.;
G4double eta3 = 0.;
G4double eta4 = 0.;
G4double eta5 = 0.;

G4double diff1 = 0.;
G4double diff2 = 0.;
G4double diff3 = 0.;
G4double diff4 = 0.;
G4double diff5 = 0.;



G4double maxstp = 1;
G4double valveheight = 0.;

// SBI v2 Sanity Checks
G4double converterRadius = 37.5;
G4double converterLength = 3000./2;
// End SBI v2
;//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExUCNDetectorConstruction::ExUCNDetectorConstruction() : fVacuum(0), fGuideMaterial(0)
{

	// material
  fMessenger = new ExUCNDetectorConstructionMessenger(this);

  DefineMaterials();
//  fMessenger = new ExUCNDetectorConstructionMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExUCNDetectorConstruction::~ExUCNDetectorConstruction()
{
  if (fField) delete fField;
  delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExUCNDetectorConstruction::DefineMaterials()
{
//# tube above chopper blades, tube to detector 
optpot1 = 183;
eta1 = 0.0001;
diff1 = 0.02;
//# converter volume (INSIDE LIQUID HELIUM, therfore 18.5 neV reduced potential
optpot2 = 96.5;
eta2 = 0.00005;
diff2 = 0.02;
//# cube inside the converter. reduce by 18.5 because submerged in lHe
optpot3 = 233.5;
eta3 = 0.00005;
diff3 = 0.02;
// #end plugs of converter, reduce by 18.5neV becuase submerged in lHe 
optpot4 = 233.5;
eta4 = 0.00005;
diff4 = 0.02; 
  
  

	G4String name, symbol ;             // a=mass of a mole;
  G4double a, z, density ;            // z=mean number of protons;  
  G4int nel;
  G4int ncomponents;
  G4double fractionmass, pressure, temperature;

//  a = 1.01*g/mole;
//  G4Element* elH  = new G4Element(name="Hydrogen",symbol="H" , z= 1., a);
  a = 2.0*g/mole;
  G4Element* elD  = new G4Element(name="Deuterium",symbol="D" , z= 1., a);
//  a = 12.01*g/mole;
//  G4Element* elC = new G4Element(name="Carbon", symbol="C", z=6., a);
  a = 14.01*g/mole;
  G4Element* elN  = new G4Element(name="Nitrogen",symbol="N" , z= 7., a);
  a = 16.00*g/mole;
  G4Element* elO  = new G4Element(name="Oxygen"  ,symbol="O" , z= 8., a);
  a = 39.948*g/mole;
  G4Element* elAr = new G4Element(name="Argon", symbol="Ar", z=18., a);
  a = 9.01*g/mole;
  G4Element* elBe = new G4Element(name="Ber", symbol="Be", z=4., a);
  a = 26.98*g/mole;
  G4Element* elAl = new G4Element(name="Aluminium", symbol="Al", z=13., a);
  a = 55.8*g/mole;
  G4Element* elFe = new G4Element(name="Iron", symbol="Fe", z=26., a);
  a = 58.7*g/mole;
  G4Element* elNi = new G4Element(name="Nickel", symbol="Ni", z=28., a);
  a = 58.*g/mole;
  G4Element* elNi58 = new G4Element(name="Nickel58", symbol="Ni58", z=28., a);
  a = 58.*g/mole;
  G4Element* elCytop = new G4Element(name="Cytopmaterial", symbol="Cy", z=15., a);


// Beryllium
  density = 1.848*g/cm3;
  G4Material* Beryllium = new G4Material(name="Beryllium", density, nel=1);
  Beryllium->AddElement(elBe,1);
  // Aluminium
  density = 2.7*g/cm3;
  G4Material* Aluminium = new G4Material(name="Aluminium", density, nel=1);
  Aluminium->AddElement(elAl,1);
  // Iron
  density = 7.874*g/cm3;
  G4Material* Iron = new G4Material(name="Iron", density, nel=1);
  Iron->AddElement(elFe,1);
  // Nickel
  density = 8.908*g/cm3;
  G4Material* Nickel = new G4Material(name="Nickel", density, nel=1);
  Nickel->AddElement(elNi,1);
  // Nickel58
  density = 8.908*g/cm3;
  G4Material* Nickel58 = new G4Material(name="Nickel58", density, nel=1);
  Nickel58->AddElement(elNi58,1);

// Cytop
  density = 1.6*g/cm;   // PROPERTIES TO BE CHECKED !!!!
  G4Material* Cytop = new G4Material(name="Cytop", density, nel=1);
  Cytop->AddElement(elCytop,1);



	G4NistManager* nistMan = G4NistManager::Instance();
  fVacuum = nistMan->FindOrBuildMaterial("G4_Galactic");
  //fGuideMaterial = nistMan->FindOrBuildMaterial("Iron");
  GuideMaterial1 = nistMan->FindOrBuildMaterial("Aluminium");
  GuideMaterial2 = nistMan->FindOrBuildMaterial("Nickel");
  GuideMaterial3 = nistMan->FindOrBuildMaterial("Cytop");
  GuideMaterial4 = nistMan->FindOrBuildMaterial("Nickel58");
  GuideMaterial5 = nistMan->FindOrBuildMaterial("Iron");
  DetMat = nistMan->FindOrBuildMaterial("G4_Ni");
/*
  G4UCNMaterialPropertiesTable* MPT = new G4UCNMaterialPropertiesTable();
 // MPT->AddConstProperty("TAB", 0);
 MPT->AddConstProperty("DIFFUSION", 0);
  MPT->AddConstProperty("FERMIPOT", -49);  // Golub, Table 2.1 in neiV
  MPT->AddConstProperty("LOSS", 1e-3);  //  Golub, Table 2.1
  MPT->AddConstProperty("LOSSCS", 0.);
  MPT->AddConstProperty("ABSCS", 4.49);  // 1/v loss cross-section  at room temp.
  MPT->AddConstProperty("SPINFLIP", 0.0);  
  MPT->AddConstProperty("SCATCS", 5); // (incoherent) "elastic" scattering cs

  fGuideMaterial->SetMaterialPropertiesTable(MPT);
*/
  G4double neV = 1e-9 * eV;

  G4UCNMaterialPropertiesTable* MPT = new G4UCNMaterialPropertiesTable();
//  MPT2->AddConstProperty("TAB", 5);
  MPT->AddConstProperty("DIFFUSION", diff1);
  MPT->AddConstProperty("FERMIPOT", optpot1);  // Golub, Table 2.1 in neV
  MPT->AddConstProperty("LOSS", eta1);  //  Golub, Table 2.1
  MPT->AddConstProperty("LOSSCS", 0.);
  MPT->AddConstProperty("SPINFLIP", 0.0);

  MPT->AddConstProperty("ABSCS", 4.49);  // 1/v loss cross-section  at room temp.
  MPT->AddConstProperty("SCATCS", 5);  // (incoherent) "elastic" scattering cs
  GuideMaterial5->SetMaterialPropertiesTable(MPT);

//  MPT->SetMicroRoughnessParameters(30 * nm, 1 * nm, 180, 1000, 0 * degree, 90 * degree, 1 * neV,
//                                   1000 * neV, 15, 15, 0.01 * degree);


  G4UCNMaterialPropertiesTable* MPT2 = new G4UCNMaterialPropertiesTable();
//  MPT2->AddConstProperty("TAB", 1);
  MPT2->AddConstProperty("DIFFUSION", diff1);
  MPT2->AddConstProperty("FERMIPOT", optpot1);  // Golub, Table 2.1 in neV
  MPT2->AddConstProperty("LOSS", eta1);  //  Golub, Table 2.1
  MPT2->AddConstProperty("LOSSCS", 0.);
  MPT2->AddConstProperty("SPINFLIP", 0.0);

  MPT2->AddConstProperty("ABSCS", 4.49);  // 1/v loss cross-section  at room temp.
  MPT2->AddConstProperty("SCATCS", 1);  // (incoherent) "elastic" scattering cs
  GuideMaterial1->SetMaterialPropertiesTable(MPT2);

 G4UCNMaterialPropertiesTable* MPT3 = new G4UCNMaterialPropertiesTable();
//  MPT3->AddConstProperty("TAB", 2);
  MPT3->AddConstProperty("DIFFUSION", diff2);
  MPT3->AddConstProperty("FERMIPOT", optpot2);  // Golub, Table 2.1 in neiV
  MPT3->AddConstProperty("LOSS", eta2);  //  Golub, Table 2.1
  MPT3->AddConstProperty("LOSSCS", 0.);
  MPT3->AddConstProperty("SPINFLIP", 0.0);
  MPT3->AddConstProperty("ABSCS", 4.49);  // 1/v loss cross-section  at room temp.
  MPT3->AddConstProperty("SCATCS", 2);  // (incoherent) "elastic" scattering cs
  GuideMaterial2->SetMaterialPropertiesTable(MPT3);
  
G4cout << "************************************optical potential 2 , wich is the converter tube " << optpot2 << G4endl;

 G4UCNMaterialPropertiesTable* MPT4 = new G4UCNMaterialPropertiesTable();
//  MPT4->AddConstProperty("TAB", 3);
  MPT4->AddConstProperty("DIFFUSION", diff3);
  MPT4->AddConstProperty("FERMIPOT", optpot3);  // Golub, Table 2.1 in neiV
  MPT4->AddConstProperty("LOSS", eta3);  //  Golub, Table 2.1
  MPT4->AddConstProperty("LOSSCS", 0.);
  MPT4->AddConstProperty("SPINFLIP", 0.0);
  MPT4->AddConstProperty("ABSCS", 4.49);  // 1/v loss cross-section  at room temp.
  MPT4->AddConstProperty("SCATCS", 3);  // (incoherent) "elastic" scattering cs 
  GuideMaterial3->SetMaterialPropertiesTable(MPT4);

  
  G4UCNMaterialPropertiesTable* MPT5 = new G4UCNMaterialPropertiesTable();
//  MPT5->AddConstProperty("TAB", 4);
  MPT5->AddConstProperty("DIFFUSION", diff4);
  MPT5->AddConstProperty("FERMIPOT", optpot4);  // Golub, Table 2.1 in neiV
  MPT5->AddConstProperty("LOSS", eta4);  //  Golub, Table 2.1
  MPT5->AddConstProperty("LOSSCS", 0.);
  MPT5->AddConstProperty("SPINFLIP", 0.0);
  MPT5->AddConstProperty("ABSCS", 4.49);  // 1/v loss cross-section  at room temp.
  MPT5->AddConstProperty("SCATCS", 4);  // (incoherent) "elastic" scattering cs

  GuideMaterial4->SetMaterialPropertiesTable(MPT5);

  //G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* ExUCNDetectorConstruction::Construct()
{

  //
  // World
  //

  G4double worldSizeX = 30. * m;
  G4double worldSizeY = 30. * m;
  G4double worldSizeZ = 30. * m;

  G4Box* solidWorld = new G4Box("World", worldSizeX / 2., worldSizeY / 2., worldSizeZ / 2.);

  G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld, fVacuum, "World");

  G4VPhysicalVolume* physiWorld =
    new G4PVPlacement(0, G4ThreeVector(), "World", logicWorld, 0, false, 0);

  G4double maxStep = maxstp * mm;
  G4double maxTime = 10000000000000 * s;
  G4UserLimits* stepLimit = new G4UserLimits(maxStep, DBL_MAX, maxTime);
  logicWorld->SetUserLimits(stepLimit);

  G4VisAttributes* guideColor = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0));
  guideColor->SetVisibility(true);
  guideColor->SetForceWireframe(true);

  G4VisAttributes* endPlateColor = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0));
  endPlateColor->SetVisibility(true);
  endPlateColor->SetForceSolid(true);

  logicWorld->SetVisAttributes(G4VisAttributes::GetInvisible());

  //// Supersun
  //
   
  // upper cube: guidematerial1
  // connecting tube: guidematerial1
  // second cube: guidematerial1
  // downward tube: guidematerial1
  // cube in converter: guidematerial3
  // convertervolume: guidenmaterial2
  // end plugs of converter: guidematerial4
  // source exit tube: guidematerial1
  // vacuumseparationfoilvolume: fVacuum
  // guide away from source: guidematerial1
  // 90 deg downwards bend: guidematerial1 
  // tube down to chopper blades: guidematerial1
  // chopperblades: fGuideMaterial
  // "cone" below: GuideMaterial1
  // tube to detector: GuideMaterial1 
  // extra tube piece before detector: GuideMaterial1
  // detector: detMat 
  //
  //
  // the upper cube
  //
  G4double cubesizehalf = 40.;
 G4double rMax = 4*cm;
  G4double rMin = 2.5*cm;
    G4double longzero = 0.;
 G4double sideoffset = 270.;

  G4Box * solidCube = new G4Box("Cube",cubesizehalf*mm,cubesizehalf*mm,cubesizehalf*mm);       //its size
  G4Tubs *solidHoleInCube = new G4Tubs("SolidHole", 0., 25., 35., 0., twopi);
  G4VSolid* subtract = new G4SubtractionSolid("box-cylinder", solidCube, solidHoleInCube,0, G4ThreeVector(0.,0.,10.));
  G4Tubs *solidHoleInCube2 = new G4Tubs("SolidHole", 0., 25., 35., 0., twopi);
G4RotationMatrix * zRot001 = new G4RotationMatrix();
 zRot001->rotateX(3.14159/2*rad);
  G4VSolid* subtract2 = new G4SubtractionSolid("box-cylinder2", subtract, solidHoleInCube2,zRot001, G4ThreeVector(0.,-10.,0.));
  G4LogicalVolume *logicCube = new G4LogicalVolume(subtract2,          //its solid
                                   GuideMaterial4,       //its material
                                   "cube");            //its name

 G4RotationMatrix * zRot00 = new G4RotationMatrix();
 zRot00->rotateZ( 3.14159/2*rad);
 zRot00->rotateY(2*3.14159/2*rad);
 zRot00->rotateZ( 3.14159*rad);

 logicCube->SetUserLimits(stepLimit);


 G4VPhysicalVolume *physiCube = new G4PVPlacement(zRot00, G4ThreeVector(-sideoffset - 2*cubesizehalf,0.,longzero+sideoffset+2*cubesizehalf), "Cube1", logicCube, physiWorld, false, 0);

G4cout << "*** Geometry infos *** " << G4endl;
G4cout << "y+ is vertical upwards, z + is the direction towards the reactor" << G4endl;
G4cout << "xyz positions are geometrical centers of highest level geometry (eg a cube with holes has the cube center at the xyz coordinates) " << G4endl;
G4cout << "***upper cube (the top of the vertical extraction from the source), a cube with two cylindrical holes " << G4endl;
G4cout << " size: " << cubesizehalf * 2 << G4endl;
G4cout << " vertical hole r = 25 mm , length = 70 mm, shifted by 10 mm off center " << G4endl;
G4cout << " horizontal hole r = 25 mm, length = 70 mm, shifted by 10 mm off center " << G4endl;
G4cout << " x-position " << -sideoffset - 2*cubesizehalf << G4endl;
G4cout << " y-position " << 0 << G4endl;
G4cout << " z-position " << longzero+sideoffset+2*cubesizehalf << G4endl;


// the connecting tube
  rMax = 4*cm;
  rMin = 2.5*cm;
  G4Tubs *solidConnectingTube = new G4Tubs("SolidTubeC", rMin, rMax, sideoffset/2., 0., twopi);
  G4LogicalVolume *logicConnectingTube = new  G4LogicalVolume(solidConnectingTube,  GuideMaterial1, "SolidTubeC");
  G4RotationMatrix * zRot5 = new G4RotationMatrix();
  G4VPhysicalVolume *physiConnectingTube = new G4PVPlacement(zRot5, G4ThreeVector(-sideoffset-2*cubesizehalf,0.,longzero+ sideoffset/2 + cubesizehalf), "ConnectingTube", logicConnectingTube, physiWorld, false, 0);

   logicConnectingTube->SetUserLimits(stepLimit);

G4cout << "*** tube connecting the cube horizontally to a second cube towards extraction " << G4endl;
G4cout << " inside radius " << rMin << G4endl;
G4cout << " length " << sideoffset << G4endl;
G4cout << " x-position " << -sideoffset-2*cubesizehalf << G4endl;
G4cout << " y-position " << 0 << G4endl;
G4cout << " z-position " << longzero+ sideoffset/2 + cubesizehalf << G4endl;

// the second cube
  G4Box * solidCube3 = new G4Box("Cube3", cubesizehalf*mm,cubesizehalf*mm,cubesizehalf*mm);       //its size
  G4Tubs *solidHoleInCube3 = new G4Tubs("SolidHole3", 0., 25., 35., 0., twopi);
  G4VSolid* subtract3 = new G4SubtractionSolid("box-cylinder1", solidCube3, solidHoleInCube3,0, G4ThreeVector(0.,0.,10.));
  G4Tubs *solidHoleInCube4 = new G4Tubs("SolidHole4", 0., 25., 35., 0., twopi);
G4RotationMatrix * zRot002 = new G4RotationMatrix();
 zRot002->rotateX(3.14159/2*rad);
  G4VSolid* subtract4 = new G4SubtractionSolid("box-cylinder3", subtract3, solidHoleInCube4,zRot002, G4ThreeVector(0.,-10.,0.));

  G4LogicalVolume *logicCube2 = new G4LogicalVolume(subtract4,          //its solid
                                   GuideMaterial4,       //its material
                                   "cube2");            //its name
 logicCube2->SetUserLimits(stepLimit);


 G4RotationMatrix * zRot003 = new G4RotationMatrix();
 zRot003->rotateZ(0.);//-3.14159/2*rad);

G4VPhysicalVolume *physiCube2 = new G4PVPlacement(zRot003, G4ThreeVector(-sideoffset - 2*cubesizehalf,0.,longzero), "Cube2", logicCube2, physiWorld, false, 0);

G4cout << "*** secound uppper cube (the top of the vertical extraction from the source), a cube with two cylindrical holes " << G4endl;
G4cout << " size: " << cubesizehalf * 2 << G4endl;
G4cout << " vertical hole r = 25 mm , length = 70 mm, shifted by 10 mm off center " << G4endl;
G4cout << " horizontal hole r = 25 mm, length = 70 mm, shifted by 10 mm off center " << G4endl;
G4cout << " x-position " << -sideoffset - 2*cubesizehalf << G4endl;
G4cout << " y-position " << 0 << G4endl;
G4cout << " z-position " << longzero << G4endl;



G4double convcubesize = 50.; // halfsize of convertercube
G4double tiefe = 286.;

// the downward tube into the converter volume
 G4double reducedheight = 0;
  rMax = 7*cm;
//rMin = 0.0*cm; // Use this for closed converter

  rMin = 2.5*cm; // Use this for vTOF / Extraction
  G4Tubs *solidDownwardTube = new G4Tubs("SolidTubeD", rMin, rMax, 196./2. - reducedheight /2., 0., twopi);
  G4LogicalVolume *logicDownwardTube = new  G4LogicalVolume(solidDownwardTube,  GuideMaterial1, "SolidTubeD");
  G4RotationMatrix * zRot7 = new G4RotationMatrix();
 // zRot7->rotateZ(3.14159/4*rad);
  zRot7->rotateX(3.14159/2*rad);
 logicDownwardTube->SetUserLimits(stepLimit);


  G4VPhysicalVolume *physiDownwardTube = new G4PVPlacement(zRot7, G4ThreeVector(-sideoffset - 2*cubesizehalf,-tiefe+196./2.+convcubesize + reducedheight/2.,longzero), "DownwardTube", logicDownwardTube, physiWorld, false, 0);

G4cout << "*** tube downward into converter " << G4endl;
G4cout << " inside radius " << rMin << G4endl;
G4cout << " length " << 196./2. - reducedheight /2. << G4endl;
G4cout << " x-position " << -sideoffset - 2*cubesizehalf << G4endl;
G4cout << " y-position " << -tiefe+196./2.+convcubesize + reducedheight/2. << G4endl;
G4cout << " z-position " << longzero  << G4endl;


// the disc valve sealing the converter
   rMax = 3.0*cm;
   rMin = 0*cm;
   G4double valveh = 0.5*cm;
  G4Tubs *solidValve = new G4Tubs("SolidValve", rMin, rMax, valveh/2, 0., twopi);
  G4LogicalVolume *logicValve = new  G4LogicalVolume(solidValve,  GuideMaterial1, "SolidValve");
 logicValve->SetUserLimits(stepLimit);


//G4VPhysicalVolume *physiValve = new G4PVPlacement(zRot7, G4ThreeVector(-sideoffset - 2*cubesizehalf,-tiefe+196./2.+convcubesize + reducedheight/2. + valveheight,longzero), "Valve", logicValve, physiWorld, false, 0);
G4VPhysicalVolume* physiValve = new G4PVPlacement(zRot7, G4ThreeVector( -sideoffset - 2*cubesizehalf, -tiefe + converterRadius - valveh/2. + valveheight, longzero), "Valve", logicValve, physiWorld, false, 0);
G4cout << "*** valve " << G4endl;
G4cout << " valve height " << valveheight <<G4endl;
G4cout << " inside radius " << rMin << G4endl;
G4cout << " length " << 196./2. - reducedheight /2. << G4endl;
G4cout << " x-position " << -sideoffset - 2*cubesizehalf << G4endl;
G4cout << " Valve y-position " << valveheight -tiefe+196./2.+convcubesize + reducedheight/2. << G4endl;
G4cout << " z-position " << longzero  << G4endl;

// Rod actuating the valve
G4double fixed_ceiling = 25.; //+ cubesizehalf;  // fixed top of rod
G4double y_valve_top   = -tiefe + converterRadius + valveheight;
G4double rodh          = fixed_ceiling - y_valve_top;  // shrinks as valveheight increases

G4Tubs *solidRod = new G4Tubs("SolidRod", 0., 0.7*mm, rodh/2., 0., twopi);
G4LogicalVolume *logicRod = new G4LogicalVolume(solidRod, GuideMaterial1, "SolidRod");
logicRod->SetUserLimits(stepLimit);

// Place rod so its bottom touches the top of the valve
G4double y_rod_centre = y_valve_top + rodh/2.;
G4VPhysicalVolume* physiRod = new G4PVPlacement(
    zRot7,
    G4ThreeVector(-sideoffset - 2*cubesizehalf, y_rod_centre, longzero),
    "Rod", logicRod, physiWorld, false, 0);

G4cout << "*** rod " << G4endl;
G4cout << " rod height " << rodh << G4endl;
G4cout << " rod y-centre " << y_rod_centre << G4endl;


/*
/////
 G4UserLimits* stepLimit2 = new G4UserLimits(0.1, DBL_MAX, maxTime);

//the (empty) volume that defines the helium volume entry
  G4double hetiefe = 280;
  G4Tubs *solidDownwardTubeempty = new G4Tubs("SolidTubeDempty", 0., 24.0, 20 /2., 0., twopi);
  G4LogicalVolume *logicDownwardTubeempty = new  G4LogicalVolume(solidDownwardTubeempty,  GuideMaterial1, "SolidTubeDempty");
 // G4LogicalVolume *logicDownwardTubeempty = new  G4LogicalVolume(solidDownwardTubeempty,  fVacuum, "SolidTubeDempty");

  //  G4RotationMatrix * zRot7 = new G4RotationMatrix();
 // zRot7->rotateZ(3.14159/4*rad);
//  zRot7->rotateX(3.14159/2*rad);
 logicDownwardTubeempty->SetUserLimits(stepLimit2);

 G4VPhysicalVolume *physiDownwardTubeempty = new G4PVPlacement(zRot7, G4ThreeVector(-sideoffset - 2*cubesizehalf,-hetiefe+20./2.+convcubesize + reducedheight/2.,longzero), "DownwardTubeempty", logicDownwardTubeempty, physiWorld, false, 0);

*/
















////  

// cube in converter volume with connection to upward tube:
  G4Box * solidCubeConverter = new G4Box("convertreCube3", convcubesize*mm,convcubesize*mm,convcubesize*mm);       //its size


  G4Tubs *solidHoleInCubeConverter = new G4Tubs("SolidHoleConverter", 0., 25., 25., 0., twopi) ;
  G4VSolid* subtract6 = new G4SubtractionSolid("box-cylinder5",
                  solidCubeConverter, solidHoleInCubeConverter,0, G4ThreeVector(0.,0.,25.));
  //G4VSolid* subtract6 = solidCubeConverter;                 
  G4Tubs *solidHoleInCubeConverter2 = new G4Tubs("SolidHoleConverter2 ", 0., 37.5, 50.1, 0., twopi);
G4RotationMatrix * zRot006 = new G4RotationMatrix();
 zRot006->rotateX(3.14159/2*rad);
  G4VSolid* subtract7 = new G4SubtractionSolid("box-cylinder6", subtract6, solidHoleInCubeConverter2,zRot006, G4ThreeVector(0.,0.,0.));

  G4LogicalVolume *logicCubeConverter = new G4LogicalVolume(subtract7,          //its solid
                                   GuideMaterial3,       //its material
                                   "cubeconverter");            //its name
 G4RotationMatrix * zRot008 = new G4RotationMatrix();
 zRot008->rotateX(1*3.14159/2*rad);
logicCubeConverter->SetUserLimits(stepLimit);
 G4VPhysicalVolume *physiCubeConverter = new G4PVPlacement(zRot008, G4ThreeVector(-sideoffset - 2*cubesizehalf,-tiefe,longzero), "CubeConvertre", logicCubeConverter, physiWorld, false, 0);




G4cout << "*** cube at the converter: a cube with two cylindrical holes, one with the diameter of the converter and one with the diameter of the extraction " << G4endl;

G4cout << "(i know this needs to be 56 m diameter for the upwards tube, tbd)" << G4endl;

G4cout << " size: " << cubesizehalf * 2 << G4endl;
G4cout << " vertical hole r = 25 mm , length = 70 mm, shifted by 10 mm off center " << G4endl;
G4cout << " horizontal hole r = 37.5 mm, length = 100 mm (through going) " << G4endl;
G4cout << " x-position " << -sideoffset - 2*cubesizehalf << G4endl;
G4cout << " y-position " << -tiefe << G4endl;
G4cout << " z-position " << longzero << G4endl;


// converter volume
  G4Tubs *solidConverter = new G4Tubs("SolidConverter", converterRadius, 150., converterLength, 0., twopi);
  G4LogicalVolume *logicConverter = new  G4LogicalVolume(solidConverter,  GuideMaterial2, "SolidConverter");
  G4RotationMatrix * zRot13 = new G4RotationMatrix();
  logicConverter->SetUserLimits(stepLimit);
  G4VPhysicalVolume *physiConvertre = new G4PVPlacement(zRot13, G4ThreeVector(-sideoffset - 2*cubesizehalf,-tiefe,longzero-converterLength-convcubesize ), "Convertre", logicConverter, physiWorld, false, 0);

G4cout << "*** converter tube" << G4endl;
G4cout << " inside radius " << converterRadius << G4endl;
G4cout << " length " << converterLength*2 << G4endl;
G4cout << " x-position " << -sideoffset - 2*cubesizehalf << G4endl;
G4cout << " y-position " << -tiefe << G4endl;
G4cout << " z-position " << longzero-1500.-convcubesize  << G4endl;
  
  // We scale up the plugs by 50/37.5 * converterRadius to keep the volume sealed
  // converter volume PLUGs
  G4Tubs *solidConverterPlug = new G4Tubs("SolidConverterPlug1", 0., 50./37.5 * converterRadius, 5., 0., twopi);
  G4LogicalVolume *logicConverterPlug = new  G4LogicalVolume(solidConverterPlug,  GuideMaterial4, "SolidConverterPlug");
  logicConverterPlug->SetUserLimits(stepLimit);
  G4VPhysicalVolume *physiConverterPlug = new G4PVPlacement(zRot13, G4ThreeVector(-sideoffset - 2*cubesizehalf,-tiefe,longzero-(converterLength*2+5.)-convcubesize), "Converterplug1", logicConverterPlug, physiWorld, false, 0);

G4cout << "*** converter tube plug1" << G4endl;
G4cout << " radius " << 50 << G4endl;
G4cout << " length " << 10 << G4endl;
G4cout << " x-position " << -sideoffset - 2*cubesizehalf << G4endl;
G4cout << " y-position " << -tiefe << G4endl;
G4cout << " z-position " << longzero-3005.-convcubesize  << G4endl;


//converter volume PLUGs
  G4Tubs *solidConverterPlug2 = new G4Tubs("SolidConverterPlg2", 0., 50./37.5 * converterRadius, 5., 0., twopi);
  G4LogicalVolume *logicConverterPlug2 = new  G4LogicalVolume(solidConverterPlug2,  GuideMaterial4, "SolidConverterPlug2");
  logicConverterPlug2->SetUserLimits(stepLimit);
  G4VPhysicalVolume *physiConverterPlug2 = new G4PVPlacement(zRot13, G4ThreeVector(-sideoffset - 2*cubesizehalf,-tiefe,longzero+55.), "Converterplug2", logicConverterPlug2, physiWorld, false, 0);

G4cout << "*** converter tube plug2" << G4endl;
G4cout << " radius " << 50 << G4endl;
G4cout << " length " << 10 << G4endl;
G4cout << " x-position " << -sideoffset - 2*cubesizehalf << G4endl;
G4cout << " y-position " << -tiefe << G4endl;
G4cout << " z-position " << longzero+55  << G4endl;


G4double hLength = 430.0*mm;

// guide from source exit into the sourc
//
//
// e
  rMax = 7.0 * cm;
  rMin = 2.5 * cm;  
  G4Tubs *solidGuideintosource = new G4Tubs("SolidGuide", rMin, rMax, hLength/2., 0., twopi);
  G4LogicalVolume *logicGuideintosource = new  G4LogicalVolume(solidGuideintosource,  GuideMaterial1, "GuideLV");
  G4RotationMatrix * zRot0 = new G4RotationMatrix();
  zRot0->rotateY(3.14159/2*rad);
    logicGuideintosource->SetUserLimits(stepLimit);

  G4VPhysicalVolume *physiGuideintosource = new G4PVPlacement(zRot0, G4ThreeVector(-sideoffset + hLength/2 - cubesizehalf,0.,sideoffset+ 2*cubesizehalf), "GuideIntoSource", logicGuideintosource, physiWorld, false, 0);

 // logicGuideintosource-> SetVisAttributes (LightBleu);

G4cout << "*** tube connecting the second upper cube horizontally to the source exit " << G4endl;
G4cout << " inside radius " << rMin << G4endl;
G4cout << " length " << hLength << G4endl;
G4cout << " x-position " << -sideoffset + hLength/2 - cubesizehalf << G4endl;
G4cout << " y-position " << 0 << G4endl;
G4cout << " z-position " << sideoffset+ 2*cubesizehalf << G4endl;




// vacuum separation foil at  source exit
  rMax = 2.5*cm;
  rMin = 0*cm;
  hLength = 20.0*mm;
  G4Tubs *solidFoil = new G4Tubs("SolidFoil", rMin, rMax, hLength/2, 0., twopi);
  G4LogicalVolume *logicFoil = new  G4LogicalVolume(solidFoil, fVacuum, "Foil");
  G4RotationMatrix * zRotf = new G4RotationMatrix();
  zRotf->rotateZ(3.14159/2*rad);
    zRotf->rotateY(3.14159/2*rad);
          zRotf->rotateX(3.14159/2*rad);
logicFoil->SetUserLimits(stepLimit);
  G4VPhysicalVolume *physiFoil = new G4PVPlacement(zRotf, G4ThreeVector( hLength +cubesizehalf/2,0.,sideoffset+ 2*cubesizehalf), "Foil", logicFoil, physiWorld, false, 0);


G4cout << "*** vacuum separation foil cylinder at  the source exit (placeholder for explicit cacluation)" << G4endl;
G4cout << " radius " << rMax << G4endl;
G4cout << " length " << hLength << G4endl;
G4cout << " x-position " << hLength +cubesizehalf/2 << G4endl;
G4cout << " y-position " << 0 << G4endl;
G4cout << " z-position " << sideoffset+ 2*cubesizehalf << G4endl;

//    logicGuide-> SetVisAttributes (LightGreen);





//// here is the TOF system

// guide from source
  rMax = 4*cm;
  rMin = 2.5*cm;
  hLength = 200.*mm;
  G4Tubs *solidGuide = new G4Tubs("SolidGuide", rMin, rMax, hLength/2, 0., twopi);
  G4LogicalVolume *logicGuide = new  G4LogicalVolume(solidGuide,  GuideMaterial1, "GuideLV");
  G4RotationMatrix * zRot = new G4RotationMatrix();
  zRot->rotateZ(3.14159/2*rad);

    zRot->rotateY(3.14159/2*rad);
          zRot->rotateX(3.14159/2*rad);
logicGuide->SetUserLimits(stepLimit);
  G4VPhysicalVolume *physiGuide = new G4PVPlacement(zRot, G4ThreeVector( hLength +cubesizehalf/2,0.,sideoffset+ 2*cubesizehalf), "GuidePV", logicGuide, physiWorld, false, 0);

G4cout << "*** guide from source exit to bend to vertical TOF" << G4endl;
G4cout << " radius " << rMin << G4endl;
G4cout << " length " << hLength << G4endl;
G4cout << " x-position " << hLength +cubesizehalf/2 << G4endl;
G4cout << " y-position " << 0 << G4endl;
G4cout << " z-position " << sideoffset+ 2*cubesizehalf << G4endl;


    //logicGuide-> SetVisAttributes (LightGreen);
hLength = 200 + 30;

// 90 degree bend downwards
   G4double trMax = 4*cm;
  G4double trMin = 2.5*cm;
  G4double sweptrad = 70*mm;
  G4Torus *solidTorus = new G4Torus("SolidTorus", trMin, trMax, sweptrad, 0., CLHEP::pi/2);
  G4LogicalVolume *logicTorus = new  G4LogicalVolume(solidTorus,  GuideMaterial1, "GuideTorus");
 G4RotationMatrix * zRot1 = new G4RotationMatrix();
 // zRot1->rotateY(3*pi/2);
 // zRot1->rotateY(pi/2);
logicTorus->SetUserLimits(stepLimit);
  G4VPhysicalVolume *physiTorus = new G4PVPlacement(zRot1, G4ThreeVector(hLength +cubesizehalf/2+sweptrad ,-sweptrad, sideoffset+ 2*cubesizehalf ), "GuideTorus", logicTorus, physiWorld, false, 0);

G4cout << "*** quarter of a torus to bend the neutrons downwards to the vertical TOF" << G4endl;
G4cout << " radius " << trMin << G4endl;
G4cout << " swept radius " << sweptrad << G4endl;
G4cout << " x-position " << hLength +cubesizehalf/2+sweptrad << G4endl;
G4cout << " y-position " << -sweptrad << G4endl;
G4cout << " z-position " << sideoffset+ 2*cubesizehalf << G4endl;


///////
// chopper

// should go from - 70 mm (lower end of torus) to -195 (entrance of chopper)

// tube from bend to chopper blades
  G4double RMaxc = 5.5*cm;
  G4double RMinc = 2.5*cm;
  G4double hh = 195. + 575.;// depth below top tube center height
 G4double distc = hh-sweptrad;

  G4Tubs *solidlargeguidec = new G4Tubs("SolidLargeGuidec", RMinc, RMaxc, distc/2, 0., twopi);
  G4LogicalVolume *logiclargeGuidec = new  G4LogicalVolume(solidlargeguidec,  GuideMaterial1, "GuideLVc");
  G4RotationMatrix * zRot3c = new G4RotationMatrix();
logiclargeGuidec->SetUserLimits(stepLimit);
  zRot3c->rotateX(pi/2);
  G4VPhysicalVolume *physiLargeGuidec = new G4PVPlacement(zRot3c, G4ThreeVector(hLength +cubesizehalf/2+2*sweptrad, -hh+distc/2, sideoffset+ 2*cubesizehalf), "GuideLargec", logiclargeGuidec, physiWorld, false, 0);

G4cout << "*** guide from 90 deg bend to chopper blades" << G4endl;
G4cout << " radius " << RMinc << G4endl;
G4cout << " length " << distc << G4endl;
G4cout << " x-position " << hLength +cubesizehalf/2+2*sweptrad << G4endl;
G4cout << " y-position " <<  -hh+distc/2 << G4endl;
G4cout << " z-position " << sideoffset+ 2*cubesizehalf << G4endl;



  G4double HLength = 300.0*mm + 30*mm;

  G4double chthickness = 5.1;
 G4double chwidth = 3.;
 G4double cwidth = 3.;
 G4double chlength = 100.;
 G4double sideoffset2 =  - HLength-2*sweptrad;

  G4Box * solidCh = new G4Box("Chopperblade",chthickness/2*mm,chwidth/2*mm,chlength/2*mm);       //its size
  G4LogicalVolume *logicCh = new  G4LogicalVolume(solidCh,  GuideMaterial5, "Chopperblade");
logicCh->SetUserLimits(stepLimit);
  G4RotationMatrix * zRot009 = new G4RotationMatrix();


 G4VPhysicalVolume *physiCh = new G4PVPlacement(zRot009, G4ThreeVector(-sideoffset2 - 2*cubesizehalf- cubesizehalf,-hh-cwidth/2,longzero+sideoffset+2*cubesizehalf), "chopperfin1", logicCh, physiWorld, false, 0);

G4cout << "***chopper blade 1" << G4endl;
G4cout << " x " << chthickness << G4endl;
G4cout << " y " << chwidth << G4endl;
G4cout << " z " << chlength << G4endl;
G4cout << " x-position " << -sideoffset2 - 2*cubesizehalf- cubesizehalf << G4endl;
G4cout << " y-position " << -hh-cwidth/2 << G4endl;
G4cout << " z-position " << longzero+sideoffset+2*cubesizehalf << G4endl;


 G4VPhysicalVolume *physiCh2 = new G4PVPlacement(zRot009, G4ThreeVector(-sideoffset2 - 2*cubesizehalf- cubesizehalf + chwidth+5.1,-hh-cwidth/2,longzero+sideoffset+2*cubesizehalf), "chopperfin2", logicCh, physiWorld, false, 0);

G4cout << "***chopper blade 2" << G4endl;
G4cout << " x " << chthickness << G4endl;
G4cout << " y " << chwidth << G4endl; 
G4cout << " z " << chlength << G4endl;
G4cout << " x-position " << -sideoffset2 - 2*cubesizehalf- cubesizehalf + chwidth+5.1 << G4endl;
G4cout << " y-position " << -hh-cwidth/2 << G4endl;
G4cout << " z-position " << longzero+sideoffset+2*cubesizehalf << G4endl;


G4VPhysicalVolume *physiCh3 = new G4PVPlacement(zRot009, G4ThreeVector(-sideoffset2 - 2*cubesizehalf- cubesizehalf + (chwidth+5.1)*2,-hh-cwidth/2,longzero+sideoffset+2*cubesizehalf), "chopperfin2", logicCh, physiWorld, false, 0);

G4cout << "***chopper blade 3" << G4endl;
G4cout << " x " << chthickness << G4endl;
G4cout << " y " << chwidth << G4endl;
G4cout << " z " << chlength << G4endl;
G4cout << " x-position " << -sideoffset2 - 2*cubesizehalf- cubesizehalf + (chwidth+5.1)*2<< G4endl;
G4cout << " y-position " << -hh-cwidth/2 << G4endl;
G4cout << " z-position " << longzero+sideoffset+2*cubesizehalf << G4endl;


G4VPhysicalVolume *physiCh4 = new G4PVPlacement(zRot009, G4ThreeVector(-sideoffset2 - 2*cubesizehalf- cubesizehalf + (chwidth+5.1)*3,-hh-cwidth/2,longzero+sideoffset+2*cubesizehalf), "chopperfin2", logicCh, physiWorld, false, 0);

G4cout << "***chopper blade 4" << G4endl;
G4cout << " x " << chthickness << G4endl;
G4cout << " y " << chwidth << G4endl;
G4cout << " z " << chlength << G4endl;
G4cout << " x-position " << -sideoffset2 - 2*cubesizehalf- cubesizehalf + (chwidth+5.1)*3<< G4endl;
G4cout << " y-position " << -hh-cwidth/2 << G4endl;
G4cout << " z-position " << longzero+sideoffset+2*cubesizehalf << G4endl;


G4VPhysicalVolume *physiCh5 = new G4PVPlacement(zRot009, G4ThreeVector(-sideoffset2 - 2*cubesizehalf- cubesizehalf + (chwidth+5.1)*4,-hh-cwidth/2,longzero+sideoffset+2*cubesizehalf), "chopperfin2", logicCh, physiWorld, false, 0);

G4cout << "***chopper blade 5" << G4endl;
G4cout << " x " << chthickness << G4endl;
G4cout << " y " << chwidth << G4endl;
G4cout << " z " << chlength << G4endl;
G4cout << " x-position " << -sideoffset2 - 2*cubesizehalf- cubesizehalf + (chwidth+5.1)*4<< G4endl;
G4cout << " y-position " << -hh-cwidth/2 << G4endl;
G4cout << " z-position " << longzero+sideoffset+2*cubesizehalf << G4endl;


G4VPhysicalVolume *physiCh6 = new G4PVPlacement(zRot009, G4ThreeVector(-sideoffset2 - 2*cubesizehalf- cubesizehalf + (chwidth+5.1)*5,-hh-cwidth/2,longzero+sideoffset+2*cubesizehalf), "chopperfin2", logicCh, physiWorld, false, 0);

G4cout << "***chopper blade 6" << G4endl;
G4cout << " x " << chthickness << G4endl;
G4cout << " y " << chwidth << G4endl;
G4cout << " z " << chlength << G4endl;
G4cout << " x-position " << -sideoffset2 - 2*cubesizehalf- cubesizehalf + (chwidth+5.1)*5<< G4endl;
G4cout << " y-position " << -hh-cwidth/2 << G4endl;
G4cout << " z-position " << longzero+sideoffset+2*cubesizehalf << G4endl;


G4VPhysicalVolume *physiCh7 = new G4PVPlacement(zRot009, G4ThreeVector(-sideoffset2 - 2*cubesizehalf- cubesizehalf + (chwidth+5.1)*6,-hh-cwidth/2,longzero+sideoffset+2*cubesizehalf), "chopperfin2", logicCh, physiWorld, false, 0);

G4cout << "***chopper blade 7" << G4endl;
G4cout << " x " << chthickness << G4endl;
G4cout << " y " << chwidth << G4endl;
G4cout << " z " << chlength << G4endl;
G4cout << " x-position " << -sideoffset2 - 2*cubesizehalf- cubesizehalf + (chwidth+5.1)*6<< G4endl;
G4cout << " y-position " << -hh-cwidth/2 << G4endl;
G4cout << " z-position " << longzero+sideoffset+2*cubesizehalf << G4endl;


G4VPhysicalVolume *physiCh8 = new G4PVPlacement(zRot009, G4ThreeVector(-sideoffset2 - 2*cubesizehalf- cubesizehalf + (chwidth+5.1)*7,-hh-cwidth/2,longzero+sideoffset+2*cubesizehalf), "chopperfin2", logicCh, physiWorld, false, 0);

G4VPhysicalVolume *physiCh9 = new G4PVPlacement(zRot009, G4ThreeVector(-sideoffset2 - 2*cubesizehalf- cubesizehalf + (chwidth+5.1)*8,-hh-cwidth/2,longzero+sideoffset+2*cubesizehalf), "chopperfin2", logicCh, physiWorld, false, 0);

G4cout << "***chopper blade 8" << G4endl;
G4cout << " x " << chthickness << G4endl;
G4cout << " y " << chwidth << G4endl;
G4cout << " z " << chlength << G4endl;
G4cout << " x-position " << -sideoffset2 - 2*cubesizehalf- cubesizehalf + (chwidth+5.1)*7<< G4endl;
G4cout << " y-position " << -hh-cwidth/2 << G4endl;
G4cout << " z-position " << longzero+sideoffset+2*cubesizehalf << G4endl;

G4VPhysicalVolume *physiCh10 = new G4PVPlacement(zRot009, G4ThreeVector(-sideoffset2 - 2*cubesizehalf- cubesizehalf + (chwidth+5.1)*9,-hh-cwidth/2,longzero+sideoffset+2*cubesizehalf), "chopperfin2", logicCh, physiWorld, false, 0);

G4cout << "***chopper blade 9" << G4endl;
G4cout << " x " << chthickness << G4endl;
G4cout << " y " << chwidth << G4endl;
G4cout << " z " << chlength << G4endl;
G4cout << " x-position " << -sideoffset2 - 2*cubesizehalf- cubesizehalf + (chwidth+5.1)*8<< G4endl;
G4cout << " y-position " << -hh-cwidth/2 << G4endl;
G4cout << " z-position " << longzero+sideoffset+2*cubesizehalf << G4endl;




// cone to make 80 mm from 50 mm guide
   G4double crMax = 4.*cm;
  G4double crMin = 4.*cm;
   G4double cRMax = 5.5*cm;
  G4double cRMin = 4*cm;
    G4double cz = 10*cm;
  G4Cons *solidCone = new G4Cons("SolidCone", crMin, crMax, cRMin, cRMax, cz, 0, 2*CLHEP::pi);
  G4LogicalVolume *logicCone = new  G4LogicalVolume(solidCone,  GuideMaterial1, "GuideCone");
 G4RotationMatrix * zRot2 = new G4RotationMatrix();
  zRot2->rotateX(3*pi/2);
  logicCone->SetUserLimits(stepLimit);
  G4VPhysicalVolume *physiCone = new G4PVPlacement(zRot2, G4ThreeVector(hLength +cubesizehalf/2+2*sweptrad,
                          -sweptrad-cz/2-2*trMin-distc-cwidth,
                          sideoffset+ 2*cubesizehalf), "GuideCone", logicCone, physiWorld, false, 0);


G4cout << "*** guide from chopper blade downwards" << G4endl;
G4cout << " radius " << crMin << G4endl;
G4cout << " length " << cz*2 << G4endl;
G4cout << " x-position " << hLength +cubesizehalf/2+2*sweptrad << G4endl;
G4cout << " y-position " << -sweptrad-cz/2-2*trMin-distc-cwidth << G4endl;
G4cout << " z-position " << sideoffset+ 2*cubesizehalf << G4endl;



//    G4VisAttributes* LightGreen = new G4VisAttributes( G4Colour(153/255. ,255/255. ,153/255. ));

  //logicCone -> SetVisAttributes (LightGreen);


// tube downwards to detector
  G4double RMax = 5.5*cm;
  G4double RMin = 4.*cm;

  G4Tubs *solidlargeguide = new G4Tubs("SolidLargeGuide", RMin, RMax, HLength/2, 0., twopi);
  G4LogicalVolume *logiclargeGuide = new  G4LogicalVolume(solidlargeguide,  GuideMaterial1, "GuideLV");
  G4RotationMatrix * zRot3 = new G4RotationMatrix();
  //G4double HLength = 300.0*mm;
logiclargeGuide->SetUserLimits(stepLimit);

  zRot3->rotateX(pi/2);
  G4VPhysicalVolume *physiLargeGuide = new G4PVPlacement(zRot3, G4ThreeVector(hLength +cubesizehalf/2+2*sweptrad,-sweptrad-cz/2-6*trMin-HLength/2-distc-cwidth, sideoffset+ 2*cubesizehalf), "GuideLarge", logiclargeGuide, physiWorld, false, 0);

G4cout << "*** guide further downwards" << G4endl;
G4cout << " radius " << RMin << G4endl;
G4cout << " length " << HLength << G4endl;
G4cout << " x-position " << hLength +cubesizehalf/2+2*sweptrad<< G4endl;
G4cout << " y-position " << -sweptrad-cz/2-6*trMin-HLength/2-distc-cwidth << G4endl;
G4cout << " z-position " << sideoffset+ 2*cubesizehalf << G4endl;



  // extra tube downwards to detector
  //G4double RMax = 5.5*cm;
  //G4double RMin = 4.*cm;

  G4double extralength = 90.*mm;
  G4Tubs *solidlargeguide2 = new G4Tubs("SolidLargeGuide2", RMin, RMax, extralength/2, 0., twopi);
  G4LogicalVolume *logiclargeGuide2 = new  G4LogicalVolume(solidlargeguide2,  GuideMaterial1, "GuideLV2");
  //G4RotationMatrix * zRot3 = new G4RotationMatrix();
  //G4double HLength = 300.0*mm;
logiclargeGuide2->SetUserLimits(stepLimit);

  //zRot3->rotateX(pi/2);
 G4VPhysicalVolume *physiLargeGuide2 = new G4PVPlacement(zRot3, G4ThreeVector(hLength +cubesizehalf/2+2*sweptrad,-sweptrad-cz/2-6*trMin-HLength-distc-cwidth - extralength/2, sideoffset+ 2*cubesizehalf), "GuideLarge2", logiclargeGuide2, physiWorld, false, 0);

G4cout << "*** one more guide further downwards" << G4endl;
G4cout << " radius " << RMin << G4endl;
G4cout << " length " << extralength << G4endl;
G4cout << " x-position " << hLength +cubesizehalf/2+2*sweptrad<< G4endl;
G4cout << " y-position " << -sweptrad-cz/2-6*trMin-HLength-distc-cwidth - extralength/2 << G4endl;
G4cout << " z-position " << sideoffset+ 2*cubesizehalf << G4endl;




//  logiclargeGuide2 -> SetVisAttributes (LightGreen);
  // detector window
  G4double Dmin = 0.*cm;
  G4double Dmax = 5.*cm;
  G4double Dheight = 10.0*mm;
  G4Tubs *solidDet = new G4Tubs("Det", Dmin, Dmax, Dheight, 0., twopi);
  G4LogicalVolume *logicDet = new  G4LogicalVolume(solidDet,  DetMat, "Det");
  G4RotationMatrix * zRot4 = new G4RotationMatrix();
  zRot4->rotateX(pi/2);
  logicDet->SetUserLimits(stepLimit);

  G4VPhysicalVolume *physiDet = new G4PVPlacement(zRot4, G4ThreeVector(hLength +cubesizehalf/2+2*sweptrad,-sweptrad-cz/2-6*trMin-2*HLength/2-Dheight/2-distc-cwidth - extralength, sideoffset+ 2*cubesizehalf), "Det", logicDet, physiWorld, false, 0);


G4cout << "*** Detector windows" << G4endl;
G4cout << " radius " << Dmin << G4endl;
G4cout << " length " << Dheight*2<< G4endl;
G4cout << " x-position " << hLength +cubesizehalf/2+2*sweptrad<< G4endl;
G4cout << " y-position " << -sweptrad-cz/2-6*trMin-2*HLength/2-Dheight/2-distc-cwidth - extralength << G4endl;
G4cout << " z-position " << sideoffset+ 2*cubesizehalf << G4endl;





  //
  // always return the physical World
  //
  return physiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal G4UniformGravityField* ExUCNDetectorConstruction::fField = 0;

void ExUCNDetectorConstruction::ConstructSDandField()
{
  if (!fField) {
    fField = new G4UniformGravityField();

    G4RepleteEofM* equation = new G4RepleteEofM(fField);
    //     G4RepleteEofM* equation = new G4RepleteEofM(fField,12);
    //     G4EqGravityField* equation = new G4EqGravityField(fField);

    G4FieldManager* fieldManager =
      G4TransportationManager::GetTransportationManager()->GetFieldManager();
    fieldManager->SetDetectorField(fField);

    //G4MagIntegratorStepper* stepper = new G4DormandPrince745(equation);
        G4MagIntegratorStepper* stepper = new G4ClassicalRK4(equation,8);
	//G4MagIntegratorStepper* stepper = new G4CashKarpRKF45(equation,8);

    G4double minStep = 0.00001 * mm;
//    G4double minStep = 0.10 * mm;

    G4ChordFinder* chordFinder = new G4ChordFinder((G4MagneticField*)fField, minStep, stepper);

    // Set accuracy parameters
    G4double deltaChord = 0.000010 * mm;
    chordFinder->SetDeltaChord(deltaChord);

    G4double deltaOneStep = 0.00001 * mm;
    fieldManager->SetAccuraciesWithDeltaOneStep(deltaOneStep);

    G4double deltaIntersection = 0.00001 * mm;
    fieldManager->SetDeltaIntersection(deltaIntersection);

    G4TransportationManager* transportManager = G4TransportationManager::GetTransportationManager();

    G4PropagatorInField* fieldPropagator = transportManager->GetPropagatorInField();

    // Dimensionless limits for relative accuracy of integration
    G4double epsMin = 2.5e-6;
    G4double epsMax = 0.01;  // Will soon be maximum without warning.
    // The relative accuracy used for a step of length 'l'
    //                    a.)  epsMin              if deltaOneStep / l < epsMin
    //    epsilon_step =  b.)  epsMax              if deltaOneStep / l > epsMax
    //                    c.)  deltaOneStep / l    otherwise

    fieldPropagator->SetMinimumEpsilonStep(epsMin);
    fieldPropagator->SetMaximumEpsilonStep(epsMax);

    fieldManager->SetChordFinder(chordFinder);
  }
}

// messenger interaction
//
void ExUCNDetectorConstruction::SetProperty1(G4double e){
optpot1 = e;
G4cout << "Detector construction, pass userprop1 from messenger: " << e << G4endl;}
void ExUCNDetectorConstruction::SetProperty2(G4double e){
optpot2 = e;
G4cout << "Detector construction, pass userprop2 from messenger: " << e << G4endl;}

void ExUCNDetectorConstruction::SetProperty3(G4double e){
optpot3 = e;
G4cout << "Detector construction, pass userprop3 from messenger: " << e << G4endl;}

void ExUCNDetectorConstruction::SetProperty4(G4double e){
optpot4 = e;
G4cout << "Detector construction, pass userprop4 from messenger: " << e << G4endl;}

void ExUCNDetectorConstruction::SetProperty5(G4double e){
eta1 = e;
G4cout << "Detector construction, pass userprop5 from messenger: " << e << G4endl;}

void ExUCNDetectorConstruction::SetProperty6(G4double e){
eta2 = e;
G4cout << "Detector construction, pass userprop6 from messenger: " << e << G4endl;}

void ExUCNDetectorConstruction::SetProperty7(G4double e){
eta3 = e;
G4cout << "Detector construction, pass userprop7 from messenger: " << e << G4endl;}
void ExUCNDetectorConstruction::SetProperty8(G4double e){
eta4 = e;
G4cout << "Detector construction, pass userprop8 from messenger: " << e << G4endl;}

void ExUCNDetectorConstruction::SetProperty9(G4double e){
diff1 = e;
G4cout << "Detector construction, pass userprop9 from messenger: " << e << G4endl;}

void ExUCNDetectorConstruction::SetProperty10(G4double e){
diff2 = e;
G4cout << "Detector construction, pass userprop10 from messenger: " << e << G4endl;}


void ExUCNDetectorConstruction::SetProperty11(G4double e){
diff3 = e;
G4cout << "Detector construction, pass userprop11 from messenger: " << e << G4endl;}


void ExUCNDetectorConstruction::SetProperty12(G4double e){
diff4 = e;
G4cout << "Detector construction, pass userprop12 from messenger: " << e << G4endl;}



void ExUCNDetectorConstruction::SetProperty13(G4double e){
valveheight  = e;
G4cout << "Detector construction, pass userprop13 from messenger: " << e << G4endl;}

void ExUCNDetectorConstruction::SetProperty14(G4double e){
maxstp  = e;
G4cout << "Detector construction, pass userprop14 from messenger: " << e << G4endl;}

// SBI v2
void ExUCNDetectorConstruction::SetProperty15(G4double e){
converterRadius  = e;
fConverterRadius = e;
G4cout << "Detector construction, pass userprop15 from messenger: " << e << G4endl;}

void ExUCNDetectorConstruction::SetProperty16(G4double e){
converterLength  = e;
fConverterLength  = e;
G4cout << "Detector construction, pass userprop16 from messenger: " << e << G4endl;}
// End SBI v2


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
