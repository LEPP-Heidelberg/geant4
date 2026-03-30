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
/// \file exoticphysics/ucn/src/ExUCNPrimaryGeneratorAction.cc
/// \brief Implementation of the ExUCNPrimaryGeneratorAction class
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "ExUCNPrimaryGeneratorAction.hh"

#include "ExUCNDetectorConstruction.hh"

#include "G4Event.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include "PrimaryGeneratorMessenger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExUCNPrimaryGeneratorAction::ExUCNPrimaryGeneratorAction(ExUCNDetectorConstruction* DC)
: fDetector(DC)
//ExUCNPrimaryGeneratorAction::ExUCNPrimaryGeneratorAction(ExUCNDetectorConstruction* DC):fDetector(DC)
{
  G4int n_particle = 1;
  fParticleGun = new G4ParticleGun(n_particle);

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  G4ParticleDefinition* particle = particleTable->FindParticle("neutron");
  fParticleGun->SetParticleDefinition(particle);

  // create a messenger for this class
  fGunMessenger = new PrimaryGeneratorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExUCNPrimaryGeneratorAction::~ExUCNPrimaryGeneratorAction()
{
  delete fParticleGun;
  delete fGunMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExUCNPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // This function is called at the begining of event
/*
  fParticleGun->SetParticlePosition(G4ThreeVector(0.0, 0.0, 0.0));
  fParticleGun->SetParticlePolarization(G4ThreeVector(0, 1, 0));

  G4double particleEnergy = 60*1e-9 * eV;// + G4UniformRand() * (1e-7 * eV - 1e-9 * eV);
  fParticleGun->SetParticleEnergy(particleEnergy);

  G4double theta = pi/2;//0;//2 * pi * G4UniformRand();
  G4double phi = 0;// std::acos(1 - 2 * G4UniformRand());
  if (phi > pi / 2 && phi < pi) phi = pi - phi;

  G4double z = 1;//std::sin(phi) * std::cos(theta);
  G4double x = 1;//std::sin(phi) * std::sin(theta);
  G4double y = 1;//std::cos(phi);

  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(x, y, z));

  fParticleGun->GeneratePrimaryVertex(anEvent);
*/
   // There variables read from the detector construction and defined the generator's properties here
   G4double R = fDetector->GetConverterRadius() - 1.5;
   G4double L = fDetector->GetConverterLength()*2;
	
   G4double px=100.0,py=0.0,pz=0.0 ;
   G4ThreeVector spin(1,1,0);
   G4ThreeVector xaxis(1,0,0);
   G4ThreeVector yaxis(0,1,0);

   G4double x1 = G4UniformRand()*72.-R;
   G4double z1  = G4UniformRand()*72.-R;
   G4double radius1 = sqrt(x1*x1+z1*z1);
   while (radius1 > R)
   {
    x1 = G4UniformRand()*72.-R;
    z1  = G4UniformRand()*72.-R;
   radius1 = sqrt(x1*x1+z1*z1);
   }

   G4double longi = G4UniformRand()*(L-10.);
   G4ThreeVector gunPositionShift(x1 - 270 - 80 , z1- 286., -L + longi);
   G4double offset = 0;
   G4double particleEnergy = (Espec(Especfile)+offset) * 1e-9*eV;

   //G4double particleEnergy = gunEnergy_neV * 1e-9 * eV;
   //gunPositionShift.setY(-286.);
  // monoenergetic spectrum
//particleEnergy = 60*1e-15;


   fParticleGun->SetParticlePosition(gunPositionShift);
   fParticleGun->SetParticleEnergy(particleEnergy);
   G4cout << "E = " << particleEnergy << G4endl;
   //G4cout << " energy " << gunEnergy_neV << G4endl;

std::ofstream myfile("start.txt", std::ofstream::app);
//std::ofstream myfile("start.txt");
   //myfile << particleEnergy*1e15 <<  " "  << gunPositionShift.getX() << " " << gunPositionShift.getY()  << " " << gunPositionShift.getZ() <<  std::endl;


  //G4cout << "Distribution of angular distribution following: " << Angspecfile << G4endl;
  /*
   G4ThreeVector startmom(0.,1.,0.);
   // the main UCN beam direction
    G4double theta = CLHEP::pi*(G4UniformRand());//Espec(Angspecfile)/90.*pi/2; // 0 degrees is forward, 90 degrees in transverse
    G4double phi = 2*CLHEP::pi*(G4UniformRand()); // rotationsymmetric distribution
    startmom.rotateX(theta);
    startmom.rotateY(phi);
*/
    G4double mx = 1;
    G4double my = 1;
    G4double mz = 1;

    G4double mradius = sqrt(mx*mx +my*my+mz*mz);
    while (mradius > 1){
      mx = 2 * G4UniformRand()-1;
      my = 2 * G4UniformRand()-1;
      mz = 2 * G4UniformRand()-1;
      mradius = sqrt(mx*mx +my*my+mz*mz);
    } 
    G4ThreeVector startmom(mx,my,mz);
	    
    
    myfile << particleEnergy*1e15 <<  " " << startmom.getX() << " " << startmom.getY() << " " << startmom.getZ() << " " << gunPositionShift.getX() << " " << gunPositionShift.getY()  << " " << gunPositionShift.getZ() <<  std::endl;  
    G4cout << "startmomem " << startmom.getX() << ";" << startmom.getY() << ";" << startmom.getZ() << G4endl;

    fParticleGun->SetParticleMomentumDirection(startmom.unit());
    
    // SBI v2 
    if (tA > 0.0){fParticleGun->SetParticleTime(tA * G4UniformRand() * s);} 
    // End SBI v2	
/*
  if (Angspeclogfile.length() > 0){
         char stri[1000]; sprintf(stri, "%f %f %f %f %f %f %f", particleEnergy*1e9*1e6 ,
         startmom.getX(), startmom.getY(),
         startmom.getZ(),gunPositionShift.getX(),gunPositionShift.getY(), gunPositionShift.getZ());
         std::ofstream examplefile (Angspeclogfile, std::ofstream::app);
         if (examplefile.is_open()){ examplefile << stri << " " ; examplefile.close();}
  }
*/
	


    fParticleGun->GeneratePrimaryVertex(anEvent);



  }

G4double ExUCNPrimaryGeneratorAction::Espec(G4String fnam){

G4double e = G4UniformRand();
   double energy = 0;
   double sum = 0;
   double ref = 0;
   double refold = 0;
   int i = 0;
   int k = 0;
   double entries[1000];
   double velo[1000];

    int n = read_two_column_file(fnam, velo, entries, 1000);

     ref = 0;
     sum = 0;
       for (int a = 0; a <n;a++){
       sum += entries[a];
       		}

     for (k =0;k<n; k++){
        refold = ref;
        ref = ref + entries[k] / sum;

        if (e < ref){
         energy = (velo[k]-velo[k-1])/(ref-refold) * (e-refold) + velo[k-1];

        return energy;
        }
     }
     //}

   //}
   return e;
}

int ExUCNPrimaryGeneratorAction::read_two_column_file(const char *filename, double x[], double y[], size_t maxSize) {
        FILE *fp = fopen(filename, "r");
    if (!fp) {
        perror("Spectrum file issue");
        return -1;
    }

    size_t count = 0;
    char line[256];

    while (fgets(line, sizeof(line), fp)) {
        // Kommentarzeilen überspringen
        if (line[0] == '#') {
            continue;
        }
        double val1, val2;
        if (sscanf(line, "%lf %lf", &val1, &val2) == 2) {
            if (count < maxSize) {
                x[count] = val1;
                y[count] = val2;
                count++;
            //G4cout << " v 1 " << val1 << ", v2 " << val2 << G4endl;
            } else {
                fprintf(stderr, "Maximale Grö?~_e erreicht (%zu Werte)\n", maxSize);
                break;
            }
        }
    }
    fclose(fp);
    return (int)count;
}


void ExUCNPrimaryGeneratorAction::SetGunEnergy_neV(G4double e)
{
  gunEnergy_neV = e ;
  G4cout << " Energy the gun = " << gunEnergy_neV << " neV" << G4endl;
}

// SBI v2
// This function takes a user input to uniformly add a number betwenn [0, tA] to the particle time
void ExUCNPrimaryGeneratorAction::SettA(G4double e)
{
	tA = e;
	G4cout << "Note: Particles will be initialised with times uniformly spread over [0, " << tA << "]." << G4endl;
} 
//

void ExUCNPrimaryGeneratorAction::SetEspec(G4String fil){
Especfile = fil;
//G4cout << " set energy spec " << fil << G4endl;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
