//////////////////////////////////////////////////////////////////////
//  Material Boundary Process for ultracold neutrons
//  9.9.04 peter fierlinger
////////////////////////////////////////////////////////////////////////

#include "G4ios.hh"
#include "UCNMaterialBoundary2.hh"
#include "UCNMaterialBoundary2Messenger.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4StepPoint.hh"
#include "G4TransportationManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#define neV (1.0e-9*eV)
#include <fstream>
#include "G4SystemOfUnits.hh"

G4bool isExcited = false;
G4double chopper_t = 0;
G4double chopper_x = 0;
G4double chopper_y = 0;
G4double chopper_z = 0;
G4double chopper_px = 0;
G4double chopper_py = 0;
G4double chopper_pz = 0;
G4double st = 0;
G4double foilct =0 ; 
G4double ezero = 0;
G4double num = 0;
G4double justprinted = 0;
G4double excitedFrac = 0;
G4double E_recoil = 0;
G4double meanHeLifetime = 500.;
G4double meanDownShiftLifetime = 0.000000000001;
G4double meanLifetime = 880.; 
G4double inside_he = 0;
G4double foilabs = 0;
G4double timesteps = 0;
G4double stepsbeforekill = 100000000000;
G4double maxGlobalTime = 0.0;
G4double downShiftLifetime = 10000000000000.;
G4double betadecaylifetime = 10000000000000.;
G4double heliumlifetime = 10000000000000.;
G4String outputfile1 = "default1.txt";
G4String outputfile2 = "default2.txt";
G4String outputfile3 = "default3.txt";
G4double remember_random = 0;
//std::ofstream myfile(outputfile2, std::ofstream::app);

/**
 * Add a state (0 - closed, 1 - open) at a given time
 * keeps the array sorted so you don't have to
 */
void UCNShutterStates::AddState(int state, float time) {
  if (n >= UCN_SHUTTERS_STATES_MAX) {
    G4cerr << "Shutter states maximum reached" << G4endl;
  }
  // keep array in order
  int i, j;
  for (i=0; i<n; i++) {
    if (time < times[i]) {
      break;
    }
  }
  // shift tail of array
  for (j=n-1; j>=i; j--) {
    states[j+1] = states[j];
    times[j+1] = times[j];
  }
  // insert
  states[i] =  state;
  times[i] = time;
  n++;
}

/**
 * Get shutter state (0 - closed, 1 - open) at a given time
 * default value if not in definite state
 */
int UCNShutterStates::GetState(float time) {
  int state = UCN_SHUTTERS_DEFAULT_STATE;
  if (loop && n > 0) {
    // in case of loop, offset time (mod)
    float last_t = times[n-1];
    if (time > last_t) {
      int m = floor(time/last_t);
      time -= m * last_t;
      // set to last state in loop
      state = states[n-1];
    }
  }
  for (int i=0; i<n; i++) {
    if (time >= times[i]) {
        state = states[i];
    }
    else {
      break;
    }

  }
  return state;
}

void UCNShutterStates::SetLoop(int loop) {
  this->loop = loop;
}

UCNMaterialBoundary2* UCNMaterialBoundary2::theInstance = 0;

UCNMaterialBoundary2::UCNMaterialBoundary2(const G4String& processName, G4ProcessType type)
         : G4VDiscreteProcess(processName, type)
{
  theInstance = this;
 just_reflected = 0;
 fermipotdiff_shutter = 1e300;
 return_it = 0;
 sh_fermipot = 0;sh_spinflip = 0;sh_loss = 0;sh_diffusion = 0;sh_reflectivity = 0;
 sh_abscs = 0;sh_losscs = 0;sh_scatcs = 0;
 return_it = 0;
 useshutters = 0;
 foil = 0;

 eta1 = 0; eta2 = 0; eta3 = 0; opt1 = 0; opt2 = 0; opt3 = 0; C = 0; V = 0; f_a = 0; f_h = 0; lambda = 0; Gamma1a = 0; Gamma1h = 0; Gamma2a = 0; Gamma2h = 0; Gamma3a = 0; Gamma3h = 0; u_min = 0; u_max = 0;

 // create a messenger for this class
  theMessenger = new UCNMaterialBoundary2Messenger(this);

}

UCNMaterialBoundary2::~UCNMaterialBoundary2(){

delete theMessenger;
} 

G4VParticleChange* UCNMaterialBoundary2::PostStepDoIt(const G4Track& aTrack, 
const G4Step& aStep)
{
if (aTrack.GetCurrentStepNumber() == 1) {
st = 0;
chopper_x = 0;
chopper_y = 0;
chopper_z = 0;
chopper_px = 0;
chopper_py = 0;
chopper_pz = 0;
chopper_t = 0;
foilct = 0;
justprinted = 0;
num++;
ezero = aTrack.GetKineticEnergy();

/*
G4cout << "Theory parameters for the source are " << G4endl;
G4cout << "opt1     " << opt1 << G4endl;
G4cout << "opt2     " << opt2 << G4endl;
G4cout << "opt3     " << opt3 << G4endl;
G4cout << "eta1     " << eta1 << G4endl;
G4cout << "eta2     " << eta2 << G4endl;
G4cout << "lambda   " << lambda << G4endl;
G4cout << "f_a      " << f_a << G4endl;
G4cout << "f_h      " << f_h << G4endl;
G4cout << "u_min    " << u_min << G4endl;
G4cout << "u_max    " << u_max << G4endl;
G4cout << "Gamma1a  " << Gamma1a << G4endl;
G4cout << "Gamma1h  " << Gamma1h << G4endl;
G4cout << "Gamma2a  " << Gamma2a << G4endl;
G4cout << "Gamma2h  " << Gamma2h << G4endl;
G4cout << "Gamma3a  " << Gamma3a << G4endl;
G4cout << "Gamma3h  " << Gamma3h << G4endl;
G4cout << "C        " << C << G4endl;
G4cout << "V        " << V << G4endl;
*/


}


// dump track and energy infos in output file in specified time steps
//

if (timesteps > 0){ // if timesteps is set via macro command
// we try milliseconds
if (aTrack.GetGlobalTime() *1e-6  > justprinted){
	G4StepPoint* pPreStepPoint  = aStep.GetPreStepPoint();
  G4StepPoint* pPostStepPoint = aStep.GetPostStepPoint();

justprinted = (int)(aTrack.GetGlobalTime()*1e-6+ timesteps);
//G4cout <<  "**t = " << aTrack.GetGlobalTime() *1e-6 << ", "  << (int)(aTrack.GetGlobalTime() *1e-6)  << G4endl;

std::ofstream myfile(outputfile2, std::ofstream::app);
  myfile << aTrack.GetKineticEnergy()/neV << " " << aTrack.GetMomentumDirection().getX() << " " <<  aTrack.GetMomentumDirection().getY() << " " <<
                  aTrack.GetMomentumDirection().getZ() << " " << pPreStepPoint->GetPosition().getX() << " " << 
		  pPreStepPoint->GetPosition().getY() << " " << pPreStepPoint->GetPosition().getZ() << " " << aTrack.GetGlobalTime()*1e-6  << std::endl;
}	

}



//G4cout << ", t = " << aTrack.GetGlobalTime() *1e-6 << ", step = " << aTrack.GetCurrentStepNumber()  << ", step length " << aTrack.GetStepLength() <<  G4endl;
/*
    if (aTrack.GetCurrentStepNumber() < stepsbeforekill) {
G4StepPoint* pPreStepPoint  = aStep.GetPreStepPoint();
  G4StepPoint* pPostStepPoint = aStep.GetPostStepPoint();

//	G4cout << "UCNMaterialBoundary2 record first steps" << aTrack.GetCurrentStepNumber() << ", " << stepsbeforekill <<  G4endl;
  std::ofstream myfile(outputfile1, std::ofstream::app);
  myfile << aTrack.GetKineticEnergy()/neV << " " << aTrack.GetMomentumDirection().getX() << " " <<  aTrack.GetMomentumDirection().getY() << " " <<
                  aTrack.GetMomentumDirection().getZ() << " " << pPreStepPoint->GetPosition().getX() << " " << pPreStepPoint->GetPosition().getY() << " " << pPreStepPoint->GetPosition().getZ() << " " << aTrack.GetGlobalTime() << " " <<
                  chopper_t << " " << chopper_x << " " << foilct << " " << chopper_z << " " << chopper_px << " " << chopper_py << " " << "10000" << std::endl;

}
*/


 //G4cout << "steps " << aTrack.GetCurrentStepNumber() << ", time "
   //      << aTrack.GetProperTime()  << G4endl;


if (aTrack.GetCurrentStepNumber() == stepsbeforekill) {
G4StepPoint* pPreStepPoint  = aStep.GetPreStepPoint();
  G4StepPoint* pPostStepPoint = aStep.GetPostStepPoint();

 G4cout << "Reached maximum number of steps " << aTrack.GetCurrentStepNumber() << ", time "
	 << aTrack.GetGlobalTime()  << G4endl;
  aParticleChange.Initialize(aTrack);
 aParticleChange.ProposeTrackStatus(fStopAndKill);
 return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);

}



if (aTrack.GetGlobalTime()*1e-9 > maxGlobalTime) {

G4StepPoint* pPreStepPoint  = aStep.GetPreStepPoint();
  G4StepPoint* pPostStepPoint = aStep.GetPostStepPoint();

 G4cout << "Reached defined max age: " << maxGlobalTime << "s, at Step Number: " << aTrack.GetCurrentStepNumber() << ", global time "
         << aTrack.GetGlobalTime() *1e-9 << "s" <<G4endl;
  aParticleChange.Initialize(aTrack);
 
      std::ofstream myfile(outputfile1, std::ofstream::app);
          myfile << aTrack.GetKineticEnergy()/neV << " " << aTrack.GetMomentumDirection().getX() << " " <<  aTrack.GetMomentumDirection().getY() << " " <<
                  aTrack.GetMomentumDirection().getZ() << " " << pPreStepPoint->GetPosition().getX() << " " << pPreStepPoint->GetPosition().getY() << " " << pPreStepPoint->GetPosition().getZ() << " " << aTrack.GetGlobalTime() << " " <<
                  chopper_t << " " << chopper_x << " " << foilct << " " << chopper_z << " " << chopper_px << " " << chopper_py << " " << "13000" << std::endl;

  
  
  aParticleChange.ProposeTrackStatus(fStopAndKill);
 return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);

}


// neutron life-time
if (aTrack.GetCurrentStepNumber() == 1) {
betadecaylifetime = -std::log(1.0 - G4UniformRand()) * meanLifetime;
isExcited = false; 
}

// excited_state
if ((aTrack.GetCurrentStepNumber() == 1) && (G4UniformRand() < excitedFrac) ) {
downShiftLifetime = -std::log(1.0 - G4UniformRand()) * meanDownShiftLifetime;
isExcited = true;
}

// energy independent loss time inside helium
if (aTrack.GetCurrentStepNumber() == 1) {
heliumlifetime = -std::log(1.0 - G4UniformRand()) * meanHeLifetime + aTrack.GetGlobalTime()*1e-9;
}


  aParticleChange.Initialize(aTrack);
  verboselevel=0;//6;
  verboseLevel=0;
  G4StepPoint* pPreStepPoint  = aStep.GetPreStepPoint();
  G4StepPoint* pPostStepPoint = aStep.GetPostStepPoint();
  G4String volnam1 = pPreStepPoint->GetPhysicalVolume()->GetName();
  G4String volnam2 = pPostStepPoint->GetPhysicalVolume()->GetName();
  G4String volnam1a = pPreStepPoint->GetPhysicalVolume()->GetName();
  G4String volnam2a = pPostStepPoint->GetPhysicalVolume()->GetName();




  // leaving and entering helium
int usehe = 1;

if (usehe == 1){


  if (volnam2 == "DownwardTubeempty" && volnam1 == "World"){
//G4cout << " entering v2post = " << volnam2 << ", leave v1pre = " << volnam1 <<  ", p = " << pPostStepPoint->GetPosition() << G4endl; 
  G4ThreeVector theGlobalPoint = pPostStepPoint->GetPosition();
  G4ThreeVector thePrePoint = pPreStepPoint->GetPosition();
   G4Navigator* theNavigator =
    G4TransportationManager::GetTransportationManager()->
    GetNavigatorForTracking();
  G4ThreeVector theLocalPoint = theNavigator->
    GetGlobalToLocalTransform().
    TransformPoint(theGlobalPoint);
  G4ThreeVector theLocalNormal; // Normal points back into volume
  G4bool valid;
  theLocalNormal = -theNavigator->GetLocalExitNormal(&valid);
  //G4cout << " loc norm " << theLocalNormal << G4endl;
  if (theLocalNormal.getZ() == -1){
     //G4cout << "leaving helium" << G4endl;	  
     // reset helium lifetime to infinity
    meanHeLifetime = 100000000000000.;
    // add helium opt pot to vertical velocity component and set new energy
    G4double vel = aTrack.GetVelocity() * 1e6;
    //G4cout << "vel = " << vel << G4endl;
    G4ThreeVector momdir = aTrack.GetMomentumDirection();
    //G4cout << "momdir " << momdir << G4endl;
    G4ThreeVector vert(0,1,0);
    G4double tfakt = momdir.dot(vert);
    //G4cout << "tfakt " << tfakt << G4endl;
    G4double velnorm = vel * tfakt;
    //G4cout << "velnorm " << velnorm << G4endl;
    G4double enormal = vel*vel/2.* 1.6749e-27; // 1.6749e-27; CLHEP::neutron_mass_c2/neV;
    // 18.5 nano electron volts (n eV) is equal to 2.9637 x 10^-27 joules.
    G4double enormal_new = enormal + 2.9639e-27; // joules
    //G4cout << "enormal " << enormal  * 6.242e27 << G4endl;
    //G4cout << "enormalnew " << enormal_new * 6.242e27<< G4endl;
    G4double velnorm_new = sqrt(2.*enormal_new/1.6749e-27);
    //G4cout << "velnorm_new " << velnorm_new << G4endl;
    G4double tfaktnew = velnorm_new / velnorm;
    G4ThreeVector mom = aTrack.GetMomentum();
    //G4cout << "tfaktnew " << tfaktnew << G4endl;
    mom.setY(mom.getY()*tfaktnew);
    //G4cout << "mom dir new " << mom << G4endl;
    mom.unit();
    G4double energy = aTrack.GetKineticEnergy()/neV;
    //G4cout << "energy " << energy << G4endl;
    aParticleChange.ProposeMomentumDirection(mom.unit());
    aParticleChange.ProposeEnergy((energy + 18.5)*neV);
  }	 
  }


  if (volnam1 == "DownwardTubeempty" && volnam2 == "World"){
  //G4cout << " entering v2post = " << volnam2 << ", leave v1pre = " << volnam1 <<  ", p = " << pPostStepPoint->GetPosition() << G4endl;

	  G4ThreeVector theGlobalPoint = pPostStepPoint->GetPosition();
  G4ThreeVector thePrePoint = pPreStepPoint->GetPosition();
   G4Navigator* theNavigator =
    G4TransportationManager::GetTransportationManager()->
    GetNavigatorForTracking();
  G4ThreeVector theLocalPoint = theNavigator->
    GetGlobalToLocalTransform().
    TransformPoint(theGlobalPoint);
  G4ThreeVector theLocalNormal; // Normal points back into volume
  G4bool valid;
  theLocalNormal = -theNavigator->GetLocalExitNormal(&valid);
  //G4cout << " loc norm " << theLocalNormal << G4endl;
  if (theLocalNormal.getY() > 0.999){
    //G4cout << "re entering helium" << G4endl;
      // new helium life time
    heliumlifetime = -std::log(1.0 - G4UniformRand()) * meanHeLifetime + aTrack.GetGlobalTime()*1e-9;
// remove helium opt pot to vertical velocity component and set new energy
    G4double vel = aTrack.GetVelocity() * 1e6;
    //G4cout << "vel = " << vel << G4endl;
    G4ThreeVector momdir = aTrack.GetMomentumDirection();
    //G4cout << "momdir " << momdir << G4endl;
    G4ThreeVector vert(0,1,0);
    G4double tfakt = momdir.dot(vert);
    //G4cout << "tfakt " << tfakt << G4endl;
    G4double velnorm = vel * tfakt; 
    //G4cout << "velnorm " << velnorm << G4endl; 
    G4double enormal = vel*vel/2.* 1.6749e-27; // CLHEP::neutron_mass_c2/neV;
    // 18.5 nano electron volts (n eV) is equal to 2.9637 x 10^-27 joules.
    G4double enormal_new = enormal - 2.9639e-27; // joules
    //G4cout << "enormal " << enormal  * 6.242e27 << G4endl;
    //G4cout << "enormalnew " << enormal_new * 6.242e27<< G4endl;
    if (enormal_new < 0){
      // reflect from surface
        G4ThreeVector mom = aTrack.GetMomentum();
        mom.setY(-mom.getY());
        aParticleChange.ProposeMomentumDirection(mom.unit());
    }
    else{
    G4double velnorm_new = sqrt(2.*enormal_new/1.6749e-27);
    //G4cout << "velnorm_new " << velnorm_new << G4endl;
    G4double tfaktnew = velnorm_new / velnorm;
    G4ThreeVector mom = aTrack.GetMomentum();
    //G4cout << "tfaktnew " << tfaktnew << G4endl;
    mom.setY(mom.getY()*tfaktnew);
    //G4cout << "mom dir new " << mom << G4endl;
    mom.unit();
    G4double energy = aTrack.GetKineticEnergy()/neV;
    //G4cout << "energy " << energy << G4endl;
    aParticleChange.ProposeMomentumDirection(mom.unit());
    aParticleChange.ProposeEnergy((energy - 18.5)*neV);
    }
  
  
  }
  }


} // use helium


/// foil
int usefoil = 0;

///*
// passage through the vacuum separation foil

if ( pPreStepPoint->GetPosition().getX() > 150 ){
         if (foil < 150){
//		 G4cout << "through foil outward" << G4endl;
		 foil = 10000; // pPreStepPoint->GetPosition().getX();
usefoil = 1;
	 }
}
else{
if (foil > 150){
  //       G4cout << "thoguht foil inward " << G4endl;
         foil = 0;
usefoil = 1;
}
}

if (usefoil == 1){
foilct++;
// empiricial polypropylene foil transmission for 4 micro m foil thickness
G4double mx = aTrack.GetMomentumDirection().getX();
G4double my = aTrack.GetMomentumDirection().getY();
G4double mz = aTrack.GetMomentumDirection().getZ();
G4double mtot = sqrt(mx*mx+my*my+mz*mz);
G4double rel = std::abs(mtot / mx);
G4double v = aTrack.GetVelocity()*1e6;
//G4cout << " velocity " << v << ", mtot " << mtot << ", mx " << mx << ", rel " << rel  << G4endl;
G4double dx = 0.002 / rel; // in mm
G4double sig = foilabs;// 0.3; //// 0.3 for pp
G4double mu = sig * 2200 / v;
G4double foiltrans = exp(-mu*dx);

//G4cout << "PP foil transmission " << foiltrans << G4endl;

G4double counted = G4UniformRand();
if (counted > foiltrans){
G4cout << "Killed pp foil" << foiltrans << ", rnd=  " << counted << G4endl;

	std::ofstream myfile(outputfile1, std::ofstream::app);
          myfile << aTrack.GetKineticEnergy()/neV << " " << aTrack.GetMomentumDirection().getX() << " " <<  aTrack.GetMomentumDirection().getY() << " " <<
                  aTrack.GetMomentumDirection().getZ() << " " << pPreStepPoint->GetPosition().getX() << " " << pPreStepPoint->GetPosition().getY() << " " << pPreStepPoint->GetPosition().getZ() << " " << aTrack.GetGlobalTime() << " " <<
                  chopper_t << " " << chopper_x << " " << foilct << " " << chopper_z << " " << chopper_px << " " << chopper_py << " " << "1000" << std::endl;

  aParticleChange.ProposeTrackStatus( fStopAndKill ) ;
 return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
}
}
//

if ( pPreStepPoint->GetPosition().getY() < -3500 ){

        // falling all the way down, passing stacked turtles 
        //
//G4StepPoint* pPreStepPoint  = aStep.GetPreStepPoint();
G4cout << "Falling down" << G4endl;

              std::ofstream myfile(outputfile1, std::ofstream::app);
          myfile << aTrack.GetKineticEnergy()/neV << " " << aTrack.GetMomentumDirection().getX() << " " <<  aTrack.GetMomentumDirection().getY() << " " << 
                  aTrack.GetMomentumDirection().getZ() << " " << pPreStepPoint->GetPosition().getX() << " " << pPreStepPoint->GetPosition().getY() << " " << pPreStepPoint->GetPosition().getZ() << " " << aTrack.GetGlobalTime() << " " <<
                  chopper_t << " " << chopper_x << " " << foilct << " " << chopper_z << " " << chopper_px << " " << chopper_py << " " << "15000" << std::endl;

             
      
  aParticleChange.ProposeTrackStatus(fStopAndKill);
 return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
    
}




//// chopper


if (  pPreStepPoint->GetPosition().getX() > 150 ){
//G4cout << " 0at TOF system start " << " t = " << chopper_t << G4endl;

        if (  pPreStepPoint->GetPosition().getY() < -195. - 575.){
//G4cout << " 1at TOF system start " << " t = " << chopper_t << G4endl;

                if ( st == 0 ){

 chopper_t =  aTrack.GetGlobalTime() ;
 chopper_x =  pPreStepPoint->GetPosition().getX() ;
 chopper_y =  pPreStepPoint->GetPosition().getY() ;
 chopper_z =  pPreStepPoint->GetPosition().getZ() ;
 chopper_px = aTrack.GetMomentumDirection().getX();
 chopper_py = aTrack.GetMomentumDirection().getY();
 chopper_pz = aTrack.GetMomentumDirection().getZ();
st = 1;
 //G4cout << " 2at TOF system start " << " t = " << chopper_t << G4endl;
        }
        }
}

if (  pPreStepPoint->GetPosition().getX() > 150 ){
//G4cout << " 0at TOF system start " << " t = " << chopper_t << G4endl;
        if (  pPreStepPoint->GetPosition().getY() > -175 && st != 0){
        st = 0;
	}
}

  // make sure to be at a geometrical boundary
  if (pPostStepPoint->GetStepStatus() != fGeomBoundary) {
    just_reflected = 0;
   return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
  }

  G4Material* Material1 = pPreStepPoint->GetPhysicalVolume()->GetLogicalVolume()->GetMaterial();
  G4Material* Material2 = pPostStepPoint->GetPhysicalVolume()->GetLogicalVolume()->GetMaterial();
 
  
//if ( aTrack.GetGlobalTime()*1e-9 > betadecaylifetime){
if ( (aTrack.GetGlobalTime() * 1e-9 > betadecaylifetime) && (isExcited == false)){ 
//G4cout << "beta decay " << G4endl;

std::ofstream myfile(outputfile1, std::ofstream::app);
G4cout << "beta decay " << chopper_t << ", " << aTrack.GetGlobalTime() << ", " << aTrack.GetKineticEnergy()/neV << ", " << aTrack.GetMomentumDirection() << ", " << pPreStepPoint->GetPosition() << G4endl;

myfile << aTrack.GetKineticEnergy()/neV << " " << aTrack.GetMomentumDirection().getX() << " " <<  aTrack.GetMomentumDirection().getY() << " " <<  aTrack.GetMomentumDirection().getZ() << " " << pPreStepPoint->GetPosition().getX() << " " << pPreStepPoint->GetPosition().getY() << " " << pPreStepPoint->GetPosition().getZ() << " " << aTrack.GetGlobalTime() <<" " <<  chopper_t << " " << chopper_x << " " << foilct << " " << chopper_z << " " << chopper_px << " " << chopper_py << " " << "11000" << std::endl;

aParticleChange.ProposeTrackStatus( fStopAndKill ) ;
return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);


}

//downshift stuff
if ( (aTrack.GetGlobalTime()*1e-9 > downShiftLifetime) && (isExcited == true) && (Material1 != Material2)){
	// Remove excited state glag 
	isExcited = false;

	// Reset the beta decay life clock to begin at de-excitation
       	betadecaylifetime += aTrack.GetGlobalTime()*1e-9;	
	
	//Add momentum stuff here
	G4double NeutronMass = G4Neutron::Neutron()->GetPDGMass();
	G4double Energy = aTrack.GetKineticEnergy();
	G4ThreeVector Momentum = aTrack.GetMomentum();
	G4double PhotonEnergy = E_recoil*eV; // This should be in a macro, 5eV -> 13.29neV

	//Recoil Energy
	G4double RecoilEnergy = (PhotonEnergy * PhotonEnergy) / (2.0 * NeutronMass); // This should be the energy imparted to the neutron
	//NewEnergy = Energy + RecoilEnergy; // I think this is wrong -- Husain 

	//Momentum Direction sampling math using random variables x1 x2, should be uniformly sampled in 3D space
	G4double x1 = G4UniformRand();
	G4double x2 = G4UniformRand(); 

	G4double cosTheta = 1.0 - 2.0 * x1;
	G4double sinTheta = std::sqrt(1.0 - cosTheta * cosTheta);
	G4double phi = CLHEP::twopi * x2;

	// Create the momentum vector in xyz from the random sampling (with length of RecoilEnergy, RecoilMomentum's magnitude is same as PhotonEnergy
	G4double px = PhotonEnergy * std::cos(phi) * sinTheta;
	G4double py = PhotonEnergy * std::sin(phi) * sinTheta;
	G4double pz = PhotonEnergy * cosTheta;
	G4ThreeVector RecoilMomentum = G4ThreeVector(px,py,pz);

	G4cout << "The recoil momentum magnitude: " << RecoilMomentum.mag() << G4endl;
	G4cout << "The recoil momentum direction: " << RecoilMomentum << G4endl;
	G4cout << "The initial momentum direction: " << Momentum << G4endl;
	// Calculate recoiled momentum
	G4ThreeVector NewMomentum = Momentum - RecoilMomentum;
	G4cout << "The new momentum direction: " << NewMomentum << G4endl;

	// Normalize NewMomentum, then attach new energy and momentum to tracking
	G4double NewEnergy = NewMomentum.mag2() / (2.0 * NeutronMass);
	G4ThreeVector NewMomentum_full = NewMomentum; // This is a copy for diagnostics
	NewMomentum /= NewMomentum.mag();
	
	// Record initial state for diagnostics
	const G4double E_before = aTrack.GetKineticEnergy();
	const G4ThreeVector p_before = aTrack.GetMomentum();
	const G4ThreeVector dir_before = aTrack.GetMomentumDirection();

	// Now Propose the new calculated unit vector momentum and kinetic energy to the track	
	aParticleChange.Initialize(aTrack);
	aParticleChange.ProposeMomentumDirection(NewMomentum);
	aParticleChange.ProposeEnergy(NewEnergy);
	
	// Calculate changes
	const G4double p_before_mag = p_before.mag();
	const G4double p_after_mag  = NewMomentum_full.mag();
	const G4double E_recoil_expected = (PhotonEnergy * PhotonEnergy) / (2.0 * NeutronMass);
	
	// Print diagnostic information and return
	G4cout
  	<< "\n[Downshift recoil applied]\n"
  	<< "  Time                : " << aTrack.GetGlobalTime()/ns << " ns\n"
  	<< "  Photon energy       : " << PhotonEnergy/eV << " eV\n"
  	<< "  Expected recoil E   : " << E_recoil_expected/neV << " neV\n"
  	<< "\n"
  	<< "  Kinetic energy (before) : " << E_before/neV << " neV\n"
  	<< "  Kinetic energy (after)  : " <<  NewEnergy/neV << " neV\n"
  	<< "  Delta E                 : " << (NewEnergy - E_before)/neV << " neV\n"
  	<< "\n"
  	<< "  |p| before : " << p_before_mag/(eV) << " eV\n"
  	<< "  |p| after  : " << p_after_mag/(eV)  << " eV\n"
  	<< "\n"
  	<< "  Direction before : " << dir_before << "\n"
  	<< "  Direction after  : " << NewMomentum << "\n"
  	<< G4endl;
	
	//G4cout
        //<< "\n[Downshift recoil applied**************************]\n" << G4endl;
	//return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
	//return &aParticleChange;
}



if ( aTrack.GetGlobalTime()*1e-9 > heliumlifetime){
        //G4cout << "helium loss " << G4endl;

std::ofstream myfile(outputfile1, std::ofstream::app);
G4cout << "Lost in helium " << chopper_t << ", " << aTrack.GetGlobalTime() << ", " << aTrack.GetKineticEnergy()/neV << ", " << aTrack.GetMomentumDirection() << ", " << pPreStepPoint->GetPosition() << G4endl;

myfile << aTrack.GetKineticEnergy()/neV << " " << aTrack.GetMomentumDirection().getX() << " " <<  aTrack.GetMomentumDirection().getY() << " " <<  aTrack.GetMomentumDirection().getZ() << " " << pPreStepPoint->GetPosition().getX() << " " << pPreStepPoint->GetPosition().getY() << " " << pPreStepPoint->GetPosition().getZ() << " " << aTrack.GetGlobalTime() <<" " <<  chopper_t << " " << chopper_x << " " << foilct << " " << chopper_z << " " << chopper_px << " " << chopper_py << " " << "12000" << std::endl;

aParticleChange.ProposeTrackStatus( fStopAndKill ) ;
return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);


}


if ( volnam1 == "Det"  ){
            G4cout << " at TOP OF TOF " << chopper_t << ", " << aTrack.GetGlobalTime() << ", " << aTrack.GetKineticEnergy()/neV << ", " << aTrack.GetMomentumDirection() << ", " << pPreStepPoint->GetPosition() << G4endl;

 G4ThreeVector theGlobalNormal(0,1,0);
  G4double vel = aTrack.GetVelocity();
  G4ThreeVector momdir = aTrack.GetMomentumDirection();
  G4ThreeVector mom = aTrack.GetMomentum();
  G4double tfakt = momdir.dot(theGlobalNormal);
  G4double momnorm = mom.dot(theGlobalNormal);
  G4double velnorm = vel * tfakt;
  G4double enormal = momnorm*momnorm/2./CLHEP::neutron_mass_c2/neV;
  G4double energy = aTrack.GetKineticEnergy()/neV;


            // calculate backreflection probability at top of detector entrance foil
   // calculate loss and reflectivity june 27 2025
   G4double W = 1e-3 * 54;
   G4double alphaW = sqrt((54-enormal)*(54-enormal)+W*W);
   G4double part1 = enormal - sqrt(enormal)*sqrt(2*alphaW - 2*(54-enormal)) + alphaW;
   G4double part2 = enormal + sqrt(enormal)*sqrt(2*alphaW - 2*(54-enormal)) + alphaW;
   G4double R_sq = (part1)/(part2);
   G4double reflection = R_sq;
   G4double rnd1 = G4UniformRand();
   if (rnd1 <  reflection){   // reflect it 
  //    G4cout << "reflect it! " << G4endl;
      G4ThreeVector ref = reflect(0, 0, momdir, theGlobalNormal);
      ref = ref.unit();
      aParticleChange.ProposeMomentumDirection(ref);
      just_reflected = 1;
}  


// AlMG3 foil transmission for 0.1 mm foil thickness
G4double my = aTrack.GetMomentumDirection().getY();
G4double v = aTrack.GetVelocity()*1e6-my;
            G4cout << "vertical velocity " << aTrack.GetVelocity()*1e6 << " my " << my << " vy " << v << G4endl;
G4double foiltrans = 0;

foiltrans = exp(-0.627043/v);
//  =-0.0002*v*v*v*v + 0.0058*v*v*v - 0.0675*v*v + 0.3474*v + 0.2076;
G4cout << "Alu foil transmission v = " << v << ", " << foiltrans << G4endl;



G4double counted0 = G4UniformRand();
if (counted0 > foiltrans){ 

std::ofstream myfile(outputfile1, std::ofstream::app);
myfile << aTrack.GetKineticEnergy()/neV << " " << aTrack.GetMomentumDirection().getX() << " " <<  aTrack.GetMomentumDirection().getY() << " " <<  aTrack.GetMomentumDirection().getZ() << " " << pPreStepPoint->GetPosition().getX() << " " << pPreStepPoint->GetPosition().getY() << " " << pPreStepPoint->GetPosition().getZ() << " " << aTrack.GetGlobalTime() <<" " <<  chopper_t << " " << chopper_x << " " << foilct << " " << chopper_z << " " << chopper_px << " " << chopper_py << " " << "8000" << std::endl;

   aParticleChange.ProposeTrackStatus( fStopAndKill ) ;
 return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
}

// Boron absorption  we assume a 100 nm thick 96% enriched boron layer
G4double layerd = 2.5e-7; // Doubled by Husain on 27/08/2025 because faster neutrons were being transmitted
G4double enr = 0.96;
G4double sigmab = (2200./ v) * 3835.;
G4double mfp = 1/(sigmab*13.2);
G4double btrans = exp(-layerd/mfp/enr );
G4double beff = 1 - btrans;
G4cout << "Boron foil efficiency" << beff << G4endl;

G4double counted = G4UniformRand();
if (counted < beff){
	std::ofstream myfile(outputfile1, std::ofstream::app);
	  myfile << aTrack.GetKineticEnergy()/neV << " " << aTrack.GetMomentumDirection().getX() << " " <<  aTrack.GetMomentumDirection().getY() << " " <<  
		  aTrack.GetMomentumDirection().getZ() << " " << pPreStepPoint->GetPosition().getX() << " " << pPreStepPoint->GetPosition().getY() << " " << pPreStepPoint->GetPosition().getZ() << " " << aTrack.GetGlobalTime() << " " <<  
		  chopper_t << " " << chopper_x << " " << foilct << " " << chopper_z << " " << chopper_px << " " << chopper_py << " " << "6000" << std::endl;

chopper_x = 0;
chopper_y = 0;
chopper_z = 0;
chopper_px = 0;
chopper_py = 0;
chopper_pz = 0;
chopper_t = 0;
ezero = 0;
}
else {
std::ofstream myfile(outputfile1, std::ofstream::app);
	myfile << aTrack.GetKineticEnergy()/neV << " " << aTrack.GetMomentumDirection().getX() << " " <<  aTrack.GetMomentumDirection().getY() << " " <<  aTrack.GetMomentumDirection().getZ() << " " << pPreStepPoint->GetPosition().getX() << " " << pPreStepPoint->GetPosition().getY() << " " << pPreStepPoint->GetPosition().getZ() << " " << aTrack.GetGlobalTime() <<" " <<  chopper_t << " " << chopper_x << " " << foilct << " " << chopper_z << " " << chopper_px << " " << chopper_py << " " << "7000" << std::endl;
}
	  aParticleChange.ProposeTrackStatus( fStopAndKill ) ;
 return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
	  
}
/*
  if ( Material2->GetName() == "Detector1"){
	  myfile1 << aTrack.GetKineticEnergy()/neV << aTrack.GetMomentumDirection() << pPreStepPoint->GetPosition()   <<  aTrack.GetGlobalTime() << aTrack.GetPolarization() << std::endl;
	
 return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
  }
  if ( Material2->GetName() == "Detector2"){
	  myfile2 << aTrack.GetKineticEnergy()/neV << aTrack.GetMomentumDirection() << pPreStepPoint->GetPosition()  <<  aTrack.GetGlobalTime() <<aTrack.GetPolarization() << std::endl;
	  aParticleChange.ProposeTrackStatus( fStopAndKill ) ;
 return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);

  }
*/

if (Material1 == Material2)  return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
  
  G4ThreeVector theGlobalPoint = pPostStepPoint->GetPosition();
  G4ThreeVector thePrePoint = pPreStepPoint->GetPosition();
  
  G4Navigator* theNavigator =
    G4TransportationManager::GetTransportationManager()->
    GetNavigatorForTracking();
  
  G4ThreeVector theLocalPoint = theNavigator->
    GetGlobalToLocalTransform().
    TransformPoint(theGlobalPoint);
  
  G4ThreeVector theLocalNormal;	// Normal points back into volume
  
  G4bool valid;
  theLocalNormal = -theNavigator->GetLocalExitNormal(&valid);
  
  if (!valid) {
    G4cout << "local normal: " << theLocalNormal <<G4endl;
  } 
  
  G4ThreeVector theGlobalNormal = theNavigator->GetLocalToGlobalTransform().
    TransformAxis(theLocalNormal);
  
  G4double vel = aTrack.GetVelocity();
  G4ThreeVector momdir = aTrack.GetMomentumDirection();
  G4ThreeVector mom = aTrack.GetMomentum();
  G4double tfakt = momdir.dot(theGlobalNormal);
  G4double momnorm = mom.dot(theGlobalNormal);
  G4double velnorm = vel * tfakt;
  G4double enormal = momnorm*momnorm/2./CLHEP::neutron_mass_c2/neV; 
  G4double energy = aTrack.GetKineticEnergy()/neV; 
   
  G4MaterialPropertiesTable* aMaterialPropertiesTable;
  aMaterialPropertiesTable = Material2->GetMaterialPropertiesTable();

  G4double fermipot = 0.;
  G4double pdiffus = 0.;
  G4double pspinflip = 0.;
  G4double pupscatter = 0.;

 /*
  // properties in the new volume
  if (aMaterialPropertiesTable) { 
    fermipot = aMaterialPropertiesTable->GetConstProperty("FERMIPOT");
  }
 if (aMaterialPropertiesTable) { 
    pdiffus = aMaterialPropertiesTable->GetConstProperty("DIFFUSION");
  }  
 if (aMaterialPropertiesTable) {
    pspinflip = aMaterialPropertiesTable->GetConstProperty("SPINFLIP");
  }
 if (aMaterialPropertiesTable) {
    pupscatter = aMaterialPropertiesTable->GetConstProperty("LOSS");
  }

  // properties of the old volume
  G4MaterialPropertiesTable* aMaterialPropertiesTable2;
  aMaterialPropertiesTable2 = Material1->GetMaterialPropertiesTable();
  G4double fermipot_previous = 0.;
  G4double pspinflip_previous = 0.;
  G4double pdiffus_previous = 0.;
  G4double pupscatter_previous = 0.;
 if (aMaterialPropertiesTable2) { 
    fermipot_previous = aMaterialPropertiesTable2->GetConstProperty("FERMIPOT");
  }
 if (aMaterialPropertiesTable2) {
    pdiffus_previous = aMaterialPropertiesTable2->GetConstProperty("DIFFUSION");
  }
 if (aMaterialPropertiesTable2) {
    pspinflip_previous = aMaterialPropertiesTable2->GetConstProperty("SPINFLIP");
  }
 if (aMaterialPropertiesTable2) {
    pupscatter_previous = aMaterialPropertiesTable2->GetConstProperty("LOSS");
  }

*/

    G4MaterialPropertiesTable* aMaterialPropertiesTabl;
  aMaterialPropertiesTabl = Material2->GetMaterialPropertiesTable();
    if (aMaterialPropertiesTabl) { 
    G4double tab = aMaterialPropertiesTabl->GetConstProperty("SCATCS");
      if (tab == 1) {fermipot = opt1; pupscatter = eta1; pdiffus = diff1;}
      if (tab == 2) {fermipot = opt2; pupscatter = eta2; pdiffus = diff2;}
      if (tab == 3) {fermipot = opt3; pupscatter = eta3; pdiffus = diff3;}
      if (tab == 4) {fermipot = opt4; pupscatter = eta4; pdiffus = diff4;}
      if (tab == 5) {fermipot = opt5; pupscatter = eta5; pdiffus = diff5;}
    }
  G4MaterialPropertiesTable* aMaterialPropertiesTabl2;
  aMaterialPropertiesTabl2 = Material1->GetMaterialPropertiesTable();
  G4double fermipot_previous = 0.;
  G4double pspinflip_previous = 0.;
  G4double pdiffus_previous = 0.;
  G4double pupscatter_previous = 0.;
 if (aMaterialPropertiesTabl2) { 
   G4double tab = aMaterialPropertiesTabl2->GetConstProperty("SCATCS");
      if (tab == 1) {fermipot_previous = opt1; pupscatter_previous = eta1; pdiffus_previous = diff1;}
      if (tab == 2) {fermipot_previous = opt2; pupscatter_previous = eta2; pdiffus_previous = diff2;}
      if (tab == 3) {fermipot_previous = opt3; pupscatter_previous = eta3; pdiffus_previous = diff3;}
      if (tab == 4) {fermipot_previous = opt4; pupscatter_previous = eta4; pdiffus_previous = diff4;}
      if (tab == 5) {fermipot_previous = opt5; pupscatter_previous = eta5; pdiffus_previous = diff5;}
 }



  G4double fermipot_diff = fermipot - fermipot_previous;
  
  if (fermipotdiff_shutter < 1e300){
      if (verboselevel > 2){ 
	G4cout << "override by shutter " << fermipotdiff_shutter << G4endl;
      }
      fermipot_diff = fermipot - fermipot_previous + fermipotdiff_shutter;
      fermipotdiff_shutter = 1e300;
  }
 
  // calculate loss and reflectivity june 27 2025
  // golub eq 2.71
   G4double W = pupscatter * fermipot;
   G4double alphaW = sqrt((fermipot-enormal)*(fermipot-enormal)+W*W);
   G4double part1 = enormal - sqrt(enormal)*sqrt(2*alphaW - 2*(fermipot-enormal)) + alphaW;
   G4double part2 = enormal + sqrt(enormal)*sqrt(2*alphaW - 2*(fermipot-enormal)) + alphaW;
   G4double R_sq = (part1)/(part2);
   G4double reflection = R_sq;
   // until here
   //G4cout << "Just reflected: " << just_reflected << "\nStep number: " << aTrack.GetCurrentStepNumber() << "\nP(up scatter): " << pupscatter << "\nFermipot: " << fermipot <<  "\nReflection: " << reflection << "\n\n" << aTrack.GetKineticEnergy()/neV << " " << aTrack.GetMomentumDirection().getX() << " " <<  aTrack.GetMomentumDirection().getY() << " " <<  aTrack.GetMomentumDirection().getZ() << " " << pPreStepPoint->GetPosition().getX() << " " << pPreStepPoint->GetPosition().getY() << " " << pPreStepPoint->GetPosition().getZ() << "\n\n" << G4endl;  

// if (true): All particle hit the wall once and die
// if (rnd1 > reflection && just_reflected == 0): Particle is stuck (possible geom overlap)
//

  if (comparepot(enormal, fermipot_diff) == 0){             // below critical velocity
      /////// losses
      G4double rnd1 = G4UniformRand();
      remember_random = rnd1;
//G4cout << " below crit, ref = " << reflection << " rnd " << rnd1 << G4endl;

      if (rnd1 > reflection){             // loss on reflection
		G4cout << " loss!" << G4endl;
        std::ofstream myfile(outputfile1, std::ofstream::app);
        myfile << aTrack.GetKineticEnergy()/neV << " " << aTrack.GetMomentumDirection().getX() << " " <<  aTrack.GetMomentumDirection().getY() << " " <<  aTrack.GetMomentumDirection().getZ() << " " << pPreStepPoint->GetPosition().getX() << " " << pPreStepPoint->GetPosition().getY() << " " << pPreStepPoint->GetPosition().getZ() << " " << aTrack.GetGlobalTime() <<" " <<  chopper_t << " " << chopper_x << " " << foilct << " " << chopper_z << " " << chopper_px << " " << chopper_py << " " << "5000" << std::endl;
        // kill it.
        aParticleChange.ProposeTrackStatus( fStopAndKill ) ;
        return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
      }
      ////// reflect it
      ///

          //    G4cout << "reflect it! " << G4endl;
      G4ThreeVector ref = reflect(0, pdiffus, momdir, theGlobalNormal);
      ref = ref.unit();
      aParticleChange.ProposeMomentumDirection(ref);
      just_reflected = 1;
  }
  else {                                      // transmit material
	  
	  if (just_reflected == 0){
          G4double rnd2 = G4UniformRand();
	  if (remember_random){ rnd2 = remember_random; remember_random = 0; }
          if (rnd2 <  reflection){
            //G4cout << "above v_c reflect " << G4endl;
                ////// reflect it
          G4ThreeVector ref = reflect(0, pdiffus, momdir, theGlobalNormal);
          aParticleChange.ProposeMomentumDirection(ref.unit());
          just_reflected = 1;
          } 
	  

        /*
	  else {

		G4cout << "KILL WHEN ENTERING CHOPPER BLADES!!!! " << G4endl;
std::ofstream myfile(outputfile1, std::ofstream::app);
        myfile << aTrack.GetKineticEnergy()/neV << " " << aTrack.GetMomentumDirection().getX() << " " <<  aTrack.GetMomentumDirection().getY() << " " <<  aTrack.GetMomentumDirection().getZ() << " " << pPreStepPoint->GetPosition().getX() << " " << pPreStepPoint->GetPosition().getY() << " " << pPreStepPoint->GetPosition().getZ() << " " << aTrack.GetGlobalTime() <<" " <<  chopper_t << " " << chopper_x << " " << foilct << " " << chopper_z << " " << chopper_px << " " << chopper_py << " " << "17000" << std::endl;
        // kill it.
        aParticleChange.ProposeTrackStatus( fStopAndKill ) ;
        return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
	}*/

	          else {
//G4cout << " pot 1 = " << fermipot << ", prev = " << fermipot_previous << G4endl;

//G4cout << " material1 " << Material1 << ", mart2 " << Material2 << ", reflection = " << reflection << G4endl;


                //if (reflection < 0.000000001) return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
	
         // If NOT in chopper
         // SANITY CHECK --------------------------------------------------- 
	 std::ofstream myfile(outputfile1, std::ofstream::app);   
	 myfile << aTrack.GetKineticEnergy()/neV << " " << aTrack.GetMomentumDirection().getX() << " " <<  aTrack.GetMomentumDirection().getY() << " " <<  aTrack.GetMomentumDirection().getZ() << " " << pPreStepPoint->GetPosition().getX() << " " << pPreStepPoint->GetPosition().getY() << " " << pPreStepPoint->GetPosition().getZ() << " " << aTrack.GetGlobalTime() <<" " <<  chopper_t << " " << chopper_x << " " << foilct << " " << chopper_z << " " << chopper_px << " " << chopper_py << " " << "5000" << std::endl;
          // kill it.
          aParticleChange.ProposeTrackStatus( fStopAndKill ) ;
          return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
          // END SANITY CHECK ----------------------------------------------- 

       /* if ( fermipot < 0  ){





                G4cout << "KILL WHEN ENTERING CHOPPER BLADES!!!! " << G4endl;
std::ofstream myfile(outputfile1, std::ofstream::app);
   return &aParticleChange;     myfile << aTrack.GetKineticEnergy()/neV << " " << aTrack.GetMomentumDirection().getX() << " " <<  aTrack.GetMomentumDirection().getY() << " " <<  aTrack.GetMomentumDirection().getZ() << " " << pPreStepPoint->GetPosition().getX() << " " << pPreStepPoint->GetPosition().getY() << " " << pPreStepPoint->GetPosition().getZ() << " " << aTrack.GetGlobalTime() <<" " <<  chopper_t << " " << chopper_x << " " << foilct << " " << chopper_z << " " << chopper_px << " " << chopper_py << " " << "17000" << std::endl;
        // kill it.
        aParticleChange.ProposeTrackStatus( fStopAndKill ) ;                                                                                 return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);                                                                      
          } */ 

	    
          G4cout << "above vc transmit" << G4endl; 
          // --- transmission because it is faster than the critical velocity
          G4double enew = transmit(fermipot_diff, energy);
          G4double m = -sqrt(momnorm*momnorm - CLHEP::neutron_mass_c2*2.*fermipot_diff*neV);
          G4ThreeVector ref = mom - (momnorm-m)*theGlobalNormal;
          aParticleChange.ProposeMomentumDirection(ref.unit());
          aParticleChange.ProposeEnergy(enew*neV); 
	  
         }    
       }
       else if (just_reflected == 1){
        just_reflected = 0 ;
      }
  }   


// this is the approximated reflection probability as it was used in the past

/*
  ///// below critical velocity
    if (comparepot(enormal, fermipot_diff) == 0){             // reflect from surface
      if (verboselevel >1) G4cout << "MATERIALBOUNDARY: reflect " << G4endl;
      /////// losses
      if (loss(pupscatter, velnorm, fermipot_diff) == 1){             // loss on reflection
	if (verboselevel >1) G4cout << "MATERIALBOUNDARY: loss on surface " << G4endl;
std::ofstream myfile("detector.txt", std::ofstream::app);
myfile << aTrack.GetKineticEnergy()/neV << " " << aTrack.GetMomentumDirection().getX() << " " <<  aTrack.GetMomentumDirection().getY() << " " <<  aTrack.GetMomentumDirection().getZ() << " " << pPreStepPoint->GetPosition().getX() << " " << pPreStepPoint->GetPosition().getY() << " " << pPreStepPoint->GetPosition().getZ() << " " << aTrack.GetGlobalTime() <<" " <<  chopper_t << " " << chopper_x << " " << foilct << " " << chopper_z << " " << chopper_px << " " << chopper_py << " " << "5000" << std::endl;



	// kill it.
	aParticleChange.ProposeTrackStatus( fStopAndKill ) ;
   return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);

      }
      /////// spinflips
      if (spinflip(pspinflip) == 1) {
	if (verboselevel >1) G4cout << "MATERIALBOUNDARY: spinflip " << G4endl;
	G4ThreeVector spin = aTrack.GetPolarization();
	spin *= -1;
	aParticleChange.ProposePolarization(spin);
	
      }
      ////// reflect it
      G4ThreeVector ref = reflect(0, pdiffus, momdir, theGlobalNormal);
      ref = ref.unit();
      aParticleChange.ProposeMomentumDirection(ref);
      just_reflected = 1;
    }
    
    /////// above critical velocity
    else {                                      // transmit material
      if (just_reflected == 0){
	// if it is faster than the crticial velocity, there is a probability to be still
	// reflected. this formula is (only) valid for low loss materials
	G4double refl2 = reflectivity(0, fermipot_diff,enormal);
	if (verboselevel > 5)G4cout << "jr " << just_reflected << " " << G4endl; 
	if (verboselevel >1) G4cout << "MATERIALBOUNDARY: reflectivity " << refl2 << G4endl;
	if (G4UniformRand() < refl2){ 
	G4cout << "above v_c reflect " << G4endl;
	      	////// reflect it
	  G4ThreeVector ref = reflect(0, pdiffus, momdir, theGlobalNormal);
	  aParticleChange.ProposeMomentumDirection(ref);
	  // set the just_reflected variable 1
	  just_reflected = 1;
	} else {

	G4cout << "above vc transmi" << G4endl;	
	  // --- transmission because it is faster than the critical velocity
	 
	  // --- kinetic energy in the new media
	  G4double enew = transmit(fermipot_diff, energy);
	  
	  // --- change of the normal momentum component
	  //     p = sqrt(2*m*Ekin)
	  G4double m = -sqrt(momnorm*momnorm - CLHEP::neutron_mass_c2*2.*fermipot_diff*neV);
	  
	  // --- momentum direction in new media
	  G4ThreeVector ref = mom - (momnorm-m)*theGlobalNormal;
	  aParticleChange.ProposeMomentumDirection(ref.unit());
	  aParticleChange.ProposeEnergy(enew*neV);

	}
      }
      else if (just_reflected == 1){
	just_reflected = 0 ;
      }
    }
*/


  return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
    }


G4double UCNMaterialBoundary2::GetMeanFreePath(const G4Track&,
 				         G4double ,
				         G4ForceCondition* condition)
{
//  SetGPILSelection(CandidateForSelection);
  *condition = Forced;
//G4cout << "MFP " << G4endl;
  return DBL_MAX;
}
/*
G4VParticleChange* UCNMaterialBoundary2::AlongStepDoIt(
                             const G4Track& aTrack,
                             const G4Step&)
{
	G4cout << "along... " << G4endl;
// 2at TOF system start  t = 3810  
aParticleChange.Initialize(aTrack);
  return &aParticleChange;
}

G4double UCNMaterialBoundary2::GetContinuousStepLimit(const G4Track &  aTrack,
                                                                G4double  previousStepSize,
                                                                G4double currentMinimumStep,
                                                                G4double& )


{

G4cout << "contstep" << G4endl;      
	if ((previousStepSize == 0 || just_reflected==1) && aTrack.GetCurrentStepNumber()-1) {
    return .1*micrometer;
  }
  return std::min(currentMinimumStep, 10*cm);
}
*/


int UCNMaterialBoundary2::comparepot(G4double energy, G4double fermipot_diff)
{
  return (energy > fermipot_diff);    
}


int UCNMaterialBoundary2::loss(G4double coefficient, G4double velnorm, G4double
fermipot){
  
  /// the surface roughness is not taken into account here, 
  // one could use e.g. ultracold neutrons, r.golub, p.35,
  // where mu is increased by roughness parameters sigma and
  // omega, which are related to the height and the distance of
  // "disturbances" on the surface 
  
  G4double v_bound = sqrt(2.*fermipot*neV/CLHEP::neutron_mass_c2*CLHEP::c_squared);
  G4double v_ratio = velnorm/v_bound;
  G4double loss_von_E = (2*coefficient*v_ratio)/(sqrt(1-(v_ratio*v_ratio)));

  return (G4UniformRand() <= fabs(loss_von_E));
}


int UCNMaterialBoundary2::spinflip(G4double coefficient){
  G4double rnd_loss = G4UniformRand();
  if (rnd_loss <= coefficient){
   //G4cout << "MATERIALBOUNDARY: spinflip " << rnd_loss << G4endl;
   return 1;
  }
  else{
   return 0;
  }
}

G4double UCNMaterialBoundary2::reflectivity(int, G4double fpot, G4double enormal){
  
  G4double r = (sqrt(enormal) - sqrt(enormal - fpot)) /
  (sqrt(enormal) + sqrt(enormal - fpot));
  
  return r*r;
} 
G4ThreeVector UCNMaterialBoundary2::reflect(int, G4double coefficient, 
					   G4ThreeVector momdir, G4ThreeVector localnormal){
  // reflect specular
  G4double tfakt = momdir.dot(localnormal);
  G4ThreeVector reflected = momdir - 2 * localnormal * tfakt;
  if (verboselevel > 5)G4cout << "specular reflected << " << reflected << G4endl;

  // reflect diffuse 
  if (reflected==momdir || G4UniformRand() < coefficient) {
    G4ThreeVector diffus = cos_diff(localnormal);
    if (verboselevel > 5)G4cout << "diffus localnormal " << localnormal << ", " << diffus << G4endl;
    return diffus;
  }
  return reflected;
}

G4double UCNMaterialBoundary2::transmit(G4double fermipot, G4double energy){
  
  return energy - fermipot;

}

G4ThreeVector UCNMaterialBoundary2::cos_diff(G4ThreeVector localnormal)
{
  G4ThreeVector momentum;
  // cosine distribution - Lambert's law
  momentum.setRThetaPhi(1., acos(sqrt(G4UniformRand())), 2.*CLHEP::pi*G4UniformRand());
  momentum.rotateUz(localnormal);
  if(momentum*localnormal<0) {momentum*=-1;G4cout << "!" << G4endl;}
  
  //if(momentum*localnormal<0.03){
  //  while(momentum*localnormal<0.03){
  //     momentum.setRThetaPhi(1., acos(sqrt(G4UniformRand())), 2.*CLHEP::pi*G4UniformRand());
  //      momentum.rotateUz(localnormal);
  //  }

  //}

  return momentum;
}

UCNMaterialBoundary2* UCNMaterialBoundary2::GetInstance()
{
	return theInstance;
}

void UCNMaterialBoundary2::SetFermiPotDiff(G4double state)
{
	fermipotdiff_shutter = state;
}

void UCNMaterialBoundary2::setVerbose(G4int level)
{
	verboselevel = level ;
}

void UCNMaterialBoundary2::SetShutterClose(G4String newval){
  int nr = 0;
  float ti = 0.;
  sscanf(newval, "%d %f", &nr, &ti); // shutternr time
  if (nr >= UCN_SHUTTERS_MAX || nr < 0) {
    G4cerr << "invalid shutter number " << nr << G4endl;
    return;
  }
  G4cout << "ucnshutter: setshutterclose (shutter " << nr
	 << ", t=" << ti << G4endl;
  shutter_states[nr].AddState(0, ti*second);
}

void UCNMaterialBoundary2::SetShutterOpen(G4String newval){
  int nr = 0;
  float ti = 0.;
  sscanf(newval, "%d %f", &nr, &ti); // shutternr time
  if (nr >= UCN_SHUTTERS_MAX || nr < 0) {
    G4cerr << "invalid shutter number " << nr << G4endl;
    return;
  }
  G4cout << "ucnshutter: setshutteropen (shutter " << nr
	 << ", t=" << ti << G4endl;
  shutter_states[nr].AddState(1, ti*second);
}

void UCNMaterialBoundary2::SetShutterLoop(G4String newval){
  int nr = 0;
  int loop = 0;
  sscanf(newval, "%d %d", &nr, &loop);
  if (nr >= UCN_SHUTTERS_MAX || nr < 0) {
    G4cerr << "invalid shutter number " << nr << G4endl;
    return;
  }
  G4cout << "ucnshutter: setshutter loop " << loop
	 << " (shutter " << nr << ")" << G4endl;
  shutter_states[nr].SetLoop(loop);
}

void UCNMaterialBoundary2::SetUseShutters(G4String newval){
if (newval == "1"){
 useshutters = 1;
}

}
void UCNMaterialBoundary2::SetShutterVerbose(G4String newval){
}

void UCNMaterialBoundary2::SetFileName1(G4String fil){
outputfile1 = fil;
G4cout << "List of Material Properties read from macro:" << G4endl;
G4cout << "opt1     " << opt1 << G4endl;
G4cout << "diff1     " << diff1 << G4endl;
G4cout << "eta1     " << eta1 << G4endl;
G4cout << "opt2     " << opt2 << G4endl;
G4cout << "diff2     " << diff2 << G4endl;
G4cout << "eta2     " << eta2 << G4endl;
G4cout << "opt3     " << opt3 << G4endl;
G4cout << "diff3     " << diff3 << G4endl;
G4cout << "eta3     " << eta3 << G4endl;



std::ofstream myfile(outputfile1, std::ofstream::app);
  myfile << "# Exit Code Key:"<<std::endl;
  myfile << "# 1000: Killed at pp foil"<<std::endl;
  myfile << "# 2000: Reached maximum defined steps"<<std::endl;
  myfile << "# 3000: Undefined"<<std::endl;
  myfile << "# 4000: Undefined"<<std::endl;
  myfile << "# 5000: Lost due to reflection"<<std::endl;
  myfile << "# 6000: Absorbed by Boron in the detector"<<std::endl;
  myfile << "# 7000: Other interaction with detector"<<std::endl;
  myfile << "# 8000: Absorbed by AlMG3 foil"<<std::endl;
  myfile << "# 9000: Undefined"<<std::endl;
  myfile << "# 10000: Undefined"<<std::endl;
  myfile << "# 11000: Beta decay"<<std::endl;
  myfile << "# 12000: Helium loss"<<std::endl;
  myfile << "# 13000: Reached maximum global time"<<std::endl;
  myfile << "# 14000: Undefined"<<std::endl;
  myfile << "# 15000: Falling down"<<std::endl;
  myfile << "# 16000: Undefined"<<std::endl;
  myfile << "# 17000: Killed in chopper blades"<<std::endl;
  myfile << "# -------------------------------------------------------------------------------------------"<<std::endl;
  myfile << "# KE[neV] px[unitless] py[unitless] pz[unitless] x[mm] y[mm] z[mm] t_global[s] t_chopper[s] x_chopper[mm] n_foil z_chopper[mm] px_chopper[unitless] py_chopper[unitless] Exit Code" << std::endl;
//G4cout << "MaterialBoundary2: set outputfile1 " << fil << G4endl;
}

void UCNMaterialBoundary2::SetFileName2(G4String fil){
outputfile2 = fil;

//G4cout << " set outputfile2 " << fil << G4endl;
}
void UCNMaterialBoundary2::SetFileName3(G4String fil){
outputfile3 = fil;

//G4cout << " set outputfile3 " << fil << G4endl;
}

void UCNMaterialBoundary2::SetSteps(G4double fil){
stepsbeforekill = fil;

//G4cout << " set steps before kill" << fil << G4endl;
}

void UCNMaterialBoundary2::SetLifeTime(G4double fil){
meanLifetime = fil;

//G4cout << " set n lifetime" << fil << G4endl;
}
void UCNMaterialBoundary2::SetHeLifeTime(G4double fil){
meanHeLifetime = fil;

//G4cout << " set he lifetime" << fil << G4endl;
}
void UCNMaterialBoundary2::SetFoil(G4double fil){
foilabs = fil;

//G4cout << " set vac sep foil abs" << fil << G4endl;
}
void UCNMaterialBoundary2::SetTime(G4double fil){
timesteps = fil;
if (timesteps > 0){ 
std::ofstream myfile(outputfile2, std::ofstream::app);
myfile << "# KE[neV] px[unitless] py[unitless] pz[unitless] x[mm] y[mm] z[mm] t_global[s]" << std::endl;
}

//G4cout << " set timesteps" << fil << G4endl;
}
void UCNMaterialBoundary2::SetOpt1(G4double fil){opt1 = fil;}
void UCNMaterialBoundary2::SetOpt2(G4double fil){opt2 = fil;}
void UCNMaterialBoundary2::SetOpt3(G4double fil){opt3 = fil;}
void UCNMaterialBoundary2::SetEta1(G4double fil){eta1 = fil;}
void UCNMaterialBoundary2::SetOpt4(G4double fil){opt4 = fil;}
void UCNMaterialBoundary2::SetOpt5(G4double fil){opt5 = fil;}
void UCNMaterialBoundary2::SetEta4(G4double fil){eta4 = fil;}
void UCNMaterialBoundary2::SetEta2(G4double fil){eta2 = fil;}
void UCNMaterialBoundary2::SetEta3(G4double fil){eta3 = fil;}
void UCNMaterialBoundary2::SetEta5(G4double fil){eta5 = fil;}
void UCNMaterialBoundary2::SetDiff1(G4double fil){diff1 = fil;}
void UCNMaterialBoundary2::SetDiff2(G4double fil){diff2 = fil;}
void UCNMaterialBoundary2::SetDiff3(G4double fil){diff3 = fil;}
void UCNMaterialBoundary2::SetDiff4(G4double fil){diff4 = fil;}
void UCNMaterialBoundary2::SetDiff5(G4double fil){diff5 = fil;}
void UCNMaterialBoundary2::SetC(G4double fil){C = fil;}
void UCNMaterialBoundary2::SetV(G4double fil){V = fil;}
void UCNMaterialBoundary2::Setf_a(G4double fil){f_a = fil;}
void UCNMaterialBoundary2::Setf_h(G4double fil){f_h = fil;}
// Excited States
void UCNMaterialBoundary2::SetexcitedFrac(G4double fil){excitedFrac = fil;}
void UCNMaterialBoundary2::SetmeanDownShiftLifetime(G4double fil){meanDownShiftLifetime = fil;}
void UCNMaterialBoundary2::SetE_recoil(G4double fil){E_recoil = fil;}
void UCNMaterialBoundary2::SetmaxGlobalTime(G4double fil){maxGlobalTime = fil;}
// End Excited States

// replace Gamma1a e.g. with meanHeLifetime
// replace Gamma2a e.g. with meanLifetime
void UCNMaterialBoundary2::SetGamma1a(G4double fil){Gamma1a = fil;}
void UCNMaterialBoundary2::SetGamma1h(G4double fil){Gamma1h = fil;}
void UCNMaterialBoundary2::SetGamma2a(G4double fil){Gamma2a = fil;}
void UCNMaterialBoundary2::SetGamma2h(G4double fil){Gamma2h = fil;}
void UCNMaterialBoundary2::SetGamma3a(G4double fil){Gamma3a = fil;}
void UCNMaterialBoundary2::SetGamma3h(G4double fil){Gamma3h = fil;}
void UCNMaterialBoundary2::Setlambda(G4double fil){lambda = fil;}
void UCNMaterialBoundary2::Setu_min(G4double fil){u_min = fil;}
void UCNMaterialBoundary2::Setu_max(G4double fil){u_max = fil;}

