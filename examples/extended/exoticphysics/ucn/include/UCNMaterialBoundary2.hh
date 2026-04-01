
////////////////////////////////////////////////////////////////////////
// 4.9.04 UCN material boundary process, peter fierlinger
// 
////////////////////////////////////////////////////////////////////////


#ifndef UCNMaterialBoundary2_h
#define UCNMaterialBoundary2_h 1

/////////////
// Includes
/////////////

#include "globals.hh"
#include "templates.hh"
#include "Randomize.hh"
#include "G4Step.hh"
#include "G4VContinuousDiscreteProcess.hh"
#include "G4VDiscreteProcess.hh"
#include "G4DynamicParticle.hh"
#include "G4Material.hh"
#include "G4Neutron.hh"

#define UCN_SHUTTERS_MAX 10
#define UCN_SHUTTERS_STATES_MAX 20
#define UCN_SHUTTERS_DEFAULT_STATE 0
#define UCN_SHUTTERS_DEFAULT_LOOP 0

class UCNMaterialBoundary2Messenger;
class UCNShutterMessenger;

class UCNShutterStates
{
 public:
  int n;
  int loop;
  float times[UCN_SHUTTERS_STATES_MAX];
  int states[UCN_SHUTTERS_STATES_MAX];
  UCNShutterStates() : n(0), loop(UCN_SHUTTERS_DEFAULT_LOOP) {};
  void AddState(int state, float time);
  int GetState(float time);
  void SetLoop(int loop);
};

//class UCNMaterialBoundary2 : public G4VContinuousDiscreteProcess
class UCNMaterialBoundary2 : public G4VDiscreteProcess
{

private:
	// Read Chopper Function for vTOF 
	// --- Jitter PDF sampler (built once in constructor) ---
    	void LoadAndBuildCDF(const G4String& filename);
        G4double SampleSmear() const;
	
	G4String fChopperFile;                  // set by mac command
	std::vector<G4double> fCDF_times;   // centered times in seconds (x-axis of inverse CDF)
	std::vector<G4double> fCDF_values;  // CDF values in [0,1]  (y-axis, used to look up)
	// Stored inverted: fCDF_values[i] -> fCDF_times[i]
	// so sampling = find where uniform random falls in fCDF_values, interpolate fCDF_times
	// End Chopper Logic


public: 

//2(const G4String& processName, G4ProcessType type)
// G4VDiscreteProcess(const G4String& aName,
    //                   G4ProcessType aType = fNotDefined );
   // G4VDiscreteProcess(G4VDiscreteProcess&);



         UCNMaterialBoundary2(const G4String& processName = 
"UCNMaterialBoundary", const G4ProcessType aType = fNotDefined);

	~UCNMaterialBoundary2();

public: 

        inline G4bool IsApplicable(const G4ParticleDefinition& aParticleType);
	//G4double GetContinuousStepLimit(const G4Track& aTrack,
        //                     G4double  previousStepSize,
        //                     G4double  currentMinimumStep,
	//		     G4double& currentSafety);

	G4double GetMeanFreePath(const G4Track& aTrack,
				 G4double ,
				 G4ForceCondition* condition);

 	G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
 				        const G4Step&  aStep);

	//G4VParticleChange* AlongStepDoIt(
        //                     const G4Track&,
        //                     const G4Step&);

	int comparepot(double energy, double fermipot);
	int loss(double coefficient, double velnorm, double fermipot);
	int spinflip(double coefficient);
	G4ThreeVector reflect(int modus, double coefficient, G4ThreeVector momdir, G4ThreeVector localnormal);
	double reflectivity(int, double fpot, double enormal);

	double transmit(double fermipot, double energy);
	G4ThreeVector cos_diff(G4ThreeVector localnormal);
        
	
        static UCNMaterialBoundary2* GetInstance();
	void SetFermiPotDiff(double state);
	void setVerbose(G4int level);
	void SetShutterClose(G4String newval);
	void SetShutterOpen(G4String newval);
	void SetShutterLoop(G4String newval);
        void SetUseShutters(G4String newval);
	void SetShutterVerbose(G4String newval);

 void SetFileName1(G4String fil);
  void SetFileName2(G4String fil);
   void SetFileName3(G4String fil);
   void SetSteps(G4double fil);
   void SetLifeTime(G4double fil);
   void SetHeLifeTime(G4double fil);
   void SetFoil(G4double fil);
   void SetTime(G4double fil);

  void SetDiff1(G4double fil);
  void SetDiff2(G4double fil);
  void SetDiff3(G4double fil);
  void SetDiff4(G4double fil);
  void SetDiff5(G4double fil);
  void SetOpt1(G4double fil);
  void SetOpt2(G4double fil);
  void SetOpt3(G4double fil);
   void SetOpt4(G4double fil);
 void SetOpt5(G4double fil);
  void SetEta1(G4double fil);
  void SetEta2(G4double fil);
  void SetEta3(G4double fil);
  void SetEta4(G4double fil);
  void SetEta5(G4double fil);

void SetDiff6(G4double val);
void SetOpt6(G4double val);
void SetEta6(G4double val);

   void SetC(G4double fil);
  void SetV(G4double fil);
  void Setf_a(G4double fil);
  void Setf_h(G4double fil);
  void Setlambda(G4double fil);
  void SetGamma1a(G4double fil);
  void SetGamma1h(G4double fil);
  void SetGamma2a(G4double fil);
  void SetGamma2h(G4double fil);
  void SetGamma3a(G4double fil);
  void SetGamma3h(G4double fil);
  void Setu_min(G4double fil);
  void Setu_max(G4double fil);
// Excited States
  void SetexcitedFrac(G4double fil);
  void SetmeanDownShiftLifetime(G4double fil);
  void SetE_recoil(G4double fil);
  void SetmaxGlobalTime(G4double fil);
// End Sxcited States
//int read_two_column_file(const char* filename, double x[], double y[], std::size_t maxSize);
//G4double Chopperfunc(G4String fnam, G4double t) ;
//void SetChopperFunction(G4String fil); // redone below as inline function on 30.03.2026 
// Called by the messenger when /mat/chopperFile is parsed
void SetChopperFunction(const G4String& filename) { 
        fChopperFile = filename; 
        LoadAndBuildCDF(filename);  // builds the CDF right then
    }
G4String GetChopperFile() const { return fChopperFile; }


private:
	UCNShutterStates shutter_states[UCN_SHUTTERS_MAX];
	UCNShutterMessenger * theShutterMessenger;
	UCNMaterialBoundary2Messenger * theMessenger;
	G4int just_reflected;
	G4int foil; 
        G4int just_he_reflected;
	static UCNMaterialBoundary2* theInstance;
        double fermipotdiff_shutter;
	G4int verboselevel;
	int return_it; // shutter 
	int useshutters; // shutter
	double sh_fermipot ;
	double sh_spinflip ;
	double sh_loss ;
	double sh_diffusion ;
	double sh_reflectivity ;
	double sh_abscs ;
	double sh_losscs ;
	double sh_scatcs ;	
        std::ofstream out;
        G4bool atReflection, transmited, stop_gravity;

G4double diff1;
G4double diff2;
G4double diff3;
G4double diff4;
G4double diff5;
G4double diff6;

G4double opt1;
G4double opt2;
G4double opt3;
G4double opt6;

G4double eta1;
G4double eta2;
G4double eta3;
G4double eta4;
G4double eta5;
G4double eta6;

G4double opt4;
G4double opt5;
G4double lambda;
G4double f_a;
G4double f_h;
G4double u_min;
G4double u_max;
G4double C;
G4double V;
G4double Gamma1a;
G4double Gamma1h;
G4double Gamma2a;
G4double Gamma2h;
G4double Gamma3a;
G4double Gamma3h;

// Excited States
G4double excitedFrac;
G4double meanDownShiftLifetime;
G4double E_recoil;
G4double setMaxGlobalTime;
// End Excited States

};

////////////////////
// Inline methods
////////////////////

inline
G4bool UCNMaterialBoundary2::IsApplicable(const G4ParticleDefinition& 
aParticleType)
{
//G4cout << "is applicable " << (&aParticleType == G4Neutron::Neutron()) << G4endl;
return ( &aParticleType == G4Neutron::Neutron() );
}

#endif /* UCNMaterialBoundary2_h */
