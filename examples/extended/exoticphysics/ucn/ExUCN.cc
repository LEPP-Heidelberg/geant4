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
/// \file exoticphysics/ucn/ExUCN.cc
/// \brief Main program of the exoticphysics/ucn example
//
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef WIN32
#  include <unistd.h>
#endif

#include "ExUCNActionInitialization.hh"
#include "ExUCNDetectorConstruction.hh"
#include "ExUCNVisOnlyWorld.hh"
#include "ExUCNPhysicsList.hh"

//#include "G4RunManagerFactory.hh"
#include "G4Types.hh"
#include "G4UIExecutive.hh"
#include "G4UIcommand.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "Randomize.hh"
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

namespace
{
void PrintUsage()
{
  G4cerr << " Usage: " << G4endl;
  G4cerr << " ExUCN [-m macro ] [-u UIsession] [-t nThreads] [-r seed] [-visbody 0|1] " << G4endl;
  G4cerr << "   note: -t option is available only for multi-threaded mode." << G4endl;
  G4cerr << "   -visbody 1: draw the converter end-unit enclosure body (vis-only parallel world)." << G4endl;
  G4cerr << "               WARNING: this switches neutron transport to CoupledTransportation, which" << G4endl;
  G4cerr << "               transports WITHOUT GRAVITY. For looking at the geometry only, never for physics." << G4endl;
  G4cerr << "               Default 0 (off)." << G4endl;
}
}  // namespace

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc, char** argv)
{
  // Evaluate arguments
  //
  if (argc > 11) {
    PrintUsage();
    return 1;
  }

  G4String macro;
  G4String session;
  G4int nThreads = 0;
  G4bool drawEnclosureBody = false;

  G4long myseed = time(NULL); // 15733297276;//time(NULL); //1234*G4UniformRand();
  for (G4int i = 1; i < argc; i = i + 2) {
    if (G4String(argv[i]) == "-m")
      macro = argv[i + 1];
    else if (G4String(argv[i]) == "-u")
      session = argv[i + 1];
    else if (G4String(argv[i]) == "-r")
      myseed = atoi(argv[i + 1]);
    else if (G4String(argv[i]) == "-t") {
      nThreads = G4UIcommand::ConvertToInt(argv[i + 1]);
    }
    else if (G4String(argv[i]) == "-visbody" && i + 1 < argc) {
      drawEnclosureBody = G4UIcommand::ConvertToBool(argv[i + 1]);
    }
    else {
      PrintUsage();
      return 1;
    }
  }


  // Instantiate G4UIExecutive if interactive mode
  // (NEW (Claude, 2026-11-22): interactive = no macro given; was argc == 1, which made
  //  "ExUCN -visbody 1" run without a UI session and crash.)
  G4UIExecutive* ui = nullptr;
  if (macro.empty()) {
    ui = new G4UIExecutive(argc, argv);
  }

  // Choose the Random engine
  //
  //G4Random::setTheEngine(new CLHEP::RanecuEngine);

  // Construct the default run manager
  //
  //auto* runManager = G4RunManagerFactory::CreateRunManager();
  auto* runManager = new G4RunManager();
   if (nThreads > 0) runManager->SetNumberOfThreads(nThreads);



  // Seed the random number generator manually
  G4Random::setTheSeed(myseed);

  // Set mandatory initialization classes
  //
  // Detector construction
  auto detector = new ExUCNDetectorConstruction();
  // NEW (Claude): vis-only parallel world with a drawable copy of the converter end-unit
  // enclosure body (the real Enclosure solid cannot be meshed by the viewer).
  // ONLY behind the -visbody switch (2026-11-22): registering ANY parallel world makes Geant4
  // use CoupledTransportation (G4RunManagerKernel::InitializePhysics), and with it neutrons in
  // this application are transported WITHOUT GRAVITY (measured: straight paths, constant
  // kinetic energy). See GDML_ENCLOSURE_INTEGRATION.md. Off by default.
  if (drawEnclosureBody) {
    detector->RegisterParallelWorld(new ExUCNVisOnlyWorld("EnclosureVisOnlyWorld", detector));
    G4cout << "\n*****************************************************************************\n"
           << "*** -visbody 1: vis-only enclosure body ON (parallel world registered).     ***\n"
           << "*** Neutron transport is now CoupledTransportation: NO GRAVITY.             ***\n"
           << "*** Use this only to look at the geometry; do NOT use the results.          ***\n"
           << "*****************************************************************************\n"
           << G4endl;
  }
  else {
    G4cout << "*** Vis-only enclosure body OFF (default; physics unaffected). "
           << "Start with -visbody 1 to draw it (disables gravity)." << G4endl;
  }
  runManager->SetUserInitialization(detector);
  // Physics list
  runManager->SetUserInitialization(new ExUCNPhysicsList());
  // User action initialization
  runManager->SetUserInitialization(new ExUCNActionInitialization(detector));

  // Initialize G4 kernel
  //
  runManager->Initialize();
  //G4Random::setTheSeed(myseed);

  // Initialize visualization
  //
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
 
 
  visManager->Initialize();

  // Get the pointer to the User Interface manager
  //
  
  G4UImanager* UImanager = G4UImanager::GetUIpointer();
 std::ofstream myfile("start.txt", std::ofstream::app);                                                                                   
 myfile << "# Energy[neV] px py pz x[mm] y[mm] z[mm]" << std::endl;

 if (macro.size()) {
    // batch mode
    G4String command = "/control/execute ";
 //UImanager->ApplyCommand("/control/execute vis.mac");
																	   UImanager->ApplyCommand(command + macro);
  }
  else {  // interactive mode : define UI session
    UImanager->ApplyCommand("/control/execute vis.mac");
    if (ui->IsGUI()) UImanager->ApplyCommand("/control/execute gui.mac");
    ui->SessionStart();
    delete ui;
}


  G4cout << "End of ExUCN.cc" << G4endl;
  delete visManager;
  G4cout << "Deleted visManager" << G4endl;
  delete UImanager;
  G4cout << "Deleted UImanager" << G4endl;
  //delete runManager;
  //G4cout << "Deleted runManager" << G4endl;
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
