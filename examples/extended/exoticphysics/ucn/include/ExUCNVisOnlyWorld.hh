// NEW (Claude, 2026-11-17): vis-only parallel world for the converter end-unit enclosure.
//
// The real Enclosure solid (enclosure_plus_window.gdml) is used for tracking in the mass world,
// but Geant4 cannot build a polyhedron for it (meshing does not finish), so it is hidden from vis.
// This parallel world holds a drawable stand-in, enclosure_body_visonly.gdml (make_gdml.py,
// enclosure(simplified=True): same body, no R10 pocket blends / small holes, split into 8 pieces).
//
// CORRECTION (2026-11-22): it is NOT physics-neutral. Registering any parallel world makes
// Geant4 use CoupledTransportation for all particles (G4RunManagerKernel::InitializePhysics),
// and with it neutrons in this application are transported without gravity (measured).
// Therefore it is only registered with the command-line switch "ExUCN -visbody 1"
// (ExUCN.cc), for looking at the geometry. Draw it with "/vis/drawVolume worlds" (vis.mac).
#ifndef ExUCNVisOnlyWorld_h
#define ExUCNVisOnlyWorld_h 1

#include "G4VUserParallelWorld.hh"

class ExUCNDetectorConstruction;

class ExUCNVisOnlyWorld : public G4VUserParallelWorld
{
  public:
    ExUCNVisOnlyWorld(const G4String& worldName, const ExUCNDetectorConstruction* det);
    ~ExUCNVisOnlyWorld() override = default;
    void Construct() override;

  private:
    const ExUCNDetectorConstruction* fDetector;
};

#endif
