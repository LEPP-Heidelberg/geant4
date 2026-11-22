// NEW (Claude, 2026-11-17): see include/ExUCNVisOnlyWorld.hh.
#include "ExUCNVisOnlyWorld.hh"

#include "ExUCNDetectorConstruction.hh"

#include "G4Colour.hh"
#include "G4GDMLParser.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4ios.hh"

ExUCNVisOnlyWorld::ExUCNVisOnlyWorld(const G4String& worldName,
                                     const ExUCNDetectorConstruction* det)
  : G4VUserParallelWorld(worldName), fDetector(det)
{}

void ExUCNVisOnlyWorld::Construct()
{
  G4LogicalVolume* worldLV = GetWorld()->GetLogicalVolume();
  worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());

  // /run/reinitializeGeometry calls this again on the SAME parallel world: drop the old pieces.
  while (worldLV->GetNoDaughters() > 0) {
    G4VPhysicalVolume* old = worldLV->GetDaughter(0);
    worldLV->RemoveDaughter(old);
    delete old;
  }

  G4GDMLParser parser;
  parser.Read("../ucn/enclosure_body_visonly.gdml", false);
  // The GDML world is named "enclosure_body_visonly_World" (not "World"), so this name lookup
  // cannot resolve to the simulation's own world volume.
  G4LogicalVolume* gdmlWorldLV = parser.GetWorldVolume()->GetLogicalVolume();

  auto* bodyVis = new G4VisAttributes(G4Colour(0.6, 0.6, 0.65, 0.35));
  // Geant4's boolean mesher (BooleanProcessor) is unreliable for these pieces at high circle
  // resolution (vis.mac sets 100 segments/circle: several pieces fail or come out partial),
  // but all mesh correctly at 24 -- force that for these vis-only volumes only.
  bodyVis->SetForceLineSegmentsPerCircle(24);
  const G4ThreeVector origin = fDetector->GetEnclosurePartOrigin();
  G4int nPlaced = 0;
  for (std::size_t i = 0; i < gdmlWorldLV->GetNoDaughters(); ++i) {
    G4VPhysicalVolume* piece = gdmlWorldLV->GetDaughter(i);
    G4LogicalVolume* lv = piece->GetLogicalVolume();
    lv->SetVisAttributes(bodyVis);
    new G4PVPlacement(nullptr, origin + piece->GetTranslation(), lv, lv->GetName(), worldLV,
                      false, 0, false);
    ++nPlaced;
  }
  G4cout << "*** Vis-only enclosure body: " << nPlaced << " pieces in parallel world \""
         << GetName() << "\" at part origin " << origin
         << " (draw with /vis/drawVolume worlds). WARNING: with this parallel world registered"
         << " (-visbody 1), transport is CoupledTransportation and neutrons have NO GRAVITY." << G4endl;
}
