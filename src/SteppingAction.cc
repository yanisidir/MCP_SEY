#include "SteppingAction.hh"

#include "DetectorConstruction.hh"
#include "EventAction.hh"
#include "GenerationInfo.hh"

#include "G4Electron.hh"
#include "G4GeometryTolerance.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VPhysicalVolume.hh"

SteppingAction::SteppingAction(const DetectorConstruction* detector, EventAction* eventAction)
  : fDetector(detector), fEventAction(eventAction)
{}

void SteppingAction::UserSteppingAction(const G4Step* step)
{
    auto* track = step->GetTrack();
    const auto* pre = step->GetPreStepPoint();
    const auto* post = step->GetPostStepPoint();
    if (track->GetParticleDefinition() != G4Electron::Definition()
        || post->GetStepStatus() != fGeomBoundary || !pre->GetPhysicalVolume())
    {
        return;
    }

    const auto tolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
    
    auto* info = static_cast<GenerationInfo*>(track->GetUserInformation());

    if (pre->GetPhysicalVolume()->GetName() == "Channel"
        && post->GetPhysicalVolume() && post->GetPhysicalVolume()->GetName() == "Drift"
        && post->GetPosition().z() >= fDetector->Thickness() - tolerance
        && post->GetMomentumDirection().z() > 0 && !info->exitRecorded)
    {
        info->exitRecorded = true;
        fEventAction->RecordExit(track, post);  // On continue la meme trajectoire.
    }

    if (pre->GetPhysicalVolume()->GetName() == "Drift") {
        if (info->exitRecorded && post->GetMomentumDirection().z() > 0
            && post->GetPosition().z() >= fDetector->Thickness() + fDetector->DriftDistance() - tolerance)
        {
            fEventAction->RecordDownstream(track, post);
        }
        // Sortie du domaine de mesure, aval ou laterale.
        track->SetTrackStatus(fStopAndKill);
    }
}
