#include "TrackingAction.hh"

#include "EventAction.hh"

TrackingAction::TrackingAction(EventAction* eventAction) : fEventAction(eventAction) {}

void TrackingAction::PreUserTrackingAction(const G4Track* track)
{
    fEventAction->StartTrack(track);
}
