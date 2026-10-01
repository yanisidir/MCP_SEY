#include "ActionInitialization.hh"

#include "DetectorConstruction.hh"
#include "EventAction.hh"
#include "GaussianSource.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "SteppingAction.hh"
#include "TrackingAction.hh"

ActionInitialization::ActionInitialization(const DetectorConstruction* detector,
                                           const GaussianSource* source)
  : fDetector(detector), fSource(source)
{}

void ActionInitialization::BuildForMaster() const
{
    SetUserAction(new RunAction(fDetector, fSource));
}

void ActionInitialization::Build() const
{
    // Instance propre au thread : TrackingAction, SteppingAction et
    // FurmanPiviProcess ecrivent tous dedans.
    auto* eventAction = new EventAction(fDetector);

    SetUserAction(new PrimaryGeneratorAction(fDetector, fSource));
    SetUserAction(new RunAction(fDetector, fSource));
    SetUserAction(eventAction);
    SetUserAction(new TrackingAction(eventAction));
    SetUserAction(new SteppingAction(fDetector, eventAction));
}
