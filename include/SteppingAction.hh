#ifndef SteppingAction_h
#define SteppingAction_h 1

#include "G4UserSteppingAction.hh"

class DetectorConstruction;
class EventAction;
class G4Step;

/// Detection des deux franchissements qui comptent : la face aval de la
/// plaque, puis le plan de mesure a +50 um ou la trajectoire est arretee.
class SteppingAction : public G4UserSteppingAction
{
  public:
    SteppingAction(const DetectorConstruction* detector, EventAction* eventAction);
    ~SteppingAction() override = default;

    void UserSteppingAction(const G4Step*) override;

  private:
    const DetectorConstruction* fDetector = nullptr;
    EventAction* fEventAction = nullptr;
};

#endif
