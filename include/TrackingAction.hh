#ifndef TrackingAction_h
#define TrackingAction_h 1

#include "G4UserTrackingAction.hh"

class EventAction;
class G4Track;

/// Point d'entree de chaque nouvelle trajectoire : c'est la que la generation
/// d'emission est posee sur le primaire et que les tracks sont comptes.
class TrackingAction : public G4UserTrackingAction
{
  public:
    explicit TrackingAction(EventAction* eventAction);
    ~TrackingAction() override = default;

    void PreUserTrackingAction(const G4Track* track) override;

  private:
    EventAction* fEventAction = nullptr;
};

#endif
