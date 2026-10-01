#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class DetectorConstruction;
class G4Event;
class G4StepPoint;
class G4Track;

/// Compteurs d'une avalanche. Alimente par TrackingAction, SteppingAction et
/// FurmanPiviProcess ; EndOfEventAction ecrit la ligne de la table "events".
class EventAction : public G4UserEventAction
{
  public:
    explicit EventAction(const DetectorConstruction* detector);
    ~EventAction() override = default;

    void BeginOfEventAction(const G4Event*) override;
    void EndOfEventAction(const G4Event*) override;

    void StartTrack(const G4Track* track);

    /// Nombre de secondaires que l'appelant peut creer sous le plafond.
    G4int ReserveSecondaries(G4int requested);

    void RecordExit(const G4Track* track, const G4StepPoint* point);

    void RecordDownstream(const G4Track* track, const G4StepPoint* point);

    /// Instrumentation seule : ne consomme aucun tirage aleatoire.
    void RecordImpact(G4double energy, G4double angle);

  private:
    const DetectorConstruction* fDetector = nullptr;

    G4int fEventID = 0;
    G4int fPrimary = 0;
    G4int fSecondary = 0;
    G4int fMaxGeneration = 0;
    G4int fAt50um = 0;
    G4long fTracks = 0;
    G4long fCreatedTracks = 1;
    G4bool fMultiplicationLimited = false;

    G4double fT0 = 0.;
    G4double fEnergy = 0.;

    G4long fImpacts = 0;
    G4double fImpactEnergy = 0.;
    G4double fImpactAngle = 0.;

    // Moments courants, algorithme en ligne.
    G4double fFirst = 0.;
    G4double fLast = 0.;
    G4double fMean = 0.;
    G4double fM2 = 0.;
};

#endif
