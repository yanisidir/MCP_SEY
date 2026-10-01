#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class DetectorConstruction;
class G4Event;
class G4StepPoint;
class G4Track;

/// Compteurs et statistiques de temps propres a un evenement, c'est-a-dire a
/// une avalanche. L'objet appartient au worker : TrackingAction, SteppingAction
/// et FurmanPiviProcess l'alimentent au fil du suivi, puis EndOfEventAction
/// ecrit la ligne resumee dans la table "events".
class EventAction : public G4UserEventAction
{
  public:
    explicit EventAction(const DetectorConstruction* detector);
    ~EventAction() override = default;

    void BeginOfEventAction(const G4Event*) override;
    void EndOfEventAction(const G4Event*) override;

    /// Nouvelle trajectoire suivie : marquage de generation et comptage.
    void StartTrack(const G4Track* track);

    /// Places encore disponibles sous le plafond de tracks de l'evenement ;
    /// retourne le nombre de secondaires que l'appelant peut reellement creer.
    G4int ReserveSecondaries(G4int requested);

    /// Premier passage d'un electron par la face aval de la plaque.
    void RecordExit(const G4Track* track, const G4StepPoint* point);

    /// Passage par le plan de mesure, a +50 um de la face aval.
    void RecordDownstream(const G4Track* track, const G4StepPoint* point);

    /// Impact d'un electron sur la paroi du canal : energie reconstruite et
    /// angle d'incidence. Instrumentation seule, sans effet sur le transport
    /// ni sur l'emission. L'energie d'impact est la grandeur qui fixe la pente
    /// du gain, k = dln(delta)/dE evaluee a cette energie.
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

    // Impacts sur la paroi : comptage et sommes pour les moyennes par evenement.
    G4long fImpacts = 0;
    G4double fImpactEnergy = 0.;
    G4double fImpactAngle = 0.;

    // Temps d'arrivee : extremes et moments courants (algorithme en ligne).
    G4double fFirst = 0.;
    G4double fLast = 0.;
    G4double fMean = 0.;
    G4double fM2 = 0.;
};

#endif
