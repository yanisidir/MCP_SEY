#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
class G4Event;
class G4ParticleGun;
class DetectorConstruction;
class GaussianSource;

// Un electron par evenement, tire selon les trois lois de GaussianSource,
// injecte sur l'axe du pore a la hauteur z tiree.
class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
    PrimaryGeneratorAction(const DetectorConstruction* detector, const GaussianSource* source);
    ~PrimaryGeneratorAction() override;
    void GeneratePrimaries(G4Event*) override;
private:
    static G4double TruncatedGauss(G4double mu, G4double sigma, G4double low, G4double high);
    G4ParticleGun* fGun;
    const DetectorConstruction* fDetector;
    const GaussianSource* fSource;
};

#endif
