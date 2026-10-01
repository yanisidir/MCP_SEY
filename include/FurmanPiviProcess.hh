#ifndef FurmanPiviProcess_h
#define FurmanPiviProcess_h 1

#include "G4VDiscreteProcess.hh"
#include "G4ParticleChange.hh"
#include "FurmanPiviModel.hh"
#include "WuModel.hh"
class DetectorConstruction;

// Une reflexion ramene aussi la navigation et le materiau dans le pore.
class ReflectionParticleChange : public G4ParticleChange {
public:
    G4bool returnToPore = false;
    G4Step* UpdateStepForPostStep(G4Step*) override;
};

class FurmanPiviProcess : public G4VDiscreteProcess {
public:
    explicit FurmanPiviProcess(const DetectorConstruction* detector);
    G4bool IsApplicable(const G4ParticleDefinition&) override;
    G4double GetMeanFreePath(const G4Track&, G4double, G4ForceCondition*) override;
    G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&) override;
private:
    G4ThreeVector InsidePore(const G4ThreeVector& impact, const G4ThreeVector& normal) const;
    const DetectorConstruction* fDetector;
    FurmanPiviParameters fParameters;
    FurmanPiviModel fModel;
    WuParameters fWuParameters;
    WuModel fWuModel;
    ReflectionParticleChange fChange;
};

#endif
