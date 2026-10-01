#ifndef FurmanPiviPhysics_h
#define FurmanPiviPhysics_h 1

#include "G4VPhysicsConstructor.hh"
class DetectorConstruction;
class FurmanPiviPhysics : public G4VPhysicsConstructor {
public:
    explicit FurmanPiviPhysics(const DetectorConstruction* detector) : fDetector(detector) {}
    void ConstructParticle() override;
    void ConstructProcess() override;
private:
    const DetectorConstruction* fDetector;
};

#endif
