#ifndef WuModel_h
#define WuModel_h 1

#include "FurmanPiviTypes.hh"
#include "WuParameters.hh"

#include "globals.hh"
#include "G4ThreeVector.hh"

#include <vector>

// Emission secondaire selon Wu et al., Rev. Sci. Instrum. 79 (2008) 073104.
// Deux composantes : reflexion elastique et vrais secondaires. Eq. (3) donne le
// rendement VRAI et la reflexion s'y ajoute par Eq. (7), au contraire de
// FurmanPiviModel ou Eq. (2) est le rendement TOTAL dont on soustrait les
// reflexions. Meme interface, pour que le processus passe de l'un a l'autre.
class WuModel
{
  public:
    explicit WuModel(const WuParameters& parameters);

    FurmanEmissionResult GenerateEmission(
        G4double incidentEnergy,
        G4double incidentAngle,
        const G4ThreeVector& incidentDirection,
        const G4ThreeVector& materialToVacuumNormal) const;

    // Scholtz et al., Philips J. Res. 50 (1996) 375, via Wu Eq. (6), x R.
    G4double ReflectedFraction(G4double incidentEnergy) const;

    G4double TrueSecondaryYield(G4double incidentEnergy, G4double incidentAngle) const;
    G4double TotalYield(G4double incidentEnergy, G4double incidentAngle) const;

  private:
    G4double SampleEmissionEnergy() const;
    G4ThreeVector SampleCosineDirection(const G4ThreeVector& materialToVacuumNormal) const;
    G4ThreeVector SpecularDirection(
        const G4ThreeVector& incidentDirection,
        const G4ThreeVector& materialToVacuumNormal) const;

    WuParameters fParameters;
};

#endif
