#ifndef FurmanPiviModel_h
#define FurmanPiviModel_h 1

#include "FurmanPiviParameters.hh"
#include "FurmanPiviTypes.hh"

#include "globals.hh"
#include "G4ThreeVector.hh"

#include <vector>


class FurmanPiviModel
{
public:
    explicit FurmanPiviModel(
        const FurmanPiviParameters& parameters);

    FurmanEmissionResult GenerateEmission(
        G4double incidentEnergy,
        G4double incidentAngle,
        const G4ThreeVector& incidentDirection,
        const G4ThreeVector& materialToVacuumNormal) const;

    FurmanYields ComputeYields(
        G4double incidentEnergy,
        G4double incidentAngle) const;

    G4double SampleElasticEnergy(G4double incidentEnergy) const;

    G4double SampleRediffusedEnergy(G4double incidentEnergy) const;

    std::vector<G4double> SampleTrueSecondaryEnergies(G4double incidentEnergy, G4int multiplicity) const;

    // Publique pour les outils de validation ; GenerateEmission l'utilise aussi.
    G4ThreeVector SampleDiffuseDirection(
        const G4ThreeVector& materialToVacuumNormal) const;

    /// Normale inversee, tangentielle conservee. La normale est perturbee au
    /// prealable si reflectionRoughness est non nul.
    G4ThreeVector SpecularDirection(
        const G4ThreeVector& incidentDirection,
        const G4ThreeVector& materialToVacuumNormal) const;

private:
    struct TrueSecondaryEnergySample
    {
        std::vector<G4double> energies;
    };

    TrueSecondaryEnergySample SampleTrueSecondaryEnergyGroup(
        G4double incidentEnergy,
        G4int multiplicity) const;

    G4double ComputeElasticYield(
        G4double incidentEnergy,
        G4double incidentAngle) const;

    G4double ComputeRediffusedYield(
        G4double incidentEnergy,
        G4double incidentAngle) const;

    /// Rendement TOTAL, Peng Eq. (2). La composante vraie s'en deduit par
    /// soustraction dans ComputeYields, interpretation retenue de Peng Eq. (7).
    G4double ComputeTotalYield(
        G4double incidentEnergy,
        G4double incidentAngle) const;



    std::vector<G4double> BuildMultiplicityProbabilities(
        const FurmanYields& yields) const;


    G4int SampleMultiplicity(
        const std::vector<G4double>& probabilities) const;


    SEEType SampleSingleElectronType(
        const FurmanYields& yields,
        G4double probabilityOneTrueSecondary) const;
    
    G4double SampleGamma(G4double shape, G4double scale) const;
    G4double GetTrueSecondaryShape(G4int multiplicity) const;
    G4double GetTrueSecondaryScale(G4int multiplicity) const;
    
        

private:
    FurmanPiviParameters fParameters;
};

#endif
