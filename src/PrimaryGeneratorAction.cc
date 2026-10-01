#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
#include "GaussianSource.hh"
#include "G4Electron.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include <algorithm>
#include <cmath>
#include <limits>

PrimaryGeneratorAction::PrimaryGeneratorAction(const DetectorConstruction* detector,
                                               const GaussianSource* source)
    : fGun(new G4ParticleGun(1)), fDetector(detector), fSource(source)
{
    fGun->SetParticleDefinition(G4Electron::Definition());
    fGun->SetParticleTime(0);
}
PrimaryGeneratorAction::~PrimaryGeneratorAction() { delete fGun; }

// Loi normale tronquee a [low, high] par rejet. Avec un sigma nul, ou si les
// tirages s'epuisent, la moyenne est rabattue dans l'intervalle.
G4double PrimaryGeneratorAction::TruncatedGauss(G4double mu, G4double sigma, G4double low, G4double high)
{
    if (sigma > 0) {
        for (int attempt = 0; attempt < 1000; ++attempt) {
            const auto x = G4RandGauss::shoot(mu, sigma);
            if (x >= low && x <= high) return x;
        }
    }
    return std::clamp(mu, low, high);
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
    const auto thickness = fDetector->Thickness();
    if (fSource->MuZ() < 0 || fSource->MuZ() > thickness) {
        G4Exception("PrimaryGeneratorAction", "Source", FatalException,
                    "La hauteur d'injection moyenne est hors de la plaque : verifier muZ et zOrigin.");
    }
    const auto energy = TruncatedGauss(fSource->MuE(), fSource->SigmaE(),
                                       0, std::numeric_limits<G4double>::infinity());
    const auto z = TruncatedGauss(fSource->MuZ(), fSource->SigmaZ(), 0, thickness);

    // pz est borne par l'impulsion totale que l'energie impose.
    const auto p = std::sqrt(energy * (energy + 2 * electron_mass_c2));
    
    const auto pz = TruncatedGauss(fSource->MuPz(), fSource->SigmaPz(), -p, p);

    const auto cosTheta = p > 0 ? pz / p : 1.0;
    const auto sinTheta = std::sqrt(std::max(0.0, 1 - cosTheta * cosTheta));
    G4ThreeVector direction(0, 0, cosTheta);
    if (sinTheta > 0) {  // Azimut uniforme ; pas de tirage pour un tir axial.
        const auto phi = twopi * G4UniformRand();
        direction.set(sinTheta * std::cos(phi), sinTheta * std::sin(phi), cosTheta);
    }

    // Sur l'axe du pore : incline vers +y, il passe par l'origine.
    fGun->SetParticlePosition({0, z * std::tan(fDetector->Angle()), z});
    fGun->SetParticleMomentumDirection(direction);
    fGun->SetParticleEnergy(energy);
    fGun->GeneratePrimaryVertex(event);
}
