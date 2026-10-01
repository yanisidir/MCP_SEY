#include "WuModel.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"
#include "G4Poisson.hh"

#include <algorithm>
#include <cmath>

WuModel::WuModel(const WuParameters& parameters) : fParameters(parameters) {}

G4double WuModel::ReflectedFraction(G4double incidentEnergy) const
{
    const G4double energy = incidentEnergy / eV;
    if (energy <= 0.0) return 0.0;

    // Wu Eq. (6), ajustement de Scholtz. Les deux ancrages cites dans l'article
    // sont restitues : f0(10 eV) = 0.84 et f0(100 eV) = 0.046, pour "~95 %" et
    // "~5 %" annonces.
    const G4double l = std::log(energy);
    const G4double f0 =
        std::exp(1.59 + 3.75 * l - 1.37 * l * l + 0.12 * l * l * l) / 100.0;

    return std::clamp(fParameters.reflectionScale * f0, 0.0, 1.0);
}

G4double WuModel::TrueSecondaryYield(G4double incidentEnergy, G4double incidentAngle) const
{
    const G4double cosTheta = std::clamp(std::cos(incidentAngle), 1.0e-6, 1.0);

    // Eq. (2) et la loi de Vm(theta) qui la precede.
    const G4double peakEnergy = fParameters.peakEnergy / std::sqrt(cosTheta);
    const G4double peakYield =
        fParameters.peakYield * std::exp(fParameters.angularYield * (1.0 - cosTheta));

    // Eq. (3) : fonction universelle s x / (s - 1 + x^s).
    const G4double x = incidentEnergy / peakEnergy;
    if (x <= 0.0) return 0.0;
    const G4double s = fParameters.yieldShape;
    return peakYield * s * x / (s - 1.0 + std::pow(x, s));
}

G4double WuModel::TotalYield(G4double incidentEnergy, G4double incidentAngle) const
{
    const G4double trueYield = TrueSecondaryYield(incidentEnergy, incidentAngle);
    const G4double f = ReflectedFraction(incidentEnergy);
    const G4double denominator = 1.0 - f + trueYield * f;
    if (denominator <= 0.0) return trueYield;
    return trueYield / denominator;                       // Eq. (7)
}

G4double WuModel::SampleEmissionEnergy() const
{
    // Eq. (4) : P(E) proportionnel a (E/E0) exp(-E/E0), soit une loi Gamma de
    // forme 2, c'est-a-dire la somme de deux exponentielles.
    const G4double u1 = std::max(G4UniformRand(), 1.0e-300);
    const G4double u2 = std::max(G4UniformRand(), 1.0e-300);
    return -fParameters.emissionEnergy * (std::log(u1) + std::log(u2));
}

G4ThreeVector WuModel::SampleCosineDirection(const G4ThreeVector& materialToVacuumNormal) const
{
    const G4double alpha = fParameters.emissionAngularExponent;
    const G4double cosTheta = std::pow(G4UniformRand(), 1.0 / (alpha + 1.0));
    const G4double sinTheta = std::sqrt(std::max(0.0, 1.0 - cosTheta * cosTheta));
    const G4double phi = twopi * G4UniformRand();

    const G4ThreeVector normal = materialToVacuumNormal.unit();
    const G4ThreeVector tangent1 = normal.orthogonal().unit();
    const G4ThreeVector tangent2 = normal.cross(tangent1).unit();

    return (sinTheta * std::cos(phi) * tangent1
            + sinTheta * std::sin(phi) * tangent2
            + cosTheta * normal).unit();
}

G4ThreeVector WuModel::SpecularDirection(
    const G4ThreeVector& incidentDirection,
    const G4ThreeVector& materialToVacuumNormal) const
{
    const G4ThreeVector normal = materialToVacuumNormal.unit();
    const G4ThreeVector incident = incidentDirection.unit();
    G4ThreeVector reflected = incident - 2.0 * incident.dot(normal) * normal;
    if (reflected.mag2() <= 0.0) return normal;
    return reflected.unit();
}

FurmanEmissionResult WuModel::GenerateEmission(
    G4double incidentEnergy,
    G4double incidentAngle,
    const G4ThreeVector& incidentDirection,
    const G4ThreeVector& materialToVacuumNormal) const
{
    FurmanEmissionResult result;
    if (incidentEnergy <= 0.0) return result;

    const G4double trueYield = TrueSecondaryYield(incidentEnergy, incidentAngle);
    const G4double totalYield = TotalYield(incidentEnergy, incidentAngle);
    const G4double f = ReflectedFraction(incidentEnergy);

    // Eq. (8) : delta_e = f delta_t est la PROBABILITE qu'un electron soit
    // reflechi elastiquement. Le complement tire une multiplicite de Poisson
    // de moyenne delta_s : la moyenne totale vaut alors exactement delta_t,
    // ce que Eq. (7) est construite pour assurer.
    const G4double reflectionProbability = std::clamp(f * totalYield, 0.0, 1.0);

    if (G4UniformRand() < reflectionProbability) {
        EmittedElectron electron;
        electron.type = SEEType::Elastic;
        electron.kineticEnergy = incidentEnergy;      // reflexion elastique
        electron.direction = SpecularDirection(incidentDirection, materialToVacuumNormal);
        result.electrons.push_back(electron);
        result.multiplicity = 1;
        return result;
    }

    G4int multiplicity = static_cast<G4int>(G4Poisson(trueYield));
    if (multiplicity <= 0) return result;

    // Conservation de l'energie : Wu retire l'ensemble du groupe tant que la
    // somme depasse l'energie d'impact. Si la multiplicite tiree est trop
    // grande pour l'energie disponible, la retirer indefiniment ne convergerait
    // pas : on la reduit d'un cran, ce qui revient a dire que l'energie
    // d'impact ne peut pas alimenter autant d'electrons.
    std::vector<G4double> energies;
    while (multiplicity > 0) {
        bool accepted = false;
        for (G4int attempt = 0; attempt < fParameters.maximumSamplingAttempts; ++attempt) {
            energies.clear();
            G4double sum = 0.0;
            for (G4int i = 0; i < multiplicity; ++i) {
                const G4double energy = SampleEmissionEnergy();
                energies.push_back(energy);
                sum += energy;
            }
            if (sum < incidentEnergy) { accepted = true; break; }
        }
        if (accepted) break;
        --multiplicity;
    }
    if (multiplicity <= 0) return result;

    for (G4int i = 0; i < multiplicity; ++i) {
        EmittedElectron electron;
        electron.type = SEEType::TrueSecondary;
        electron.kineticEnergy = energies[i];
        electron.direction = SampleCosineDirection(materialToVacuumNormal);
        result.electrons.push_back(electron);
    }
    result.multiplicity = multiplicity;
    return result;
}
