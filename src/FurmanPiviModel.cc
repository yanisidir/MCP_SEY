#include "FurmanPiviModel.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <vector>

FurmanPiviModel::FurmanPiviModel(const FurmanPiviParameters& parameters)
    : fParameters(parameters) {}

FurmanYields FurmanPiviModel::ComputeYields(G4double energy, G4double angle) const
{
    const G4double elastic = std::max(0.0, ComputeElasticYield(energy, angle));
    const G4double rediffused = std::max(0.0, ComputeRediffusedYield(energy, angle));

    // Interpretation retenue de Peng Eq. (7) : Eq. (2) donne le rendement TOTAL
    // et la composante vraie s'obtient par soustraction des reflexions. La
    // troncature a zero porte le seuil d'emission E'. Cette lecture n'a pas pu
    // etre confirmee : avec les parametres de Peng dans sa geometrie, le gain
    // obtenu reste bien inferieur a celui qu'il publie (cf. WuModel.hh).
    const G4double total = std::max(0.0, ComputeTotalYield(energy, angle));
    const G4double trueSecondary = std::max(0.0, total - elastic - rediffused);

    return {elastic, rediffused, trueSecondary};
}

G4double FurmanPiviModel::ComputeElasticYield(
    G4double incidentEnergy,
    G4double incidentAngle) const
{
    if (fParameters.elasticYieldWidth <= 0.0 ||
        fParameters.elasticYieldShape <= 0.0) {
        throw std::runtime_error("Furman-Pivi: elastic yield width and shape must be positive.");
    }

    const G4double scaledEnergy =
        std::abs(
            incidentEnergy
            - fParameters.elasticPeakEnergy)
        / fParameters.elasticYieldWidth;

    const G4double deltaNormal =
        fParameters.elasticHighEnergyYield
        + (
            fParameters.elasticPeakYield
            - fParameters.elasticHighEnergyYield
          )
        * std::exp(
            -std::pow(
                scaledEnergy,
                fParameters.elasticYieldShape)
            / fParameters.elasticYieldShape);

    const G4double cosTheta =
        std::clamp(
            std::cos(incidentAngle),
            0.0,
            1.0);

    const G4double angularFactor =
        1.0
        + fParameters.elasticAngular1
        * (
            1.0
            - std::pow(
                cosTheta,
                fParameters.elasticAngular2)
          );

    return deltaNormal * angularFactor;
}

G4double FurmanPiviModel::ComputeRediffusedYield(
    G4double incidentEnergy,
    G4double incidentAngle) const
{
    const G4double energy = std::max(0.0, incidentEnergy);

    if (fParameters.rediffusedYieldEnergyScale <= 0.0 ||
        fParameters.rediffusedYieldShape <= 0.0) {
        throw std::runtime_error("Furman-Pivi: rediffused yield energy scale and shape must be positive.");
    }

    const G4double deltaNormal =
        fParameters.rediffusedHighEnergyYield
        * (
            1.0
            - std::exp(
                -std::pow(
                    energy
                    / fParameters.rediffusedYieldEnergyScale,
                    fParameters.rediffusedYieldShape))
          );

    const G4double cosTheta =
        std::clamp(
            std::cos(incidentAngle),
            0.0,
            1.0);

    const G4double angularFactor =
        1.0
        + fParameters.rediffusedAngular1
        * (
            1.0
            - std::pow(
                cosTheta,
                fParameters.rediffusedAngular2)
          );

    return deltaNormal * angularFactor;
}

G4double FurmanPiviModel::ComputeTotalYield(
    G4double incidentEnergy,
    G4double incidentAngle) const
{
    if (incidentEnergy <= 0.0) {
        return 0.0;
    }

    const G4double cosTheta = std::clamp(std::cos(incidentAngle), 0.0, 1.0);

    const G4double deltaPeak =
        fParameters.totalPeakYield
        * (
            1.0
            + fParameters.totalAngular1
            * (
                1.0
                - std::pow(
                    cosTheta,
                    fParameters.totalAngular2)
              )
          );

    const G4double peakEnergy =
        fParameters.totalPeakIncidentEnergy
        * (
            1.0
            + fParameters.totalAngular3
            * (
                1.0
                - std::pow(
                    cosTheta,
                    fParameters.totalAngular4)
              )
          );

    if (peakEnergy <= 0.0) {return 0.0;}

    const G4double x = incidentEnergy / peakEnergy;

    const G4double s = fParameters.totalYieldShape;

    if (s <= 1.0) {
        throw std::runtime_error(
            "Furman-Pivi: total yield "
            "shape must be greater than 1.");
    }
    const G4double denominator = s - 1.0 + std::pow(x, s);

    if (denominator <= 0.0) {return 0.0;}

    const G4double universalFunction = s * x / denominator;

    return deltaPeak * universalFunction;
}

std::vector<G4double>
FurmanPiviModel::BuildMultiplicityProbabilities(
    const FurmanYields& yields) const
{
    const G4int maximumMultiplicity = fParameters.maximumMultiplicity;

    if (maximumMultiplicity <= 0) {throw std::runtime_error("maximumMultiplicityy must be positive.");}

    std::vector<G4double> probabilities(maximumMultiplicity + 1, 0.0); // M + proba = 0 

    const G4double reflectionProbability = yields.elastic + yields.rediffused;

    if (reflectionProbability < 0.0 ||
        reflectionProbability > 1.0) {
        throw std::runtime_error("Invalid elastic + rediffused probability.");
    }

    const G4double penetrationProbability = 1.0 - reflectionProbability;

    if (penetrationProbability <= 0.0) {
        probabilities[1] = 1.0;
        return probabilities;
    }

    const G4double conditionalMean = yields.trueSecondary / penetrationProbability;

    if (conditionalMean < 0.0 ||
        conditionalMean > static_cast<G4double>(maximumMultiplicity)) {
        throw std::runtime_error(
            "True-secondary conditional mean "
            "outside binomial range.");
    }

    std::vector<G4double> trueSecondaryProbabilities(maximumMultiplicity + 1, 0.0);

    /*
    * Binomial Law
    */

    const G4double p = conditionalMean / static_cast<G4double>(maximumMultiplicity);

    if (p <= 0.0) {trueSecondaryProbabilities[0] = 1.0;}
    else if (p >= 1.0) {trueSecondaryProbabilities[maximumMultiplicity] = 1.0;}
    else {
        trueSecondaryProbabilities[0] = std::pow(1.0 - p, maximumMultiplicity);

        for (G4int n = 1; n <= maximumMultiplicity; ++n) {
            trueSecondaryProbabilities[n] = trueSecondaryProbabilities[n-1] * static_cast<G4double>(maximumMultiplicity - n + 1)
                                                                            / static_cast<G4double>(n) * p / (1.0 - p);
        }
    }

    probabilities[0] = penetrationProbability * trueSecondaryProbabilities[0];

    probabilities[1] = reflectionProbability + penetrationProbability * trueSecondaryProbabilities[1];

    for (G4int n = 2; n <= maximumMultiplicity; ++n) {
        probabilities[n] = penetrationProbability * trueSecondaryProbabilities[n];
    }

    const G4double probabilitySum = std::accumulate(probabilities.begin(), probabilities.end(), static_cast<G4double>(0.0));

    if (std::abs(probabilitySum -1.0) > 1.0e-12 )  {
        throw std::runtime_error("Multiplicity probabilities do not sum to one.");
    }

    return probabilities;
}

G4int FurmanPiviModel::SampleMultiplicity(
    const std::vector<G4double>& probabilities) const
{
    const G4double randomNumber = G4UniformRand();

    G4double cumulativeProbability = 0.0;

    for (std::size_t n = 0; n < probabilities.size(); ++n) {
        cumulativeProbability += probabilities[n];

        if (randomNumber < cumulativeProbability) {return static_cast<G4int>(n);}
    }

    return static_cast<G4int>(probabilities.size() - 1);
}

SEEType FurmanPiviModel::SampleSingleElectronType(
    const FurmanYields& yields,
    G4double probabilityOneTrueSecondary) const
{
    const G4double totalWeight = yields.elastic + yields.rediffused + probabilityOneTrueSecondary;

    if (totalWeight <= 0.0) {return SEEType::TrueSecondary;}

    const G4double randomNumber = G4UniformRand() * totalWeight;

    if (randomNumber < yields.elastic) {return SEEType::Elastic;}

    if (randomNumber < yields.elastic + yields.rediffused) {return SEEType::Rediffused;}

    return SEEType::TrueSecondary;
}

FurmanEmissionResult FurmanPiviModel::GenerateEmission(
    G4double incidentEnergy,
    G4double incidentAngle,
    const G4ThreeVector& incidentDirection,
    const G4ThreeVector& materialToVacuumNormal) const
{
    FurmanEmissionResult result;

    const FurmanYields yields = ComputeYields(incidentEnergy, incidentAngle);
    
    if (yields.elastic < 0.0 ||
        yields.rediffused < 0.0 ||
        yields.trueSecondary < 0.0) {
        throw std::runtime_error("Negative Furman-Pivi yield.");
    }
    const std::vector<G4double> probabilities = BuildMultiplicityProbabilities(yields);

    result.multiplicity = SampleMultiplicity(probabilities);

    if (result.multiplicity == 0) {return result;}

    const G4double penetrationProbability = std::max( 0.0, 1.0 - yields.elastic - yields.rediffused);

    const G4double conditionalMean =
        penetrationProbability > 0.0
        ? yields.trueSecondary / penetrationProbability
        : 0.0;

    const G4int M = fParameters.maximumMultiplicity;

    const G4double p = conditionalMean / M;

    G4double probabilityOneTrueSecondary = 0.0;

    if (p > 0.0 && p < 1.0) {
        probabilityOneTrueSecondary = 
            penetrationProbability
            * static_cast<G4double>(M)
            * p 
            * std::pow(1.0 - p, M - 1);
    }
    else if ( p >= 1.0 && M == 1) {
        probabilityOneTrueSecondary = penetrationProbability;
    }

    if (result.multiplicity == 1) {
        const SEEType type = SampleSingleElectronType(yields, probabilityOneTrueSecondary);
        EmittedElectron electron;
        electron.type = type;

        if (type == SEEType::Elastic) {
            electron.kineticEnergy = SampleElasticEnergy(incidentEnergy);
            // Wu et al., Rev. Sci. Instrum. 79 (2008) 073104 : "the axial and
            // angular components of the electrons' velocity are left unchanged,
            // while the radial component is reversed".
            electron.direction = SpecularDirection(incidentDirection, materialToVacuumNormal);
        }
        else if (type == SEEType::Rediffused) {
            electron.kineticEnergy = SampleRediffusedEnergy(incidentEnergy);
            electron.direction = SampleDiffuseDirection(materialToVacuumNormal);
        }
        else {
            const auto sample =
                SampleTrueSecondaryEnergyGroup(incidentEnergy, 1);
            electron.kineticEnergy = sample.energies.front();
            electron.direction = SampleDiffuseDirection(materialToVacuumNormal);
        }

        result.electrons.push_back(electron);

        const G4double totalEmittedEnergy = 
            std::accumulate(result.electrons.begin(), result.electrons.end(), 0.0,
                            [](G4double sum, 
                               const EmittedElectron& electron) {
                                return sum + electron.kineticEnergy;
                               });
        const G4double tolerance = 
            std::max(0.0, fParameters.energyConservationRelativeTolerance) * std::max(incidentEnergy, 1.0 * eV);
        
        if (totalEmittedEnergy > incidentEnergy + tolerance) {
            throw std::runtime_error("Furman-Pivi energy conservation Violation");
        }
        return result;
    }

    
    const auto sample =
        SampleTrueSecondaryEnergyGroup(incidentEnergy, result.multiplicity);
    const auto& energies = sample.energies;

    result.electrons.clear();
    result.electrons.reserve(energies.size());

    for (const G4double energy : energies) {
        EmittedElectron electron;

        electron.type = SEEType::TrueSecondary;
        electron.kineticEnergy = energy;
        electron.direction = SampleDiffuseDirection(materialToVacuumNormal);
        result.electrons.push_back(electron);
    }
    
    return result;
}

G4ThreeVector FurmanPiviModel::SpecularDirection(
    const G4ThreeVector& incidentDirection,
    const G4ThreeVector& materialToVacuumNormal) const
{
    const G4ThreeVector normal = materialToVacuumNormal.unit();
    const G4ThreeVector incident = incidentDirection.unit();

    // Microfacettes : reflexion autour d'une normale perturbee. A rugosite
    // nulle aucun tirage n'est consomme.
    G4ThreeVector axis = normal;
    const G4double roughness = fParameters.reflectionRoughness;
    if (roughness > 0.0) {
        const G4double theta = std::abs(G4RandGauss::shoot(0.0, roughness));
        const G4double phi = twopi * G4UniformRand();
        const G4ThreeVector tangent1 = normal.orthogonal().unit();
        const G4ThreeVector tangent2 = normal.cross(tangent1).unit();
        axis = (std::cos(theta) * normal
                + std::sin(theta) * (std::cos(phi) * tangent1 + std::sin(phi) * tangent2)).unit();
    }

    G4ThreeVector reflected = incident - 2.0 * incident.dot(axis) * axis;
    if (reflected.mag2() <= 0.0) {return normal;}
    reflected = reflected.unit();

    // Une normale inclinee peut renvoyer la direction dans le materiau.
    const G4double along = reflected.dot(normal);
    if (along <= 0.0) {reflected = (reflected - 2.0 * along * normal).unit();}
    return reflected;
}

G4double FurmanPiviModel::SampleElasticEnergy(
    G4double incidentEnergy) const
{
    if (incidentEnergy <= 0.0) {return 0.0;}

    if (fParameters.elasticEnergySigma <= 0.0) {return incidentEnergy;}

    const G4int maximumAttempts = std::max(1, fParameters.maximumSamplingAttempts);

    for (G4int attempt = 0; attempt < maximumAttempts; ++attempt) {

        const G4double gaussian = G4RandGauss::shoot(0.0, 1.0);

        const G4double energy = 
            incidentEnergy - fParameters.elasticEnergySigma * std::abs(gaussian);
        
        if (energy >= 0.0) {return energy;}

    }

    throw std::runtime_error("SampleElasticEnergy: rejection sampler failed");

}

G4double FurmanPiviModel::SampleRediffusedEnergy(G4double incidentEnergy) const
{
    if (incidentEnergy <= 0.0) {return 0.0;}

    const G4double exponent = fParameters.rediffusedEnergyExponent;

    if (exponent <= -1.0) {
        throw std::runtime_error("SampleRediffusedEnergy: q must be greater than -1.");
    }

    const G4double u = std::max(G4UniformRand(), std::numeric_limits<G4double>::min());

    return incidentEnergy * std::pow(u, 1.0 / (exponent + 1.0));

}

std::vector<G4double> 
FurmanPiviModel::SampleTrueSecondaryEnergies(G4double incidentEnergy,
                                             G4int multiplicity) const
{
    return SampleTrueSecondaryEnergyGroup(
        incidentEnergy, multiplicity).energies;
}

FurmanPiviModel::TrueSecondaryEnergySample
FurmanPiviModel::SampleTrueSecondaryEnergyGroup(
    G4double incidentEnergy,
    G4int multiplicity) const
{
    if (incidentEnergy <= 0.0) {throw std::invalid_argument("Incident energy must be positive.");}

    if (multiplicity <= 0) {throw std::invalid_argument("Multiplicity must be positive.");}

    const G4double shape = GetTrueSecondaryShape(multiplicity);
    const G4double scale = GetTrueSecondaryScale(multiplicity);

    const G4int maximumAttempts =
        std::max(1, fParameters.maximumSamplingAttempts);
    const G4int directAttempts = std::clamp(
        fParameters.directSamplingAttempts,
        0,
        maximumAttempts);
    
    TrueSecondaryEnergySample sample;
    auto& energies = sample.energies;
    energies.resize(static_cast<std::size_t>(multiplicity));

    // Chemin rapide : tout groupe accepte suit la loi voulue, produit de
    // Gamma conditionne a sum(E_i) <= E0.
    for (G4int attempt = 0; attempt < directAttempts; ++attempt)
    {
        G4double totalEnergy = 0.0;

        for (auto& energy : energies) {
            energy = SampleGamma(shape,scale);
            totalEnergy += energy;
        }

        /*
        * Si la somme est autorisée, chaque énergie 
        * individuelle est automatiquement <= E0.
        */
       if (totalEnergy <= incidentEnergy) {

           return sample;
       }

    }

    // Le rejet direct devient inefficace a haute multiplicite et basse energie.
    // La somme suit Gamma(m * shape, scale) et les fractions normalisees une
    // Dirichlet independante : on tire cette loi conditionnelle par une
    // proposition en puissance adaptative. Le rejet reste exact.
    const G4double totalShape =
        static_cast<G4double>(multiplicity) * shape;
    const G4double matchedExponent =
        totalShape - incidentEnergy / scale;
    const G4bool useMatchedProposal = matchedExponent > 0.0;
    const G4double proposalExponent =
        useMatchedProposal ? matchedExponent : totalShape;

    G4double totalEnergy = 0.0;
    G4bool acceptedTotal = false;
    for (G4int attempt = 0; attempt < maximumAttempts; ++attempt) {
        const G4double uniform = std::max(
            G4UniformRand(),
            std::numeric_limits<G4double>::min());
        totalEnergy = incidentEnergy *
            std::pow(uniform, 1.0 / proposalExponent);

        G4double acceptance = std::exp(-totalEnergy / scale);
        if (useMatchedProposal) {
            acceptance =
                std::pow(
                    totalEnergy / incidentEnergy,
                    totalShape - proposalExponent)
                * std::exp((incidentEnergy - totalEnergy) / scale);
        }
        if (G4UniformRand() <= acceptance) {
            acceptedTotal = true;
            break;
        }
    }
    if (!acceptedTotal) {
        throw std::runtime_error(
            "SampleTrueSecondaryEnergies: truncated gamma sampler failed.");
    }

    G4double fractionSum = 0.0;
    for (auto& energy : energies) {
        energy = SampleGamma(shape, 1.0);
        fractionSum += energy;
    }
    if (fractionSum <= 0.0) {
        throw std::runtime_error(
            "SampleTrueSecondaryEnergies: invalid Dirichlet weights.");
    }

    for (auto& energy : energies) {
        energy = totalEnergy * energy / fractionSum;
    }

    return sample;
}

G4double FurmanPiviModel::SampleGamma(
    G4double shape,
    G4double scale) const
{
    if (shape <= 0.0 || scale <= 0.0) {
        throw std::invalid_argument(
            "SampleGamma: shape and scale must be positive."
        );
    }

    /*
     * Pour 0 < shape < 1 :
     *
     * Gamma(shape) peut être obtenu à partir de
     * Gamma(shape + 1) et d'une puissance uniforme.
    */

    if (shape < 1.0) {
        const G4double u = std::max(G4UniformRand(), std::numeric_limits<G4double>::min());

        return SampleGamma(shape + 1.0, scale) * std::pow(u, 1.0 / shape);
    } 

    /*
     * Algorithme de Marsaglia et Tsang pour shape >= 1.
    */

    const G4double d = shape - 1.0 / 3.0;
    const G4double c = 1.0 / std::sqrt(9.0 * d);

    while (true) {
        const G4double x = G4RandGauss::shoot(0.0, 1.0);
        const G4double factor = 1.0 + c * x;
        
        if (factor <= 0.0) {continue;}

        const G4double v = factor * factor * factor;

        const G4double u = G4UniformRand();

        if (u < 1.0 - 0.0331 * x * x * x* x) {return scale * d * v;}

        if (std::log(u) < 0.5 * x * x + d * (1.0 - v + std::log(v))) {
            return scale * d * v;
        }
    }
}

G4double FurmanPiviModel::GetTrueSecondaryShape(
    G4int multiplicity) const
{
    if (multiplicity <= 0) {
        throw std::invalid_argument("Multiplicity must be positive.");
    }

    const std::size_t index = std::min<std::size_t>(
        static_cast<std::size_t>(multiplicity - 1), fParameters.trueSecondaryEnergyShape.size() - 1);
    
    return fParameters.trueSecondaryEnergyShape[index];
}

G4double FurmanPiviModel::GetTrueSecondaryScale(
    G4int multiplicity) const
{
    if (multiplicity <= 0) {
        throw std::invalid_argument("Multiplicity must be positive.");
    }

    const std::size_t index = std::min<std::size_t>(
        static_cast<std::size_t>(multiplicity - 1), fParameters.trueSecondaryEnergyScale.size() - 1);
    
    return fParameters.trueSecondaryEnergyScale[index];
    
}

G4ThreeVector FurmanPiviModel::SampleDiffuseDirection(
    const G4ThreeVector& materialToVacuumNormal) const
{
    G4double cosTheta = 1.0;
    G4double sinTheta = 0.0;
    G4double phi = 0.0;

    // Loi cosinus : deux tirages, comme dans la simulation d’origine.
    const G4double alpha = fParameters.emissionAngularExponent;
    if (alpha <= -1.0) throw std::runtime_error("Angular exponent must exceed -1");
    const G4double u1 = G4UniformRand();
    const G4double u2 = G4UniformRand();
    cosTheta = std::pow(u1, 1.0 / (alpha + 1.0));
    sinTheta = std::sqrt(std::max(0.0, 1.0 - cosTheta * cosTheta));
    phi = twopi * u2;

    const G4ThreeVector normal =
        materialToVacuumNormal.unit();

    const G4ThreeVector tangent1 =
        normal.orthogonal().unit();

    const G4ThreeVector tangent2 =
        normal.cross(tangent1).unit();

    return (
        sinTheta * std::cos(phi) * tangent1
        + sinTheta * std::sin(phi) * tangent2
        + cosTheta * normal
    ).unit();
}

