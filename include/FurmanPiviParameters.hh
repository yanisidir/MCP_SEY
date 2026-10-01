#ifndef FurmanPiviParameters_h
#define FurmanPiviParameters_h 1

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"

#include <array>

struct FurmanPiviParameters
{
    // Rendements : valeurs fournies le 16/09/2026.
    // Elastic yield component
    // ========================================================

    // Ajustement sur les courbes digitalisees de Peng et al. Fig. 4 (peng-sey-fit).
    // Ces valeurs coincident avec les parametres cuivre de Furman-Pivi 2002
    // (0.02 / 0.496 / 0 / 60.86 eV / 1.0), Peng et al. n'ayant pas mesure les
    // composantes de reflexion du verre au plomb separement.
    G4double elasticHighEnergyYield = 0.0197991598770481;
    G4double elasticPeakYield = 0.500809176899711;

    G4double elasticPeakEnergy = 0.0 * eV;
    G4double elasticYieldWidth = 60.7290930799462 * eV;

    G4double elasticYieldShape = 0.985795034334466;

    G4double elasticAngular1 = 0.26;
    G4double elasticAngular2 = 2.0;

    // ========================================================
    // Elastic outgoing-energy distribution
    // ========================================================

    // Furman-Pivi copper fallback; no complete lead-glass value is available
    // in Peng et al. or Li et al.
    G4double elasticEnergySigma = 1.9 * eV;

    // ========================================================
    // Rediffused yield component
    // ========================================================

    // Asymptote : ajustement sur Peng et al. Fig. 4 (peng-sey-fit), contrainte
    // par les points au-dessus de 100 eV.
    G4double rediffusedHighEnergyYield = 0.193031961294623;
    // Montee a basse energie : l'ajustement restituait E_r = 0.119 eV et r = 0.145,
    // valeurs cuivre de Furman-Pivi, qui saturent des 1 eV. Un electron rediffuse
    // doit penetrer le materiau, y perdre de l'energie et ressortir : le rendement
    // doit tendre vers 0 a basse energie. Forme physique reprise de Li et al.,
    // Photonics 9 (2022) 978. Verification : delta_r(1 eV) passe de 0.144 a 0.005,
    // et la fraction elastique a 10 eV de 72 % a 91 % (Scholtz, via Wu et al.,
    // Rev. Sci. Instrum. 79 (2008) 073104 : ~95 % a 10 eV, ~5 % a 100 eV).
    G4double rediffusedYieldEnergyScale = 40.0 * eV;
    G4double rediffusedYieldShape = 1.0;

    G4double rediffusedAngular1 = 0.26;
    G4double rediffusedAngular2 = 2.0;

    // ========================================================
    // Rediffused outgoing-energy distribution
    // ========================================================

    // Furman-Pivi copper fallback.
    G4double rediffusedEnergyExponent = 1.0;

    // ========================================================
    // Total yield curve
    // ========================================================

    // Peng et al. (2024) Eq. (8) : ajustement du rendement TOTAL mesure sur
    // 100-1500 eV, par la forme fonctionnelle de Furman Eq. (2). Ces
    // parametres decrivent delta_total, pas la seule composante vraie : celle-ci
    // s'en deduit par soustraction, Peng Eq. (7), dans ComputeYields.

    // Amplitude : 2.40 est la valeur mesuree par Peng et al. sur du verre au plomb
    // NU. La paroi active d'un MCP est un verre au plomb REDUIT sous hydrogene, dont
    // le rendement maximal est situe entre 2.5 et 4 par la litterature (Fijol et al.
    // 2.5 ; Then & Pantano ~3.5 ; Wu et al. 3-4 ; Li et al. 3.46). Valeur centrale
    // retenue ici ; c'est le parametre a caler sur une mesure de gain.
    // La forme (E_pic et s ci-dessous) reste celle de l'ajustement de Peng.
    G4double totalPeakYield = 3.00;

    G4double totalPeakIncidentEnergy =
        243.71 * eV;

    G4double totalYieldShape = 1.39;

    // Li et al. (2022), because Peng et al. do not publish t1, t2, t3 and t4.
    G4double totalAngular1 = 0.66;
    G4double totalAngular2 = 0.80;
    G4double totalAngular3 = 0.70;
    G4double totalAngular4 = 1.00;

    // ========================================================
    // True-secondary outgoing-energy distributions
    // ========================================================

    // Furman-Pivi copper fallback; Peng et al. and Li et al. do not provide a
    // complete outgoing-energy spectrum parameterization for lead glass.

    static constexpr G4int numberOfMultiplicityBins = 10;

    // p_n for n = 1, ..., 10.
    std::array<G4double, numberOfMultiplicityBins>
        trueSecondaryEnergyShape {
            1.6,
            2.0,
            1.8,
            4.7,
            1.8,
            2.4,
            1.8,
            1.8,
            2.3,
            1.8
        };

    // epsilon_n for n = 1, ..., 10.
    std::array<G4double, numberOfMultiplicityBins>
        trueSecondaryEnergyScale {
             3.90 * eV,
             6.20 * eV,
            13.00 * eV,
             8.80 * eV,
             6.25 * eV,
             2.25 * eV,
             9.20 * eV,
             5.30 * eV,
            17.80 * eV,
            10.00 * eV
        };

    // ========================================================
    // Multiplicity
    // ========================================================

    // The 8-degree channel produces grazing impacts. With the lead-glass
    // angular correction their conditional mean can exceed the old Cu
    // validation limit M=10, so the binomial support must be wider.
    G4int maximumMultiplicity = 50;

    // ========================================================
    // Emitted-direction distribution
    // ========================================================


    G4double emissionAngularExponent = 1.0;

    // ========================================================
    // Rugosite de paroi
    // ========================================================

    // Ecart-type, en radians, de l'inclinaison de la normale locale par rapport
    // a la normale moyenne : modele de microfacettes. La reflexion elastique est
    // speculaire autour de cette normale perturbee.
    //
    // 0 redonne la reflexion speculaire exacte (Wu et al. : "the axial and
    // angular components of the electrons' velocity are left unchanged, while
    // the radial component is reversed"). Une valeur croissante ramene
    // progressivement vers le tirage diffus de Furman-Pivi.
    //
    // Parametre libre, a caler sur une courbe gain-tension. C'est le seul levier
    // connu qui corrige ENSEMBLE le niveau du gain et la pente : allonger les
    // sauts monte l'energie d'impact, donc le gain, et k vaut environ 1/E.
    // Les poignees de rendement (s, delta_max) corrigent l'un en degradant
    // l'autre.
    G4double reflectionRoughness = 0.0;

    // ========================================================
    // Numerical and geometrical parameters
    // ========================================================

    G4double surfaceOffset = 1.0 * nm;


    // Fast path for the direct rejection sampler. Difficult low-energy,
    // high-multiplicity cases then switch to the exact truncated-Gamma plus
    // Dirichlet sampler instead of wasting up to maximumSamplingAttempts
    // complete groups of Gamma draws.
    G4int directSamplingAttempts = 32;

    // Safety limit for rejection samplers, including the adaptive exact
    // sampler used after the direct fast path.
    G4int maximumSamplingAttempts = 100000;

    G4double energyConservationRelativeTolerance =
        1.0e-10;
};

#endif
