#ifndef FurmanPiviParameters_h
#define FurmanPiviParameters_h 1

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"

#include <array>

// Modele de Furman-Pivi, Phys. Rev. ST Accel. Beams 5 (2002) 124404, avec les
// rendements ajustes sur Peng et al., Nucl. Instrum. Methods A 1062 (2024)
// 169163, et les lois angulaires de Li et al., Photonics 9 (2022) 978.
struct FurmanPiviParameters
{
    // ========================================================
    // Composante elastique
    // ========================================================

    // Ajustement sur Peng et al. Fig. 4 (dossier peng-sey-fit). Les valeurs
    // coincident avec le cuivre de Furman-Pivi : Peng n'a pas mesure les
    // composantes de reflexion du verre au plomb separement.
    G4double elasticHighEnergyYield = 0.0197991598770481;
    G4double elasticPeakYield = 0.500809176899711;
    G4double elasticPeakEnergy = 0.0 * eV;
    G4double elasticYieldWidth = 60.7290930799462 * eV;
    G4double elasticYieldShape = 0.985795034334466;

    G4double elasticAngular1 = 0.26;
    G4double elasticAngular2 = 2.0;

    G4double elasticEnergySigma = 1.9 * eV;     // cuivre, faute de valeur publiee

    // ========================================================
    // Composante rediffusee
    // ========================================================

    G4double rediffusedHighEnergyYield = 0.193031961294623;   // Peng Fig. 4

    // L'ajustement restituait les valeurs cuivre (Er = 0.119 eV, r = 0.145),
    // qui saturent des 1 eV. Un electron rediffuse ayant penetre le materiau,
    // on attend un rendement decroissant vers 0 a basse energie ; la forme a ete
    // remplacee sur cet argument, non sur une mesure du verre au plomb.
    // Forme reprise de Li et al. Verification : delta_r(1 eV) passe de 0.144 a
    // 0.005, et la fraction elastique a 10 eV de 72 % a 91 % (Scholtz et al.,
    // Philips J. Res. 50 (1996) 375 : ~95 % a 10 eV, ~5 % a 100 eV).
    G4double rediffusedYieldEnergyScale = 40.0 * eV;
    G4double rediffusedYieldShape = 1.0;

    G4double rediffusedAngular1 = 0.26;
    G4double rediffusedAngular2 = 2.0;
    G4double rediffusedEnergyExponent = 1.0;    // cuivre

    // ========================================================
    // Rendement total
    // ========================================================

    // Peng Eq. (8) ajuste le rendement TOTAL mesure sur 100-1500 eV. La
    // composante vraie s'en deduit par soustraction, Peng Eq. (7), dans
    // ComputeYields.

    // Peng mesure 2.40 sur du verre au plomb NU. La paroi d'un MCP est du verre
    // REDUIT sous hydrogene, que la litterature situe entre 2.5 et 4 (Fijol 2.5 ;
    // Then & Pantano ~3.5 ; Wu 3-4 ; Li 3.46). Valeur centrale de cette plage,
    // non ajustee sur une mesure de gain : elle reste le parametre le plus
    // naturel a caler si une telle mesure est disponible.
    G4double totalPeakYield = 3.00;

    G4double totalPeakIncidentEnergy = 243.71 * eV;
    G4double totalYieldShape = 1.39;

    // Li et al. : Peng ne publie pas t1 a t4.
    G4double totalAngular1 = 0.66;
    G4double totalAngular2 = 0.80;
    G4double totalAngular3 = 0.70;
    G4double totalAngular4 = 1.00;

    // ========================================================
    // Energies des vrais secondaires
    // ========================================================

    // Cuivre de Furman-Pivi : ni Peng ni Li ne publient de spectre complet
    // pour le verre au plomb.

    static constexpr G4int numberOfMultiplicityBins = 10;

    std::array<G4double, numberOfMultiplicityBins>
        trueSecondaryEnergyShape {1.6, 2.0, 1.8, 4.7, 1.8, 2.4, 1.8, 1.8, 2.3, 1.8};

    std::array<G4double, numberOfMultiplicityBins>
        trueSecondaryEnergyScale {
             3.90 * eV,  6.20 * eV, 13.00 * eV,  8.80 * eV,  6.25 * eV,
             2.25 * eV,  9.20 * eV,  5.30 * eV, 17.80 * eV, 10.00 * eV};

    // Les impacts rasants du canal incline peuvent depasser la limite M = 10
    // des validations cuivre.
    G4int maximumMultiplicity = 50;

    G4double emissionAngularExponent = 1.0;     // loi cosinus

    // ========================================================
    // Rugosite de paroi
    // ========================================================

    // Modele de microfacettes : ecart-type, en radians, de l'inclinaison de la
    // normale locale. 0 redonne la reflexion speculaire exacte sans consommer
    // de tirage. Parametre libre, non calibre. Ignore par WuModel.
    G4double reflectionRoughness = 0.0;

    // ========================================================
    // Parametres numeriques
    // ========================================================

    G4double surfaceOffset = 1.0 * nm;

    // Chemin rapide du tirage par rejet ; au-dela, le tirage exact
    // (Gamma tronquee + Dirichlet) prend le relais.
    G4int directSamplingAttempts = 32;
    G4int maximumSamplingAttempts = 100000;

    G4double energyConservationRelativeTolerance = 1.0e-10;
};

#endif
