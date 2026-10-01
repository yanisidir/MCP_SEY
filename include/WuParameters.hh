#ifndef WuParameters_h
#define WuParameters_h 1

#include "globals.hh"
#include "G4SystemOfUnits.hh"

/// Parametres du modele de Wu et al., Rev. Sci. Instrum. 79 (2008) 073104.
/// Toutes les valeurs proviennent de sa Table 1 : il s'agit d'une tentative de
/// REPRODUCTION, pas d'un ajustement. Aucune de ces valeurs n'est libre.
struct WuParameters
{
    // Table 1. Rendement maximal a incidence normale et energie du pic.
    G4double peakYield = 4.0;
    G4double peakEnergy = 260.0 * eV;

    // Eq. (3), parametre de forme de la fonction universelle. Wu : "typically
    // around 1.3 (the value used in this work)".
    G4double yieldShape = 1.3;

    // Eq. (2), variation du rendement maximal avec l'angle d'incidence :
    // delta_m(theta) = delta_m(0) exp[alpha (1 - cos theta)].
    G4double angularYield = 0.6;

    // Energie la plus probable des secondaires, Eq. (4). Table 1 : 3 eV.
    G4double emissionEnergy = 3.0 * eV;

    // Facteur d'echelle sur la fraction reflechie de Scholtz, Eq. (6).
    // Table 1 : R = 0.85.
    G4double reflectionScale = 0.85;

    // Conservation de l'energie : Wu retire l'ensemble des energies tant que
    // leur somme depasse l'energie d'impact. Plafond pour eviter une boucle
    // infinie quand l'energie d'impact est trop faible pour la multiplicite.
    G4int maximumSamplingAttempts = 100;

    // Loi cosinus d'emission, Eq. (5) : exposant 1 = loi de Lambert.
    G4double emissionAngularExponent = 1.0;
};

#endif
