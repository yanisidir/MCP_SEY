#ifndef WuParameters_h
#define WuParameters_h 1

#include "globals.hh"
#include "G4SystemOfUnits.hh"

// Wu et al., Rev. Sci. Instrum. 79 (2008) 073104, Table 1. Aucune valeur libre.
struct WuParameters
{
    G4double peakYield = 4.0;                  // delta_m(0)
    G4double peakEnergy = 260.0 * eV;          // V_m(0)
    G4double yieldShape = 1.3;                 // s, Eq. (3)
    G4double angularYield = 0.6;               // alpha, Eq. (2)
    G4double emissionEnergy = 3.0 * eV;        // E_0, Eq. (4)
    G4double reflectionScale = 0.85;           // R, facteur d'echelle sur Eq. (6)
    G4double emissionAngularExponent = 1.0;    // loi cosinus, Eq. (5)

    // Borne le rejet sur la conservation de l'energie, qui ne converge pas si
    // la multiplicite tiree est trop grande pour l'energie d'impact.
    G4int maximumSamplingAttempts = 100;
};

#endif
