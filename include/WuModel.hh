#ifndef WuModel_h
#define WuModel_h 1

#include "FurmanPiviTypes.hh"
#include "WuParameters.hh"

#include "globals.hh"
#include "G4ThreeVector.hh"

#include <vector>

/// Modele d'emission secondaire de Wu et al., Rev. Sci. Instrum. 79 (2008)
/// 073104. Il differe de Furman-Pivi sur cinq points :
///
///   1. deux composantes seulement, reflexion elastique et vrais secondaires,
///      sans composante rediffusee ;
///   2. Eq. (3) donne le rendement VRAI, la reflexion s'y ajoute par Eq. (7),
///      alors que Peng Eq. (7) traite Eq. (2) comme le rendement TOTAL et
///      soustrait les composantes de reflexion ;
///   3. multiplicite de Poisson, non binomiale ;
///   4. energies d'emission en Maxwell-Boltzmann Eq. (4), retirees en bloc
///      tant que leur somme depasse l'energie d'impact ;
///   5. loi angulaire Vm(theta) = Vm(0)/sqrt(cos theta).
///
/// L'interface est celle de FurmanPiviModel pour que le processus puisse
/// passer de l'un a l'autre sans autre changement.
class WuModel
{
  public:
    explicit WuModel(const WuParameters& parameters);

    FurmanEmissionResult GenerateEmission(
        G4double incidentEnergy,
        G4double incidentAngle,
        const G4ThreeVector& incidentDirection,
        const G4ThreeVector& materialToVacuumNormal) const;

    /// Fraction reflechie de Scholtz et al., Philips J. Res. 50 (1996) 375,
    /// ajustee par Wu Eq. (6), multipliee par le facteur d'echelle R.
    G4double ReflectedFraction(G4double incidentEnergy) const;

    /// Eq. (3) : rendement des VRAIS secondaires.
    G4double TrueSecondaryYield(G4double incidentEnergy, G4double incidentAngle) const;

    /// Eq. (7) : rendement total, reflechis compris.
    G4double TotalYield(G4double incidentEnergy, G4double incidentAngle) const;

  private:
    /// Eq. (4) : Maxwell-Boltzmann d'energie la plus probable E0, soit une
    /// loi Gamma de forme 2 et d'echelle E0.
    G4double SampleEmissionEnergy() const;

    /// Eq. (5) : loi cosinus autour de la normale sortante.
    G4ThreeVector SampleCosineDirection(const G4ThreeVector& materialToVacuumNormal) const;

    /// Reflexion speculaire : "the axial and angular components of the
    /// electrons' velocity are left unchanged, while the radial component is
    /// reversed".
    G4ThreeVector SpecularDirection(
        const G4ThreeVector& incidentDirection,
        const G4ThreeVector& materialToVacuumNormal) const;

    WuParameters fParameters;
};

#endif
