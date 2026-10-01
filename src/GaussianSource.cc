#include "GaussianSource.hh"
#include "G4GenericMessenger.hh"

GaussianSource::GaussianSource()
{
    fMessenger = new G4GenericMessenger(this, "/source/gauss/", "Source gaussienne d'electrons");
    fMessenger->DeclarePropertyWithUnit("muE", "keV", fMuE, "Energie cinetique moyenne")
        .SetParameterName("e", false).SetRange("e>0").SetStates(G4State_PreInit, G4State_Idle);
    fMessenger->DeclarePropertyWithUnit("sigmaE", "keV", fSigmaE, "Ecart-type de l'energie ; 0 = fixe")
        .SetParameterName("s", false).SetRange("s>=0").SetStates(G4State_PreInit, G4State_Idle);
    fMessenger->DeclarePropertyWithUnit("muZ", "mm", fMuZ, "Hauteur d'injection moyenne, sur l'axe du pore")
        .SetParameterName("z", false).SetStates(G4State_PreInit, G4State_Idle);
    fMessenger->DeclarePropertyWithUnit("sigmaZ", "mm", fSigmaZ, "Ecart-type de la hauteur ; 0 = fixe")
        .SetParameterName("s", false).SetRange("s>=0").SetStates(G4State_PreInit, G4State_Idle);
    fMessenger->DeclarePropertyWithUnit("zOrigin", "mm", fZOrigin,
        "z de la face d'entree dans le repere des valeurs de muZ ; 0 = repere de la plaque")
        .SetParameterName("z", false).SetStates(G4State_PreInit, G4State_Idle);
    // Impulsion en unites d'energie (c = 1) : 100 keV signifie 100 keV/c.
    fMessenger->DeclarePropertyWithUnit("muPz", "keV", fMuPz,
        "Composante pz moyenne de l'impulsion ; au-dela de p, tir axial +z")
        .SetParameterName("p", false).SetStates(G4State_PreInit, G4State_Idle);
    fMessenger->DeclarePropertyWithUnit("sigmaPz", "keV", fSigmaPz, "Ecart-type de pz ; 0 = fixe")
        .SetParameterName("s", false).SetRange("s>=0").SetStates(G4State_PreInit, G4State_Idle);
}
GaussianSource::~GaussianSource() { delete fMessenger; }
