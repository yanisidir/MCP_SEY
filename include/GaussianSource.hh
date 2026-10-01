#ifndef GaussianSource_h
#define GaussianSource_h 1

#include "G4SystemOfUnits.hh"
#include "globals.hh"
class G4GenericMessenger;

// Parametres de la source, commandes /source/gauss/. Trois lois normales
// independantes : energie cinetique, hauteur z de l'injection et composante
// pz de l'impulsion. L'objet est cree dans le master, qui recoit les
// commandes, et lu par les workers a chaque evenement ; il ne change
// qu'entre deux runs, comme DetectorConstruction.
class GaussianSource {
public:
    GaussianSource();
    ~GaussianSource();
    G4double MuE() const { return fMuE; }
    G4double SigmaE() const { return fSigmaE; }
    G4double MuZ() const { return fMuZ - fZOrigin; }  // Dans le repere de la plaque.
    G4double SigmaZ() const { return fSigmaZ; }
    G4double ZOrigin() const { return fZOrigin; }
    G4double MuPz() const { return fMuPz; }
    G4double SigmaPz() const { return fSigmaPz; }
private:
    G4GenericMessenger* fMessenger = nullptr;
    G4double fMuE = 300 * eV;
    G4double fSigmaE = 0;
    G4double fMuZ = 0;
    G4double fSigmaZ = 0;
    G4double fZOrigin = 0;
    G4double fMuPz = 1 * MeV;  // Au-dela de p : tir axial +z.
    G4double fSigmaPz = 0;
};

#endif
