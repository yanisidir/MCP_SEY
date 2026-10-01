#ifndef FurmanPiviTypes_h
#define FurmanPiviTypes_h 1

#include "globals.hh"
#include "G4ThreeVector.hh"

#include <vector>

enum class SEEType
{
    Elastic,
    Rediffused,
    TrueSecondary
};

struct FurmanYields
{
    G4double elastic = 0.0;
    G4double rediffused = 0.0;
    G4double trueSecondary = 0.0;

    G4double Total() const
    {
        return elastic + rediffused + trueSecondary;
    }
};

struct EmittedElectron
{
    SEEType type;

    G4double kineticEnergy = 0.0;
    
    G4ThreeVector direction;
};

struct FurmanEmissionResult
{
    G4int multiplicity = 0;

    std::vector<EmittedElectron> electrons;

};

#endif
