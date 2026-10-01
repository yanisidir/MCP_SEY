#ifndef GenerationInfo_h
#define GenerationInfo_h 1

#include "G4VUserTrackInformation.hh"
#include "G4Track.hh"

// Generation d'emission : 0 pour le primaire, +1 a chaque vrai secondaire.
// Une reflexion garde la generation incidente.
class GenerationInfo : public G4VUserTrackInformation {
public:
    explicit GenerationInfo(G4int generation) : fGeneration(generation) {}
    G4bool exitRecorded = false; // Premier passage a la face aval du MCP.
    G4int Generation() const { return fGeneration; }
    static G4int Of(const G4Track& track)
    {
        const auto* info = dynamic_cast<const GenerationInfo*>(track.GetUserInformation());
        return info ? info->Generation() : 0;  // Sans marque : electron primaire.
    }
private:
    G4int fGeneration;
};

#endif
