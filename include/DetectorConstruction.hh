#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
class G4GenericMessenger;
class G4LogicalVolume;

class DetectorConstruction : public G4VUserDetectorConstruction {
public:
    DetectorConstruction();
    ~DetectorConstruction() override;
    G4VPhysicalVolume* Construct() override;
    void ConstructSDandField() override;
    void SetPlateDiameter(G4double v) { fPlateDiameter = v; }
    void SetThickness(G4double v);
    void SetChannelDiameter(G4double v) { fChannelDiameter = v; }
    void SetAngle(G4double v) { fAngle = v; }
    void SetVoltage(G4double v);
    G4ThreeVector Field() const { return {0, 0, -fVoltage / fThickness}; }
    G4double DriftDistance() const { return 50 * um; } // Plan fixe apres la plaque.
    G4double Thickness() const { return fThickness; }
    G4double Diameter() const { return fChannelDiameter; }
    G4double Angle() const { return fAngle; }
    G4double Voltage() const { return fVoltage; }
    G4int MaxTracks() const { return fMaxTracks; }
    G4int MaxGenerations() const { return fMaxGenerations; }
    G4bool SaveExits() const { return fSaveExits; }

    /// Modele d'emission : true = Wu et al. 2008 (defaut), false =
    /// Furman-Pivi avec les parametres de Peng. Le modele retenu par defaut est
    /// celui de Wu parce qu'il reproduit les trois configurations publiees a un
    /// facteur ~2, la ou l'autre s'en ecarte d'un facteur 4 a 41.
    G4bool UseWuModel() const { return fUseWuModel; }
    const G4String& Output() const { return fOutput; }
private:
    G4GenericMessenger* fMessenger = nullptr;
    G4LogicalVolume* fPoreLogic = nullptr;
    G4double fPlateDiameter = 0.5 * mm;
    G4double fThickness = 1.0 * mm;
    G4double fChannelDiameter = 15.0 * um;
    G4double fAngle = 8.0 * degree;
    G4double fVoltage = 2000.0 * volt;
    G4int fMaxTracks = 100000;
    G4int fMaxGenerations = 0;
    G4bool fSaveExits = false;
    G4bool fUseWuModel = true;
    G4String fOutput = "mcp.root";
};

#endif
