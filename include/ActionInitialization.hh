#ifndef ActionInitialization_h
#define ActionInitialization_h 1

#include "G4VUserActionInitialization.hh"

class DetectorConstruction;
class GaussianSource;

/// Construit les User Actions. Le master ne recoit que le RunAction, qui
/// ferme le fichier fusionne ; chaque worker recoit sa propre chaine
/// complete, reliee au meme EventAction.
class ActionInitialization : public G4VUserActionInitialization
{
  public:
    ActionInitialization(const DetectorConstruction* detector, const GaussianSource* source);
    ~ActionInitialization() override = default;

    void BuildForMaster() const override;
    void Build() const override;

  private:
    const DetectorConstruction* fDetector = nullptr;
    const GaussianSource* fSource = nullptr;
};

#endif
