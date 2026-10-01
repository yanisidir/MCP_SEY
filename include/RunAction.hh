#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class DetectorConstruction;
class GaussianSource;
class G4Run;

/// Ouvre le fichier ROOT du run, definit les trois tables d'analyse et
/// ecrit la ligne de configuration en fin de run.
class RunAction : public G4UserRunAction
{
  public:
    RunAction(const DetectorConstruction* detector, const GaussianSource* source);
    ~RunAction() override = default;

    void BeginOfRunAction(const G4Run*) override;
    void EndOfRunAction(const G4Run*) override;

  private:
    const DetectorConstruction* fDetector = nullptr;
    const GaussianSource* fSource = nullptr;
};

#endif
