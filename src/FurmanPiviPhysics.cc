#include "FurmanPiviPhysics.hh"
#include "FurmanPiviProcess.hh"
#include "G4Electron.hh"
#include "G4ProcessManager.hh"
void FurmanPiviPhysics::ConstructParticle() { G4Electron::Definition(); }
void FurmanPiviPhysics::ConstructProcess()
{
    G4Electron::Definition()->GetProcessManager()->AddDiscreteProcess(new FurmanPiviProcess(fDetector));
}
