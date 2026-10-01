#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "FurmanPiviPhysics.hh"
#include "GaussianSource.hh"

#include "FTFP_BERT.hh"
#include "G4EmParameters.hh"
#include "G4RunManagerFactory.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"

#include <iostream>
#include <string>

int main(int argc, char** argv)
{
    if (argc > 2 || (argc == 2 && std::string(argv[1]) == "--help")) {
        std::cout << "Usage: MCP_SEY [fichier.mac]\n"
                  << "Sans macro: visualisation interactive. Reglages dans les commandes /mcp/, /gps/, /run/.\n";
        return argc > 2 ? 1 : 0;
    }

    G4UIExecutive* ui = nullptr;
    if (argc == 1) ui = new G4UIExecutive(argc, argv);

    // La visualisation interactive reste en serie.
    auto* runManager = G4RunManagerFactory::CreateRunManager(
        ui ? G4RunManagerType::SerialOnly : G4RunManagerType::Default, 1);

    auto* detector = new DetectorConstruction;
    runManager->SetUserInitialization(detector);

    auto* physicsList = new FTFP_BERT;
    physicsList->SetVerboseLevel(0);
    // Ne pas arreter artificiellement les electrons lents dans le vide.
    G4EmParameters::Instance()->SetLowestElectronEnergy(0.0);
    physicsList->RegisterPhysics(new FurmanPiviPhysics(detector));
    runManager->SetUserInitialization(physicsList);

    // La source vit dans le master : les workers la lisent a chaque evenement.
    auto* source = new GaussianSource;
    runManager->SetUserInitialization(new ActionInitialization(detector, source));

    auto* visManager = new G4VisExecutive;
    visManager->Initialize();

    auto* uiManager = G4UImanager::GetUIpointer();
    int status = 0;
    if (!ui) {
        status = uiManager->ApplyCommand(G4String("/control/execute ") + argv[1]);
    }
    else {
        // vis.mac initialise le run et ouvre la vue.
        status = uiManager->ApplyCommand("/control/execute macros/vis.mac");
        if (!status) ui->SessionStart();
        delete ui;
    }

    delete visManager;
    delete runManager;
    return status ? 1 : 0;
}
