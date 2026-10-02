#include "RunAction.hh"

#include "DetectorConstruction.hh"
#include "GaussianSource.hh"

#include "G4AnalysisManager.hh"
#include "G4Run.hh"
#include "G4SystemOfUnits.hh"
#include "G4Threading.hh"

#include <filesystem>
#include <system_error>
#include <unistd.h>

RunAction::RunAction(const DetectorConstruction* detector, const GaussianSource* source)
  : fDetector(detector), fSource(source)
{
    auto* analysisManager = G4AnalysisManager::Instance();
    analysisManager->SetDefaultFileType("root");
    analysisManager->SetNtupleMerging(true);

    // Table 0 : une ligne par evenement.
    analysisManager->CreateNtuple("events", "Un electron injecte par evenement");
    analysisManager->CreateNtupleIColumn("RunID");
    analysisManager->CreateNtupleIColumn("EventID");
    analysisManager->CreateNtupleIColumn("NDetectedPrimary");
    analysisManager->CreateNtupleIColumn("NDetectedSecondary");
    analysisManager->CreateNtupleIColumn("Gain");
    analysisManager->CreateNtupleIColumn("IsMultiplicationLimited");
    analysisManager->CreateNtupleDColumn("NTracks");
    analysisManager->CreateNtupleDColumn("IncidentEnergy_eV");
    analysisManager->CreateNtupleDColumn("FirstTime_ns");
    analysisManager->CreateNtupleDColumn("MeanTime_ns");
    analysisManager->CreateNtupleDColumn("LastTime_ns");
    analysisManager->CreateNtupleDColumn("TimeSpread_ns");
    analysisManager->CreateNtupleIColumn("MaxGeneration");
    analysisManager->CreateNtupleDColumn("NCreatedTracks");
    analysisManager->CreateNtupleIColumn("NAt50um");
    analysisManager->CreateNtupleDColumn("NImpacts");
    analysisManager->CreateNtupleDColumn("MeanImpactEnergy_eV");
    analysisManager->CreateNtupleDColumn("MeanImpactAngle_deg");
    analysisManager->FinishNtuple();

    // Table 1 : une ligne par electron atteignant le plan de mesure.
    analysisManager->CreateNtuple("downstream", "Position, temps et energie au plan de mesure");
    analysisManager->CreateNtupleIColumn("RunID");
    analysisManager->CreateNtupleIColumn("EventID");
    analysisManager->CreateNtupleIColumn("TrackID");
    analysisManager->CreateNtupleIColumn("ParentID");
    analysisManager->CreateNtupleDColumn("Time_ns");
    analysisManager->CreateNtupleDColumn("KineticEnergy_eV");
    analysisManager->CreateNtupleDColumn("X_mm");
    analysisManager->CreateNtupleDColumn("Y_mm");
    analysisManager->CreateNtupleDColumn("Z_mm");
    analysisManager->CreateNtupleIColumn("Generation");
    analysisManager->CreateNtupleDColumn("DirectionX");
    analysisManager->CreateNtupleDColumn("DirectionY");
    analysisManager->CreateNtupleDColumn("DirectionZ");
    analysisManager->FinishNtuple();

    // Table 2 : une seule ligne, la configuration du run.
    analysisManager->CreateNtuple("configuration", "Geometrie et limites du run");
    analysisManager->CreateNtupleIColumn("RunID");
    analysisManager->CreateNtupleDColumn("Thickness_mm");
    analysisManager->CreateNtupleDColumn("PoreDiameter_um");
    analysisManager->CreateNtupleDColumn("Angle_deg");
    analysisManager->CreateNtupleDColumn("Voltage_V");
    analysisManager->CreateNtupleIColumn("MaxTracksPerEvent");
    analysisManager->CreateNtupleIColumn("MaxGenerations");
    analysisManager->CreateNtupleIColumn("TrackLimitMode"); // 1 : plafond de creation, transport poursuivi.
    analysisManager->CreateNtupleDColumn("DriftDistance_um");
    analysisManager->CreateNtupleIColumn("FieldOnlyInPore");
    analysisManager->CreateNtupleIColumn("ScoringAtPlateFace");
    // 1 = Wu, 0 = Furman-Pivi / Peng. Absente des fichiers anterieurs a
    // octobre 2026, produits avec Furman-Pivi / Peng.
    analysisManager->CreateNtupleIColumn("EmissionModel");
    analysisManager->CreateNtupleDColumn("SourceMuE_keV");
    analysisManager->CreateNtupleDColumn("SourceSigmaE_keV");
    analysisManager->CreateNtupleDColumn("SourceMuZ_mm");
    analysisManager->CreateNtupleDColumn("SourceSigmaZ_mm");
    analysisManager->CreateNtupleDColumn("SourceMuPz_keV");
    analysisManager->CreateNtupleDColumn("SourceSigmaPz_keV");
    analysisManager->FinishNtuple();

    // Spectres d'impact, cumules sur le run et fusionnes entre threads.
    analysisManager->CreateH1("impactEnergy",
        "Energie d'impact sur la paroi (eV)", 400, 0., 400.);
    analysisManager->CreateH1("impactAngle",
        "Angle d'incidence sur la paroi (deg)", 90, 0., 90.);
}

namespace {
// G4AnalysisManager::OpenFile renvoie true meme quand l'ouverture echoue ; le
// programme mourait alors bien plus loin, par une faute de segmentation. Le
// chemin est donc verifie avant. Appelee par chaque thread : create_directories
// est idempotente et sa variante a code d'erreur tolere la course.
void PreparerDossierDeSortie(const G4String& sortie)
{
    const std::filesystem::path chemin(sortie.c_str());
    if (chemin.filename().empty()) {
        G4ExceptionDescription detail;
        detail << "Le chemin de sortie \"" << sortie
               << "\" ne designe pas un fichier.";
        G4Exception("RunAction::BeginOfRunAction", "Output", FatalException, detail);
    }

    auto dossier = chemin.parent_path();
    if (dossier.empty()) {dossier = ".";}

    std::error_code erreur;
    if (!std::filesystem::is_directory(dossier, erreur)) {
        std::filesystem::create_directories(dossier, erreur);
        if (!std::filesystem::is_directory(dossier)) {
            G4ExceptionDescription detail;
            detail << "Impossible de creer le dossier de sortie \"" << dossier.string()
                   << "\" : " << erreur.message() << ".";
            G4Exception("RunAction::BeginOfRunAction", "Output", FatalException, detail);
        }
        if (G4Threading::G4GetThreadId() <= 0) {
            G4cout << "Dossier de sortie cree : " << dossier.string() << G4endl;
        }
    }

    if (::access(dossier.c_str(), W_OK) != 0) {
        G4ExceptionDescription detail;
        detail << "Le dossier de sortie \"" << dossier.string()
               << "\" n'est pas accessible en ecriture.";
        G4Exception("RunAction::BeginOfRunAction", "Output", FatalException, detail);
    }
}
}

void RunAction::BeginOfRunAction(const G4Run*)
{
    if (G4Threading::G4GetThreadId() == 0 || !G4Threading::IsMultithreadedApplication()) {
        G4cout << "Modele d'emission secondaire : "
               << (fDetector->UseWuModel() ? "Wu et al. 2008"
                                           : "Furman-Pivi / Peng")
               << G4endl;
    }

    PreparerDossierDeSortie(fDetector->Output());

    if (!G4AnalysisManager::Instance()->OpenFile(fDetector->Output())) {
        G4Exception("RunAction::BeginOfRunAction", "Output", FatalException,
                    "Impossible d'ouvrir le fichier ROOT.");
    }
}

void RunAction::EndOfRunAction(const G4Run* run)
{
    auto* analysisManager = G4AnalysisManager::Instance();

    // Une seule ligne, ecrite par le worker 0 (ou en serie).
    if (G4Threading::G4GetThreadId() == 0 || !G4Threading::IsMultithreadedApplication()) {
        constexpr G4int table = 2;
        analysisManager->FillNtupleIColumn(table, 0, run->GetRunID());
        analysisManager->FillNtupleDColumn(table, 1, fDetector->Thickness() / mm);
        analysisManager->FillNtupleDColumn(table, 2, fDetector->Diameter() / um);
        analysisManager->FillNtupleDColumn(table, 3, fDetector->Angle() / degree);
        analysisManager->FillNtupleDColumn(table, 4, fDetector->Voltage() / volt);
        analysisManager->FillNtupleIColumn(table, 5, fDetector->MaxTracks());
        analysisManager->FillNtupleIColumn(table, 6, fDetector->MaxGenerations());
        analysisManager->FillNtupleIColumn(table, 7, 1);
        analysisManager->FillNtupleDColumn(table, 8, fDetector->DriftDistance() / um);
        analysisManager->FillNtupleIColumn(table, 9, 1);
        analysisManager->FillNtupleIColumn(table, 10, 1);
        analysisManager->FillNtupleIColumn(table, 11, fDetector->UseWuModel() ? 1 : 0);
        analysisManager->FillNtupleDColumn(table, 12, fSource->MuE() / keV);
        analysisManager->FillNtupleDColumn(table, 13, fSource->SigmaE() / keV);
        analysisManager->FillNtupleDColumn(table, 14, fSource->MuZ() / mm);
        analysisManager->FillNtupleDColumn(table, 15, fSource->SigmaZ() / mm);
        analysisManager->FillNtupleDColumn(table, 16, fSource->MuPz() / keV);
        analysisManager->FillNtupleDColumn(table, 17, fSource->SigmaPz() / keV);
        analysisManager->AddNtupleRow(table);
    }

    analysisManager->Write();
    analysisManager->CloseFile();
}
