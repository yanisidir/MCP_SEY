#include "EventAction.hh"

#include "DetectorConstruction.hh"
#include "GenerationInfo.hh"

#include "G4AnalysisManager.hh"
#include "G4Event.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"

#include <algorithm>
#include <cmath>
#include <limits>

EventAction::EventAction(const DetectorConstruction* detector) : fDetector(detector) {}

void EventAction::BeginOfEventAction(const G4Event* event)
{
    const auto* vertex = event->GetPrimaryVertex();
    if (event->GetNumberOfPrimaryVertex() != 1 || !vertex || vertex->GetNumberOfParticle() != 1
        || vertex->GetPrimary()->GetPDGcode() != 11)
    {
        G4Exception("EventAction::BeginOfEventAction", "Source", FatalException,
                    "Injecter exactement un electron par evenement.");
    }

    fEventID = event->GetEventID();
    fT0 = vertex->GetT0();
    fEnergy = vertex->GetPrimary()->GetKineticEnergy();
    fPrimary = fSecondary = 0;
    fMaxGeneration = 0;
    fTracks = 0;
    fAt50um = 0;
    fCreatedTracks = 1;  // Le primaire est deja cree.
    fImpacts = 0;
    fImpactEnergy = fImpactAngle = 0.;
    fMultiplicationLimited = fDetector->MaxTracks() == 1;
    fFirst = std::numeric_limits<G4double>::infinity();
    fLast = fMean = fM2 = 0.;
}

void EventAction::StartTrack(const G4Track* track)
{
    ++fTracks;
    if (!track->GetUserInformation()) {
        const_cast<G4Track*>(track)->SetUserInformation(new GenerationInfo(0));
    }
    fMaxGeneration = std::max(fMaxGeneration, GenerationInfo::Of(*track));
}

G4int EventAction::ReserveSecondaries(G4int requested)
{
    // Compter avant la mise en attente : les descendants non encore suivis
    // occupent deja leur place. Etat propre a cet evenement et a son worker.
    const auto limit = fDetector->MaxTracks();
    const auto accepted =
        limit > 0
            ? static_cast<G4int>(std::min<G4long>(requested, std::max<G4long>(0, limit - fCreatedTracks)))
            : requested;
    fCreatedTracks += accepted;
    if (limit > 0 && fCreatedTracks >= limit) fMultiplicationLimited = true;
    return accepted;
}

void EventAction::RecordExit(const G4Track* track, const G4StepPoint* point)
{
    if (track->GetParentID() == 0) ++fPrimary;
    else ++fSecondary;

    const G4double t = point->GetGlobalTime() - fT0;
    fFirst = std::min(fFirst, t);
    fLast = std::max(fLast, t);
    // Moyenne et variance en ligne : pas de tableau geant par avalanche.
    const auto delta = t - fMean;
    fMean += delta / (fPrimary + fSecondary);
    fM2 += delta * (t - fMean);
}

void EventAction::RecordDownstream(const G4Track* track, const G4StepPoint* point)
{
    ++fAt50um;
    // Table 1 : impacts individuels uniquement a +50 um.
    constexpr G4int table = 1;
    const auto t = point->GetGlobalTime() - fT0;
    if (!fDetector->SaveExits()) return;

    auto* analysisManager = G4AnalysisManager::Instance();
    analysisManager->FillNtupleIColumn(table, 0, G4RunManager::GetRunManager()->GetCurrentRun()->GetRunID());
    analysisManager->FillNtupleIColumn(table, 1, fEventID);
    analysisManager->FillNtupleIColumn(table, 2, track->GetTrackID());
    analysisManager->FillNtupleIColumn(table, 3, track->GetParentID());
    analysisManager->FillNtupleDColumn(table, 4, t / ns);
    analysisManager->FillNtupleDColumn(table, 5, point->GetKineticEnergy() / eV);
    for (G4int i = 0; i < 3; ++i) {
        analysisManager->FillNtupleDColumn(table, 6 + i, point->GetPosition()[i] / mm);
    }
    analysisManager->FillNtupleIColumn(table, 9, GenerationInfo::Of(*track));
    for (G4int i = 0; i < 3; ++i) {
        analysisManager->FillNtupleDColumn(table, 10 + i, point->GetMomentumDirection()[i]);
    }
    analysisManager->AddNtupleRow(table);
}

void EventAction::RecordImpact(G4double energy, G4double angle)
{
    ++fImpacts;
    fImpactEnergy += energy;
    fImpactAngle += angle;

    // Spectres cumules sur tout le run. Les histogrammes sont fusionnes entre
    // threads par G4AnalysisManager, contrairement aux sommes ci-dessus qui
    // restent locales a l'evenement.
    auto* analysisManager = G4AnalysisManager::Instance();
    analysisManager->FillH1(0, energy / eV);
    analysisManager->FillH1(1, angle / degree);
}

void EventAction::EndOfEventAction(const G4Event*)
{
    auto* analysisManager = G4AnalysisManager::Instance();
    constexpr G4int table = 0;
    const G4int gain = fPrimary + fSecondary;
    const G4double nan = std::numeric_limits<G4double>::quiet_NaN();

    analysisManager->FillNtupleIColumn(table, 0, G4RunManager::GetRunManager()->GetCurrentRun()->GetRunID());
    analysisManager->FillNtupleIColumn(table, 1, fEventID);
    analysisManager->FillNtupleIColumn(table, 2, fPrimary);
    analysisManager->FillNtupleIColumn(table, 3, fSecondary);
    analysisManager->FillNtupleIColumn(table, 4, gain);
    analysisManager->FillNtupleIColumn(table, 5, fMultiplicationLimited);
    analysisManager->FillNtupleDColumn(table, 6, fTracks);
    analysisManager->FillNtupleDColumn(table, 7, fEnergy / eV);
    analysisManager->FillNtupleDColumn(table, 8, gain ? fFirst / ns : nan);
    analysisManager->FillNtupleDColumn(table, 9, gain ? fMean / ns : nan);
    analysisManager->FillNtupleDColumn(table, 10, gain ? fLast / ns : nan);
    analysisManager->FillNtupleDColumn(table, 11, gain ? std::sqrt(std::max(0., fM2 / gain)) / ns : nan);
    analysisManager->FillNtupleIColumn(table, 12, fMaxGeneration);
    analysisManager->FillNtupleDColumn(table, 13, fCreatedTracks);
    analysisManager->FillNtupleIColumn(table, 14, fAt50um);
    analysisManager->FillNtupleDColumn(table, 15, fImpacts);
    analysisManager->FillNtupleDColumn(
        table, 16, fImpacts ? fImpactEnergy / fImpacts / eV : nan);
    analysisManager->FillNtupleDColumn(
        table, 17, fImpacts ? fImpactAngle / fImpacts / degree : nan);
    analysisManager->AddNtupleRow(table);
}
