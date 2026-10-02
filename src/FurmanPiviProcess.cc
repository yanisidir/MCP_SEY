#include "FurmanPiviProcess.hh"
#include "DetectorConstruction.hh"
#include "GenerationInfo.hh"
#include "EventAction.hh"
#include "G4RunManager.hh"
#include "Randomize.hh"
#include "G4Electron.hh"
#include "G4DynamicParticle.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4VPhysicalVolume.hh"
#include "G4TransportationManager.hh"
#include "G4Navigator.hh"
#include "G4LogicalVolume.hh"
#include "G4GeometryTolerance.hh"
#include "G4VSolid.hh"
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <limits>

namespace {
// Un electron emis exactement sur la paroi serait aussitot compte dans le verre.
// On le decale vers l'interieur du pore, en verifiant le resultat plutot qu'en
// le supposant : le pore est coupe aux faces, et pres d'un bord un decalage vers
// l'axe peut sortir par z=0 ou z=L.
//
// Deux replis successifs si le point obtenu n'est pas interieur : reduire le
// decalage par dichotomie, puis reculer vers le point precedent. L'exception
// finale est volontaire -- mieux vaut s'arreter que d'emettre hors du volume.
G4ThreeVector EmissionPosition(const G4Step& step, const G4ThreeVector& normal, G4double offset)
{
    const auto* pre = step.GetPreStepPoint();
    const auto impact = step.GetPostStepPoint()->GetPosition();
    const auto& transform = pre->GetTouchableHandle()->GetHistory()->GetTopTransform();
    const auto* solid = pre->GetPhysicalVolume()->GetLogicalVolume()->GetSolid();
    auto inside = [&](const G4ThreeVector& p) {
        return solid->Inside(transform.TransformPoint(p)) == kInside;
    };
    const auto tolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
    for (auto distance=offset; distance>=4*tolerance; distance*=0.5) {
        const auto position = impact + distance*normal;
        if (inside(position)) return position;
    }
    // Au bord d'une face, reduire le decalage ne suffit pas toujours : reculer
    // aussi vers le point precedent.
    const auto displacement = pre->GetPosition()-impact;
    const auto back = displacement.unit();
    for (int factor : {1,2,4,8}) {
        const auto distance = std::min(factor*offset, 0.5*displacement.mag());
        const auto position = impact + distance*back;
        if (inside(position)) return position;
    }
    G4ExceptionDescription detail;
    detail << std::setprecision(17) << "Aucun point interieur apres repositionnement. Impact (mm): "
           << impact/mm << "; pre (mm): " << pre->GetPosition()/mm << "; normale: " << normal;
    G4Exception("EmissionPosition", "EmissionOutsidePore", FatalException, detail);
    return impact; // L'exception fatale termine normalement l'execution.
}
}

// Une reflexion doit ANNULER un franchissement que Geant4 a deja opere : quand
// ce code s'execute, Transportation a place la trajectoire dans le verre. On
// surcharge donc la mise a jour du pas pour la ramener dans le pore.
//
// Les champs du point final sont alors incoherents (ils decrivent le verre) et
// doivent tous etre reecrits a la main. En oublier un ne provoque aucune erreur
// visible : la trajectoire continue avec le mauvais materiau ou le mauvais
// volume, et le resultat est silencieusement faux.
G4Step* ReflectionParticleChange::UpdateStepForPostStep(G4Step* step)
{
    G4ParticleChange::UpdateStepForPostStep(step);
    if (!returnToPore) return step;

    auto* post = step->GetPostStepPoint();
    const auto* pre = step->GetPreStepPoint();
    auto* navigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
    const auto direction = post->GetMomentumDirection();
    // Recherche ABSOLUE (et non relative au volume courant) : c'est elle qui
    // annule le franchissement et retrouve le pore au point decale.
    const auto* volume = navigator->LocateGlobalPointAndSetup(
        post->GetPosition(), &direction, false, false);
    if (volume != pre->GetPhysicalVolume()) {
        G4ExceptionDescription detail;
        detail << std::setprecision(17) << "Point (mm): " << post->GetPosition()/mm
               << "; volume trouve: " << (volume ? volume->GetName() : "hors monde")
               << "; volume attendu: " << pre->GetPhysicalVolume()->GetName();
        G4Exception("ReflectionParticleChange", "ReflectionOutsidePore", FatalException, detail);
    }
    post->SetTouchableHandle(navigator->CreateTouchableHistoryHandle());
    post->SetMaterial(pre->GetMaterial());
    post->SetMaterialCutsCouple(pre->GetMaterialCutsCouple());
    post->SetSensitiveDetector(pre->GetSensitiveDetector());
    post->SetSafety(0.0);
    post->SetStepStatus(fPostStepDoItProc);
    return step;
}

FurmanPiviProcess::FurmanPiviProcess(const DetectorConstruction* detector)
    : G4VDiscreteProcess("FurmanPiviSEE"), fDetector(detector), fModel(fParameters), fWuModel(fWuParameters) {}

G4bool FurmanPiviProcess::IsApplicable(const G4ParticleDefinition& p)
{ return &p == G4Electron::Definition(); }

G4double FurmanPiviProcess::GetMeanFreePath(const G4Track&, G4double, G4ForceCondition* condition)
{
    *condition = StronglyForced;
    return std::numeric_limits<G4double>::max();
}

G4VParticleChange* FurmanPiviProcess::PostStepDoIt(const G4Track& track, const G4Step& step)
{
    fChange.Initialize(track);
    fChange.returnToPore = false;
    const auto* pre = step.GetPreStepPoint();
    const auto* post = step.GetPostStepPoint();
    const auto* before = pre->GetPhysicalVolume();
    const auto* after = post->GetPhysicalVolume();
    if (post->GetStepStatus() != fGeomBoundary || !before || !after
        || before->GetName() != "Channel"
        || (after->GetName() != "MCP" && before != after)) return &fChange;

    // Channel -> Drift ou World est deja exclu : seule une paroi peut emettre.

    // Le plafond bloque seulement les vrais secondaires, pas les reflexions.
    const auto generation = GenerationInfo::Of(track);
    const auto maxGenerations = fDetector->MaxGenerations();
    const bool atGenerationLimit = maxGenerations > 0 && generation >= maxGenerations;

    // La normale sortante du vide pointe dans le verre : on l’inverse.
    G4bool valid = false;
    auto* navigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
    const auto normalOut = navigator->GetGlobalExitNormal(post->GetPosition(), &valid);
    if (!valid || normalOut.mag2() == 0) return &fChange;
    // Une eventuelle etape de navigation sur une face ouverte ne produit pas
    // d'emission. Pres du bord, une vraie normale laterale reste une paroi.
    const auto tolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
    if (std::abs(normalOut.z()) > 1-1e-12
        && (std::abs(post->GetPosition().z()) <= tolerance
            || std::abs(post->GetPosition().z()-fDetector->Thickness()) <= tolerance)) return &fChange;
    const auto normal = -normalOut.unit();

    // reconstruction d’énergie : K_pre + q E·dx.
    const auto displacement = post->GetPosition() - pre->GetPosition();
    const auto work = track.GetDynamicParticle()->GetCharge() * fDetector->Field().dot(displacement);
    const auto energy = std::max(0.0, pre->GetKineticEnergy() + work);
    const auto direction = post->GetMomentumDirection();
    const auto angle = std::acos(std::clamp(-direction.dot(normal), 0.0, 1.0));

    auto* event = const_cast<EventAction*>(static_cast<const EventAction*>(
        G4RunManager::GetRunManager()->GetUserEventAction()));
    // Compte aussi les impacts qui n'emettront rien. Aucun tirage consomme.
    event->RecordImpact(energy, angle);

    const auto emission = fDetector->UseWuModel()
        ? fWuModel.GenerateEmission(energy, angle, direction, normal)
        : fModel.GenerateEmission(energy, angle, direction, normal);

    // Une reflexion est la continuation de l'incident, pas un nouvel electron.
    if (emission.electrons.size() == 1) {
        const auto& electron = emission.electrons.front();
        if (electron.type != SEEType::TrueSecondary && electron.kineticEnergy > 0
            && electron.direction.dot(normal) > 0) {
            fChange.ProposeTrackStatus(fAlive);
            fChange.ProposeEnergy(electron.kineticEnergy);
            fChange.ProposeMomentumDirection(electron.direction);
            fChange.ProposePosition(EmissionPosition(step, normal, fParameters.surfaceOffset));
            fChange.ProposeLastStepInVolume(false);
            fChange.SetNumberOfSecondaries(0);
            fChange.returnToPore = true;
            return &fChange;
        }
    }

    fChange.ProposeTrackStatus(fStopAndKill);
    fChange.ProposeEnergy(0);
    std::vector<const EmittedElectron*> candidates;
    if (!atGenerationLimit) {
        for (const auto& electron : emission.electrons) {
            if (electron.type == SEEType::TrueSecondary && electron.kineticEnergy > 0
                && electron.direction.dot(normal) > 0) candidates.push_back(&electron);
        }
    }
    const auto keep = event->ReserveSecondaries(static_cast<G4int>(candidates.size()));
    // Sous-ensemble uniforme sans remise si l'emission depasse les places.
    if (keep > 0 && keep < static_cast<G4int>(candidates.size())) {
        for (G4int i = 0; i < keep; ++i) {
            const auto j = i + static_cast<G4int>(G4UniformRand() * (candidates.size() - i));
            std::swap(candidates[i], candidates[j]);
        }
    }
    std::vector<G4Track*> children;
    const auto position = keep ? EmissionPosition(step, normal, fParameters.surfaceOffset) : post->GetPosition();
    for (G4int i = 0; i < keep; ++i) {
        const auto& electron = *candidates[i];
        auto* particle = new G4DynamicParticle(G4Electron::Definition(), electron.direction, electron.kineticEnergy);
        auto* child = new G4Track(particle, post->GetGlobalTime(),
            position);
        child->SetUserInformation(new GenerationInfo(generation + 1));
        children.push_back(child);
    }
    fChange.SetNumberOfSecondaries(static_cast<G4int>(children.size()));
    for (auto* child : children) fChange.AddSecondary(child);
    return &fChange;
}
