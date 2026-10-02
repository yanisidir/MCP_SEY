#include "DetectorConstruction.hh"
#include "G4GenericMessenger.hh"
#include "G4RunManager.hh"
#include "G4StateManager.hh"
#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4PhysicalConstants.hh"
#include "G4UniformElectricField.hh"
#include "G4EqMagElectricField.hh"
#include "G4DormandPrince745.hh"
#include "G4IntegrationDriver.hh"
#include "G4ChordFinder.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4AutoDelete.hh"
#include "G4VisAttributes.hh"
#include <algorithm>
#include <cmath>

DetectorConstruction::DetectorConstruction()
{
    fMessenger = new G4GenericMessenger(this, "/mcp/", "Reglages du monocanal");
    fMessenger->DeclareMethodWithUnit("setPlateDiameter", "mm", &DetectorConstruction::SetPlateDiameter)
        .SetParameterName("d", false).SetRange("d>0").SetStates(G4State_PreInit);
    fMessenger->DeclareMethodWithUnit("setThickness", "mm", &DetectorConstruction::SetThickness)
        .SetParameterName("l", false).SetRange("l>0").SetStates(G4State_PreInit, G4State_Idle)
        .SetToBeBroadcasted(false); // Parametre partage, modifie seulement par le master.
    fMessenger->DeclareMethodWithUnit("setChannelDiameter", "um", &DetectorConstruction::SetChannelDiameter)
        .SetParameterName("d", false).SetRange("d>0").SetStates(G4State_PreInit);
    fMessenger->DeclareMethodWithUnit("setChannelAngle", "deg", &DetectorConstruction::SetAngle)
        .SetParameterName("a", false).SetRange("a>=0 && a<90").SetStates(G4State_PreInit);
    fMessenger->DeclareMethodWithUnit("setVoltage", "volt", &DetectorConstruction::SetVoltage)
        .SetParameterName("v", false).SetRange("v>=0").SetStates(G4State_PreInit, G4State_Idle)
        .SetToBeBroadcasted(false); // Parametre partage, modifie seulement par le master.
    fMessenger->DeclareProperty("setOutputFile", fOutput).SetStates(G4State_PreInit, G4State_Idle);
    fMessenger->DeclareProperty("setTransitOutput", fSaveExits).SetStates(G4State_PreInit);
    fMessenger->DeclareProperty("setWuModel", fUseWuModel)
        .SetGuidance("true : modele de Wu et al. 2008 (defaut) ; false : Furman-Pivi/Peng.")
        .SetStates(G4State_PreInit, G4State_Idle);
    fMessenger->DeclareProperty("setMaxTracksPerEvent", fMaxTracks)
        .SetParameterName("n", false).SetRange("n>=0").SetStates(G4State_PreInit);
    // Plafonne les vrais secondaires ; les reflexions restent autorisees.
    fMessenger->DeclareProperty("setMaxGenerations", fMaxGenerations)
        .SetParameterName("n", false).SetRange("n>=0").SetStates(G4State_PreInit, G4State_Idle);
}
DetectorConstruction::~DetectorConstruction() { delete fMessenger; }

void DetectorConstruction::SetThickness(G4double thickness)
{
    if (thickness == fThickness) return;
    fThickness = thickness;
    // Changer un parametre ne suffit pas : la geometrie doit etre reconstruite,
    // sinon chaque thread garde l'ancienne dans son gestionnaire de transport.
    // Ici l'argument destroyFirst vide les stores de solides -- necessaire car
    // l'epaisseur redimensionne TOUS les volumes, et sans cela ceux des runs
    // precedents s'accumulent en memoire.
    if (G4StateManager::GetStateManager()->GetCurrentState() == G4State_Idle)
        G4RunManager::GetRunManager()->ReinitializeGeometry(true);
}

void DetectorConstruction::SetVoltage(G4double voltage)
{
    if (voltage == fVoltage) return;
    fVoltage = voltage;
    // Meme reconstruction, mais sans destroyFirst : la tension ne change que
    // l'intensite du champ, les solides restent valables.
    if (G4StateManager::GetStateManager()->GetCurrentState() == G4State_Idle)
        G4RunManager::GetRunManager()->ReinitializeGeometry();
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
    // Le pore doit rester dans le disque de verre, même avec son inclinaison.
    if (fThickness * std::tan(fAngle) + fChannelDiameter / (2 * std::cos(fAngle))
        >= fPlateDiameter / 2) {
        G4Exception("DetectorConstruction", "Geometry", FatalException,
                    "Le pore sort lateralement de la plaque : augmenter setPlateDiameter.");
    }
    auto* nist = G4NistManager::Instance();
    auto* vacuum = nist->FindOrBuildMaterial("G4_Galactic");
    auto* glass = G4Material::GetMaterial("LeadGlassMCP", false);
    if (!glass) {
        glass = new G4Material("LeadGlassMCP", 4.93 * g/cm3, 3);
        glass->AddElement(nist->FindOrBuildElement("Pb"), 0.60);
        glass->AddElement(nist->FindOrBuildElement("Si"), 0.16);
        glass->AddElement(nist->FindOrBuildElement("O"), 0.24);
    }
    auto* worldSolid = new G4Box("World", std::max(2.5*mm, fPlateDiameter),

        std::max(2.5*mm, fPlateDiameter), std::max(2.5*mm, fThickness+1*mm));

    auto* worldLogic = new G4LogicalVolume(worldSolid, vacuum, "World");

    auto* world = new G4PVPlacement(nullptr, {}, worldLogic, "World", nullptr, false, 0, true);

    auto* blank = new G4Tubs("Plate", 0, fPlateDiameter/2, fThickness/2, 0, twopi);
    // Outil de soustraction etendu pour traverser les deux faces du verre.
    auto* pore = new G4Tubs("Pore", 0, fChannelDiameter/2,
        fThickness/(2*std::cos(fAngle)) + 10*um, 0, twopi);

    auto* rotation = new G4RotationMatrix;
    rotation->rotateX(fAngle);
    G4ThreeVector offset(0, fThickness/2*std::tan(fAngle), 0);

    auto* solid = new G4SubtractionSolid("GlassWithPore", blank, pore, rotation, offset);
    auto* glassLogic = new G4LogicalVolume(solid, glass, "MCP");

    new G4PVPlacement(nullptr, {0,0,fThickness/2}, glassLogic, "MCP", worldLogic, false, 0, true);

    // Le cylindre etendu sert seulement a percer le verre. Le volume de vide
    // portant le champ s'arrete exactement aux faces z=0 et z=L.
    auto* channelSolid = new G4IntersectionSolid("PoreInPlate", blank, pore, rotation, offset);
    fPoreLogic = new G4LogicalVolume(channelSolid, vacuum, "Channel");
    new G4PVPlacement(nullptr, {0,0,fThickness/2}, fPoreLogic, "Channel", worldLogic, false, 0, true);

    // Vide sans champ. La face aval de ce volume est notre plan de mesure.
    auto* driftSolid = new G4Box("Drift", std::max(2.5*mm, fPlateDiameter),
        std::max(2.5*mm, fPlateDiameter), DriftDistance()/2);
    auto* driftLogic = new G4LogicalVolume(driftSolid, vacuum, "Drift");
    new G4PVPlacement(nullptr, {0,0,fThickness+DriftDistance()/2},
        driftLogic, "Drift", worldLogic, false, 0, true);
    driftLogic->SetVisAttributes(G4VisAttributes::GetInvisible());

    worldLogic->SetVisAttributes(G4VisAttributes::GetInvisible());
    glassLogic->SetVisAttributes(G4VisAttributes(G4Colour(0.2,0.4,0.9,0.25)));

    G4VisAttributes poreView(G4Colour(1,0.3,0));

    poreView.SetForceWireframe(true);
    fPoreLogic->SetVisAttributes(poreView);

    G4cout << "L=" << fThickness/mm << " mm, D=" << fChannelDiameter/um
           << " um, angle=" << fAngle/degree << " deg, V=" << fVoltage/volt << " V" << G4endl;
    return world;
}

void DetectorConstruction::ConstructSDandField()
{
    // Champ unique : uniforme, suivant -z. L’électron est accéléré vers +z.
    auto* field = new G4UniformElectricField(Field());
    G4AutoDelete::Register(field);
    auto* equation = new G4EqMagElectricField(field);
    auto* stepper = new G4DormandPrince745(equation, 8);
    auto* driver = new G4IntegrationDriver<G4DormandPrince745>(1*nm, stepper, 8);
    auto* chord = new G4ChordFinder(driver);
    G4AutoDelete::Register(chord);
    auto* manager = new G4FieldManager;
    // G4FieldManagerStore gere la destruction du gestionnaire.
    manager->SetDetectorField(field);
    manager->SetFieldChangesEnergy(true);
    manager->SetChordFinder(chord);
    fPoreLogic->SetFieldManager(manager, true); // Aucun champ dans World ou Drift.
    manager->SetDeltaOneStep(0.01*nm);
    manager->SetDeltaIntersection(0.01*nm);
    chord->SetDeltaChord(0.1*nm);
    manager->SetMinimumEpsilonStep(1e-9);
    manager->SetMaximumEpsilonStep(1e-8);
}
