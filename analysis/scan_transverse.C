// Reutilise Spot pour garder exactement les memes definitions de largeur.
#include "transverse.C"

#include <TGraphErrors.h>
#include <TLine.h>
#include <TError.h>
#include <TRegexp.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>

#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace TransverseScan {

// ============================================================================
// 1. CONSTANTES ET TYPES
// ============================================================================

// Ordre commun aux calculs, graphiques et colonnes CSV.
constexpr int kNumObservables = 12;
constexpr int kNumConfigValues = 7;

const char* kObservableNames[kNumObservables] = {
    "SigmaXExit_um",
    "SigmaYExit_um",
    "SigmaRExit_um",
    "EnergyExit_eV",
    "SigmaX50_um",
    "SigmaY50_um",
    "SigmaR50_um",
    "Energy50_eV",
    "DeltaSigmaR_um",
    "R80Exit_um",
    "R8050_um",
    "DeltaR80_um"
};

// Indices des differences par avalanche : axe vertical signe.
constexpr int kDeltaSigmaR = 8;
constexpr int kDeltaR80 = 11;

const char* kConfigNames[] = {
    "Thickness_mm",
    "PoreDiameter_um",
    "Angle_deg",
    "Voltage_V",
    "MaxTracksPerEvent",
    "DriftDistance_um"
};


// ============================================================================
// 2. STATISTIQUES
// ============================================================================

// Statistique ENTRE avalanches : chaque avalanche a le meme poids.
struct Average {
    long n = 0;
    double mean = 0.0;
    double m2 = 0.0;

    void Add(double x)
    {
        ++n;

        const double delta = x - mean;
        mean += delta / n;
        m2 += delta * (x - mean);
    }

    double Mean() const
    {
        return n
            ? mean
            : std::numeric_limits<double>::quiet_NaN();
    }

    double SD() const
    {
        return n > 1
            ? std::sqrt(std::max(0.0, m2 / (n - 1)))
            : std::numeric_limits<double>::quiet_NaN();
    }

    double SEM() const
    {
        return SD() / std::sqrt(n);
    }
};


// ============================================================================
// 3. DONNEES D'UN POINT DE SCAN
// ============================================================================

struct Point {
    int generation = 0;

    long events = 0;
    long limited = 0;
    long generationReached = 0;
    long selected = 0;

    long electronsExit = 0;
    long electrons50 = 0;

    double reached = 0.0;

    std::array<Average, kNumObservables> values;
    std::array<double, kNumConfigValues> config;
};


// ============================================================================
// 4. OUTILS DE LECTURE ROOT
// ============================================================================

static double Value(TTree* tree, const char* name)
{
    auto* leaf = tree->GetLeaf(name);

    if (!leaf) {
        throw std::runtime_error(
            std::string("Colonne absente : ") + name
        );
    }

    return leaf->GetValue();
}


static void ValidateConfiguration(
    TTree* configuration,
    const std::string& path)
{
    if (configuration->GetEntries() != 1) {
        throw std::runtime_error(
            "Un seul run attendu par fichier : " + path
        );
    }

    configuration->GetEntry(0);

    const bool fieldOnlyInPore =
        Value(configuration, "FieldOnlyInPore") == 1;

    const bool scoringAtPlateFace =
        Value(configuration, "ScoringAtPlateFace") == 1;

    const bool correctDriftDistance =
        std::abs(
            Value(configuration, "DriftDistance_um") - 50.0
        ) <= 1e-9;

    if (!fieldOnlyInPore ||
        !scoringAtPlateFace ||
        !correctDriftDistance) {

        throw std::runtime_error(
            "Mesures attendues a la face MCP et +50 um sans champ : "
            + path
        );
    }
}


static Point Read(const std::string& path)
{
    // ------------------------------------------------------------------------
    // 4.1 Ouvrir le fichier et recuperer les arbres
    // ------------------------------------------------------------------------

    TFile file(path.c_str());

    auto* configuration = file.Get<TTree>("configuration");
    auto* events        = file.Get<TTree>("events");
    auto* exits         = file.Get<TTree>("exits");
    auto* downstream    = file.Get<TTree>("downstream");

    if (file.IsZombie() ||
        !configuration ||
        !events ||
        !exits ||
        !downstream) {

        throw std::runtime_error(
            "Tables manquantes : " + path
        );
    }


    // ------------------------------------------------------------------------
    // 4.2 Lire et verifier la configuration
    // ------------------------------------------------------------------------

    ValidateConfiguration(configuration, path);

    Point point;

    point.generation =
        static_cast<int>(
            Value(configuration, "MaxGenerations")
        );

    for (int i = 0; i < 6; ++i) {
        point.config[i] =
            Value(configuration, kConfigNames[i]);
    }

    // L'energie incidente est lue ensuite dans l'arbre events.
    point.config[6] = 0.0;


    // ------------------------------------------------------------------------
    // 4.3 Recuperer les compteurs par evenement
    // ------------------------------------------------------------------------

    using EventKey = std::pair<int, int>;

    // spots[event][0] = face MCP
    // spots[event][1] = plan +50 um
    std::map<EventKey, std::array<Spot, 2>> spots;

    // counts[event] = {Gain, NAt50um}
    std::map<EventKey, std::pair<long, long>> counts;

    point.events = events->GetEntries();

    for (Long64_t i = 0; i < events->GetEntries(); ++i) {
        events->GetEntry(i);

        const EventKey key = {
            static_cast<int>(Value(events, "RunID")),
            static_cast<int>(Value(events, "EventID"))
        };

        if (counts.count(key)) {
            throw std::runtime_error(
                "Evenement duplique : " + path
            );
        }

        counts[key] = {
            static_cast<long>(Value(events, "Gain")),
            static_cast<long>(Value(events, "NAt50um"))
        };

        point.limited +=
            Value(events, "IsMultiplicationLimited") != 0;
        point.generationReached +=
            point.generation > 0 && Value(events, "MaxGeneration") >= point.generation;

        point.reached +=
            Value(events, "MaxGeneration");

        const double energy =
            Value(events, "IncidentEnergy_eV");

        if (i > 0 &&
            std::abs(energy - point.config[6]) > 1e-8) {

            throw std::runtime_error(
                "Energie source variable : " + path
            );
        }

        point.config[6] = energy;
    }

    if (point.events > 0) {
        point.reached /= point.events;
    }


    // ------------------------------------------------------------------------
    // 4.4 Lire les passages a la face MCP et a +50 um
    // ------------------------------------------------------------------------

    const std::array<TTree*, 2> tables = {
        exits,
        downstream
    };

    for (int plane = 0; plane < 2; ++plane) {
        TTree* tree = tables[plane];

        for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
            tree->GetEntry(i);

            const EventKey key = {
                static_cast<int>(Value(tree, "RunID")),
                static_cast<int>(Value(tree, "EventID"))
            };

            if (!counts.count(key)) {
                throw std::runtime_error(
                    "Passage sans evenement : " + path
                );
            }

            const double x_um =
                Value(tree, "X_mm") * 1000.0;

            const double y_um =
                Value(tree, "Y_mm") * 1000.0;

            const double energy_eV =
                Value(tree, "KineticEnergy_eV");

            if (!std::isfinite(x_um) ||
                !std::isfinite(y_um) ||
                !std::isfinite(energy_eV)) {

                throw std::runtime_error(
                    "Mesure non finie : " + path
                );
            }

            spots[key][plane].Add(
                x_um,
                y_um,
                energy_eV
            );
        }
    }

    point.electronsExit = exits->GetEntries();
    point.electrons50   = downstream->GetEntries();


    // ------------------------------------------------------------------------
    // 4.5 Construire les observables avalanche par avalanche
    // ------------------------------------------------------------------------

    for (const auto& entry : counts) {
        const EventKey& key = entry.first;

        const auto& expectedCounts = entry.second;
        const auto& eventSpots = spots[key];

        const Spot& exitSpot = eventSpots[0];
        const Spot& spot50   = eventSpots[1];

        // Verification de coherence entre les arbres.
        if (exitSpot.n != expectedCounts.first ||
            spot50.n   != expectedCounts.second) {

            throw std::runtime_error(
                "Passages incomplets : activer setTransitOutput true. "
                "Fichier : " + path
            );
        }

        // Il faut au moins 2 electrons dans chaque plan
        // pour definir une largeur.
        if (exitSpot.n < 2 || spot50.n < 2) {
            continue;
        }

        ++point.selected;

        // Calcule une fois : chaque appel trie les distances de l'avalanche.
        const double exitR80 = exitSpot.R80();
        const double r80At50 = spot50.R80();

        const double observables[kNumObservables] = {
            exitSpot.Sx(),
            exitSpot.Sy(),
            exitSpot.Sr(),
            exitSpot.energy,

            spot50.Sx(),
            spot50.Sy(),
            spot50.Sr(),
            spot50.energy,

            spot50.Sr() - exitSpot.Sr(),

            exitR80,
            r80At50,
            r80At50 - exitR80
        };

        for (int i = 0; i < kNumObservables; ++i) {
            point.values[i].Add(observables[i]);
        }
    }

    return point;
}


// ============================================================================
// 5. RECHERCHE DES FICHIERS DU SCAN
// ============================================================================

static std::vector<Point> LoadScan(const char* motif)
{
    const TString folder =
        gSystem->DirName(motif);

    const TString pattern =
        gSystem->BaseName(motif);

    TSystemDirectory directory(folder, folder);

    auto* files = directory.GetListOfFiles();

    if (!files) {
        throw std::runtime_error(
            "Repertoire introuvable"
        );
    }

    TRegexp filter(pattern, kTRUE);
    TIter next(files);

    std::vector<Point> points;
    std::set<int> generations;

    while (auto* entry =
        static_cast<TSystemFile*>(next())) {

        const TString name = entry->GetName();
        Ssiz_t length = 0;

        const bool matchesPattern =
            name.Index(filter, &length) == 0 &&
            length == name.Length();

        if (entry->IsDirectory() ||
            !matchesPattern) {
            continue;
        }

        Point point =
            Read(
                std::string(folder.Data()) +
                "/" +
                name.Data()
            );

        if (!point.events) {
            std::cout
                << "Fichier vide ignore : "
                << name
                << '\n';

            continue;
        }

        if (!generations.insert(point.generation).second) {
            throw std::runtime_error(
                "Deux fichiers pour le meme plafond : "
                "preciser le motif."
            );
        }

        // Toutes les configurations doivent etre identiques,
        // sauf MaxGenerations qui est justement la variable du scan.
        if (!points.empty()) {
            for (int i = 0; i < kNumConfigValues; ++i) {
                if (std::abs(
                        point.config[i] -
                        points.front().config[i]
                    ) > 1e-8) {

                    throw std::runtime_error(
                        "Geometrie, tension, energie ou plafond "
                        "de tracks differents entre fichiers."
                    );
                }
            }
        }

        points.push_back(point);
    }

    if (points.empty()) {
        throw std::runtime_error(
            "Aucun run exploitable pour ce motif"
        );
    }

    std::sort(
        points.begin(),
        points.end(),
        [](const Point& a, const Point& b) {
            return a.generation < b.generation;
        }
    );

    return points;
}


// ============================================================================
// 6. EXPORT CSV ET RESUME CONSOLE
// ============================================================================

static void WriteCSV(
    const std::vector<Point>& points,
    const std::string& figure)
{
    std::ofstream csv(figure + ".csv");

    if (!csv) {
        throw std::runtime_error(
            "CSV non accessible"
        );
    }

    csv
        << std::setprecision(10)
        << "MaxGenerations,"
        << "Events,"
        << "Limited,GenerationReached,"
        << "Selected,"
        << "ElectronsExit,"
        << "Electrons50,"
        << "MeanMaxGeneration";

    for (const char* name : kObservableNames) {
        csv
            << ',' << name << "_Mean"
            << ',' << name << "_SD"
            << ',' << name << "_SEM";
    }

    csv << '\n';

    for (const Point& point : points) {
        csv
            << point.generation << ','
            << point.events << ','
            << point.limited << ','
            << point.generationReached << ','
            << point.selected << ','
            << point.electronsExit << ','
            << point.electrons50 << ','
            << point.reached;

        for (const Average& value : point.values) {
            csv
                << ',' << value.Mean()
                << ',' << value.SD()
                << ',' << value.SEM();
        }

        csv << '\n';
    }
}


static void PrintSummary(
    const std::vector<Point>& points)
{
    std::cout << "lim_gen : evenements ayant atteint le plafond de generations (0 si sans plafond), toutes injections incluses.\n";
    std::cout << "gen.max : moyenne de MaxGeneration sur toutes les injections.\n";
    printf(
        "%8s %7s %7s %7s %9s %7s %12s %12s %12s %12s %12s %12s\n",
        "plafond",
        "evts",
        "limites",
        "lim_gen",
        "gen.max",
        "retenus",
        "sigmaR_sortie",
        "sigmaR_50um",
        "R80_sortie",
        "R80_50um",
        "E_sortie",
        "E_50um"
    );

    for (const Point& point : points) {
        const std::string label =
            point.generation
                ? std::to_string(point.generation)
                : "libre";

        printf(
            "%8s %7ld %7ld %7ld %9.3f %7ld %12.4g %12.4g %12.4g %12.4g %12.4g %12.4g\n",
            label.c_str(),
            point.events,
            point.limited,
            point.generationReached,
            point.reached,
            point.selected,
            point.values[2].Mean(),
            point.values[6].Mean(),
            point.values[9].Mean(),
            point.values[10].Mean(),
            point.values[3].Mean(),
            point.values[7].Mean()
        );

        if (point.selected < 20) {
            std::cout
                << "Plafond "
                << label
                << " : seulement "
                << point.selected
                << " avalanches retenues.\n";
        }
    }

    std::cout << "um, eV ; poids egal par avalanche ; R80 = rayon a "
              << 100.0 * kContainmentFraction << " % des electrons.\n";
}


// ============================================================================
// 7. GRAPHIQUES
// ============================================================================

static bool ComputeGenerationRange(
    const std::vector<Point>& points,
    double& xmin,
    double& xmax)
{
    xmin = 1e30;
    xmax = -1e30;

    for (const Point& point : points) {
        if (!point.generation) {
            continue;
        }

        xmin = std::min(
            xmin,
            static_cast<double>(point.generation)
        );

        xmax = std::max(
            xmax,
            static_cast<double>(point.generation)
        );
    }

    if (xmin > xmax) {
        return false;
    }

    const double margin =
        std::max(
            1.0,
            0.05 * (xmax - xmin)
        );

    xmin -= margin;
    xmax += margin;

    return true;
}


static void DrawPanel(
    TCanvas& canvas,
    const std::vector<Point>& points,
    int pad,
    int first,
    int second,
    const char* title,
    double xmin,
    double xmax,
    bool showSD = false)
{
    canvas.cd(pad);

    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.06);

    double ymin = 0.0;
    double ymax = 0.0;

    // Chercher l'echelle verticale.
    for (const Point& point : points) {
        for (int index : {first, second}) {
            if (index < 0) {
                continue;
            }

            const double y =
                showSD
                    ? point.values[index].SD()
                    : point.values[index].Mean();

            const double error =
                (!showSD && point.selected > 1)
                    ? point.values[index].SEM()
                    : 0.0;

            if (std::isfinite(y)) {
                ymin = std::min(ymin, y - error);
                ymax = std::max(ymax, y + error);
            }
        }
    }

    const double gap =
        std::max(
            1e-3,
            0.35 * (ymax - ymin)
        );

    const double frameYMin =
        (first == kDeltaSigmaR || first == kDeltaR80)
            ? ymin - gap
            : std::max(0.0, ymin - gap);

    gPad->DrawFrame(
        xmin,
        frameYMin,
        xmax,
        ymax + gap,
        title
    );

    auto* legend =
        new TLegend(0.17, 0.80, 0.90, 0.90);

    legend->SetNColumns(
        second >= 0 ? 2 : 1
    );

    legend->SetTextSize(0.033);

    int curve = 0;

    for (int index : {first, second}) {
        if (index < 0) {
            continue;
        }

        auto* graph = new TGraphErrors;

        graph->SetMarkerStyle(
            curve ? 24 : 20
        );

        graph->SetMarkerColor(
            curve ? kRed : kBlue
        );

        graph->SetLineColor(
            curve ? kRed : kBlue
        );

        for (const Point& point : points) {
            const double y =
                showSD
                    ? point.values[index].SD()
                    : point.values[index].Mean();

            if (!std::isfinite(y)) {
                continue;
            }

            // generation = 0 : reference sans plafond.
            // On l'affiche comme une ligne horizontale.
            if (!point.generation) {
                auto* line =
                    new TLine(
                        xmin,
                        y,
                        xmax,
                        y
                    );

                line->SetLineColor(
                    curve ? kRed : kBlue
                );

                line->SetLineStyle(2);
                line->Draw();

                continue;
            }

            const int n =
                graph->GetN();

            graph->SetPoint(
                n,
                point.generation,
                y
            );

            const double error =
                (!showSD && point.selected > 1)
                    ? point.values[index].SEM()
                    : 0.0;

            graph->SetPointError(
                n,
                0.0,
                error
            );
        }

        if (graph->GetN()) {
            graph->Draw("PL SAME");
        }

        const char* label =
            second >= 0
                ? (curve ? "+50 um" : "Face MCP")
                : "Difference par avalanche";

        legend->AddEntry(
            graph,
            label,
            "pl"
        );

        ++curve;
    }

    legend->Draw();
}


static void DrawPlots(
    const std::vector<Point>& points,
    const std::string& figure)
{
    TCanvas canvas(
        "transverse_scan",
        "Dispersion et energie selon les generations",
        1800,
        850
    );

    canvas.Divide(4, 2);

    double xmin = 0.0;
    double xmax = 0.0;

    if (!ComputeGenerationRange(
            points,
            xmin,
            xmax)) {

        std::cout
            << "Seulement une reference libre : "
            << "CSV produit, pas de courbe en generations.\n";

        return;
    }

    DrawPanel(
        canvas,
        points,
        1,
        0,
        4,
        "Moyenne de sigma_x;"
        "Plafond de generations;"
        "sigma_x (um)",
        xmin,
        xmax
    );

    DrawPanel(
        canvas,
        points,
        2,
        1,
        5,
        "Moyenne de sigma_y;"
        "Plafond de generations;"
        "sigma_y (um)",
        xmin,
        xmax
    );

    DrawPanel(
        canvas,
        points,
        3,
        2,
        6,
        "Moyenne de sigma_r;"
        "Plafond de generations;"
        "sigma_r (um)",
        xmin,
        xmax
    );

    DrawPanel(
        canvas,
        points,
        4,
        9,
        10,
        "Moyenne de R80;"
        "Plafond de generations;"
        "R80 (um)",
        xmin,
        xmax
    );

    DrawPanel(
        canvas,
        points,
        5,
        kDeltaSigmaR,
        -1,
        "Elargissement radial moyen;"
        "Plafond de generations;"
        "Delta sigma_r (um)",
        xmin,
        xmax
    );

    DrawPanel(
        canvas,
        points,
        6,
        kDeltaR80,
        -1,
        "Elargissement R80 moyen;"
        "Plafond de generations;"
        "Delta R80 (um)",
        xmin,
        xmax
    );

    DrawPanel(
        canvas,
        points,
        7,
        3,
        7,
        "Moyenne de l'energie moyenne par avalanche;"
        "Plafond de generations;"
        "Energie (eV)",
        xmin,
        xmax
    );

    DrawPanel(
        canvas,
        points,
        8,
        2,
        6,
        "Dispersion inter-avalanches de sigma_r;"
        "Plafond de generations;"
        "SD(sigma_r) (um)",
        xmin,
        xmax,
        true
    );

    canvas.SaveAs(
        (figure + ".png").c_str()
    );

    canvas.SaveAs(
        (figure + ".pdf").c_str()
    );
}

} // namespace TransverseScan


// ============================================================================
// 8. POINT D'ENTREE
// ============================================================================
//
// Exemple :
// root -l -b -q 'analysis/scan_transverse.C("transverse_gen*.root")'
//
void scan_transverse(
    const char* motif = "transverse_gen*.root",
    const char* figure = "scan_transverse")
{
    using namespace TransverseScan;
    gErrorIgnoreLevel = kWarning;

    const std::vector<Point> points =
        LoadScan(motif);

    WriteCSV(
        points,
        figure
    );

    PrintSummary(
        points
    );

    DrawPlots(
        points,
        figure
    );
}
