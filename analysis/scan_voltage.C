#include <TCanvas.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TLeaf.h>
#include <TRegexp.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TTree.h>
#include <TVirtualPad.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace VoltageScan {

using EventKey = std::pair<int, int>;
constexpr int kNumConfigValues = 7;

// Moyenne et variance en ligne : pas besoin de stocker toutes les valeurs.
struct Stats {
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
        return n ? mean : std::numeric_limits<double>::quiet_NaN();
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

// Statistiques spatiales d'un paquet d'electrons.
struct Packet {
    Stats x;
    Stats y;
};

// Resultats associes a une tension.
struct Point {
    double voltage = 0.0;
    long events = 0;
    long limited = 0;
    long generationReached = 0;

    Stats transit;
    Stats radial;
    Stats gain;
    Stats maxGeneration;

    // L, D, angle, cap tracks, cap generations, distance, energie injectee.
    std::array<double, kNumConfigValues> config;
};

static double ReadValue(TTree* tree, const char* name)
{
    TLeaf* leaf = tree->GetLeaf(name);

    if (!leaf) {
        throw std::runtime_error(std::string("Colonne absente : ") + name);
    }

    return leaf->GetValue();
}

static EventKey ReadEventKey(TTree* tree)
{
    return {
        static_cast<int>(ReadValue(tree, "RunID")),
        static_cast<int>(ReadValue(tree, "EventID"))
    };
}

static void ReadConfiguration(
    TTree* configuration,
    Point& point,
    const std::string& path)
{
    if (configuration->GetEntries() != 1) {
        throw std::runtime_error("Un seul run attendu par fichier : " + path);
    }

    configuration->GetEntry(0);

    point.voltage = ReadValue(configuration, "Voltage_V");

    const char* fields[] = {
        "Thickness_mm",
        "PoreDiameter_um",
        "Angle_deg",
        "MaxTracksPerEvent",
        "MaxGenerations",
        "DriftDistance_um"
    };

    for (int i = 0; i < 6; ++i) {
        point.config[i] = ReadValue(configuration, fields[i]);
    }

    // Energie incidente : elle sera lue ensuite dans l'arbre events.
    point.config[6] = 0.0;

    const bool correctDistance =
        std::abs(point.config[5] - 50.0) <= 1e-8;

    const bool fieldOnlyInPore =
        ReadValue(configuration, "FieldOnlyInPore") == 1;

    const bool scoringAtPlateFace =
        ReadValue(configuration, "ScoringAtPlateFace") == 1;

    if (!correctDistance || !fieldOnlyInPore || !scoringAtPlateFace) {
        throw std::runtime_error(
            "Mesures attendues a la face MCP et a +50 um sans champ : " + path
        );
    }
}

static Point Read(const std::string& path)
{
    TFile file(path.c_str());

    TTree* configuration = file.Get<TTree>("configuration");
    TTree* events = file.Get<TTree>("events");
    TTree* hits = file.Get<TTree>("downstream");

    if (file.IsZombie() || !configuration || !events || !hits) {
        throw std::runtime_error(
            "Un run avec configuration, events et downstream attendu : " + path
        );
    }

    Point point;
    ReadConfiguration(configuration, point, path);

    std::map<EventKey, long> expectedHits;
    std::map<EventKey, Packet> packets;

    point.events = events->GetEntries();

    for (Long64_t i = 0; i < events->GetEntries(); ++i) {
        events->GetEntry(i);

        const EventKey key = ReadEventKey(events);

        if (expectedHits.count(key)) {
            throw std::runtime_error("Evenement duplique : " + path);
        }

        expectedHits[key] =
            static_cast<long>(ReadValue(events, "NAt50um"));

        point.limited +=
            ReadValue(events, "IsMultiplicationLimited") != 0;
        point.maxGeneration.Add(ReadValue(events, "MaxGeneration"));
        point.generationReached +=
            point.config[4] > 0 && ReadValue(events, "MaxGeneration") >= point.config[4];

        const double incidentEnergy =
            ReadValue(events, "IncidentEnergy_eV");

        if (i > 0 && std::abs(incidentEnergy - point.config[6]) > 1e-8) {
            throw std::runtime_error("Energie injectee variable : " + path);
        }

        point.config[6] = incidentEnergy;

        const double gain = ReadValue(events, "Gain");

        // Gain moyen sur TOUTES les injections, gain nul compris : c'est la
        // convention sortie/entree, et non une moyenne conditionnelle.
        point.gain.Add(gain);

        const double meanTimePs =
            ReadValue(events, "MeanTime_ns") * 1000.0;

        // Un barycentre temporel n'a de sens qu'avec au moins deux electrons.
        if (gain >= 2.0) {
            if (!std::isfinite(meanTimePs)) {
                throw std::runtime_error("Temps non fini : " + path);
            }

            point.transit.Add(meanTimePs);
        }
    }

    for (Long64_t i = 0; i < hits->GetEntries(); ++i) {
        hits->GetEntry(i);

        const EventKey key = ReadEventKey(hits);

        if (!expectedHits.count(key)) {
            throw std::runtime_error("Impact sans evenement : " + path);
        }

        const double xUm = ReadValue(hits, "X_mm") * 1000.0;
        const double yUm = ReadValue(hits, "Y_mm") * 1000.0;

        if (!std::isfinite(xUm) || !std::isfinite(yUm)) {
            throw std::runtime_error("Position non finie : " + path);
        }

        packets[key].x.Add(xUm);
        packets[key].y.Add(yUm);
    }

    for (const auto& item : expectedHits) {
        const EventKey& key = item.first;
        const long expected = item.second;
        const Packet& packet = packets[key];
        const long n = packet.x.n;

        if (n != expected) {
            throw std::runtime_error(
                "Impacts incomplets : activer setTransitOutput true. " + path
            );
        }

        // Largeur INTERNE au paquet d'un evenement (denominateur N, barycentre
        // propre au paquet), distincte de la dispersion ENTRE evenements.
        if (n >= 2) {
            const double sigmaR =
                std::sqrt(
                    std::max(
                        0.0,
                        (packet.x.m2 + packet.y.m2) / n
                    )
                );

            point.radial.Add(sigmaR);
        }
    }

    return point;
}

static std::vector<Point> LoadScan(const char* motif)
{
    const TString folder = gSystem->DirName(motif);
    const TString pattern = gSystem->BaseName(motif);

    TSystemDirectory directory(folder, folder);
    auto* files = directory.GetListOfFiles();

    if (!files) {
        throw std::runtime_error("Repertoire introuvable");
    }

    TRegexp filter(pattern, kTRUE);
    TIter next(files);

    std::vector<Point> points;

    while (auto* entry = static_cast<TSystemFile*>(next())) {
        const TString name = entry->GetName();
        Ssiz_t length = 0;

        const bool matchesPattern =
            name.Index(filter, &length) == 0 &&
            length == name.Length();

        if (entry->IsDirectory() || !matchesPattern) {
            continue;
        }

        Point point = Read(
            std::string(folder.Data()) + "/" + name.Data()
        );

        if (!point.events) {
            std::cout << "Run vide ignore : " << name << '\n';
            continue;
        }

        // Toutes les configurations doivent etre identiques,
        // sauf la tension, qui est la variable du scan.
        if (!points.empty()) {
            for (int i = 0; i < kNumConfigValues; ++i) {
                if (std::abs(
                        point.config[i] - points.front().config[i]
                    ) > 1e-8) {

                    throw std::runtime_error(
                        "Geometrie, energie source ou plafonds differents entre fichiers."
                    );
                }
            }
        }

        points.push_back(point);
    }

    if (points.empty()) {
        throw std::runtime_error("Aucun run exploitable pour ce motif");
    }

    std::sort(
        points.begin(),
        points.end(),
        [](const Point& a, const Point& b) {
            return a.voltage < b.voltage;
        }
    );

    for (size_t i = 1; i < points.size(); ++i) {
        if (std::abs(points[i].voltage - points[i - 1].voltage) < 1e-8) {
            throw std::runtime_error(
                "Plusieurs fichiers a la meme tension : preciser le motif."
            );
        }
    }

    return points;
}

static void WriteCSV(
    const std::vector<Point>& points,
    const std::string& figure)
{
    std::ofstream csv(figure + ".csv");

    if (!csv) {
        throw std::runtime_error("CSV non accessible");
    }

    csv
        << std::setprecision(10)
        << "Voltage_V,Events,Limited,GenerationReached,MeanMaxGeneration,"
        << "NTransit,MeanTransit_ps,SigmaTransit_ps,SEMTransit_ps,"
        << "NRadial,MeanSigmaR50_um,SDSigmaR50_um,SEMSigmaR50_um,"
        << "MeanGain,SDGain,SEMGain\n";

    for (const Point& point : points) {
        csv
            << point.voltage << ','
            << point.events << ','
            << point.limited << ','
            << point.generationReached << ','
            << point.maxGeneration.Mean() << ','
            << point.transit.n << ','
            << point.transit.Mean() << ','
            << point.transit.SD() << ','
            << point.transit.SEM() << ','
            << point.radial.n << ','
            << point.radial.Mean() << ','
            << point.radial.SD() << ','
            << point.radial.SEM() << ','
            << point.gain.Mean() << ','
            << point.gain.SD() << ','
            << point.gain.SEM() << '\n';
    }
}

static void PrintSummary(const std::vector<Point>& points)
{
    std::cout << "lim_gen : evenements ayant atteint le plafond de generations (0 si sans plafond), toutes injections incluses.\n";
    std::cout << "gen.max : moyenne de MaxGeneration sur toutes les injections.\n";
    printf(
        "%7s %7s %7s %7s %9s %7s %12s %12s %7s %12s %12s\n",
        "V",
        "evts",
        "limites",
        "lim_gen",
        "gen.max",
        "Ntemps",
        "temps(ps)",
        "sigmaT(ps)",
        "Nradial",
        "sigmaR(um)",
        "gain"
    );

    for (const Point& point : points) {
        printf(
            "%7.0f %7ld %7ld %7ld %9.3f %7ld %12.4g %12.4g %7ld %12.4g %12.4g\n",
            point.voltage,
            point.events,
            point.limited,
            point.generationReached,
            point.maxGeneration.Mean(),
            point.transit.n,
            point.transit.Mean(),
            point.transit.SD(),
            point.radial.n,
            point.radial.Mean(),
            point.gain.Mean()
        );

        if (point.transit.n < 20 || point.radial.n < 20) {
            std::cout
                << "Attention : moins de 20 avalanches utilisables a "
                << point.voltage
                << " V.\n";
        }
    }

    std::cout
        << "Gain : electrons sortant du MCP par injection, toutes les injections incluses (meme gain nul).\n"
        << "Temps : barycentre a la FACE MCP, depuis l'injection, gain>=2. SigmaT : SD entre avalanches.\n"
        << "Rayon : sqrt(<(x-<x>)^2+(y-<y>)^2>) a +50 um, pour chaque paquet avec >=2 impacts.\n"
        << "Poids egal par avalanche ; tous les evenements limites sont inclus. Barres des moyennes : SEM.\n"
        << "nan = indefini. SigmaT n'est pas la largeur temporelle interne d'un paquet.\n";
}

static void DrawPanel(
    TCanvas& canvas,
    const std::vector<Point>& points,
    int pad,
    Stats Point::* quantity,
    bool showSD,
    const char* title)
{
    canvas.cd(pad);
    gPad->SetLeftMargin(0.15);

    auto* graph = new TGraphErrors;
    double ymax = 0.0;

    for (const Point& point : points) {
        const Stats& stats = point.*quantity;
        const double y = showSD ? stats.SD() : stats.Mean();

        if (!std::isfinite(y)) {
            continue;
        }

        const double error =
            (!showSD && stats.n > 1)
                ? stats.SEM()
                : 0.0;

        const int n = graph->GetN();

        graph->SetPoint(
            n,
            point.voltage / 1000.0,
            y
        );

        graph->SetPointError(
            n,
            0.0,
            error
        );

        ymax = std::max(ymax, y + error);
    }

    const double xmin = points.front().voltage / 1000.0;
    const double xmax = points.back().voltage / 1000.0;

    const double margin =
        std::max(
            0.05,
            0.05 * (xmax - xmin)
        );

    gPad->DrawFrame(
        xmin - margin,
        0.0,
        xmax + margin,
        ymax > 0.0 ? 1.2 * ymax : 1.0,
        title
    );

    graph->SetMarkerStyle(20);
    graph->SetMarkerColor(kBlue + 1);
    graph->SetLineColor(kBlue + 1);

    if (graph->GetN()) {
        graph->Draw("PL SAME");
    }
}

static void DrawPlots(
    const std::vector<Point>& points,
    const std::string& figure)
{
    TCanvas canvas(
        "voltage_analysis",
        "Temps et dispersion selon la tension",
        1100,
        1100
    );

    canvas.Divide(2, 3);

    DrawPanel(
        canvas,
        points,
        1,
        &Point::transit,
        false,
        "Temps de transit moyen dans le MCP;Tension (kV);Temps (ps)"
    );

    DrawPanel(
        canvas,
        points,
        2,
        &Point::transit,
        true,
        "Dispersion des temps entre avalanches;Tension (kV);sigma_T (ps)"
    );

    DrawPanel(
        canvas,
        points,
        3,
        &Point::radial,
        false,
        "Largeur radiale moyenne a +50 um;Tension (kV);Moyenne de sigma_r (um)"
    );

    DrawPanel(
        canvas,
        points,
        4,
        &Point::radial,
        true,
        "Dispersion des largeurs entre avalanches;Tension (kV);SD(sigma_r) (um)"
    );

    DrawPanel(
        canvas,
        points,
        5,
        &Point::gain,
        false,
        "Gain moyen a la face MCP;Tension (kV);Gain moyen"
    );

    DrawPanel(
        canvas,
        points,
        6,
        &Point::gain,
        true,
        "Dispersion du gain entre injections;Tension (kV);SD(gain)"
    );

    canvas.SaveAs((figure + ".png").c_str());
    canvas.SaveAs((figure + ".pdf").c_str());
}

} // namespace VoltageScan

// Usage :
// root -l -b -q 'analysis/scan_voltage.C("voltage_*.root")'
void scan_voltage(
    const char* motif = "voltage_*200events.root",
    const char* figure = "scan_voltage")
{
    using namespace VoltageScan;

    const std::vector<Point> points = LoadScan(motif);

    WriteCSV(points, figure);
    PrintSummary(points);
    DrawPlots(points, figure);
}
