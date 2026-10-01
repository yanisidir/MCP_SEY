#ifndef MCP_TRANSVERSE_ANALYSIS_C
#define MCP_TRANSVERSE_ANALYSIS_C

#include <TCanvas.h>
#include <TError.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLeaf.h>
#include <TLegend.h>
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


// ============================================================================
// 1. TYPES ET CONSTANTES
// ============================================================================

// Toutes les positions de cette analyse sont en um,
// les energies sont en eV.

using EventKey = std::pair<int, int>;

constexpr double kMapRadiusUm = 150.0;

// Fraction des electrons d'une avalanche contenue dans le rayon R80.
constexpr double kContainmentFraction = 0.80;


// Statistiques d'une avalanche sur un plan donne.
struct Spot {
    long n = 0;

    double x = 0.0;
    double y = 0.0;

    double m2x = 0.0;
    double m2y = 0.0;

    double energy = 0.0;

    // Positions conservees pour le rayon de confinement, qui exige
    // le centroide final : il n'a pas de forme incrementale.
    std::vector<std::pair<double, double>> points;


    void Add(
        double px,
        double py,
        double e)
    {
        ++n;
        points.emplace_back(px, py);

        const double dx = px - x;
        const double dy = py - y;

        x += dx / n;
        y += dy / n;

        m2x += dx * (px - x);
        m2y += dy * (py - y);

        energy += (e - energy) / n;
    }


    double Sx() const
    {
        return n
            ? std::sqrt(
                std::max(
                    0.0,
                    m2x / n
                )
            )
            : 0.0;
    }


    double Sy() const
    {
        return n
            ? std::sqrt(
                std::max(
                    0.0,
                    m2y / n
                )
            )
            : 0.0;
    }


    double Sr() const
    {
        return std::hypot(
            Sx(),
            Sy()
        );
    }


    // Rayon du cercle centre sur le centroide qui contient la fraction
    // demandee des electrons : distance du k-ieme plus proche, k = ceil(f*n).
    // Insensible a la queue de la distribution, contrairement a sigma_r.
    double ContainmentRadius(
        double fraction) const
    {
        if (points.empty()) {
            return 0.0;
        }

        std::vector<double> distances;
        distances.reserve(points.size());

        for (const auto& p : points) {
            distances.push_back(
                std::hypot(
                    p.first - x,
                    p.second - y
                )
            );
        }

        const long k =
            std::clamp<long>(
                static_cast<long>(
                    std::ceil(fraction * distances.size())
                ),
                1,
                static_cast<long>(distances.size())
            );

        std::nth_element(
            distances.begin(),
            distances.begin() + (k - 1),
            distances.end()
        );

        return distances[k - 1];
    }


    double R80() const
    {
        return ContainmentRadius(kContainmentFraction);
    }
};


// ============================================================================
// 2. OUTILS DE LECTURE ROOT
// ============================================================================

static double ReadValue(
    TTree* tree,
    const char* name)
{
    TLeaf* leaf = tree->GetLeaf(name);

    if (!leaf) {
        throw std::runtime_error(
            std::string("Colonne absente : ") + name
        );
    }

    return leaf->GetValue();
}


static EventKey ReadEventKey(
    TTree* tree)
{
    return {
        static_cast<int>(
            ReadValue(tree, "RunID")
        ),
        static_cast<int>(
            ReadValue(tree, "EventID")
        )
    };
}


// ============================================================================
// 3. LECTURE ET AGRÉGATION DES AVALANCHES
// ============================================================================

static std::map<EventKey, std::array<Spot, 2>>
BuildSpots(
    TTree* exit,
    TTree* downstream)
{
    std::map<EventKey, std::array<Spot, 2>> spots;

    const std::array<TTree*, 2> tables = {
        exit,
        downstream
    };

    // plane = 0 : face MCP
    // plane = 1 : plan a +50 um
    for (int plane = 0; plane < 2; ++plane) {
        TTree* tree = tables[plane];

        for (Long64_t i = 0;
             i < tree->GetEntries();
             ++i) {

            tree->GetEntry(i);

            const EventKey key =
                ReadEventKey(tree);

            const double x_um =
                ReadValue(tree, "X_mm") * 1000.0;

            const double y_um =
                ReadValue(tree, "Y_mm") * 1000.0;

            const double energy_eV =
                ReadValue(
                    tree,
                    "KineticEnergy_eV"
                );

            spots[key][plane].Add(
                x_um,
                y_um,
                energy_eV
            );
        }
    }

    return spots;
}


// ============================================================================
// 4. CALCUL DES BORNES ET SELECTION DES AVALANCHES
// ============================================================================

struct PlotRanges {
    long selected = 0;

    double exitRadius = 1.0;
    double energyMax = 1.0;
    double sigmaMax = 1.0;
    double expansionMax = 1.0;

    double xmin = 1e30;
    double xmax = -1e30;
    double ymin = 1e30;
    double ymax = -1e30;

    long outsideMap = 0;
};


static PlotRanges ComputeRanges(
    const std::map<EventKey, std::array<Spot, 2>>& spots,
    TTree* exit,
    TTree* downstream)
{
    PlotRanges ranges;

    // ------------------------------------------------------------------------
    // 4.1 Selection des avalanches et bornes sur les barycentres / largeurs
    // ------------------------------------------------------------------------

    for (const auto& entry : spots) {
        const Spot& faceMCP = entry.second[0];
        const Spot& plane50 = entry.second[1];

        // On ne garde que les avalanches ayant au moins
        // 2 electrons dans chacun des deux plans.
        if (faceMCP.n < 2 || plane50.n < 2) {
            continue;
        }

        ++ranges.selected;

        ranges.sigmaMax =
            std::max(
                {
                    ranges.sigmaMax,
                    faceMCP.Sr(),
                    plane50.Sr(),
                    faceMCP.R80(),
                    plane50.R80()
                }
            );

        ranges.expansionMax =
            std::max(
                ranges.expansionMax,
                std::abs(
                    plane50.Sr() -
                    faceMCP.Sr()
                )
            );

        ranges.xmin =
            std::min(
                ranges.xmin,
                plane50.x
            );

        ranges.xmax =
            std::max(
                ranges.xmax,
                plane50.x
            );

        ranges.ymin =
            std::min(
                ranges.ymin,
                plane50.y
            );

        ranges.ymax =
            std::max(
                ranges.ymax,
                plane50.y
            );
    }

    if (!ranges.selected) {
        throw std::runtime_error(
            "Aucune avalanche avec au moins "
            "deux electrons aux deux plans."
        );
    }


    // ------------------------------------------------------------------------
    // 4.2 Deuxieme lecture des electrons :
    //     - rayon necessaire pour la carte de sortie
    //     - energie maximale
    //     - nombre d'impacts hors du zoom a +50 um
    // ------------------------------------------------------------------------

    const std::array<TTree*, 2> tables = {
        exit,
        downstream
    };

    for (int plane = 0; plane < 2; ++plane) {
        TTree* tree = tables[plane];

        for (Long64_t i = 0;
             i < tree->GetEntries();
             ++i) {

            tree->GetEntry(i);

            const EventKey key =
                ReadEventKey(tree);

            const auto& pair =
                spots.at(key);

            if (pair[0].n < 2 ||
                pair[1].n < 2) {
                continue;
            }

            const double x_um =
                ReadValue(tree, "X_mm") * 1000.0;

            const double y_um =
                ReadValue(tree, "Y_mm") * 1000.0;

            const double dx =
                std::abs(
                    x_um -
                    pair[plane].x
                );

            const double dy =
                std::abs(
                    y_um -
                    pair[plane].y
                );

            if (plane == 0) {
                ranges.exitRadius =
                    std::max(
                        {
                            ranges.exitRadius,
                            dx,
                            dy
                        }
                    );
            }
            else {
                if (dx > kMapRadiusUm ||
                    dy > kMapRadiusUm) {

                    ++ranges.outsideMap;
                }
            }

            ranges.energyMax =
                std::max(
                    ranges.energyMax,
                    ReadValue(
                        tree,
                        "KineticEnergy_eV"
                    )
                );
        }
    }


    // Petite marge visuelle.
    ranges.exitRadius   *= 1.05;
    ranges.energyMax    *= 1.05;
    ranges.sigmaMax     *= 1.05;
    ranges.expansionMax *= 1.05;

    return ranges;
}


// ============================================================================
// 5. EXPORT CSV
// ============================================================================

static void WriteCSV(
    const std::string& figure,
    const std::map<EventKey, std::array<Spot, 2>>& spots)
{
    std::ofstream csv(
        figure + ".csv"
    );

    if (!csv) {
        throw std::runtime_error(
            "Impossible d'ecrire le CSV."
        );
    }

    csv << std::setprecision(10);

    csv
        << "RunID,"
        << "EventID,"
        << "NExit,"
        << "N50,"
        << "MeanXExit_um,"
        << "MeanYExit_um,"
        << "SigmaXExit_um,"
        << "SigmaYExit_um,"
        << "SigmaRExit_um,"
        << "R80Exit_um,"
        << "MeanEnergyExit_eV,"
        << "MeanX50_um,"
        << "MeanY50_um,"
        << "SigmaX50_um,"
        << "SigmaY50_um,"
        << "SigmaR50_um,"
        << "R8050_um,"
        << "MeanEnergy50_eV\n";

    for (const auto& entry : spots) {
        const Spot& faceMCP = entry.second[0];
        const Spot& plane50 = entry.second[1];

        if (faceMCP.n < 2 ||
            plane50.n < 2) {
            continue;
        }

        csv
            << entry.first.first
            << ','
            << entry.first.second
            << ','
            << faceMCP.n
            << ','
            << plane50.n;

        for (const Spot& spot : entry.second) {
            csv
                << ',' << spot.x
                << ',' << spot.y
                << ',' << spot.Sx()
                << ',' << spot.Sy()
                << ',' << spot.Sr()
                << ',' << spot.R80()
                << ',' << spot.energy;
        }

        csv << '\n';
    }
}


// ============================================================================
// 6. HISTOGRAMMES
// ============================================================================

struct Histograms {
    TH2D xy0;
    TH2D xy50;
    TH2D bary;

    TH1D rms0;
    TH1D rms50;
    TH1D r80_0;
    TH1D r80_50;

    TH1D e0;
    TH1D e50;

    TH1D growth;


    Histograms(
        const PlotRanges& ranges)
        :
        xy0(
            "xy0",
            "Face MCP : impacts recentres;"
            "x-<x> (um);"
            "y-<y> (um)",
            100,
            -ranges.exitRadius,
            ranges.exitRadius,
            100,
            -ranges.exitRadius,
            ranges.exitRadius
        ),

        xy50(
            "xy50",
            "A +50 um : zoom recentre +/-150 um;"
            "x-<x> (um);"
            "y-<y> (um)",
            100,
            -kMapRadiusUm,
            kMapRadiusUm,
            100,
            -kMapRadiusUm,
            kMapRadiusUm
        ),

        bary(
            "bary",
            "Barycentres a +50 um;"
            "x (um);"
            "y (um)",
            80,
            ranges.xmin - 1.0,
            ranges.xmax + 1.0,
            80,
            ranges.ymin - 1.0,
            ranges.ymax + 1.0
        ),

        rms0(
            "rms0",
            "Largeur radiale par avalanche;"
            "sigma_r (um);"
            "Avalanches",
            60,
            0.0,
            ranges.sigmaMax
        ),

        rms50(
            "rms50",
            "",
            60,
            0.0,
            ranges.sigmaMax
        ),

        r80_0(
            "r80_0",
            "",
            60,
            0.0,
            ranges.sigmaMax
        ),

        r80_50(
            "r80_50",
            "",
            60,
            0.0,
            ranges.sigmaMax
        ),

        e0(
            "e0",
            "Energie cinetique;"
            "eV;"
            "Electrons",
            100,
            0.0,
            ranges.energyMax
        ),

        e50(
            "e50",
            "",
            100,
            0.0,
            ranges.energyMax
        ),

        growth(
            "growth",
            "Variation de largeur;"
            " sigma_r(+50 um)-sigma_r(sortie) (um);"
            "Avalanches",
            60,
            -ranges.expansionMax,
            ranges.expansionMax
        )
    {}
};


static void DisableStats(
    Histograms& h)
{
    for (TH1* hist : {
            static_cast<TH1*>(&h.xy0),
            static_cast<TH1*>(&h.xy50),
            static_cast<TH1*>(&h.bary),
            static_cast<TH1*>(&h.rms0),
            static_cast<TH1*>(&h.rms50),
            static_cast<TH1*>(&h.r80_0),
            static_cast<TH1*>(&h.r80_50),
            static_cast<TH1*>(&h.e0),
            static_cast<TH1*>(&h.e50),
            static_cast<TH1*>(&h.growth)
        }) {

        hist->SetStats(false);
    }
}


// ============================================================================
// 7. REMPLISSAGE DES HISTOGRAMMES
// ============================================================================

struct MeanWidths {
    double sumExit = 0.0;
    double sum50 = 0.0;
    double sumR80Exit = 0.0;
    double sumR8050 = 0.0;
};


static MeanWidths FillAvalancheHistograms(
    Histograms& h,
    const std::map<EventKey, std::array<Spot, 2>>& spots)
{
    MeanWidths sums;

    for (const auto& entry : spots) {
        const Spot& faceMCP = entry.second[0];
        const Spot& plane50 = entry.second[1];

        if (faceMCP.n < 2 ||
            plane50.n < 2) {
            continue;
        }

        h.rms0.Fill(
            faceMCP.Sr()
        );

        h.rms50.Fill(
            plane50.Sr()
        );

        h.r80_0.Fill(
            faceMCP.R80()
        );

        h.r80_50.Fill(
            plane50.R80()
        );

        h.growth.Fill(
            plane50.Sr() -
            faceMCP.Sr()
        );

        h.bary.Fill(
            plane50.x,
            plane50.y
        );

        sums.sumExit +=
            faceMCP.Sr();

        sums.sum50 +=
            plane50.Sr();

        sums.sumR80Exit +=
            faceMCP.R80();

        sums.sumR8050 +=
            plane50.R80();
    }

    return sums;
}


static void FillElectronHistograms(
    Histograms& h,
    const std::map<EventKey, std::array<Spot, 2>>& spots,
    TTree* exit,
    TTree* downstream)
{
    const std::array<TTree*, 2> tables = {
        exit,
        downstream
    };

    for (int plane = 0; plane < 2; ++plane) {
        TTree* tree = tables[plane];

        for (Long64_t i = 0;
             i < tree->GetEntries();
             ++i) {

            tree->GetEntry(i);

            const EventKey key =
                ReadEventKey(tree);

            const auto& pair =
                spots.at(key);

            if (pair[0].n < 2 ||
                pair[1].n < 2) {
                continue;
            }

            const double x_um =
                ReadValue(tree, "X_mm") * 1000.0;

            const double y_um =
                ReadValue(tree, "Y_mm") * 1000.0;

            const double dx =
                x_um -
                pair[plane].x;

            const double dy =
                y_um -
                pair[plane].y;

            const double energy_eV =
                ReadValue(
                    tree,
                    "KineticEnergy_eV"
                );

            if (plane == 0) {
                h.xy0.Fill(
                    dx,
                    dy
                );

                h.e0.Fill(
                    energy_eV
                );
            }
            else {
                h.xy50.Fill(
                    dx,
                    dy
                );

                h.e50.Fill(
                    energy_eV
                );
            }
        }
    }
}


// ============================================================================
// 8. RESUME CONSOLE
// ============================================================================

static void PrintSummary(
    TTree* exit,
    TTree* downstream,
    const PlotRanges& ranges,
    const MeanWidths& sums)
{
    const double n = ranges.selected;
    printf("electrons : sortie=%lld, +50um=%lld ; avalanches retenues=%ld\n",
           exit->GetEntries(), downstream->GetEntries(), ranges.selected);
    printf("%-8s %10s %10s\n", "um", "sortie", "+50um");
    printf("%-8s %10.4g %10.4g\n", "sigma_r", sums.sumExit / n, sums.sum50 / n);
    printf("%-8s %10.4g %10.4g\n", "R80", sums.sumR80Exit / n, sums.sumR8050 / n);
}


// ============================================================================
// 9. AFFICHAGE
// ============================================================================

static void DrawCanvas(
    Histograms& h,
    const std::string& figure)
{
    TCanvas canvas(
        "transverse_canvas",
        "Dispersion transverse",
        1400,
        850
    );

    canvas.Divide(
        3,
        2
    );


    // ------------------------------------------------------------------------
    // 9.1 Mise en page commune
    // ------------------------------------------------------------------------

    for (int pad = 1;
         pad <= 6;
         ++pad) {

        canvas.cd(pad);

        gPad->SetLeftMargin(0.14);
        gPad->SetRightMargin(0.14);
    }


    // ------------------------------------------------------------------------
    // 9.2 Carte des impacts a la face MCP
    // ------------------------------------------------------------------------

    canvas.cd(1);

    h.xy0.Draw(
        "COLZ"
    );


    // ------------------------------------------------------------------------
    // 9.3 Carte des impacts a +50 um
    // ------------------------------------------------------------------------

    canvas.cd(2);

    gPad->SetLogz();

    h.xy50.Draw(
        "COLZ"
    );


    // ------------------------------------------------------------------------
    // 9.4 Distribution des largeurs radiales
    // ------------------------------------------------------------------------

    canvas.cd(3);

    h.rms0.SetLineColor(
        kBlue
    );

    h.rms50.SetLineColor(
        kRed
    );

    h.r80_0.SetLineColor(
        kBlue
    );

    h.r80_50.SetLineColor(
        kRed
    );

    h.r80_0.SetLineStyle(
        2
    );

    h.r80_50.SetLineStyle(
        2
    );

    h.rms0.SetTitle(
        "Largeur radiale par avalanche;"
        "sigma_r (trait plein), R80 (tirets) (um);"
        "Avalanches"
    );

    h.rms0.SetMaximum(
        1.15 *
        std::max(
            {
                h.rms0.GetMaximum(),
                h.rms50.GetMaximum(),
                h.r80_0.GetMaximum(),
                h.r80_50.GetMaximum()
            }
        )
    );

    h.rms0.Draw();
    h.rms50.Draw("SAME");
    h.r80_0.Draw("SAME");
    h.r80_50.Draw("SAME");

    TLegend widthLegend(
        0.55,
        0.72,
        0.88,
        0.88
    );

    widthLegend.AddEntry(
        &h.rms0,
        "Face MCP",
        "l"
    );

    widthLegend.AddEntry(
        &h.rms50,
        "+50 um",
        "l"
    );

    widthLegend.AddEntry(
        &h.r80_0,
        "R80 face",
        "l"
    );

    widthLegend.AddEntry(
        &h.r80_50,
        "R80 +50 um",
        "l"
    );

    widthLegend.Draw();


    // ------------------------------------------------------------------------
    // 9.5 Carte des barycentres a +50 um
    // ------------------------------------------------------------------------

    canvas.cd(4);

    h.bary.Draw(
        "COLZ"
    );


    // ------------------------------------------------------------------------
    // 9.6 Spectres d'energie
    // ------------------------------------------------------------------------

    canvas.cd(5);

    h.e0.SetLineColor(
        kBlue
    );

    h.e50.SetLineColor(
        kRed
    );

    h.e50.SetLineStyle(
        2
    );

    h.e0.SetMaximum(
        1.15 *
        std::max(
            h.e0.GetMaximum(),
            h.e50.GetMaximum()
        )
    );

    h.e0.Draw();
    h.e50.Draw("SAME");

    TLegend energyLegend(
        0.55,
        0.72,
        0.88,
        0.88
    );

    energyLegend.AddEntry(
        &h.e0,
        "Face MCP",
        "l"
    );

    energyLegend.AddEntry(
        &h.e50,
        "+50 um",
        "l"
    );

    energyLegend.Draw();


    // ------------------------------------------------------------------------
    // 9.7 Variation de largeur
    // ------------------------------------------------------------------------

    canvas.cd(6);

    h.growth.Draw();


    // ------------------------------------------------------------------------
    // 9.8 Sauvegarde
    // ------------------------------------------------------------------------

    canvas.SaveAs(
        (figure + ".png").c_str()
    );

    canvas.SaveAs(
        (figure + ".pdf").c_str()
    );
}


// ============================================================================
// 10. POINT D'ENTREE
// ============================================================================
//
// Exemple :
// root -l -b -q 'analysis/transverse.C("transverse_gen50.root")'
//
void transverse(
    const char* fichier = "transverse_gen50.root",
    const char* figure = "transverse")
{
    gErrorIgnoreLevel = kWarning;

    // ------------------------------------------------------------------------
    // 10.1 Ouvrir le fichier ROOT
    // ------------------------------------------------------------------------

    TFile file(
        fichier
    );

    auto* events = file.Get<TTree>("events");
    if (!events || !events->GetLeaf("MaxGeneration")
        || !events->GetLeaf("IsMultiplicationLimited"))
        throw std::runtime_error("Table events, MaxGeneration ou IsMultiplicationLimited absente");
    auto* configuration = file.Get<TTree>("configuration");
    if (!configuration || configuration->GetEntries() != 1
        || !configuration->GetLeaf("MaxGenerations"))
        throw std::runtime_error("Configuration ou MaxGenerations absente : un run attendu");
    configuration->GetEntry(0);
    const int plafond =
        static_cast<int>(configuration->GetLeaf("MaxGenerations")->GetValue());
    double meanGeneration = 0;
    long limites = 0, limGen = 0;
    for (Long64_t i=0; i<events->GetEntries(); ++i) {
        events->GetEntry(i);
        const double generation = events->GetLeaf("MaxGeneration")->GetValue();
        meanGeneration += generation;
        if (events->GetLeaf("IsMultiplicationLimited")->GetValue() != 0) ++limites;
        if (plafond > 0 && generation >= plafond) ++limGen;
    }
    meanGeneration = events->GetEntries() ? meanGeneration/events->GetEntries()
        : std::numeric_limits<double>::quiet_NaN();
    printf("%10s %10s %10s %10s\n%10lld %10ld %10ld %10.3f\n",
           "evts", "limites", "lim_gen", "gen.max",
           events->GetEntries(), limites, limGen, meanGeneration);
    std::cout << "limites : evenements ayant atteint le plafond de creation de tracks (0 si sans plafond).\n"
              << "lim_gen : evenements ayant atteint le plafond de generations (plafond="
              << plafond << " ; 0 = sans plafond).\n"
              << "gen.max : moyenne de MaxGeneration sur toutes les injections.\n";
    TTree* exit =
        file.Get<TTree>(
            "exits"
        );

    TTree* downstream =
        file.Get<TTree>(
            "downstream"
        );

    if (!exit ||
        !downstream ||
        !exit->GetEntries() ||
        !downstream->GetEntries()) {

        throw std::runtime_error(
            "Mesures absentes : activer "
            "/mcp/setTransitOutput true "
            "et obtenir des sorties."
        );
    }


    // ------------------------------------------------------------------------
    // 10.2 Construire les statistiques par avalanche
    // ------------------------------------------------------------------------

    const auto spots =
        BuildSpots(
            exit,
            downstream
        );


    // ------------------------------------------------------------------------
    // 10.3 Calculer les bornes necessaires aux histogrammes
    // ------------------------------------------------------------------------

    const PlotRanges ranges =
        ComputeRanges(
            spots,
            exit,
            downstream
        );


    // ------------------------------------------------------------------------
    // 10.4 Creer les histogrammes
    // ------------------------------------------------------------------------

    Histograms histograms(
        ranges
    );

    DisableStats(
        histograms
    );


    // ------------------------------------------------------------------------
    // 10.5 Exporter les statistiques par avalanche
    // ------------------------------------------------------------------------

    WriteCSV(
        figure,
        spots
    );


    // ------------------------------------------------------------------------
    // 10.6 Remplir les histogrammes
    // ------------------------------------------------------------------------

    const MeanWidths meanWidths =
        FillAvalancheHistograms(
            histograms,
            spots
        );

    FillElectronHistograms(
        histograms,
        spots,
        exit,
        downstream
    );


    // ------------------------------------------------------------------------
    // 10.7 Afficher le resume
    // ------------------------------------------------------------------------

    PrintSummary(
        exit,
        downstream,
        ranges,
        meanWidths
    );


    // ------------------------------------------------------------------------
    // 10.8 Dessiner et sauvegarder les figures
    // ------------------------------------------------------------------------

    DrawCanvas(
        histograms,
        figure
    );
}


#endif
