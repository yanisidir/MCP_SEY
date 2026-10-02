#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TLine.h>
#include <TRegexp.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TVirtualPad.h>
#include <TTree.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

// Un point du balayage : un fichier ROOT, donc un plafond de generations.
struct PointBalayage {
    int plafond = 0;
    long evenements = 0, avecSortie = 0, utilises = 0, limites = 0, limGen = 0;
    double gain = 0, generationAtteinte = 0;
    double sigmaPremier = 0, sigmaMoyen = 0;
    double premier = 0, dPremier = 0, moyen = 0, dMoyen = 0, largeur = 0, dLargeur = 0;
};

// Moyenne et erreur standard sur la moyenne, a partir des sommes accumulees.
// L'erreur suppose des evenements independants.
static void Resume(double somme, double somme2, long n, double& moyenne, double& erreur)
{
    if (n <= 0) { moyenne = erreur = 0; return; }
    
    moyenne = somme / n;
    const double variance = std::max(0.0, somme2/n - moyenne*moyenne);
    erreur = n > 1 ? std::sqrt(variance/(n-1)) : 0;
}

static PointBalayage LirePoint(const std::string& chemin)
{
    TFile fichier(chemin.c_str());
    if (fichier.IsZombie()) throw std::runtime_error("Fichier illisible : " + chemin);

    auto* configuration = fichier.Get<TTree>("configuration");
    auto* evenements = fichier.Get<TTree>("events");

    if (!configuration || !evenements) throw std::runtime_error("Tables absentes : " + chemin);

    PointBalayage point;
    configuration->SetBranchAddress("MaxGenerations", &point.plafond);
    configuration->GetEntry(0);

    if (!evenements->GetBranch("IsMultiplicationLimited"))
        throw std::runtime_error("Colonne IsMultiplicationLimited absente : " + chemin);
    int limited = 0;
    evenements->SetBranchAddress("IsMultiplicationLimited", &limited);
    int gain = 0, generation = 0;
    double premier = 0, moyen = 0, largeur = 0;

    evenements->SetBranchAddress("Gain", &gain);
    evenements->SetBranchAddress("MaxGeneration", &generation);
    evenements->SetBranchAddress("FirstTime_ns", &premier);
    evenements->SetBranchAddress("MeanTime_ns", &moyen);
    evenements->SetBranchAddress("TimeSpread_ns", &largeur);

    double sommeGain = 0, sommeGeneration = 0;
    double s1 = 0, s1b = 0, s2 = 0, s2b = 0, s3 = 0, s3b = 0;
    long acheves = 0;

    for (long i = 0; i < evenements->GetEntries(); ++i) {
        evenements->GetEntry(i);
        if (limited) ++point.limites;
        if (point.plafond > 0 && generation >= point.plafond) ++point.limGen;
        ++point.evenements;
        if (gain >= 1) ++point.avecSortie;
        ++acheves;
        sommeGain += gain;
        sommeGeneration += generation;
        // Temps conditionnes a au moins deux electrons en sortie : en deca, un
        // barycentre temporel n'est pas defini.
        if (gain < 2 || !std::isfinite(premier) || !std::isfinite(moyen)
            || !std::isfinite(largeur)) continue;
        ++point.utilises;
        s1 += premier*1000;  s1b += premier*premier*1e6;
        s2 += moyen*1000;    s2b += moyen*moyen*1e6;
        s3 += largeur*1000;  s3b += largeur*largeur*1e6;
    }

    if (acheves) { point.gain = sommeGain/acheves; point.generationAtteinte = sommeGeneration/acheves; }
    Resume(s1, s1b, point.utilises, point.premier, point.dPremier);
    Resume(s2, s2b, point.utilises, point.moyen, point.dMoyen);
    Resume(s3, s3b, point.utilises, point.largeur, point.dLargeur);
    // Ecart-type ENTRE evenements (denominateur N-1), a ne pas confondre avec
    // l'erreur sur la moyenne que Resume renvoie (SEM = s / sqrt(N)).
    const double nan = std::numeric_limits<double>::quiet_NaN();
    point.sigmaPremier = point.utilises > 1 ? point.dPremier * std::sqrt(point.utilises) : nan;
    point.sigmaMoyen = point.utilises > 1 ? point.dMoyen * std::sqrt(point.utilises) : nan;
    return point;
}

// Depuis la racine : root -l -b -q 'analysis/scan_generations.C("scan_gen*.root")'
// Cree le dossier de sortie s'il manque : une copie fraiche du depot n'a pas de
// dossier figures/, et ROOT echoue alors sans ecrire ni figure ni CSV.
#ifndef MCP_PREPARER_SORTIE
#define MCP_PREPARER_SORTIE
static void PreparerSortie(const char* chemin)
{
    const TString dossier = gSystem->DirName(chemin);
    if (dossier != "." && dossier != "" && gSystem->AccessPathName(dossier))
        gSystem->mkdir(dossier, kTRUE);
}
#endif

void scan_generations(const char* motif = "scan_gen*.root", const char* figure = "scan_generations")
{
    PreparerSortie(figure);
    const TString dossier = gSystem->DirName(motif);
    const TString patron = gSystem->BaseName(motif);
    TSystemDirectory repertoire(dossier, dossier);
    auto* contenu = repertoire.GetListOfFiles();
    if (!contenu) throw std::runtime_error("Repertoire introuvable : " + std::string(dossier.Data()));

    std::vector<PointBalayage> points;
    PointBalayage libre;
    bool aLibre = false;
    TRegexp filtre(patron, kTRUE);
    TIter suivant(contenu);
    bool found = false;
    while (auto* entree = static_cast<TSystemFile*>(suivant())) {
        const TString nom = entree->GetName();
        Ssiz_t longueur = 0;
        if (entree->IsDirectory() || nom.Index(filtre, &longueur) != 0 || longueur != nom.Length()) continue;
        const auto point = LirePoint(std::string(dossier.Data()) + "/" + nom.Data());
        found = true;
        if (!point.evenements) continue;
        if (point.plafond == 0) { libre = point; aLibre = true; }
        else points.push_back(point);
    }
    if (!found) throw std::runtime_error("Aucun fichier ne correspond a " + std::string(motif));
    if (points.empty() && !aLibre) {
        std::cout << "Aucun evenement dans les fichiers.\n";
        return;
    }
    std::sort(points.begin(), points.end(),
              [](const PointBalayage& a, const PointBalayage& b) { return a.plafond < b.plafond; });

    printf("%8s %7s %9s %9s %9s %9s %10s %10s %10s %10s %10s\n", "plafond", "evts", "limites", "lim_gen", "sortie", "temps",
           "gain", "gen.max", "premier", "barycentre", "largeur");
    // Un tiret plutot qu'un zero la ou aucun evenement ne renseigne la colonne.
    auto ligne = [](const PointBalayage& p, const char* etiquette) {
        const bool acheve = p.evenements > 0;
        auto mesure = [](bool defini, double v, double d) {
            static char texte[4][24];
            static int tour = 0;
            tour = (tour + 1) % 4;
            if (defini) snprintf(texte[tour], sizeof(texte[0]), "%6.1f+-%-4.1f", v, d);
            else snprintf(texte[tour], sizeof(texte[0]), "%11s", "-");
            return texte[tour];
        };
        printf("%8s %7ld %9ld %9ld %9ld %9ld %10s %10s %s %s %s\n", etiquette, p.evenements, p.limites, p.limGen,
               p.avecSortie, p.utilises,
               acheve ? Form("%.1f", p.gain) : "-", acheve ? Form("%.1f", p.generationAtteinte) : "-",
               mesure(p.utilises, p.premier, p.dPremier), mesure(p.utilises, p.moyen, p.dMoyen),
               mesure(p.utilises, p.largeur, p.dLargeur));
    };
    for (const auto& point : points) ligne(point, Form("%d", point.plafond));
    if (aLibre) ligne(libre, "libre");
    std::cout << "Temps en ps, sur tous les evenements avec gain>=2. Limites : plafond de creation atteint. Lim_gen : plafond de generations atteint (0 si sans plafond).\n";

    for (const auto& point : points) {
        if (!point.utilises)
            printf("Plafond %d : aucun evenement avec gain>=2, temps indefini.\n", point.plafond);
        else if (point.utilises < 20)
            printf("Plafond %d : seulement %ld evenements dans les statistiques de temps.\n",
                   point.plafond, point.utilises);
    }
    if (points.empty()) {
        std::cout << "Un seul point sans limite de generations : rien a tracer contre le plafond.\n";
        return;
    }

    TCanvas canvas("scan_generations", "Balayage du plafond de generations", 1000, 750);
    canvas.Divide(2, 2);
    const int n = static_cast<int>(points.size());

    auto construire = [&](double (*valeur)(const PointBalayage&), double (*erreur)(const PointBalayage&)) {
        auto* graphe = new TGraphErrors();
        int j = 0;
        for (int i = 0; i < n; ++i) {
            if ((erreur && !points[i].utilises) || !std::isfinite(valeur(points[i]))) continue;
            graphe->SetPoint(j, points[i].plafond, valeur(points[i]));
            graphe->SetPointError(j++, 0, erreur ? erreur(points[i]) : 0);
        }
        graphe->SetMarkerStyle(20);
        return graphe;
    };
    auto reference = [&](double v) {
        if (!aLibre || n == 0) return;
        auto* trait = new TLine(points.front().plafond, v, points.back().plafond, v);
        trait->SetLineStyle(2);
        trait->Draw();
    };

    canvas.cd(1);
    auto* sigmaPremiers = construire([](const PointBalayage& p) { return p.sigmaPremier; }, nullptr);
    auto* sigmaMoyens = construire([](const PointBalayage& p) { return p.sigmaMoyen; }, nullptr);
    sigmaMoyens->SetTitle("Dispersion entre avalanches;Plafond de generations;Ecart-type entre evenements (ps)");
    sigmaMoyens->SetLineColor(kBlue+1);
    sigmaMoyens->SetMarkerColor(kBlue+1);
    sigmaPremiers->SetLineColor(kRed+1);
    sigmaPremiers->SetMarkerColor(kRed+1);
    sigmaPremiers->SetMarkerStyle(24);
    sigmaPremiers->SetLineStyle(2);
    if (sigmaMoyens->GetN()) {
        double maximum = 0;
        for (const auto& point : points)
            if (point.utilises > 1) maximum = std::max({maximum, point.sigmaPremier, point.sigmaMoyen});
        sigmaMoyens->SetMinimum(0);
        sigmaMoyens->SetMaximum(maximum > 0 ? 1.25*maximum : 1);
        sigmaMoyens->Draw("APL");
        sigmaPremiers->Draw("PL same");
    } else {
        gPad->DrawFrame(points.front().plafond-1, 0, points.back().plafond+1, 1,
            "Dispersion indefinie (moins de 2 avalanches);Plafond de generations;Ecart-type entre evenements (ps)");
    }
    auto* legende = new TLegend(0.15, 0.75, 0.65, 0.88);
    legende->AddEntry(sigmaMoyens, "barycentre du paquet", "pl");
    legende->AddEntry(sigmaPremiers, "premiere arrivee", "pl");
    legende->Draw();

    canvas.cd(2);
    gPad->SetLogy();
    auto* gains = construire([](const PointBalayage& p) { return std::max(1e-3, p.gain); }, nullptr);
    gains->SetTitle("Gain des evenements acheves;Plafond de generations;Electrons sortants");
    gains->Draw("APL");
    if (aLibre && libre.gain > 0) reference(libre.gain);

    canvas.cd(3);
    auto* premiers = construire([](const PointBalayage& p) { return p.premier; },
                                [](const PointBalayage& p) { return p.dPremier; });
    auto* moyens = construire([](const PointBalayage& p) { return p.moyen; },
                              [](const PointBalayage& p) { return p.dMoyen; });
    moyens->SetTitle("Temps de transit;Plafond de generations;Temps (ps)");
    if (moyens->GetN()) moyens->Draw("APL");
    else gPad->DrawFrame(points.front().plafond-1, 0, points.back().plafond+1, 1, "Temps indefinis;Plafond de generations;Temps (ps)");
    premiers->SetMarkerStyle(24);
    premiers->SetLineStyle(2);
    if (premiers->GetN()) premiers->Draw("PL same");
    auto* legendeTemps = new TLegend(0.15, 0.80, 0.35, 0.90);
    legendeTemps->AddEntry(moyens, "barycentre du paquet", "pl");
    legendeTemps->AddEntry(premiers, "premiere arrivee", "pl");
    legendeTemps->Draw();

    canvas.cd(4);
    auto* largeurs = construire([](const PointBalayage& p) { return p.largeur; },
                                [](const PointBalayage& p) { return p.dLargeur; });
    largeurs->SetTitle("Largeur interne du paquet;Plafond de generations;Ecart-type (ps)");
    if (largeurs->GetN()) largeurs->Draw("APL");
    else gPad->DrawFrame(points.front().plafond-1, 0, points.back().plafond+1, 1, "Largeur indefinie;Plafond de generations;Ecart-type (ps)");
    if (aLibre && libre.utilises) reference(libre.largeur);

    canvas.SaveAs((std::string(figure)+".png").c_str());
    canvas.SaveAs((std::string(figure)+".pdf").c_str());
}
