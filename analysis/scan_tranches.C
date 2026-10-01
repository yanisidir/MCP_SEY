#include <TCanvas.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TLeaf.h>
#include <TLegend.h>
#include <TRegexp.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TTree.h>
#include <TVirtualPad.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

// Une tranche = un fichier tranche_<i>.root, identifie par sa source.
struct Tranche {
    int index = 0;
    int plafond = 0;                                  // MaxGenerations du run ; 0 = sans plafond
    double muE = 0, sigmaE = 0, muZ = 0, sigmaZ = 0, muPz = 0, sigmaPz = 0;
    long evenements = 0, gainNul = 0, limites = 0, limGen = 0, temps = 0;
    double gain = 0, dGain = 0;                       // sur tous les evenements, gain nul compris
    double transit = 0, dTransit = 0;                 // barycentre moyen
    double sigmaBary = 0, dSigmaBary = 0, sigmaBaryRobuste = 0;   // SD entre evenements du barycentre
    double sigmaPremier = 0, dSigmaPremier = 0, sigmaPremierRobuste = 0; // idem, premiere arrivee
    double largeur = 0;                               // largeur interne moyenne d'un paquet
};

// Moyenne, ecart-type d'echantillon (N-1) et leurs erreurs standard. Les
// erreurs supposent des tirages independants et N suffisamment grand.
struct Echantillon {
    std::vector<double> v;
    long N() const { return static_cast<long>(v.size()); }
    double Moyenne() const {
        if (v.empty()) return std::nan("");
        double s = 0; for (double x : v) s += x; return s / v.size();
    }
    double SD() const {
        if (v.size() < 2) return std::nan("");
        const double m = Moyenne(); double s = 0;
        for (double x : v) s += (x - m) * (x - m);
        return std::sqrt(s / (v.size() - 1));
    }
    double ErreurMoyenne() const { return v.size() > 1 ? SD() / std::sqrt(v.size()) : std::nan(""); }
    double ErreurSD() const { return v.size() > 2 ? SD() / std::sqrt(2.0 * (v.size() - 1)) : std::nan(""); }
    // Demi-largeur inter-quantiles 16-84 : robuste aux queues et independante
    // de tout binning. Egale a sigma pour une distribution gaussienne seulement.
    double SDRobuste() const {
        if (v.size() < 10) return std::nan("");
        auto t = v; std::sort(t.begin(), t.end());
        const auto q = [&](double f) { return t[static_cast<size_t>(f * (t.size() - 1))]; };
        return (q(0.84) - q(0.16)) / 2;
    }
};

static double Valeur(TTree* t, const char* nom)
{
    auto* feuille = t->GetLeaf(nom);
    if (!feuille) throw std::runtime_error(std::string("Colonne absente : ") + nom);
    return feuille->GetValue();
}

static Tranche LireTranche(const std::string& chemin, int index)
{
    TFile fichier(chemin.c_str());
    if (fichier.IsZombie()) throw std::runtime_error("Fichier illisible : " + chemin);
    auto* configuration = fichier.Get<TTree>("configuration");
    auto* evenements = fichier.Get<TTree>("events");
    if (!configuration || !evenements) throw std::runtime_error("Tables absentes : " + chemin);
    if (configuration->GetEntries() != 1) throw std::runtime_error("Un seul run attendu par fichier : " + chemin);

    Tranche tranche;
    tranche.index = index;
    configuration->GetEntry(0);
    tranche.muE = Valeur(configuration, "SourceMuE_keV");
    tranche.sigmaE = Valeur(configuration, "SourceSigmaE_keV");
    tranche.muZ = Valeur(configuration, "SourceMuZ_mm");
    tranche.sigmaZ = Valeur(configuration, "SourceSigmaZ_mm");
    tranche.muPz = Valeur(configuration, "SourceMuPz_keV");
    tranche.sigmaPz = Valeur(configuration, "SourceSigmaPz_keV");
    tranche.plafond = static_cast<int>(Valeur(configuration, "MaxGenerations"));

    Echantillon gains, barycentres, premiers, largeurs;
    for (Long64_t i = 0; i < evenements->GetEntries(); ++i) {
        evenements->GetEntry(i);
        ++tranche.evenements;
        const auto gain = Valeur(evenements, "Gain");
        gains.v.push_back(gain);
        if (gain == 0) ++tranche.gainNul;
        if (Valeur(evenements, "IsMultiplicationLimited") != 0) ++tranche.limites;
        if (tranche.plafond > 0 && Valeur(evenements, "MaxGeneration") >= tranche.plafond) ++tranche.limGen;
        // Un barycentre temporel n'a de sens qu'avec au moins deux electrons en
        // sortie ; les gains 0 et 1 sont exclus des statistiques de temps.
        const auto bary = Valeur(evenements, "MeanTime_ns") * 1000;
        const auto premier = Valeur(evenements, "FirstTime_ns") * 1000;
        const auto largeur = Valeur(evenements, "TimeSpread_ns") * 1000;
        if (gain < 2 || !std::isfinite(bary) || !std::isfinite(premier) || !std::isfinite(largeur)) continue;
        ++tranche.temps;
        barycentres.v.push_back(bary);
        premiers.v.push_back(premier);
        largeurs.v.push_back(largeur);
    }
    tranche.gain = gains.Moyenne();          tranche.dGain = gains.ErreurMoyenne();
    tranche.transit = barycentres.Moyenne(); tranche.dTransit = barycentres.ErreurMoyenne();
    tranche.sigmaBary = barycentres.SD();      tranche.dSigmaBary = barycentres.ErreurSD();
    tranche.sigmaBaryRobuste = barycentres.SDRobuste();
    tranche.sigmaPremier = premiers.SD();      tranche.dSigmaPremier = premiers.ErreurSD();
    tranche.sigmaPremierRobuste = premiers.SDRobuste();
    tranche.largeur = largeurs.Moyenne();
    return tranche;
}

// Depuis la racine : root -l -b -q 'analysis/scan_tranches.C("tranche_*.root")'
void scan_tranches(const char* motif = "tranche_*.root", const char* figure = "scan_tranches")
{
    const TString dossier = gSystem->DirName(motif), patron = gSystem->BaseName(motif);
    TSystemDirectory repertoire(dossier, dossier);
    auto* contenu = repertoire.GetListOfFiles();
    if (!contenu) throw std::runtime_error("Repertoire introuvable : " + std::string(dossier.Data()));

    std::vector<Tranche> tranches;
    TRegexp filtre(patron, kTRUE), numero("[0-9]+");
    TIter suivant(contenu);
    while (auto* entree = static_cast<TSystemFile*>(suivant())) {
        const TString nom = entree->GetName();
        Ssiz_t longueur = 0;
        if (entree->IsDirectory() || nom.Index(filtre, &longueur) != 0 || longueur != nom.Length()) continue;
        Ssiz_t lnum = 0;
        const auto pos = nom.Index(numero, &lnum);
        const int index = pos == kNPOS ? static_cast<int>(tranches.size()) : std::stoi(nom(pos, lnum).Data());
        tranches.push_back(LireTranche(std::string(dossier.Data()) + "/" + nom.Data(), index));
    }
    if (tranches.empty()) throw std::runtime_error("Aucun fichier ne correspond a " + std::string(motif));
    std::sort(tranches.begin(), tranches.end(),
              [](const Tranche& a, const Tranche& b) { return a.muZ < b.muZ; });

    printf("%7s %8s %8s %8s %7s %7s %7s %7s %7s %12s %12s %12s %12s %8s\n", "tranche", "muZ(mm)", "E(keV)", "pz(keV)",
           "evts", "gain0", "limites", "lim_gen", "temps", "gain", "transit(ps)", "sigma bary", "sigma prem.", "largeur");
    for (const auto& t : tranches) {
        printf("%7d %8.4f %8.1f %8.1f %7ld %7ld %7ld %7ld %7ld %6.1f+-%-5.1f %6.1f+-%-5.1f %5.1f+-%-5.1f %5.1f+-%-5.1f %8.1f\n",
               t.index, t.muZ, t.muE, t.muPz, t.evenements, t.gainNul, t.limites, t.limGen, t.temps,
               t.gain, t.dGain, t.transit, t.dTransit, t.sigmaBary, t.dSigmaBary, t.sigmaPremier, t.dSigmaPremier, t.largeur);
    }
    std::cout << "limites : evenements ayant atteint le plafond de creation de tracks (0 si sans plafond).\n"
              << "lim_gen : evenements ayant atteint le plafond de generations (0 si sans plafond).\n"
              << "Gain : moyenne sur tous les evenements, gain nul compris ; erreur sur la moyenne.\n"
              << "Temps en ps depuis l'injection, evenements avec gain>=2. sigma : ecart-type (N-1) ENTRE evenements\n"
              << "du barycentre du paquet et de sa premiere arrivee, avec l'erreur sur l'ecart-type.\n"
              << "largeur : largeur interne moyenne d'un paquet, a ne pas confondre avec la dispersion\n"
              << "entre evenements ci-dessus.\n"
              << "Le CSV ajoute la demi-largeur 16-84 % de chaque sigma, robuste aux queues.\n"
              << "ATTENTION : cette macro ne donne que des ecarts-types. Le TTS est generalement\n"
              << "publie en largeur a mi-hauteur, qui ne vaut 2.355 x sigma que pour une\n"
              << "distribution gaussienne -- hypothese non verifiee ici. Pour une estimation en\n"
              << "largeur a mi-hauteur, voir gain_voltage.C.\n";
    for (const auto& t : tranches) {
        if (t.temps < 20) printf("Tranche %d : seulement %ld evenements dans les temps.\n", t.index, t.temps);
        if (t.limites > t.evenements / 2)
            printf("Tranche %d : %ld evenements sur %ld au plafond de creation, gain plafonne.\n",
                   t.index, t.limites, t.evenements);
        if (t.limGen > t.evenements / 2)
            printf("Tranche %d : %ld evenements sur %ld au plafond de generations, gain plafonne.\n",
                   t.index, t.limGen, t.evenements);
    }

    std::ofstream csv(std::string(figure) + ".csv");
    csv << std::setprecision(10)
        << "Tranche,MuZ_mm,SigmaZ_mm,MuE_keV,SigmaE_keV,MuPz_keV,SigmaPz_keV,Events,GainZero,Limited,GenerationReached,TimedEvents,"
           "Gain,Gain_SEM,Transit_ps,Transit_SEM,SigmaBary_ps,SigmaBary_Err,SigmaBary_Robust_ps,"
           "SigmaFirst_ps,SigmaFirst_Err,SigmaFirst_Robust_ps,InternalWidth_ps\n";
    for (const auto& t : tranches)
        csv << t.index << ',' << t.muZ << ',' << t.sigmaZ << ',' << t.muE << ',' << t.sigmaE << ',' << t.muPz << ','
            << t.sigmaPz << ',' << t.evenements << ',' << t.gainNul << ',' << t.limites << ',' << t.limGen << ','
            << t.temps << ','
            << t.gain << ',' << t.dGain << ',' << t.transit << ',' << t.dTransit << ',' << t.sigmaBary << ','
            << t.dSigmaBary << ',' << t.sigmaBaryRobuste << ',' << t.sigmaPremier << ',' << t.dSigmaPremier << ','
            << t.sigmaPremierRobuste << ',' << t.largeur << '\n';

    TCanvas canvas("scan_tranches", "Gain et dispersion par tranche", 1000, 750);
    canvas.Divide(2, 2);
    const int n = static_cast<int>(tranches.size());
    auto graphe = [&](double (*valeur)(const Tranche&), double (*erreur)(const Tranche&), int style) {
        auto* g = new TGraphErrors;
        for (const auto& t : tranches) {
            const double y = valeur(t);
            if (!std::isfinite(y)) continue;
            const int k = g->GetN();
            g->SetPoint(k, t.muZ, y);
            g->SetPointError(k, t.sigmaZ, erreur && std::isfinite(erreur(t)) ? erreur(t) : 0);
        }
        g->SetMarkerStyle(style);
        return g;
    };
    auto dessiner = [&](int pad, TGraphErrors* a, TGraphErrors* b, const char* titre,
                        const char* la, const char* lb, bool logy,
                        TGraphErrors* c = nullptr, const char* lc = "") {
        canvas.cd(pad);
        if (logy) gPad->SetLogy();
        // Cadre couvrant les deux courbes, barres comprises.
        double ymin = 1e300, ymax = -1e300;
        for (auto* g : {a, b, c}) if (g) for (int i = 0; i < g->GetN(); ++i) {
            ymin = std::min(ymin, g->GetY()[i] - g->GetEY()[i]);
            ymax = std::max(ymax, g->GetY()[i] + g->GetEY()[i]);
        }
        const double marge = logy ? 0 : std::max(1e-3, 0.3 * (ymax - ymin));
        a->SetTitle(titre);
        a->Draw(a->GetN() > 1 ? "APL" : "AP");
        a->GetHistogram()->SetMinimum(logy ? std::max(1e-3, ymin / 3) : std::max(0.0, ymin - marge));
        a->GetHistogram()->SetMaximum(logy ? ymax * 3 : ymax + marge);
        if (b) {
            b->SetLineStyle(2);
            b->Draw(b->GetN() > 1 ? "PL same" : "P same");
            if (c) {
                c->SetLineStyle(3);
                c->Draw(c->GetN() > 1 ? "PL same" : "P same");
            }
            auto* legende = new TLegend(0.15, c ? 0.70 : 0.75, 0.6, 0.88);
            legende->AddEntry(a, la, "pl");
            legende->AddEntry(b, lb, "pl");
            if (c) legende->AddEntry(c, lc, "pl");
            legende->Draw();
        }
    };
    dessiner(1, graphe([](const Tranche& t) { return std::max(1e-3, t.gain); },
                       [](const Tranche& t) { return t.dGain; }, 20),
             nullptr, "Gain moyen;Profondeur d'injection muZ (mm);Electrons sortants par injection", "", "", n > 0);
    dessiner(2, graphe([](const Tranche& t) { return 100.0 * t.gainNul / std::max(1L, t.evenements); }, nullptr, 20),
             graphe([](const Tranche& t) { return 100.0 * t.limites / std::max(1L, t.evenements); }, nullptr, 24),
             "Etat des evenements;Profondeur d'injection muZ (mm);Part des evenements (%)",
             "gain nul", "plafond de creation", false,
             graphe([](const Tranche& t) { return 100.0 * t.limGen / std::max(1L, t.evenements); }, nullptr, 25),
             "plafond de generations");
    dessiner(3, graphe([](const Tranche& t) { return t.transit; }, [](const Tranche& t) { return t.dTransit; }, 20),
             nullptr, "Temps de transit moyen (barycentre);Profondeur d'injection muZ (mm);Temps (ps)", "", "", false);
    dessiner(4, graphe([](const Tranche& t) { return t.sigmaBary; }, [](const Tranche& t) { return t.dSigmaBary; }, 20),
             graphe([](const Tranche& t) { return t.sigmaPremier; }, [](const Tranche& t) { return t.dSigmaPremier; }, 24),
             "Dispersion entre evenements (ecart-type);Profondeur d'injection muZ (mm);sigma (ps)",
             "barycentre du paquet", "premiere arrivee", false);
    canvas.SaveAs((std::string(figure) + ".png").c_str());
    canvas.SaveAs((std::string(figure) + ".pdf").c_str());
}
