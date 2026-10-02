#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TLeaf.h>
#include <TLegend.h>
#include <TRandom3.h>
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

// Largeur a mi-hauteur mesuree directement sur l'histogramme d'un echantillon
// TRIE, par interpolation lineaire de part et d'autre du maximum. La valeur
// depend du binning (choisi ici en n/8, borne a [10, 60]) et de la position du
// mode ; elle n'est pas definie pour une distribution multimodale.
static double FwhmDirecte(const std::vector<double>& trie)
{
    const long n = static_cast<long>(trie.size());
    if (n < 50 || trie.back() <= trie.front()) return std::nan("");
    const int bins = std::max(10, std::min(60, static_cast<int>(n / 8)));
    const double lo = trie.front(), large = (trie.back() - lo) / bins;
    std::vector<double> h(bins, 0.0);
    for (double x : trie) {
        const int k = static_cast<int>((x - lo) / large);
        ++h[std::max(0, std::min(bins - 1, k))];
    }
    const int pic = static_cast<int>(std::max_element(h.begin(), h.end()) - h.begin());
    const double demi = h[pic] / 2.0;
    const auto croisement = [&](int pas) {
        for (int k = pic + pas; k >= 0 && k < bins; k += pas) {
            if (h[k] < demi) {
                const double f = (h[k - pas] - demi) / std::max(1e-9, h[k - pas] - h[k]);
                return lo + large * (k - pas + 0.5 + pas * f);
            }
        }
        return std::nan("");
    };
    const double gauche = croisement(-1), droite = croisement(+1);
    return (std::isfinite(gauche) && std::isfinite(droite)) ? droite - gauche : std::nan("");
}

// Dispersion ENTRE evenements du barycentre temporel, a distinguer de la
// dispersion INTERNE a un evenement (TimeSpread_ns dans les tables ROOT).
// Le TTS est ici exprime en largeur a mi-hauteur, convention retenue dans cette
// analyse parce que c'est la forme la plus souvent publiee. Trois estimateurs
// en sont fournis ; ils ne coincident que si la distribution est gaussienne.
struct Dispersion {
    long n = 0;
    double sigma = std::nan(""), dSigma = std::nan("");
    double fwhmGauss = std::nan("");    // FWHM equivalente gaussienne, 2.3548 x sigma
    double fwhmRobuste = std::nan("");  // largeur inter-quantiles 16-84 convertie en
                                        // FWHM equivalente gaussienne, 1.1774 x (q84-q16)
    double fwhmDirect = std::nan(""), dFwhmDirect = std::nan("");
};

static Dispersion Disperser(std::vector<double> v)
{
    Dispersion d;
    d.n = static_cast<long>(v.size());
    if (d.n < 2) return d;

    double s = 0;
    for (double x : v) s += x;
    const double moyenne = s / d.n;
    double s2 = 0;
    for (double x : v) s2 += (x - moyenne) * (x - moyenne);
    d.sigma = std::sqrt(s2 / (d.n - 1));
    d.dSigma = d.n > 2 ? d.sigma / std::sqrt(2.0 * (d.n - 1)) : std::nan("");
    d.fwhmGauss = 2.3548 * d.sigma;

    std::sort(v.begin(), v.end());
    const auto quantile = [&v](double f) {
        const double x = f * (v.size() - 1);
        const size_t i = static_cast<size_t>(x);
        const double r = x - i;
        return i + 1 < v.size() ? v[i] * (1 - r) + v[i + 1] * r : v[i];
    };
    if (d.n >= 10) d.fwhmRobuste = 1.1774 * (quantile(0.84) - quantile(0.16));

    d.fwhmDirect = FwhmDirecte(v);
    // Pas d'erreur analytique sur cet estimateur : incertitude estimee par
    // bootstrap (200 retirages avec remise, graine fixe pour reproductibilite).
    // Elle ne couvre que la fluctuation d'echantillonnage, pas le biais de
    // binning.
    if (std::isfinite(d.fwhmDirect)) {
        TRandom3 alea(12345);
        std::vector<double> tirage(v.size()), mesures;
        for (int essai = 0; essai < 200; ++essai) {
            for (size_t i = 0; i < v.size(); ++i) tirage[i] = v[alea.Integer(v.size())];
            std::sort(tirage.begin(), tirage.end());
            const double f = FwhmDirecte(tirage);
            if (std::isfinite(f)) mesures.push_back(f);
        }
        if (mesures.size() > 2) {
            double sb = 0;
            for (double x : mesures) sb += x;
            const double mb = sb / mesures.size();
            double sb2 = 0;
            for (double x : mesures) sb2 += (x - mb) * (x - mb);
            d.dFwhmDirect = std::sqrt(sb2 / (mesures.size() - 1));
        }
    }
    return d;
}

// Un point du balayage : un fichier ROOT, donc une tension.
struct Point {
    double voltage = 0;
    long evenements = 0, gainNul = 0, limites = 0, temps = 0;
    double gain = 0, dGain = 0;       // moyenne sur tous les evenements, gain nul compris
    double transit = 0, dTransit = 0; // barycentre a la face du MCP
    double largeur = 0, dLargeur = 0; // largeur interne moyenne d'un paquet
    Dispersion dispersion;            // dispersion ENTRE evenements du barycentre :
                                      // sigma et les trois estimateurs de largeur

    double dFwhmDirect() const { return dispersion.dFwhmDirect; }
};

static double Valeur(TTree* t, const char* nom)
{
    auto* feuille = t->GetLeaf(nom);
    if (!feuille) throw std::runtime_error(std::string("Colonne absente : ") + nom);
    return feuille->GetValue();
}

static Point LirePoint(const std::string& chemin)
{
    TFile fichier(chemin.c_str());
    if (fichier.IsZombie()) throw std::runtime_error("Fichier illisible : " + chemin);
    auto* configuration = fichier.Get<TTree>("configuration");
    auto* evenements = fichier.Get<TTree>("events");
    if (!configuration || !evenements) throw std::runtime_error("Tables absentes : " + chemin);
    if (configuration->GetEntries() != 1)
        throw std::runtime_error("Un seul run attendu par fichier : " + chemin);

    Point point;
    configuration->GetEntry(0);
    point.voltage = Valeur(configuration, "Voltage_V");

    std::vector<double> barycentres;
    double sg = 0, sg2 = 0, st = 0, st2 = 0, sl = 0, sl2 = 0;
    for (Long64_t i = 0; i < evenements->GetEntries(); ++i) {
        evenements->GetEntry(i);
        ++point.evenements;
        const double gain = Valeur(evenements, "Gain");
        sg += gain;
        sg2 += gain * gain;
        if (gain == 0) ++point.gainNul;
        if (Valeur(evenements, "IsMultiplicationLimited") != 0) ++point.limites;
        // Un barycentre temporel n'a de sens qu'avec au moins deux electrons
        // en sortie ; les gains 0 et 1 sont donc exclus des statistiques de temps.
        const double transit = Valeur(evenements, "MeanTime_ns") * 1000;
        const double largeur = Valeur(evenements, "TimeSpread_ns") * 1000;
        if (gain < 2 || !std::isfinite(transit) || !std::isfinite(largeur)) continue;
        ++point.temps;
        barycentres.push_back(transit);
        st += transit;
        st2 += transit * transit;
        sl += largeur;
        sl2 += largeur * largeur;
    }

    const auto resume = [](double s, double s2, long n, double& m, double& e) {
        if (n <= 0) { m = e = std::nan(""); return; }
        m = s / n;
        e = n > 1 ? std::sqrt(std::max(0.0, s2 / n - m * m) / (n - 1)) : 0;
    };
    resume(sg, sg2, point.evenements, point.gain, point.dGain);
    resume(st, st2, point.temps, point.transit, point.dTransit);
    resume(sl, sl2, point.temps, point.largeur, point.dLargeur);
    point.dispersion = Disperser(barycentres);
    return point;
}

// Depuis la racine : root -l -b -q 'analysis/gain_voltage.C("root_files/scanV_*.root")'
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

void gain_voltage(const char* motif = "root_files/scanV_*.root",
                  const char* figure = "figures/gain_voltage")
{
    PreparerSortie(figure);
    const TString dossier = gSystem->DirName(motif), patron = gSystem->BaseName(motif);
    TSystemDirectory repertoire(dossier, dossier);
    auto* contenu = repertoire.GetListOfFiles();
    if (!contenu) throw std::runtime_error("Repertoire introuvable : " + std::string(dossier.Data()));

    std::vector<Point> points;
    TRegexp filtre(patron, kTRUE);
    TIter suivant(contenu);
    while (auto* entree = static_cast<TSystemFile*>(suivant())) {
        const TString nom = entree->GetName();
        Ssiz_t longueur = 0;
        if (entree->IsDirectory() || nom.Index(filtre, &longueur) != 0 || longueur != nom.Length())
            continue;
        points.push_back(LirePoint(std::string(dossier.Data()) + "/" + nom.Data()));
    }
    if (points.empty()) throw std::runtime_error("Aucun fichier ne correspond a " + std::string(motif));
    std::sort(points.begin(), points.end(),
              [](const Point& a, const Point& b) { return a.voltage < b.voltage; });

    printf("%8s %7s %9s %10s %14s %14s %11s %11s %11s %11s\n", "V", "evts", "gain0(%)",
           "limites(%)", "gain", "transit(ps)", "sigma(ps)", "TTS 2.35s", "TTS q16-84",
           "TTS directe");
    for (const auto& p : points) {
        printf("%8.0f %7ld %9.1f %10.1f %9.4g+-%-8.3g %6.1f+-%-6.1f %5.1f+-%-5.1f %11.1f %11.1f %11.1f\n",
               p.voltage, p.evenements, 100.0 * p.gainNul / p.evenements,
               100.0 * p.limites / p.evenements, p.gain, p.dGain,
               p.transit, p.dTransit, p.dispersion.sigma, p.dispersion.dSigma,
               p.dispersion.fwhmGauss, p.dispersion.fwhmRobuste, p.dispersion.fwhmDirect);
    }
    for (const auto& p : points) {
        if (p.limites > p.evenements / 20)
            printf("Attention : a %.0f V, %ld evenements sur %ld au plafond de tracks, "
                   "gain minore.\n", p.voltage, p.limites, p.evenements);
    }

    std::ofstream csv(std::string(figure) + ".csv");
    csv << std::setprecision(10)
        << "Voltage_V,Events,GainZero,Limited,TimedEvents,Gain,Gain_SEM,"
           "Transit_ps,Transit_SEM,Sigma_ps,Sigma_Err,TTS_FWHM_gauss_ps,"
           "TTS_FWHM_robust_ps,TTS_FWHM_direct_ps,InternalWidth_ps,InternalWidth_SEM\n";
    for (const auto& p : points)
        csv << p.voltage << ',' << p.evenements << ',' << p.gainNul << ',' << p.limites << ','
            << p.temps << ',' << p.gain << ',' << p.dGain << ',' << p.transit << ','
            << p.dTransit << ',' << p.dispersion.sigma << ',' << p.dispersion.dSigma << ','
            << p.dispersion.fwhmGauss << ',' << p.dispersion.fwhmRobuste << ',' << p.dispersion.fwhmDirect << ','
            << p.largeur << ',' << p.dLargeur << '\n';

    // Courbes contre la tension appliquee.
    TCanvas canvas("gain_voltage", "Gain et temps contre la tension", 1000, 750);
    canvas.Divide(2, 2);
    const auto graphe = [&points](double (*valeur)(const Point&), double (*erreur)(const Point&),
                                  int style) {
        auto* g = new TGraphErrors;
        for (const auto& p : points) {
            const double y = valeur(p);
            if (!std::isfinite(y)) continue;
            const int k = g->GetN();
            g->SetPoint(k, p.voltage, y);
            g->SetPointError(k, 0, erreur && std::isfinite(erreur(p)) ? erreur(p) : 0);
        }
        g->SetMarkerStyle(style);
        return g;
    };

    // Ajustement exponentiel empirique sur la plage simulee : il resume la
    // tendance, il ne decoule pas d'une loi du modele. Sans charge d'espace,
    // rien ne limite la croissance a haute tension.
    canvas.cd(1);
    gPad->SetLogy();
    auto* gGain = graphe([](const Point& p) { return std::max(1e-3, p.gain); },
                         [](const Point& p) { return p.dGain; }, 20);
    gGain->SetTitle("Gain moyen;Tension appliquee (V);Electrons sortants par injection");
    gGain->Draw("AP");
    double pente = 0;
    if (gGain->GetN() > 2) {
        auto* fit = new TF1("fit", "expo", points.front().voltage, points.back().voltage);
        gGain->Fit(fit, "Q");
        pente = fit->GetParameter(1);
        TLatex texte;
        texte.SetNDC();
        texte.SetTextSize(0.045);
        texte.DrawLatex(0.18, 0.84, Form("G #propto e^{kV}, k = %.4f V^{-1}", pente));
    }

    // 2 : etat des evenements.
    canvas.cd(2);
    auto* gZero = graphe([](const Point& p) { return 100.0 * p.gainNul / p.evenements; }, nullptr, 20);
    auto* gLim = graphe([](const Point& p) { return 100.0 * p.limites / p.evenements; }, nullptr, 24);
    gZero->SetTitle("Etat des evenements;Tension appliquee (V);Part des evenements (%)");
    gZero->GetHistogram()->SetMinimum(-5);
    gZero->GetHistogram()->SetMaximum(105);
    gZero->Draw("APL");
    gLim->SetLineStyle(2);
    gLim->Draw("PL same");
    auto* legende = new TLegend(0.35, 0.72, 0.88, 0.88);
    legende->AddEntry(gZero, "gain nul", "pl");
    legende->AddEntry(gLim, "plafond de tracks atteint", "pl");
    legende->Draw();

    // 3 : temps de transit moyen a la face du MCP.
    canvas.cd(3);
    auto* gT = graphe([](const Point& p) { return p.transit; },
                      [](const Point& p) { return p.dTransit; }, 20);
    gT->SetTitle("Temps de transit (barycentre);Tension appliquee (V);Temps (ps)");
    gT->Draw("APL");

    // TTS, estime ici par la largeur a mi-hauteur mesuree sur l'histogramme.
    // A ne pas confondre avec la largeur interne d'un paquet, qui porte sur une
    // autre population (les electrons d'un meme evenement).
    canvas.cd(4);
    auto* gF = graphe([](const Point& p) { return p.dispersion.fwhmDirect; },
                      [](const Point& p) { return p.dFwhmDirect(); }, 20);
    gF->SetTitle("TTS (largeur a mi-hauteur);Tension appliquee (V);FWHM (ps)");
    gF->Draw("APL");

    canvas.SaveAs((std::string(figure) + ".png").c_str());
    canvas.SaveAs((std::string(figure) + ".pdf").c_str());

    std::cout << "Gain : moyenne sur tous les evenements, gain nul compris ; erreur sur la moyenne.\n"
              << "sigma : ecart-type ENTRE evenements du barycentre temporel, a la face du MCP.\n"
              << "TTS : exprime ici en largeur a mi-hauteur, convention retenue pour cette analyse\n"
              << "car c'est la forme la plus souvent publiee. Trois estimateurs, non equivalents :\n"
              << "  2.355 x sigma      FWHM equivalente gaussienne, exacte seulement si la\n"
              << "                     distribution l'est ;\n"
              << "  quantiles 16-84    largeur robuste aux queues, convertie en FWHM equivalente\n"
              << "                     gaussienne par le meme facteur ;\n"
              << "  directe            mesuree sur l'histogramme (>= 50 evenements) ; depend du\n"
              << "                     binning et suppose une distribution unimodale.\n"
              << "Un ecart entre ces trois estimateurs suggere que la distribution s'ecarte de la\n"
              << "forme gaussienne, sans en preciser la nature.\n"
              << "Temps en ps depuis l'injection, evenements avec gain>=2, a la face du MCP.\n"
              << "Aux points ou le plafond de tracks est atteint, le gain moyen est minore.\n";
    if (pente > 0)
        std::cout << "Pente k = " << pente << " /V (ajustement empirique). Extrapolee, "
                     "elle place un gain 10^3 vers "
                  << Form("%.0f", points.back().voltage
                                      - std::log(points.back().gain / 1000.) / pente)
                  << " V.\n";
}
