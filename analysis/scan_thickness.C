#include <TCanvas.h>
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

// Un point du balayage : un fichier ROOT, donc une longueur de canal.
struct Point {
    double epaisseur = 0;  // mm
    double diametre = 0;   // um
    double tension = 0;    // V
    int plafond = 0;       // MaxGenerations du run ; 0 = sans plafond
    long evenements = 0, gainNul = 0, limites = 0, limGen = 0, temps = 0;
    double gain = 0, dGain = 0;       // moyenne sur tous les evenements, gain nul compris
    double transit = 0, dTransit = 0; // barycentre a la face du MCP
    double largeur = 0, dLargeur = 0; // largeur interne moyenne d'un paquet
    Dispersion dispersion;            // dispersion ENTRE evenements du barycentre :
                                      // sigma et les trois estimateurs de largeur

    double dFwhmDirect() const { return dispersion.dFwhmDirect; }

    double Alpha() const { return diametre > 0 ? epaisseur * 1000.0 / diametre : 0; }
    double Normalisee() const { return Alpha() > 0 ? tension / Alpha() : 0; }
    double Champ() const { return epaisseur > 0 ? tension / epaisseur / 1000.0 : 0; }  // kV/mm
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
    point.epaisseur = Valeur(configuration, "Thickness_mm");
    point.diametre = Valeur(configuration, "PoreDiameter_um");
    point.tension = Valeur(configuration, "Voltage_V");
    point.plafond = static_cast<int>(Valeur(configuration, "MaxGenerations"));

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
        if (point.plafond > 0 && Valeur(evenements, "MaxGeneration") >= point.plafond) ++point.limGen;
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

    // Moyenne et erreur sur la moyenne d'un echantillon deja accumule.
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

// Depuis la racine : root -l -b -q 'analysis/scan_thickness.C()'
// void scan_thickness(const char* motif = "root_files/scanThickness_*.root",
//                     const char* figure = "figures/scan_Thickness")
// {

void scan_thickness(const char* motif = "root_files/scanField_*.root",
                    const char* figure = "figures/scan_Field")
{
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
              [](const Point& a, const Point& b) { return a.epaisseur < b.epaisseur; });

    // Deux balayages ont un sens : a tension constante (optimum) ou a champ
    // constant (croissance exponentielle). On detecte lequel.
    const auto constant = [&points](double (*v)(const Point&)) {
        return std::all_of(points.begin(), points.end(), [&](const Point& p) {
            return std::abs(v(p) - v(points.front())) <= 1e-6 * std::max(1.0, std::abs(v(points.front())));
        });
    };
    const bool memeTension = constant([](const Point& p) { return p.tension; });
    const bool memeChamp = constant([](const Point& p) { return p.Champ(); });
    if (!constant([](const Point& p) { return p.diametre; }))
        printf("Attention : le diametre de pore n'est pas constant sur le balayage.\n");
    if (memeChamp && !memeTension)
        printf("Balayage a CHAMP constant (%.2f kV/mm) : le gain croit exponentiellement "
               "avec L, sans optimum.\n", points.front().Champ());
    else if (memeTension)
        printf("Balayage a TENSION constante (%.0f V) : le champ V/L baisse avec L, "
               "le gain passe par un optimum.\n", points.front().tension);
    else
        printf("Attention : ni la tension ni le champ ne sont constants sur le balayage, "
               "la comparaison melange deux effets.\n");

    printf("%8s %8s %9s %7s %8s %9s %7s %13s %13s %10s %11s %11s %13s\n", "L(mm)", "L/D",
           "V/(L/D)", "evts", "gain0(%)", "limites(%)", "temps", "gain", "transit(ps)",
           "sigma(ps)", "TTS 2.35s", "TTS q16-84", "TTS directe");
    for (const auto& p : points) {
        printf("%8.2f %8.1f %9.1f %7ld %8.1f %9.1f %7ld %8.4g+-%-6.3g %6.1f+-%-5.1f %10.1f %11.1f %11.1f %8.1f+-%-4.1f\n",
               p.epaisseur, p.Alpha(), p.Normalisee(), p.evenements,
               100.0 * p.gainNul / p.evenements, 100.0 * p.limites / p.evenements, p.temps,
               p.gain, p.dGain, p.transit, p.dTransit, p.dispersion.sigma,
               p.dispersion.fwhmGauss, p.dispersion.fwhmRobuste, p.dispersion.fwhmDirect, p.dispersion.dFwhmDirect);
    }
    for (const auto& p : points) {
        if (p.limites > p.evenements / 20)
            printf("Attention : a L = %.2f mm, %ld evenements sur %ld au plafond de tracks, "
                   "gain minore.\n", p.epaisseur, p.limites, p.evenements);
        if (p.temps < 20)
            printf("Attention : a L = %.2f mm, seulement %ld evenements dans les temps.\n",
                   p.epaisseur, p.temps);
    }

    // Optimum : la plus grande valeur mesuree, a la granularite du balayage pres.
    const auto meilleur = std::max_element(points.begin(), points.end(),
        [](const Point& a, const Point& b) { return a.gain < b.gain; });
    printf("\nGain maximal a L = %.2f mm (L/D = %.1f, V/(L/D) = %.1f V) : %.4g\n",
           meilleur->epaisseur, meilleur->Alpha(), meilleur->Normalisee(), meilleur->gain);
    if (memeTension && (meilleur == points.begin() || meilleur + 1 == points.end()))
        printf("Cet optimum est au bord du balayage : elargir la liste des longueurs.\n");

    std::ofstream csv(std::string(figure) + ".csv");
    csv << std::setprecision(10)
        << "Thickness_mm,PoreDiameter_um,Alpha,Voltage_V,NormalizedVoltage_V,Events,GainZero,"
           "Limited,GenerationReached,TimedEvents,Gain,Gain_SEM,Transit_ps,Transit_SEM,"
           "Sigma_ps,Sigma_Err,TTS_FWHM_gauss_ps,TTS_FWHM_robust_ps,"
           "TTS_FWHM_direct_ps,InternalWidth_ps,InternalWidth_SEM\n";
    for (const auto& p : points)
        csv << p.epaisseur << ',' << p.diametre << ',' << p.Alpha() << ',' << p.tension << ','
            << p.Normalisee() << ',' << p.evenements << ',' << p.gainNul << ',' << p.limites << ','
            << p.limGen << ',' << p.temps << ',' << p.gain << ',' << p.dGain << ',' << p.transit
            << ',' << p.dTransit << ',' << p.dispersion.sigma << ',' << p.dispersion.dSigma << ','
            << p.dispersion.fwhmGauss << ',' << p.dispersion.fwhmRobuste << ',' << p.dispersion.fwhmDirect << ','
            << p.largeur << ',' << p.dLargeur << '\n';

    TCanvas canvas("scan_thickness", "Gain et temps contre la longueur du canal", 1000, 750);
    canvas.Divide(2, 2);
    const auto graphe = [&points](double (*valeur)(const Point&), double (*erreur)(const Point&),
                                  int style) {
        auto* g = new TGraphErrors;
        for (const auto& p : points) {
            const double y = valeur(p);
            if (!std::isfinite(y)) continue;
            const int k = g->GetN();
            g->SetPoint(k, p.epaisseur, y);
            g->SetPointError(k, 0, erreur && std::isfinite(erreur(p)) ? erreur(p) : 0);
        }
        g->SetMarkerStyle(style);
        return g;
    };

    // Un optimum est attendu : allonger le canal augmente le nombre de
    // collisions mais diminue le champ V/L, donc l'energie par saut.
    canvas.cd(1);
    gPad->SetLogy();
    auto* gGain = graphe([](const Point& p) { return std::max(1e-3, p.gain); },
                         [](const Point& p) { return p.dGain; }, 20);
    gGain->SetTitle("Gain moyen;Longueur du canal L (mm);Electrons sortants par injection");
    gGain->Draw("APL");
    TLatex texte;
    texte.SetNDC();
    texte.SetTextSize(0.04);
    texte.DrawLatex(0.55, 0.86, memeChamp && !memeTension
                        ? Form("E = %.2f kV/mm, D = %.1f #mum", points.front().Champ(),
                               points.front().diametre)
                        : Form("V = %.0f V, D = %.1f #mum", points.front().tension,
                               points.front().diametre));

    canvas.cd(2);
    auto* gZero = graphe([](const Point& p) { return 100.0 * p.gainNul / p.evenements; }, nullptr, 20);
    auto* gLim = graphe([](const Point& p) { return 100.0 * p.limites / p.evenements; }, nullptr, 24);
    gZero->SetTitle("Etat des evenements;Longueur du canal L (mm);Part des evenements (%)");
    gZero->GetHistogram()->SetMinimum(-5);
    gZero->GetHistogram()->SetMaximum(115);
    gZero->Draw("APL");
    gLim->SetLineStyle(2);
    gLim->Draw("PL same");
    auto* legende = new TLegend(0.35, 0.74, 0.88, 0.88);
    legende->AddEntry(gZero, "gain nul", "pl");
    legende->AddEntry(gLim, "plafond de tracks atteint", "pl");
    legende->Draw();

    canvas.cd(3);
    auto* gT = graphe([](const Point& p) { return p.transit; },
                      [](const Point& p) { return p.dTransit; }, 20);
    gT->SetTitle("Temps de transit (barycentre);Longueur du canal L (mm);Temps (ps)");
    gT->Draw("APL");

    // TTS, estime par la largeur a mi-hauteur mesuree sur l'histogramme, et
    // largeur interne du paquet : deux populations differentes.
    canvas.cd(4);
    auto* gTTS = graphe([](const Point& p) { return p.dispersion.fwhmDirect; },
                        [](const Point& p) { return p.dFwhmDirect(); }, 20);
    gTTS->SetTitle("TTS (largeur a mi-hauteur);Longueur du canal L (mm);FWHM (ps)");
    gTTS->Draw("APL");

    canvas.SaveAs((std::string(figure) + ".png").c_str());
    canvas.SaveAs((std::string(figure) + ".pdf").c_str());

    std::cout << "Gain : moyenne sur tous les evenements, gain nul compris ; erreur sur la moyenne.\n"
              << "L/D : rapport longueur sur diametre. V/(L/D) : tension normalisee, le parametre\n"
              << "qui fixe l'energie acquise entre deux collisions.\n"
              << "Temps en ps depuis l'injection, evenements avec gain>=2, mesures a la face du MCP.\n"
              << "sigma : ecart-type ENTRE evenements du barycentre temporel, a ne pas confondre\n"
              << "avec la largeur interne d'un paquet, qui porte sur les electrons d'un meme\n"
              << "evenement.\n"
              << "TTS : exprime ici en largeur a mi-hauteur. Trois estimateurs non equivalents\n"
              << "(FWHM equivalente gaussienne, inter-quantiles 16-84 convertie, mesure directe sur\n"
              << "l'histogramme) ; un ecart entre eux suggere une distribution non gaussienne.\n";
}
