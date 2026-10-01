#include <TFile.h>
#include <TTree.h>
#include <TLeaf.h>
#include <TCanvas.h>
#include <TH1.h>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

// Depuis la racine : root -l -b -q 'analysis/analyse.C("mcp.root")'
void analyse(const char* fichier="tranche_0.root", const char* figure="mcp")
{
    TFile input(fichier);

    auto* events = input.Get<TTree>("events");
    if (!events) throw std::runtime_error("Table events absente");

    if (!events->GetBranch("IsMultiplicationLimited"))
        throw std::runtime_error("Colonne IsMultiplicationLimited absente : regenerer le fichier ROOT.");
    auto* config = input.Get<TTree>("configuration");
    if (!config || config->GetEntries()!=1 || !config->GetLeaf("MaxGenerations")
        || !events->GetBranch("MaxGeneration"))
        throw std::runtime_error("Configuration ou MaxGeneration absente : un run attendu.");
    config->GetEntry(0);
    const int plafond = config->GetLeaf("MaxGenerations")->GetValue();
    const auto limGen = plafond > 0
        ? events->GetEntries(("MaxGeneration>=" + std::to_string(plafond)).c_str()) : 0;
    double meanGeneration = 0;
    for (Long64_t i=0; i<events->GetEntries(); ++i) {
        events->GetEntry(i);
        meanGeneration += events->GetLeaf("MaxGeneration")->GetValue();
    }
    meanGeneration = events->GetEntries() ? meanGeneration/events->GetEntries()
        : std::numeric_limits<double>::quiet_NaN();
    std::cout << "Evenements : " << events->GetEntries()
              << ", plafond de creation atteint : " << events->GetEntries("IsMultiplicationLimited!=0")
              << ", lim_gen : " << limGen << " (plafond=" << plafond << "; 0=sans plafond)"
              << ", gen.max : " << meanGeneration
              << ", gain>=2 : " << events->GetEntries("Gain>=2") << '\n';
    TCanvas canvas("mcp_analysis", "Monocanal MCP", 1000, 750);
    canvas.Divide(2,2);
    const char* expressions[] = {"Gain", "FirstTime_ns*1000", "MeanTime_ns*1000", "TimeSpread_ns*1000"};
    const char* titles[] = {"Gain;Electrons sortants;Evenements",
        "Premiere arrivee;Temps (ps);Evenements gain>=2",
        "Barycentre temporel;Temps (ps);Evenements gain>=2",
        "Largeur RMS du paquet;Ecart-type interne (ps);Evenements gain>=2"};
    for (int i=0; i<4; ++i) {
        canvas.cd(i+1);
        const auto n = events->Draw(expressions[i], i==0 ? "" : "Gain>=2");
        if (auto* histogram = events->GetHistogram()) {
            histogram->SetTitle(titles[i]);
            if (n>0) std::cout << titles[i] << " : moyenne=" << histogram->GetMean()
                              << ", ecart-type=" << histogram->GetStdDev() << '\n';
        }
    }
    canvas.SaveAs((std::string(figure)+".png").c_str());
    canvas.SaveAs((std::string(figure)+".pdf").c_str());
}
