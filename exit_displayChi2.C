#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include "TH1F.h"
#include "TFile.h"

void exit_displayChi2(const std::string& inputTxt = "exit.txt",
              const std::string& outputRoot = "exit_displayChi2.root") {

    // Quatre histogrammes, un par categorie
    TH1F* h_ip      = new TH1F("h_1oP",     "1/p;#chi^{2}/NDF;Entries",      100, 0, 10);
    TH1F* h_ih      = new TH1F("h_Ih",      "Ih;#chi^{2}/NDF;Entries",       100, 0, 10);
    TH1F* h_bad_ip  = new TH1F("h_bad_1oP", "bad 1/p;#chi^{2}/NDF;Entries",  100, 0, 10);
    TH1F* h_bad_ih  = new TH1F("h_bad_Ih",  "bad Ih;#chi^{2}/NDF;Entries",   100, 0, 10);

    std::ifstream fin(inputTxt.c_str());
    if (!fin.is_open()) {
        std::cerr << "Erreur : impossible d'ouvrir " << inputTxt << std::endl;
        return;
    }

    std::string line;
    while (std::getline(fin, line)) {

        // On cherche la position du ':'
        size_t colon = line.find(':');
        if (colon == std::string::npos) continue;

        // Label = ce qui est avant le ':' (sans espaces de fin)
        std::string label = line.substr(0, colon);
        // trim des espaces en debut/fin du label
        size_t a = label.find_first_not_of(" \t");
        size_t b = label.find_last_not_of(" \t");
        if (a == std::string::npos) continue;
        label = label.substr(a, b - a + 1);

        // Valeur = ce qui est apres le ':'
        std::string valStr = line.substr(colon + 1);
        std::istringstream iss(valStr);
        double val;
        if (!(iss >> val)) continue;   // pas un nombre -> on saute (filtre les autres lignes)

        // Aiguillage exact sur le label
        if      (label == "1/p")     h_ip->Fill(val);
        else if (label == "Ih")      h_ih->Fill(val);
        else if (label == "bad 1/p") h_bad_ip->Fill(val);
        else if (label == "bad Ih")  h_bad_ih->Fill(val);
        // toute autre ligne (Mean histo..., Before overflow...) est ignoree
    }
    fin.close();

    // Sauvegarde
    TFile* fout = new TFile(outputRoot.c_str(), "RECREATE");
    h_ip->Write();
    h_ih->Write();
    h_bad_ip->Write();
    h_bad_ih->Write();
    fout->Close();

    std::cout << "Entrees -> 1/p: " << h_ip->GetEntries()
              << " | Ih: "      << h_ih->GetEntries()
              << " | bad 1/p: " << h_bad_ip->GetEntries()
              << " | bad Ih: "  << h_bad_ih->GetEntries() << std::endl;
    std::cout << "Histos sauves dans " << outputRoot << std::endl;
}