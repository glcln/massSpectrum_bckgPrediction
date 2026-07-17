#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TColor.h"
#include <fstream>
#include <iostream>

void th2_sfp_eta(const char* infile = "../test.txt") {

    gStyle->SetOptStat(0);
    gStyle->SetPalette(kRainBow);
    gStyle->SetNumberContours(255);

    // eta binning: 0.25-wide bins centered on the grid (-2.275, -2.05, ..., 2.275)
    // => 20 bins from -2.4 to 2.6 covers all centers with edges at center +/- 0.125
    const int    nEta   = 20;
    const double etaLow = -2.4, etaHigh = 2.6;

    // SFp (= intP/intFp) axis
    const int    nX   = 70;
    const double xLow = 0.8, xHigh = 2;

    TH2F* h = new TH2F("h_sfp_eta",
                       ";SF_{p} = #int p / #int f_{p};#eta;Entries",
                       nX, xLow, xHigh, nEta, etaLow, etaHigh);

    std::ifstream in(infile);
    if (!in.is_open()) { std::cerr << "Cannot open " << infile << std::endl; return; }

    double sfp, eta;
    long n = 0;
    while (in >> sfp >> eta) { h->Fill(sfp, eta); ++n; }
    in.close();
    std::cout << "Filled " << n << " entries" << std::endl;

    TCanvas* c = new TCanvas("c", "SFp vs eta", 800, 700);
    c->SetRightMargin(0.15);
    c->SetLeftMargin(0.12);
    h->GetXaxis()->SetTitleOffset(1.1);
    h->GetYaxis()->SetTitleOffset(1.1);
    h->Draw("COLZ");

    // "Private work" label
    TLatex tl;
    tl.SetNDC();
    tl.SetTextFont(42);
    tl.SetTextSize(0.035);
    tl.DrawLatex(0.55, 0.92, "#it{Private work}");

    c->SaveAs("th2_sfp_eta.pdf");
}
