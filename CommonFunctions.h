#pragma once

#include <TCanvas.h>
#include <TLegend.h>
#include "TFile.h"
#include "TH1.h"
#include "TDirectory.h"
#include <TRatioPlot.h>
#include <THStack.h>
#include <TROOT.h>
#include <TChain.h>
#include "TRandom3.h"
#include <TH2.h>
#include <TH3.h>
#include "TF1.h"
#include "TFitResult.h"
#include "ROOT/TProcessExecutor.hxx"
#include "ROOT/TSeq.hxx"
#include <numeric>
#include <TStyle.h>
#include <TGraphErrors.h>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <memory>

#include "Regions.h"

using namespace std::placeholders;


constexpr UInt_t kSeedBase = 1;

static const std::vector<double> RebinEta_Down_Eta1_AND_Eta1_2p4 = {-2.4, -2.0, -1.5, -1., -0.65, -0.30, 0, 0.30, 0.65, 1., 1.5, 2.0, 2.4};
static const std::vector<double> RebinEta_Down_Eta2p4            = {-2.4, -2.0, -1.6, -1.2, -0.8, -0.4, 0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4};
static const std::vector<double> RebinEta_Nom_ALLeta             = {-2.4, -2.15, -1.95, -1.75, -1.5, -1.25, -1., -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1., 1.25, 1.5, 1.75, 1.95, 2.15, 2.4};
static const std::vector<double> RebinEta_Up_ALLeta              = {-2.4, -2.2, -2.0, -1.8, -1.6, -1.4, -1.2, -1., -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1., 1.2, 1.4, 1.6, 1.8, 2., 2.2, 2.4};
static const std::vector<double> RebinEta_Down_Eta1p2_2p2        = {-2.4, -2.2, -1.85, -1.55, -1.2, +1.2, +1.55, +1.85, +2.2, +2.4};
static const std::vector<double> RebinEta_Nom_Eta1p2_2p2         = {-2.4, -2.2, -1.95, -1.70, -1.45, -1.2, +1.2, +1.45, +1.70, +1.95, +2.2, +2.4};
static const std::vector<double> RebinEta_Up_Eta1p2_2p2          = {-2.4, -2.2, -2.0, -1.8, -1.6, -1.4, -1.2, +1.2, +1.4, +1.6, +1.8, +2.0, +2.2, +2.4};
static const std::vector<double> RebinEta_Down_Eta1p2_2p4        = {-2.4, -2.0, -1.6, -1.2, +1.2, +1.6, +2.0, +2.4};
static const std::vector<double> RebinEta_Nom_Eta1p2_2p4         = {-2.4, -2.1, -1.8, -1.5, -1.2, +1.2, +1.5, +1.8, +2.1, +2.4};
static const std::vector<double> RebinEta_Up_Eta1p2_2p4          = {-2.4, -2.15, -1.90, -1.70, -1.45, -1.2, +1.2, +1.45, +1.70, +1.90, +2.15, +2.4};

const std::vector<double>* PickEtaBinning(const std::string& regionName, int rebineta) {
    if (regionName.find("Eta1p2_2p2") != std::string::npos)
        return (rebineta==8) ? &RebinEta_Down_Eta1p2_2p2
             : (rebineta==4) ? &RebinEta_Nom_Eta1p2_2p2
                             : &RebinEta_Up_Eta1p2_2p2;
    if (regionName.find("Eta1p2_2p4") != std::string::npos)
        return (rebineta==8) ? &RebinEta_Down_Eta1p2_2p4
             : (rebineta==4) ? &RebinEta_Nom_Eta1p2_2p4
                             : &RebinEta_Up_Eta1p2_2p4;
    if (regionName.find("Eta1_2p4")   != std::string::npos)
        return (rebineta==8) ? &RebinEta_Down_Eta1_AND_Eta1_2p4
             : (rebineta==4) ? &RebinEta_Nom_ALLeta
                             : &RebinEta_Up_ALLeta;
    if (regionName.find("Eta2p4")     != std::string::npos)
        return (rebineta==8) ? &RebinEta_Down_Eta2p4
             : (rebineta==4) ? &RebinEta_Nom_ALLeta
                             : &RebinEta_Up_ALLeta;
    if (regionName.find("Eta1")       != std::string::npos) 
        return (rebineta==8) ? &RebinEta_Down_Eta1_AND_Eta1_2p4
             : (rebineta==4) ? &RebinEta_Nom_ALLeta
                             : &RebinEta_Up_ALLeta;
    else {
        std::cerr << "Error: region name does not contain expected eta range for rebinning" << std::endl;
        return nullptr;
    }

    return nullptr;
}

// Les histos du step1 peuvent etre des TH1F : conversion explicite obligatoire,
// un cast direct TH1F* -> TH1D* est un comportement indefini.
TH1D* GetAsTH1D(TFile* f, const std::string& name) {
    TH1* h = dynamic_cast<TH1*>(f->Get(name.c_str()));
    if (!h || h->InheritsFrom(TH2::Class())) {
        std::cerr << "GetAsTH1D: '" << name << "' absent du fichier ou n'est pas un TH1" << std::endl;
        return nullptr;
    }
    if (h->InheritsFrom(TH1D::Class())) {
        TH1D* copy = static_cast<TH1D*>(h->Clone((name + "_copy").c_str()));
        copy->SetDirectory(nullptr);
        return copy;
    }

    const TAxis* ax = h->GetXaxis();
    TH1D* out = nullptr;
    if (ax->GetXbins()->GetSize() > 0)
        out = new TH1D((name + "_copy").c_str(), h->GetTitle(), ax->GetNbins(), ax->GetXbins()->GetArray());
    else
        out = new TH1D((name + "_copy").c_str(), h->GetTitle(), ax->GetNbins(), ax->GetXmin(), ax->GetXmax());

    out->SetDirectory(nullptr);
    out->Sumw2();
    for (int b = 0; b < h->GetNcells(); ++b) {
        out->SetBinContent(b, h->GetBinContent(b));
        out->SetBinError  (b, h->GetBinError(b));
    }
    out->SetEntries(h->GetEntries());
    return out;
}

// Les histos du step1 peuvent etre des TH2F : conversion explicite obligatoire,
// un cast direct TH2F* -> TH2D* est un comportement indefini.
TH2D* GetAsTH2D(TFile* f, const std::string& name) {
    TH2* h = dynamic_cast<TH2*>(f->Get(name.c_str()));
    if (!h) {
        std::cerr << "GetAsTH2D: '" << name << "' absent du fichier ou n'est pas un TH2" << std::endl;
        return nullptr;
    }
    if (h->InheritsFrom(TH2D::Class())) {
        TH2D* copy = static_cast<TH2D*>(h->Clone((name + "_copy").c_str()));
        copy->SetDirectory(nullptr);
        return copy;
    }

    const TAxis* ax = h->GetXaxis();
    const TAxis* ay = h->GetYaxis();
    if (ax->GetXbins()->GetSize() > 0 || ay->GetXbins()->GetSize() > 0) {
        std::cerr << "GetAsTH2D: binning variable non gere pour '" << name << "'" << std::endl;
        return nullptr;   // le step1 n'ecrit que du binning uniforme
    }
    TH2D* out = new TH2D((name + "_copy").c_str(), h->GetTitle(),
                         ax->GetNbins(), ax->GetXmin(), ax->GetXmax(),
                         ay->GetNbins(), ay->GetXmin(), ay->GetXmax());
    out->SetDirectory(nullptr);
    out->Sumw2();
    for (int b = 0; b < h->GetNcells(); ++b) {
        out->SetBinContent(b, h->GetBinContent(b));
        out->SetBinError  (b, h->GetBinError(b));
    }
    out->SetEntries(h->GetEntries());
    return out;
}


TH2D* RebinTH2Y_varBins(TH2D* h, int nEta, const double* eEta) {
    int nX = h->GetNbinsX();
    const TArrayD* xArr = h->GetXaxis()->GetXbins();
    TH2D* hNew;
    if (xArr->GetSize() > 0)
        hNew = new TH2D(Form("%s_rebinEta", h->GetName()), h->GetTitle(),
                        nX, xArr->GetArray(),
                        nEta, eEta);
    else
        hNew = new TH2D(Form("%s_rebinEta", h->GetName()), h->GetTitle(),
                        nX, h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax(),
                        nEta, eEta);
    hNew->SetDirectory(nullptr);
    hNew->Sumw2();

    for (int ix = 0; ix <= nX+1; ix++)
        for (int iy = 0; iy <= h->GetNbinsY()+1; iy++) {
            double eta = h->GetYaxis()->GetBinCenter(iy);
            int newBin = hNew->GetYaxis()->FindBin(eta);
            hNew->SetBinContent(ix, newBin,
                hNew->GetBinContent(ix, newBin) + h->GetBinContent(ix, iy));
            hNew->SetBinError(ix, newBin,
                std::sqrt(std::pow(hNew->GetBinError(ix, newBin), 2)
                        + std::pow(h->GetBinError(ix, iy), 2)));
        }
    return hNew;
}

TH2D* RebinTH2X_varBins(TH2D* h, int nEta, const double* eEta) {
    int nY = h->GetNbinsY();
    const TArrayD* yArr = h->GetYaxis()->GetXbins();
    TH2D* hNew;
    if (yArr->GetSize() > 0)
        hNew = new TH2D(Form("%s_rebinEta", h->GetName()), h->GetTitle(),
                        nEta, eEta,
                        nY, yArr->GetArray());
    else
        hNew = new TH2D(Form("%s_rebinEta", h->GetName()), h->GetTitle(),
                        nEta, eEta,
                        nY, h->GetYaxis()->GetXmin(), h->GetYaxis()->GetXmax());
    hNew->SetDirectory(nullptr);
    hNew->Sumw2();

    for (int ix = 0; ix <= h->GetNbinsX()+1; ix++)
        for (int iy = 0; iy <= nY+1; iy++) {
            double eta = h->GetXaxis()->GetBinCenter(ix);
            int newBin = hNew->GetXaxis()->FindBin(eta);
            hNew->SetBinContent(newBin, iy,
                hNew->GetBinContent(newBin, iy) + h->GetBinContent(ix, iy));
            hNew->SetBinError(newBin, iy,
                std::sqrt(std::pow(hNew->GetBinError(newBin, iy), 2)
                        + std::pow(h->GetBinError(ix, iy), 2)));
        }
    return hNew;
}

TH2D* FoldAbsTH2X(TH2D* h, const std::string& newName) {
    int nx = h->GetNbinsX();
    int ny = h->GetNbinsY();
    const TAxis* ax = h->GetXaxis();

    // Index du premier bin dont la borne basse >= 0 (frontiere a eta=0)
    int izero = ax->FindBin(0.0 + 1e-9);
    // Nombre de bins positifs
    int nxPos = nx - izero + 1;

    // Nouveau binning positif
    std::vector<double> xedges;
    for (int i = izero; i <= nx + 1; ++i) xedges.push_back(ax->GetBinLowEdge(i));
    std::vector<double> yedges;
    for (int j = 1; j <= ny + 1; ++j) yedges.push_back(h->GetYaxis()->GetBinLowEdge(j));

    TH2D* hf = new TH2D(newName.c_str(), h->GetTitle(),
                        nxPos, xedges.data(), ny, yedges.data());
    hf->SetDirectory(nullptr);
    hf->Sumw2();

    for (int j = 1; j <= ny; ++j) {
        for (int i = 1; i <= nx; ++i) {
            double xc = ax->GetBinCenter(i);
            int io = hf->GetXaxis()->FindBin(std::fabs(xc)); // bin positif cible
            double c = hf->GetBinContent(io, j) + h->GetBinContent(i, j);
            double e = std::hypot(hf->GetBinError(io, j), h->GetBinError(i, j));
            hf->SetBinContent(io, j, c);
            hf->SetBinError(io, j, e);
        }
    }
    return hf;
}

TH2D* FoldAbsTH2Y(TH2D* h, const std::string& newName) {
    int nx = h->GetNbinsX();
    int ny = h->GetNbinsY();
    const TAxis* ay = h->GetYaxis();

    int jzero = ay->FindBin(0.0 + 1e-9);
    int nyPos = ny - jzero + 1;

    std::vector<double> xedges;
    for (int i = 1; i <= nx + 1; ++i) xedges.push_back(h->GetXaxis()->GetBinLowEdge(i));
    std::vector<double> yedges;
    for (int j = jzero; j <= ny + 1; ++j) yedges.push_back(ay->GetBinLowEdge(j));

    TH2D* hf = new TH2D(newName.c_str(), h->GetTitle(),
                        nx, xedges.data(), nyPos, yedges.data());
    hf->SetDirectory(nullptr);
    hf->Sumw2();

    for (int j = 1; j <= ny; ++j) {
        double yc = ay->GetBinCenter(j);
        int jo = hf->GetYaxis()->FindBin(std::fabs(yc));
        for (int i = 1; i <= nx; ++i) {
            double c = hf->GetBinContent(i, jo) + h->GetBinContent(i, j);
            double e = std::hypot(hf->GetBinError(i, jo), h->GetBinError(i, j));
            hf->SetBinContent(i, jo, c);
            hf->SetBinError(i, jo, e);
        }
    }
    return hf;
}

template <typename TH>
TH MeanOfToys(const std::vector<TH>& toys, const char* name, bool useSEM = false) {
    if (toys.empty()) throw std::runtime_error("MeanOfToys: vecteur vide");
    const double N = static_cast<double>(toys.size());

    TH hMean(toys.front());
    hMean.SetName(name);
    hMean.SetTitle(name);
    hMean.Reset("ICESM");
    hMean.SetBinErrorOption(TH1::EBinErrorOpt::kNormal);
    hMean.SetDirectory(nullptr);
    hMean.Sumw2();

    const int nTot = hMean.GetNcells();          // couvre 1D et 2D, under/overflow inclus
    for (int b = 0; b < nTot; ++b) {
        double sum = 0., sum2 = 0.;
        for (const auto& h : toys) {
            const double v = h.GetBinContent(b);
            sum += v; sum2 += v*v;
        }
        const double mean = sum / N;
        double var = (N > 1) ? (sum2 - N*mean*mean) / (N - 1.) : 0.;
        if (var < 0.) var = 0.;
        hMean.SetBinContent(b, mean);
        hMean.SetBinError(b, std::sqrt(var) / (useSEM ? std::sqrt(N) : 1.));
    }
    return hMean;
}


bool loadHistograms(Region& r, 
                    TFile* f,
                    const std::string& regionName,
                    bool bool_rebin = true,
                    int rebineta = 1,
                    int rebinp = 1,
                    int rebinih = 1,
                    bool TakeAbsEta = false) {

    std::cout << "loading region " << regionName << "    rebineta=" << rebineta << ", rebinp=" << rebinp << ", rebinih=" << rebinih << std::endl;

    if (rebinp==4) rebinp = 8;
    if (rebinp==2) rebinp = 6;
    if (rebinp==1) rebinp = 4;

    r.eta_p    = GetAsTH2D(f, "eta_1oP_"  + regionName);
    r.ih_eta   = GetAsTH2D(f, "ih_eta_"   + regionName);
    r.mass     = GetAsTH1D(f, "mass_"     + regionName);
    r.mass_eta = GetAsTH2D(f, "mass_eta_" + regionName);

    if (!r.eta_p || !r.ih_eta || !r.mass || !r.mass_eta) {
        std::cerr << "loadHistograms: region '" << regionName << "' incomplete -> abandon" << std::endl;
        return false;
    }

    r.pred_mass     = (TH1D*) r.mass->Clone();
    r.pred_mass->SetDirectory(nullptr);
    r.pred_mass->SetName(("pred_mass_"+regionName).c_str());
    r.pred_mass->Reset();

    r.pred_mass_eta = (TH2D*) r.mass_eta->Clone();
    r.pred_mass_eta->SetDirectory(nullptr);
    r.pred_mass_eta->SetName(("pred_mass_eta_"+regionName).c_str());
    r.pred_mass_eta->Reset();

    if (bool_rebin) {

        if (rebineta==2 || rebineta==4 || rebineta==8) {

            const std::vector<double>* RebinEtaVecPtr = PickEtaBinning(regionName, rebineta);
            if (!RebinEtaVecPtr) {
                std::cerr << "loadHistograms: pas de binning eta pour '" << regionName
                        << "' -> region non chargee" << std::endl;
                return false;
            }
            const int     nEta = static_cast<int>(RebinEtaVecPtr->size()) - 1;
            const double* eEta = RebinEtaVecPtr->data();

            // eta_p : X=p (uniform), Y=eta (moving)
            r.eta_p->RebinX(rebinp);
            TH2D* tmp = RebinTH2Y_varBins(r.eta_p, nEta, eEta);
            delete r.eta_p; r.eta_p = tmp;

            // ih_eta : X=eta (moving), Y=ih (uniform)
            r.ih_eta->RebinY(rebinih);
            tmp = RebinTH2X_varBins(r.ih_eta, nEta, eEta);
            delete r.ih_eta; r.ih_eta = tmp;

            // mass_eta : X=mass (uniform), Y=eta (moving)
            r.mass_eta->RebinX(1);
            tmp = RebinTH2Y_varBins(r.mass_eta, nEta, eEta);
            delete r.mass_eta; r.mass_eta = tmp;

            // pred_mass_eta : X=mass (uniform), Y=eta (moving)
            r.pred_mass_eta->RebinX(1);
            tmp = RebinTH2Y_varBins(r.pred_mass_eta, nEta, eEta);
            delete r.pred_mass_eta; r.pred_mass_eta = tmp;
        }
        else {
            r.eta_p->Rebin2D(rebinp, rebineta);
            r.ih_eta->Rebin2D(rebineta, rebinih);
            r.mass_eta->Rebin2D(1, rebineta);
            r.pred_mass_eta->Rebin2D(1, rebineta);
        }
    }

    // ---- Repliement |eta| ----
    // eta_p        : eta sur l'axe Y  -> FoldAbsTH2Y
    // ih_eta       : eta sur l'axe X  -> FoldAbsTH2X
    // mass_eta     : eta sur l'axe Y  -> FoldAbsTH2Y
    // pred_mass_eta: eta sur l'axe Y  -> FoldAbsTH2Y
    if (TakeAbsEta) {
        TH2D* tmp;

        tmp = FoldAbsTH2Y(r.eta_p, ("eta_1oP_"+regionName+"_absEta").c_str());
        delete r.eta_p; r.eta_p = tmp;

        tmp = FoldAbsTH2X(r.ih_eta, ("ih_eta_"+regionName+"_absEta").c_str());
        delete r.ih_eta; r.ih_eta = tmp;

        tmp = FoldAbsTH2Y(r.mass_eta, ("mass_eta_"+regionName+"_absEta").c_str());
        delete r.mass_eta; r.mass_eta = tmp;

        tmp = FoldAbsTH2Y(r.pred_mass_eta, ("pred_mass_eta_"+regionName+"_absEta").c_str());
        delete r.pred_mass_eta; r.pred_mass_eta = tmp;
    }

    std::cout << "1/p bin width = " << r.eta_p->GetXaxis()->GetBinWidth(1) << " GeV" << std::endl;

    return true;
}


template <typename TH>
TH* poissonHisto(const TH& h, TRandom3* RNG) {
    TH* hres = static_cast<TH*>(h.Clone());
    hres->SetDirectory(nullptr);                       // n'encombre pas le fichier de sortie
    const int n = hres->GetNcells();                   // under/overflow inclus, 1D comme 2D
    for (int b = 0; b < n; ++b) {
        const double mu = hres->GetBinContent(b);
        const double v  = (mu > 0.) ? RNG->Poisson(mu) : 0.;
        hres->SetBinContent(b, v);
        hres->SetBinError(b, std::sqrt(v));            // sinon l'erreur reste celle du clone
    }
    return hres;
}


// Function doing the eta reweighing between two 2D-histograms as done in the Hscp background estimate method,
// because of the correlation between variables (momentum & transverse momentum). 
// The first given 2D-histogram is weighted in respect to the 1D-histogram 
void etaReweighingP_Y(TH2D* eta_p_1, const TH1D* eta2_) {
    std::unique_ptr<TH1D> eta1(eta_p_1->ProjectionY());
    std::unique_ptr<TH1D> eta2(static_cast<TH1D*>(eta2_->Clone()));
    eta1->SetDirectory(nullptr);
    eta2->SetDirectory(nullptr);

    eta1->Scale(1./eta1->Integral(0,eta1->GetNbinsX()+1));
    eta2->Scale(1./eta2->Integral(0,eta2->GetNbinsX()+1));
    eta2->Divide(eta1.get());
    for(int i=0;i<eta_p_1->GetNbinsX()+2;i++)
    {
        for(int j=0;j<eta_p_1->GetNbinsY()+2;j++)
        {
            float val_ij = eta_p_1->GetBinContent(i,j);
            float err_ij = eta_p_1->GetBinError(i,j);
            
            eta_p_1->SetBinContent(i,j,val_ij*eta2->GetBinContent(j));
            eta_p_1->SetBinError(i,j,err_ij*eta2->GetBinContent(j));
        }
    }
}


// Same but for matching D -> reweighting = B*C/A
void etaReweighingP_X(TH2D* ih_eta_C, const TH1D* eta_B_, const TH1D* eta_A_) {
    std::unique_ptr<TH1D> eta_B(static_cast<TH1D*>(eta_B_->Clone()));
    std::unique_ptr<TH1D> eta_A(static_cast<TH1D*>(eta_A_->Clone()));
    eta_B->SetDirectory(nullptr); eta_A->SetDirectory(nullptr);

    eta_B->Scale(1./eta_B->Integral(0,eta_B->GetNbinsX()+1));
    eta_A->Scale(1./eta_A->Integral(0,eta_A->GetNbinsX()+1));
    eta_B->Divide(eta_A.get());

    for(int i=0; i<ih_eta_C->GetNbinsY()+2; i++)  // ih bins
    {
        for(int j=0; j<ih_eta_C->GetNbinsX()+2; j++)  // eta bins
        {
            ih_eta_C->SetBinContent(j, i, ih_eta_C->GetBinContent(j,i)*eta_B->GetBinContent(j));
            ih_eta_C->SetBinError(j, i, ih_eta_C->GetBinError(j,i)*eta_B->GetBinContent(j));
        }
    }
}

// Reponderation B/A appliquee a un TH2 dont eta est sur l'axe Y
void etaReweighingP_Y(TH2D* h, const TH1D* eta_B_, const TH1D* eta_A_) {
    std::unique_ptr<TH1D> eta_B(static_cast<TH1D*>(eta_B_->Clone()));
    std::unique_ptr<TH1D> eta_A(static_cast<TH1D*>(eta_A_->Clone()));
    eta_B->SetDirectory(nullptr);
    eta_A->SetDirectory(nullptr);

    eta_B->Scale(1./eta_B->Integral(0, eta_B->GetNbinsX()+1));
    eta_A->Scale(1./eta_A->Integral(0, eta_A->GetNbinsX()+1));
    eta_B->Divide(eta_A.get());

    for (int j = 0; j <= h->GetNbinsY()+1; ++j) {
        const double w = eta_B->GetBinContent(j);
        for (int i = 0; i <= h->GetNbinsX()+1; ++i) {
            h->SetBinContent(i, j, h->GetBinContent(i, j) * w);
            h->SetBinError  (i, j, h->GetBinError  (i, j) * w);
        }
    }
}

// add the overflow bin to the last one
void overflowLastBin(TH1D* h) {
    h->SetBinContent(h->GetNbinsX(),h->GetBinContent(h->GetNbinsX())+h->GetBinContent(h->GetNbinsX()+1));
    h->SetBinContent(h->GetNbinsX()+1,0);
}


// rebinning histogram according to an array of bins
TH1D* rebinHisto(TH1D* h) {
    static const double xbins[36]={0.,20.,40.,60.,80.,100.,120.,140.,160.,180.,200.,220.,240.,260.,280.,300.,
                      320.,340.,360.,380.,410.,440.,480.,530.,590.,660.,760.,880.,1030.,1210.,1440.,
                      1730.,2000.,2500.,3200.,4000.};
    constexpr int nb = static_cast<int>(std::size(xbins)) - 1;   // 35
    return (TH1D*) h->Rebin(nb, (std::string(h->GetName())+"_rebinned").c_str(), xbins);
}


// Function returning the ratio of right integer (from x to infty) for two 1D-histograms
// This function is used in the Hscp data-driven background estimate to test the mass shape prediction
// The argument to use this type of ratio is that we're in case of cut & count experiment 
TH1D* ratioIntegral(TH1D* h1, TH1D* h2) {    
    float SystError = systErr_;
    TH1D* res = (TH1D*) h1->Clone(); res->Reset();
    for(int i=1;i<h1->GetNbinsX()+1;i++)
    {   
        double Perr=0, Derr=0;
        double P=h1->IntegralAndError(i,h1->GetNbinsX()+1,Perr); if(P<=0) continue;
        double D=h2->IntegralAndError(i,h2->GetNbinsX()+1,Derr);
        Perr = sqrt(Perr*Perr + pow(P*SystError,2));
        res->SetBinContent(i,D/P);
        res->SetBinError(i,sqrt(pow(Derr*P,2)+pow(Perr*D,2))/pow(P,2));
    }
    return res;
}



void saveHistoRatio(TH1D* h1,TH1D* h2,std::string st1,std::string st2,std::string st3,bool rebin=false) {
    h1->SetName(st1.c_str());
    h2->SetName(st2.c_str());
    if(rebin){
        h1 = rebinHisto(h1);
        h2 = rebinHisto(h2);
    }
    h1->Write();
    h2->Write();
    std::unique_ptr<TH1D> R(ratioIntegral(h2, h1));
    R->SetDirectory(nullptr);
    if (rebin) st3 += "_rebinned";
    R->SetName(st3.c_str());
    R->Write();
}

struct BckgOptions {
    int    nPE             = 200;
    bool   useFit          = true;
    bool   useOldIhFit     = false;
    bool   useOld1oPFit    = true;
    bool   corrTemplateIh  = false;
    bool   corrTemplate1oP = false;
    std::string etaName    = "";
    bool   saveFits        = false;
    int    fitIh           = 1;
    int    fitP            = 1;
    bool   blind           = false;
    unsigned nWorkers      = 25;
};

bool bckgEstimate(const std::string& filename, 
                  const DeDxCalib& calib,
                  const Region& B,
                  const Region& C,
                  const Region& BC,
                  const Region& A,
                  const Region& D,
                  bool ifIhpSAME,
                  const Region& B_ifIhpSAME,
                  const std::string& st,
                  const BckgOptions& o) {

    // Deballage local : le corps et le std::bind restent inchanges.
    const int         nPE             = o.nPE;
    const bool        useFit          = o.useFit;
    const bool        useOldIhFit     = o.useOldIhFit;
    const bool        useOld1oPFit    = o.useOld1oPFit;
    const bool        corrTemplateIh  = o.corrTemplateIh;
    const bool        corrTemplate1oP = o.corrTemplate1oP;
    const std::string etaName         = o.etaName;
    const bool        saveFits        = o.saveFits;
    const int         fitIh           = o.fitIh;
    const int         fitP            = o.fitP;
    const bool        blind           = o.blind;
    const unsigned      nWorkers      = o.nWorkers;

    ROOT::TProcessExecutor workers(nWorkers);

    Region c = C;
    Region bc = BC;
    Region d = D;

    TH2D c_eta_p_base(*c.eta_p);

    // Pre-fit of 1/p to get the parameters for the next fits in the toys
    TH1D* p_base = (TH1D*)c_eta_p_base.ProjectionX();
    p_base->SetDirectory(nullptr);

    float rangemax_p = 30;
    if (p_base->GetBinCenter(p_base->GetMaximumBin()) < rangemax_p) rangemax_p = 0.8 * p_base->GetBinCenter(p_base->GetMaximumBin());
    
    TF1 f_p_base("f_p_base","[0]*([1]+erf((log(x)-[2])/[3]))",0,rangemax_p);
    f_p_base.SetParameter(0,560);
    f_p_base.FixParameter(1,1.0);
    f_p_base.SetParameter(2,3.50116e+00);
    f_p_base.SetParameter(3,0.60152e+00);

    p_base->Fit(&f_p_base, "QRSL", "", 0, rangemax_p);

    double par_p2 = f_p_base.GetParameter(2);
    double par_p3 = f_p_base.GetParameter(3);

    delete p_base;


    // Toys lambda function
    auto workItem = [] (UInt_t workerID, const std::string& filename, 
                        const Region& B, const Region& C, const Region& BC, 
                        const Region& A, bool ifIhpSAME, 
                        const Region& B_ifIhpSAME, const std::string& st,
                        const DeDxCalib calib,
                        const bool useFit = true, const bool useOldIhFit = false, const bool useOld1oPFit = false,
                        const bool corrTemplateIh = false, const bool corrTemplate1oP = false, const std::string& etaName = "", const bool& saveFits = false, 
                        const int& fitIh = 1, const int& fitP = 1,
                        const double& par_p2 = 4.70839, const double& par_p3 = 1.05005)
                        -> std::tuple<TH1D, TH2D, double, TH2D, TH2D> {
        gROOT->cd();
        
        // Setup
        Region a = A;
        Region b = B;
        Region b_ifIhpSAME = B_ifIhpSAME;
        Region c = C;
        Region bc = BC;

        TH2D a_ih_eta_base(*a.ih_eta);
        TH1D* a_eta_base = (TH1D*)a_ih_eta_base.ProjectionX();
        TH2D b_ih_eta_base(*b.ih_eta);
        TH2D b_ifIhpSAME_ih_eta_base(*b_ifIhpSAME.ih_eta);
        TH2D b_ifIhpSAME_eta_p_base(*b_ifIhpSAME.eta_p);
        TH1D* b_ifIhpSAME_eta_base = (TH1D*)b_ifIhpSAME_eta_p_base.ProjectionY();
        TH2D b_eta_p_base(*b.eta_p);
        TH1D* b_eta_base = (TH1D*)b_eta_p_base.ProjectionY();
        TH2D c_eta_p_base(*c.eta_p);

        if (corrTemplateIh) corrIh(&b_ih_eta_base, etaName);
        if (corrTemplate1oP) corr1oP(&c_eta_p_base, etaName);
        

        // 1/p fit
        TH1D* p_base = (TH1D*)c_eta_p_base.ProjectionX();
        float rangemax_p = 30;
        if (p_base->GetBinCenter(p_base->GetMaximumBin()) < rangemax_p) rangemax_p = 0.8 * p_base->GetBinCenter(p_base->GetMaximumBin());
        
        TF1 f_p("f_p", useOld1oPFit ? "[0]*([1]+erf((log(x)-[2])/[3]))" : "0.5*(exp([0]*x*x+[1]*x)+exp(-[0]*x*x-[1]*x))-1", 0, rangemax_p);
        if (useOld1oPFit) {
            f_p.SetParLimits(0, 0, 1);
            f_p.FixParameter(1, 1.0);
            f_p.SetParameter(2, par_p2);
            f_p.SetParameter(3, par_p3);
        }
        else {
            f_p.SetParLimits(0, 0, 1e-3);
            f_p.SetParLimits(1, -1e-2, 1e-2);
        }

        // Ih fit
        TH1D* ih_base = (TH1D*)b_ih_eta_base.ProjectionX();
        float max_ih = ih_base->GetBinCenter(ih_base->GetMaximumBin());

        TF1 f_ihg("f_ihg", "gaus", max_ih, 8);
        f_ihg.SetParameter(0, 0.01*ih_base->Integral());
        f_ihg.SetParameter(1, max_ih);
        f_ihg.SetParameter(2, ih_base->GetStdDev());


        // Toys
        TRandom3* RNG = new TRandom3(kSeedBase + workerID);
        bc.pred_mass->Reset();
        bc.pred_mass_eta->Reset();

        TH2D* a_ih_eta = poissonHisto(a_ih_eta_base, RNG);
        TH2D* b_ih_eta = poissonHisto(b_ih_eta_base, RNG);
        TH2D* b_ifIhpSAME_ih_eta = poissonHisto(b_ifIhpSAME_ih_eta_base, RNG);
        TH2D* b_eta_p = poissonHisto(b_eta_p_base, RNG);
        TH2D* c_eta_p = poissonHisto(c_eta_p_base, RNG);
        
        TH1D* b_eta = poissonHisto(*b_eta_base, RNG);
        TH1D* b_ifIhpSAME_eta = poissonHisto(*b_ifIhpSAME_eta_base, RNG);
        TH1D* a_eta = poissonHisto(*a_eta_base, RNG);
        
        if(ifIhpSAME) etaReweighingP_X(b_ih_eta, b_ifIhpSAME_eta, a_eta);
        else etaReweighingP_Y(c_eta_p,b_eta);

        etaReweighingP_Y(b_eta_p, b_ifIhpSAME_eta, a_eta);


        // Mass prediction in the BC region
        bc.eta_p = c_eta_p;
        bc.ih_eta = b_ih_eta;
        bc.fillPredMass(filename, st, calib, f_p, f_ihg, useFit, fitIh, fitP, 
                        useOldIhFit, useOld1oPFit, etaName, saveFits,
                        par_p2, par_p3, workerID);
        
        double normA = a_ih_eta->Integral(0, a_ih_eta->GetNbinsX()+1, 0, a_ih_eta->GetNbinsY()+1);
        double normB = b_ih_eta->Integral(0, b_ih_eta->GetNbinsX()+1, 0, b_ih_eta->GetNbinsY()+1);
        double normC = c_eta_p->Integral(0, c_eta_p->GetNbinsX()+1, 0, c_eta_p->GetNbinsY()+1);
        double normB_ifIhpSAME = b_ifIhpSAME_ih_eta->Integral(0, b_ifIhpSAME_ih_eta->GetNbinsX()+1, 0, b_ifIhpSAME_ih_eta->GetNbinsY()+1);

        double normalisationABC = normB * normC / normA;
        if (ifIhpSAME) normalisationABC = normB_ifIhpSAME * normC / normA;


        // TEMPORARY
        //normalisationABC = bc.pred_mass->Integral();
        // TEMPORARY
        

        const double itg = bc.pred_mass->Integral();
        const double itg2D = bc.pred_mass_eta->Integral();
        if (itg > 0) bc.pred_mass->Scale(normalisationABC/itg);
        else std::cerr << "toy " << workerID << " : empty prediction" << std::endl;
        if (itg2D > 0) bc.pred_mass_eta->Scale(normalisationABC/itg2D);
        else std::cerr << "toy " << workerID << " : empty prediction" << std::endl;

        TH1D out_pred_mass    (*bc.pred_mass);
        TH2D out_pred_mass_eta(*bc.pred_mass_eta);
        TH2D out_ih_eta       (*bc.ih_eta);
        TH2D out_p_eta    (*b_eta_p);
        out_pred_mass.SetDirectory(nullptr);
        out_pred_mass_eta.SetDirectory(nullptr);
        out_ih_eta.SetDirectory(nullptr);
        out_p_eta.SetDirectory(nullptr);

        // End
        delete a_eta;              delete a_ih_eta;           delete a_eta_base;
        delete b_ih_eta;           delete b_eta_p;            delete b_eta;
        delete b_ifIhpSAME_eta;    delete b_ifIhpSAME_ih_eta; delete b_eta_base;
        delete b_ifIhpSAME_eta_base;
        delete c_eta_p;
        delete p_base;             delete ih_base;
        bc.ih_eta = nullptr;
        bc.eta_p  = nullptr;
        delete RNG;

        return {out_pred_mass, out_pred_mass_eta, normalisationABC, out_ih_eta, out_p_eta};
    };

    
    // Loop on the toys
    auto workItemToRun = std::bind (workItem, _1, filename, B, C, BC, A, ifIhpSAME, B_ifIhpSAME, st, calib,
                                    useFit, useOldIhFit, useOld1oPFit, corrTemplateIh, corrTemplate1oP, etaName, saveFits, fitIh, fitP, 
                                    par_p2, par_p3);
    
    auto vPE = workers.Map(workItemToRun, ROOT::TSeqI(nPE));
    if (vPE.size() != (size_t)nPE) {
        std::cerr << "ERREUR: " << vPE.size() << "/" << nPE
                  << " toys revenus — des workers ont crashe" << std::endl;
        return false;
    }


    // Get the results
    std::vector<TH1D> histo_pred_mass;
    std::vector<TH2D> histo_pred_mass_eta;
    std::vector<double> normalisations;
    std::vector<TH1D> Ih_eta, oP_eta;

    for (const auto& result : vPE) {
        histo_pred_mass.push_back(std::get<0>(result));   // Récupère *bc.pred_mass
        histo_pred_mass_eta.push_back(std::get<1>(result)); // Récupère *bc.pred_mass_eta
        normalisations.push_back(std::get<2>(result)); // Récupère normalisationABC
        
        const TH2D& h2 = std::get<3>(result);
        std::unique_ptr<TH1D> proj(h2.ProjectionY(Form("ih_eta_py_%zu", Ih_eta.size())));
        proj->SetDirectory(nullptr);
        Ih_eta.push_back(*proj);
        Ih_eta.back().SetDirectory(nullptr);

        const TH2D& h3 = std::get<4>(result);
        std::unique_ptr<TH1D> proj2(h3.ProjectionX(Form("p_eta_px_%zu", oP_eta.size())));
        proj2->SetDirectory(nullptr);
        oP_eta.push_back(*proj2);
    }

    TH1D h_ih_eta_mean = MeanOfToys(Ih_eta, ("ih_eta_mean_"+st).c_str());
    TH1D h_oP_eta_mean = MeanOfToys(oP_eta, ("oP_eta_mean_"+st).c_str());
    TH1D h_temp        = MeanOfToys(histo_pred_mass,     ("pred_mass_mean_"+st).c_str());
    TH2D h_temp_eta    = MeanOfToys(histo_pred_mass_eta, ("pred_mass_eta_mean_"+st).c_str());
        
    bc.pred_mass = &h_temp;
    bc.pred_mass_eta = &h_temp_eta;
    double avgNormalisation = std::accumulate(normalisations.begin(), normalisations.end(), 0.0) / normalisations.size();


    // Case Ih cut !!!
    TH1D* d_mass = (TH1D*) d.mass->Clone(("mass_obs_"+st).c_str());
    d_mass->SetDirectory(nullptr);
    
    if(blind) blindMass(d_mass,300);
    overflowLastBin(d_mass);
    overflowLastBin(bc.pred_mass);


    // Saving histograms
    saveHistoRatio(d_mass, bc.pred_mass, ("mass_obs_"+st).c_str(), ("mass_predBC_"+st).c_str(), ("mass_predBCR_"+st).c_str());

    bc.pred_mass_eta->Write();
    bc.pred_mass_eta->ProjectionY()->Write();

    d.mass_eta->SetName(("mass_eta_D_"+st).c_str());        d.mass_eta->Write();
    d.mass_eta->ProjectionY(("mass_eta_D_"+st+"_py").c_str())->Write();
    c.mass_eta->SetName(("mass_eta_C_"+st).c_str());        c.mass_eta->Write();
    c.mass_eta->ProjectionY(("mass_eta_C_"+st+"_py").c_str())->Write();

    A.ih_eta->Write();
    A.eta_p->Write();
    B.ih_eta->Write();
    C.eta_p->Write();
    D.eta_p->Write();
    D.ih_eta->Write();
    
    D.ih_eta->ProjectionY()->Write();
    B.ih_eta->ProjectionY()->Write();
    D.ih_eta->ProjectionX()->Write();
    B.ih_eta->ProjectionX()->Write();
    D.eta_p->ProjectionX()->Write();
    C.eta_p->ProjectionX()->Write();   
    C.eta_p->ProjectionY()->Write();

    std::unique_ptr<TH1D> cMassScaled(static_cast<TH1D*>(C.mass->Clone(("mass_C_"+st).c_str())));
    cMassScaled->SetDirectory(nullptr);
    const double itgC = cMassScaled->Integral();
    if (itgC > 0) cMassScaled->Scale(avgNormalisation/itgC);
    cMassScaled->Write();

    // Draw on one canvas all the histos of the vector histo_pred_mass
    TCanvas* c1  = new TCanvas(("toys_"+st).c_str(), "toys", 800, 800);
    c1->cd();
    for (size_t i=0; i<histo_pred_mass.size(); i++) histo_pred_mass[i].Rebin(10);
    histo_pred_mass[0].Draw("hist");
    for (size_t i=1; i<histo_pred_mass.size(); i++) histo_pred_mass[i].Draw("hist same");
    c1->Write();

    // fill a new histo with norm/d.mass->integral:
    TH1D* h_norm = new TH1D(("h_norm_"+st).c_str(), "h_norm", 40, 0.9, 1.1);
    h_norm->GetXaxis()->SetTitle("BC/AD");
    h_norm->GetYaxis()->SetTitle("Entries");
    for(size_t i=0; i<normalisations.size(); i++) h_norm->Fill(normalisations[i]/d.mass->Integral());
    h_norm->Write();

    h_ih_eta_mean.Write();
    TH1D* ih_eta_VR = D.ih_eta->ProjectionY();
    ih_eta_VR->SetName(("ih_VR_"+st).c_str());
    ih_eta_VR->Write();
    h_oP_eta_mean.Write();
    TH1D* eta_p_VR = D.eta_p->ProjectionX();
    eta_p_VR->SetName(("eta_VR_"+st).c_str());
    eta_p_VR->Write();

    return true;
}