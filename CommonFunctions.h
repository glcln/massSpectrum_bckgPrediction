#pragma once

// =============================================================================
//  CommonFunctions.h
// -----------------------------------------------------------------------------
//  Everything that sits between the step1 histograms and the final background
//  prediction:
//
//    * histogram I/O helpers (safe TH1F/TH2F -> TH1D/TH2D conversion);
//    * eta rebinning tables and variable-binning rebinners;
//    * |eta| folding;
//    * loadHistograms(), which reads one ABCD region from a step1 file and
//      applies the rebinning / folding requested by the configuration;
//    * the pseudo-experiment machinery (Poisson smearing, toy averaging);
//    * the eta reweighting used to correct for the p / pt correlation;
//    * bckgEstimate(), the top-level driver that runs the toys, averages them,
//      normalises the prediction with the ABCD factor and writes everything out.
//
//  Axis conventions (must match step1 and Regions.h):
//    eta_p        : X = 1/p (scaled by kInvPScale), Y = eta
//    ih_eta       : X = eta,                        Y = Ih
//    mass_eta     : X = mass,                       Y = eta
//    pred_mass_eta: X = mass,                       Y = eta
// =============================================================================

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


// Base of the random seeds: toy n uses seed kSeedBase + n. Deriving the seed from
// the worker index (rather than from the clock or from gRandom) is what makes a
// run bit-for-bit reproducible, which matters when comparing systematic variations.
constexpr UInt_t kSeedBase = 1;

// ---------------------------------------------------------------------------
// Eta rebinning tables.
//
// Three granularities are provided per eta range, used as the "eta binning"
// systematic variation:
//    Down (rebineta == 8): coarse
//    Nom  (rebineta == 4): nominal
//    Up   (else)         : fine
// The tables are edge arrays, so the bins do not have to be uniform: they are
// deliberately narrower at low |eta|, where the dE/dx response varies fastest.
// The "Eta1p2_*" tables have a single wide central bin covering |eta| < 1.2,
// because those regions only use the endcap part of the detector.
// ---------------------------------------------------------------------------
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

// ---------------------------------------------------------------------------
// Selects the eta binning table matching a region name and a granularity code.
//
// CRITICAL: the order of the tests matters, because they use substring matching.
// "Eta1p2_2p2" and "Eta1_2p4" both contain "Eta1", so the most specific names
// must be tested first; reordering these blocks silently makes some branches
// unreachable.
//
// Returns nullptr when no pattern matches, which the caller must treat as a hard
// failure (the region is then not loaded at all).
// ---------------------------------------------------------------------------
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

// ---------------------------------------------------------------------------
// Reads a 1D histogram from a step1 file and returns it as a TH1D.
//
// step1 may write TH1F: a direct TH1F* -> TH1D* cast is undefined behaviour and
// produces garbage silently, hence the explicit element-by-element conversion.
// The returned histogram is always a fresh object detached from any directory,
// so the caller owns it and must delete it.
//
// The conversion preserves the binning (uniform or variable), the bin errors,
// the under/overflow bins (loop over GetNcells) and the entry count.
// Returns nullptr when the object is missing or is not a 1D histogram.
// ---------------------------------------------------------------------------
TH1D* GetAsTH1D(TFile* f, const std::string& name) {
    TH1* h = dynamic_cast<TH1*>(f->Get(name.c_str()));
    // TH2 inherits from TH1, so an explicit exclusion is required here.
    if (!h || h->InheritsFrom(TH2::Class())) {
        std::cerr << "GetAsTH1D: '" << name << "' not in file, or is not a TH1" << std::endl;
        return nullptr;
    }
    // Already a TH1D: a plain clone is enough.
    if (h->InheritsFrom(TH1D::Class())) {
        TH1D* copy = static_cast<TH1D*>(h->Clone((name + "_copy").c_str()));
        copy->SetDirectory(nullptr);
        return copy;
    }

    // Rebuild the axis, keeping a variable binning if there is one.
    const TAxis* ax = h->GetXaxis();
    TH1D* out = nullptr;
    if (ax->GetXbins()->GetSize() > 0)
        out = new TH1D((name + "_copy").c_str(), h->GetTitle(), ax->GetNbins(), ax->GetXbins()->GetArray());
    else
        out = new TH1D((name + "_copy").c_str(), h->GetTitle(), ax->GetNbins(), ax->GetXmin(), ax->GetXmax());

    out->SetDirectory(nullptr);
    out->Sumw2();
    // GetNcells covers bins 0..N+1, i.e. under/overflow are copied too.
    for (int b = 0; b < h->GetNcells(); ++b) {
        out->SetBinContent(b, h->GetBinContent(b));
        out->SetBinError  (b, h->GetBinError(b));
    }
    out->SetEntries(h->GetEntries());
    return out;
}

// ---------------------------------------------------------------------------
// 2D counterpart of GetAsTH1D. Same rationale: step1 may write TH2F and a direct
// TH2F* -> TH2D* cast is undefined behaviour.
//
// Variable binning is explicitly rejected here rather than silently mishandled,
// because step1 only ever writes uniform 2D binnings; a variable-binned input
// would signal that something upstream changed.
// ---------------------------------------------------------------------------
TH2D* GetAsTH2D(TFile* f, const std::string& name) {
    TH2* h = dynamic_cast<TH2*>(f->Get(name.c_str()));
    if (!h) {
        std::cerr << "GetAsTH2D: '" << name << "' not in file, or is not a TH2" << std::endl;
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
        std::cerr << "GetAsTH2D: variable binning not supported for '" << name << "'" << std::endl;
        return nullptr;   // step1 only writes uniform binnings
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


// ---------------------------------------------------------------------------
// Rebins the Y axis of a TH2D onto an arbitrary (variable) binning.
// ROOT's Rebin2D only supports integer merging factors, which cannot produce the
// eta tables above, hence this manual implementation.
//
// The X axis is copied as is (variable binning preserved if present). Contents
// are summed and errors added in quadrature. Under/overflow bins on both axes are
// carried along: the loops run 0..N+1 and FindBin maps them onto the new axis.
//
// The returned histogram is detached from any directory; the caller owns it.
// ---------------------------------------------------------------------------
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
            // The old bin is assigned to the new bin containing its centre.
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

// X-axis counterpart of RebinTH2Y_varBins, used for ih_eta where eta is on X.
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

// ---------------------------------------------------------------------------
// Folds a TH2D whose X axis is eta onto |eta|: bins at -eta and +eta are merged.
//
// This doubles the statistics per eta bin at the price of assuming a symmetric
// detector response. It requires the binning to be symmetric around 0 and to have
// a bin boundary exactly at eta = 0 (all the tables above satisfy this).
//
// The +1e-9 offset in FindBin(0.0 + 1e-9) makes sure the bin *starting* at 0 is
// picked, not the one ending there.
//
// NB: under/overflow of the X axis are not transferred to the folded histogram
// (the loop runs 1..nx), which is intended since they lie outside |eta| < 2.4.
// ---------------------------------------------------------------------------
TH2D* FoldAbsTH2X(TH2D* h, const std::string& newName) {
    int nx = h->GetNbinsX();
    int ny = h->GetNbinsY();
    const TAxis* ax = h->GetXaxis();

    // Index of the first bin whose lower edge is >= 0 (the eta = 0 boundary)
    int izero = ax->FindBin(0.0 + 1e-9);
    // Number of positive-side bins
    int nxPos = nx - izero + 1;

    // New, positive-only X binning; the Y binning is copied unchanged.
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
            int io = hf->GetXaxis()->FindBin(std::fabs(xc)); // target positive bin
            double c = hf->GetBinContent(io, j) + h->GetBinContent(i, j);
            // hypot = quadratic sum, i.e. errors added in quadrature.
            double e = std::hypot(hf->GetBinError(io, j), h->GetBinError(i, j));
            hf->SetBinContent(io, j, c);
            hf->SetBinError(io, j, e);
        }
    }
    return hf;
}

// Y-axis counterpart of FoldAbsTH2X, used for eta_p / mass_eta / pred_mass_eta
// where eta sits on the Y axis.
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

// ---------------------------------------------------------------------------
// MeanOfToys: bin-by-bin mean over a set of pseudo-experiments.
//
// Templated on the histogram type so that the same code serves TH1D and TH2D
// (GetNcells covers both, under/overflow included).
//
// The bin error is the *spread* of the toys, not a Poisson error:
//   - useSEM == false: sample standard deviation, i.e. the statistical
//     uncertainty of the prediction as estimated by the toys;
//   - useSEM == true : standard error on the mean (sigma / sqrt(N)), i.e. how
//     well the toys determine the central value.
//
// SetBinErrorOption(kNormal) is essential: if the input histograms carry Poisson
// (Garwood) errors, ROOT would recompute the errors from the contents on read-out
// and contaminate the RMS-based uncertainty.
//
// Reset("ICESM") clears contents, errors, statistics and the sum of weights while
// keeping the axes, so the mean histogram inherits the exact binning of the toys.
// ---------------------------------------------------------------------------
template <typename TH>
TH MeanOfToys(const std::vector<TH>& toys, const char* name, bool useSEM = false) {
    if (toys.empty()) throw std::runtime_error("MeanOfToys: empty vector");
    const double N = static_cast<double>(toys.size());

    TH hMean(toys.front());
    hMean.SetName(name);
    hMean.SetTitle(name);
    hMean.Reset("ICESM");
    hMean.SetBinErrorOption(TH1::EBinErrorOpt::kNormal);
    hMean.SetDirectory(nullptr);
    hMean.Sumw2();

    const int nTot = hMean.GetNcells();          // covers 1D and 2D, under/overflow included
    for (int b = 0; b < nTot; ++b) {
        double sum = 0., sum2 = 0.;
        for (const auto& h: toys) {
            const double v = h.GetBinContent(b);
            sum += v; sum2 += v*v;
        }
        const double mean = sum / N;
        // Unbiased variance; a single toy gives no spread information.
        double var = (N > 1) ? (sum2 - N*mean*mean) / (N - 1.): 0.;
        // Guard against tiny negative values from floating-point cancellation.
        if (var < 0.) var = 0.;
        hMean.SetBinContent(b, mean);
        hMean.SetBinError(b, std::sqrt(var) / (useSEM ? std::sqrt(N): 1.));
    }
    return hMean;
}


// ---------------------------------------------------------------------------
// loadHistograms: reads one ABCD region from a step1 output file and prepares it.
//
// Steps:
//   1. read eta_1oP, ih_eta, mass and mass_eta with the safe TH*F -> TH*D helpers;
//   2. book the (empty) prediction histograms by cloning the observed ones, so
//      they automatically share the same binning;
//   3. apply the requested rebinning, either with the variable eta tables
//      (rebineta in {2, 4, 8}) or with a plain integer Rebin2D;
//   4. optionally fold onto |eta|.
//
// Returns false when anything is missing, so the caller can abort cleanly instead
// of dereferencing a null histogram later on.
//
// Arguments:
//   bool_rebin: master switch for step 3
//   rebineta  : eta granularity code (see PickEtaBinning)
//   rebinp    : 1/p granularity code, remapped just below
//   rebinih   : integer merging factor on the Ih axis
//   TakeAbsEta: fold onto |eta|
// ---------------------------------------------------------------------------
bool loadHistograms(Region& r, 
                    TFile* f,
                    const std::string& regionName,
                    bool bool_rebin = true,
                    int rebineta = 1,
                    int rebinp = 1,
                    int rebinih = 1,
                    bool TakeAbsEta = false) {

    std::cout << "loading region " << regionName << "    rebineta=" << rebineta << ", rebinp=" << rebinp << ", rebinih=" << rebinih << std::endl;

    // The configuration exposes rebinp as a Down/Nom/Up code (4/2/1); it is
    // translated here into the actual merging factor applied to the 1/p axis.
    // The order matters: 4 must be handled before it would be re-matched.
    if (rebinp==4) rebinp = 8;
    if (rebinp==2) rebinp = 6;
    if (rebinp==1) rebinp = 4;

    r.eta_p    = GetAsTH2D(f, "eta_1oP_"  + regionName);
    r.ih_eta   = GetAsTH2D(f, "ih_eta_"   + regionName);
    r.mass     = GetAsTH1D(f, "mass_"     + regionName);
    r.mass_eta = GetAsTH2D(f, "mass_eta_" + regionName);

    if (!r.eta_p || !r.ih_eta || !r.mass || !r.mass_eta) {
        std::cerr << "loadHistograms: region '" << regionName << "' incomplete -> skipping region" << std::endl;
        return false;
    }

    // Prediction histograms: same binning as the observed ones, emptied.
    r.pred_mass     = (TH1D*) r.mass->Clone();
    r.pred_mass->SetDirectory(nullptr);
    r.pred_mass->SetName(("pred_mass_"+regionName).c_str());
    r.pred_mass->Reset();

    r.pred_mass_eta = (TH2D*) r.mass_eta->Clone();
    r.pred_mass_eta->SetDirectory(nullptr);
    r.pred_mass_eta->SetName(("pred_mass_eta_"+regionName).c_str());
    r.pred_mass_eta->Reset();

    if (bool_rebin) {

        // --- Variable eta binning ------------------------------------------
        if (rebineta==2 || rebineta==4 || rebineta==8) {

            const std::vector<double>* RebinEtaVecPtr = PickEtaBinning(regionName, rebineta);
            if (!RebinEtaVecPtr) {
                std::cerr << "loadHistograms: no eta binning eta for '" << regionName
                        << "' -> region not loaded" << std::endl;
                return false;
            }
            const int     nEta = static_cast<int>(RebinEtaVecPtr->size()) - 1;
            const double* eEta = RebinEtaVecPtr->data();

            // For each histogram: first the integer rebinning of the non-eta axis,
            // then the variable rebinning of the eta axis. The temporary is swapped
            // in and the old histogram deleted, since the rebinners return new
            // objects rather than modifying in place.

            // eta_p: X=p (uniform), Y=eta (moving)
            r.eta_p->RebinX(rebinp);
            TH2D* tmp = RebinTH2Y_varBins(r.eta_p, nEta, eEta);
            delete r.eta_p; r.eta_p = tmp;

            // ih_eta: X=eta (moving), Y=ih (uniform)
            r.ih_eta->RebinY(rebinih);
            tmp = RebinTH2X_varBins(r.ih_eta, nEta, eEta);
            delete r.ih_eta; r.ih_eta = tmp;

            // mass_eta: X=mass (uniform), Y=eta (moving)
            // RebinX(1) is a no-op kept for symmetry: the mass binning is never
            // merged here, it is handled at plotting time.
            r.mass_eta->RebinX(1);
            tmp = RebinTH2Y_varBins(r.mass_eta, nEta, eEta);
            delete r.mass_eta; r.mass_eta = tmp;

            // pred_mass_eta: X=mass (uniform), Y=eta (moving)
            // Must follow mass_eta exactly, otherwise the prediction and the
            // observation would no longer be comparable bin by bin.
            r.pred_mass_eta->RebinX(1);
            tmp = RebinTH2Y_varBins(r.pred_mass_eta, nEta, eEta);
            delete r.pred_mass_eta; r.pred_mass_eta = tmp;
        }
        // --- Plain uniform rebinning ---------------------------------------
        // Note the transposed argument order for ih_eta, whose eta axis is X.
        else {
            r.eta_p->Rebin2D(rebinp, rebineta);
            r.ih_eta->Rebin2D(rebineta, rebinih);
            r.mass_eta->Rebin2D(1, rebineta);
            r.pred_mass_eta->Rebin2D(1, rebineta);
        }
    }

    // ---- |eta| folding ----
    // eta_p        : eta on the Y axis -> FoldAbsTH2Y
    // ih_eta       : eta on the X axis -> FoldAbsTH2X
    // mass_eta     : eta on the Y axis -> FoldAbsTH2Y
    // pred_mass_eta: eta on the Y axis -> FoldAbsTH2Y
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

    // Useful sanity print: the 1/p bin width directly drives the mass resolution
    // of the convolution.
    std::cout << "1/p bin width = " << r.eta_p->GetXaxis()->GetBinWidth(1) << " (10^4/GeV units)" << std::endl;

    return true;
}


// ---------------------------------------------------------------------------
// poissonHisto: one pseudo-experiment.
//
// Replaces every bin content mu by a Poisson random draw of mean mu, and sets the
// error to sqrt(v) accordingly (the clone would otherwise keep the errors of the
// original histogram, which no longer describe the smeared content).
//
// GetNcells is used so that under/overflow are smeared too. Empty or negative
// bins are set to 0 rather than passed to TRandom3::Poisson.
//
// Templated on the histogram type: the same function serves TH1D and TH2D.
// The returned histogram is detached and owned by the caller.
// ---------------------------------------------------------------------------
template <typename TH>
TH* poissonHisto(const TH& h, TRandom3* RNG) {
    TH* hres = static_cast<TH*>(h.Clone());
    hres->SetDirectory(nullptr);                       // keeps the output file clean
    const int n = hres->GetNcells();                   // under/overflow included, 1D and 2D alike
    for (int b = 0; b < n; ++b) {
        const double mu = hres->GetBinContent(b);
        const double v  = (mu > 0.) ? RNG->Poisson(mu): 0.;
        hres->SetBinContent(b, v);
        hres->SetBinError(b, std::sqrt(v));            // otherwise the clone's error would survive
    }
    return hres;
}


// Function doing the eta reweighing between two 2D-histograms as done in the Hscp background estimate method,
// because of the correlation between variables (momentum & transverse momentum). 
// The first given 2D-histogram is weighted in respect to the 1D-histogram 
//
// Overload 1 (target eta distribution given directly): the eta profile of h is
// reweighted so that it matches eta2_. Both distributions are normalised to unit
// area first, so the weight is a pure shape correction and the overall yield of
// eta_p_1 is preserved up to shape effects.
//
// NOTE: as used in this package, this overload only acts on a debug control plot,
// not on the prediction template itself.
void etaReweighingP_Y(TH2D* eta_p_1, const TH1D* eta2_) {
    std::unique_ptr<TH1D> eta1(eta_p_1->ProjectionY());
    std::unique_ptr<TH1D> eta2(static_cast<TH1D*>(eta2_->Clone()));
    eta1->SetDirectory(nullptr);
    eta2->SetDirectory(nullptr);

    // Normalise both eta profiles, then form the ratio target/current.
    eta1->Scale(1./eta1->Integral(0,eta1->GetNbinsX()+1));
    eta2->Scale(1./eta2->Integral(0,eta2->GetNbinsX()+1));
    eta2->Divide(eta1.get());
    // Apply one weight per eta row (index j), under/overflow included.
    for(int i=0;i<eta_p_1->GetNbinsX()+2;i++)
    {
        for(int j=0;j<eta_p_1->GetNbinsY()+2;j++)
        {
            float val_ij = eta_p_1->GetBinContent(i,j);
            float err_ij = eta_p_1->GetBinError(i,j);
            
            // Content and error scaled by the same factor: the relative
            // statistical uncertainty of each bin is left unchanged.
            eta_p_1->SetBinContent(i,j,val_ij*eta2->GetBinContent(j));
            eta_p_1->SetBinError(i,j,err_ij*eta2->GetBinContent(j));
        }
    }
}


// Same but for matching D -> reweighting = B*C/A
//
// Overload for a TH2 with eta on the X axis (ih_eta). The weight is the ratio of
// the normalised eta profiles B/A, which is what makes the Ih template of the
// prediction match the eta composition expected in the signal region.
void etaReweighingP_X(TH2D* ih_eta_C, const TH1D* eta_B_, const TH1D* eta_A_) {
    std::unique_ptr<TH1D> eta_B(static_cast<TH1D*>(eta_B_->Clone()));
    std::unique_ptr<TH1D> eta_A(static_cast<TH1D*>(eta_A_->Clone()));
    eta_B->SetDirectory(nullptr); eta_A->SetDirectory(nullptr);

    eta_B->Scale(1./eta_B->Integral(0,eta_B->GetNbinsX()+1));
    eta_A->Scale(1./eta_A->Integral(0,eta_A->GetNbinsX()+1));
    eta_B->Divide(eta_A.get());   // eta_B now holds the B/A weight per eta bin

    // Outer loop over Ih (Y), inner loop over eta (X): the weight index is j,
    // the eta bin. Under/overflow included on both axes.
    for(int i=0; i<ih_eta_C->GetNbinsY()+2; i++)  // ih bins
    {
        for(int j=0; j<ih_eta_C->GetNbinsX()+2; j++)  // eta bins
        {
            ih_eta_C->SetBinContent(j, i, ih_eta_C->GetBinContent(j,i)*eta_B->GetBinContent(j));
            ih_eta_C->SetBinError(j, i, ih_eta_C->GetBinError(j,i)*eta_B->GetBinContent(j));
        }
    }
}

// B/A reweighting applied to a TH2 whose eta axis is Y (eta_p).
// Same weight definition as etaReweighingP_X, transposed. The weight is computed
// once per eta row and reused across the whole row.
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

// Moves the overflow content into the last visible bin and empties the overflow.
// Required before the integral-ratio test below, which only sums bins 1..N+1 of
// the *displayed* range; without this, the highest-mass events would be lost.
// add the overflow bin to the last one
void overflowLastBin(TH1D* h) {
    h->SetBinContent(h->GetNbinsX(),h->GetBinContent(h->GetNbinsX())+h->GetBinContent(h->GetNbinsX()+1));
    h->SetBinContent(h->GetNbinsX()+1,0);
}


// Rebins a mass spectrum onto the analysis binning: fine at low mass, growing
// progressively towards the tail so that every bin keeps a usable population.
// The returned histogram is a new object (TH1D::Rebin allocates when a name and
// an edge array are given).
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
//
// For every bin i, integrates both histograms from i to the end and returns the
// ratio h2/h1. This mirrors the way the analysis is actually used: a mass
// threshold is applied and everything above it is counted, so the relevant
// comparison is between cumulative yields, not between individual bins.
//
// h1 is the denominator (the prediction, when called from saveHistoRatio) and
// receives the flat systErr_ relative uncertainty in addition to its statistical
// error. The final error is the standard ratio propagation, assuming h1 and h2
// are uncorrelated.
//
// Bins where the denominator integral is non-positive are left empty.
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



// Renames, optionally rebins, and writes the observed and predicted spectra plus
// their cumulative ratio into the current ROOT directory.
//
// Note the argument order in the ratioIntegral call: it is invoked as
// ratioIntegral(prediction, observation), so the stored ratio is obs/pred.
//
// h1 and h2 are expected to be detached clones owned by the caller; when rebin is
// true the local pointers are replaced by the rebinned copies, which are the ones
// actually written.
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

// ---------------------------------------------------------------------------
// All the knobs of bckgEstimate(), grouped in a struct.
//
// Passing them as a single object instead of a 20+ argument list removes an
// entire class of bugs (silent argument transposition) and lets the caller set
// only what differs from the defaults.
// ---------------------------------------------------------------------------
struct BckgOptions {
    int    nPE             = 200;    // number of pseudo-experiments
    bool   useFit          = true;   // master switch for the template tail fits
    bool   useOldIhFit     = false;  // legacy Ih fit shape
    bool   useOld1oPFit    = true;   // legacy 1/p fit shape (erf of log)
    bool   corrTemplateIh  = false;  // apply the Fpixel correlation correction to Ih
    bool   corrTemplate1oP = false;  // apply the Fpixel correlation correction to 1/p
    std::string etaName    = "";     // eta range tag, selects the correction parameters
    bool   saveFits        = false;  // dump every fit into DebugFit/
    int    fitIh           = 1;      // 1 = nominal, 2 = +1 sigma, else -1 sigma
    int    fitP            = 1;      // idem for the 1/p template
    bool   blind           = false;  // zero the observed spectrum above 300 GeV
    unsigned nWorkers      = 25;     // parallel processes for the toys
};

// =============================================================================
// bckgEstimate - top-level driver of the background prediction.
//
// Overall flow:
//   1. pre-fit the 1/p spectrum of the whole C region to obtain stable seeds
//      (par_p2, par_p3) for the per-toy fits;
//   2. run nPE pseudo-experiments in parallel; each one Poisson-smears the input
//      templates, applies the eta reweighting, runs the mass convolution and
//      returns its prediction plus the ABCD normalisation factor;
//   3. average the toys bin by bin (mean = prediction, RMS = uncertainty);
//   4. write the prediction, the observation, their cumulative ratio and a set of
//      control plots into the current ROOT directory.
//
// Region arguments:
//   B : source of the Ih template
//   C : source of the 1/p template
//   BC: region in which the prediction is built
//   A : normalisation region, and reference eta profile for the reweighting
//   D : signal/validation region, source of the observed spectrum
//   ifIhpSAME / B_ifIhpSAME: alternative configuration in which the Ih and 1/p
//   cuts coincide, requiring a different B region for the normalisation and a
//   different reweighting path (see the two branches below).
//
// IMPORTANT: Region copies are shallow (see Regions.h). The local copies c, bc, d
// share their histograms with C, BC, D; only their *pointers* are ever reassigned,
// never the pointed-to objects, so the originals stay intact.
//
// Returns false if any toy failed to come back, so the caller can mark the whole
// configuration as failed instead of writing a partial result.
// =============================================================================
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

    // Local unpacking: keeps the body and the std::bind below unchanged.
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

    // Multi-process executor: each toy runs in a forked process.
    // WARNING: EnableImplicitMT must NOT be active at the same time, otherwise the
    // CPU is heavily over-subscribed (nWorkers x nThreads).
    // WARNING: shell-level stdout redirection (> log 2>&1) breaks the forked
    // workers, which is why the Python launcher captures the output through
    // subprocess instead.
    ROOT::TProcessExecutor workers(nWorkers);

    // Shallow copies: only their histogram pointers are reassigned below.
    Region c = C;
    Region bc = BC;
    Region d = D;

    TH2D c_eta_p_base(*c.eta_p);

    // Pre-fit of 1/p to get the parameters for the next fits in the toys
    // Fitting the unsmeared, full-statistics spectrum once gives seed values that
    // make the per-toy fits converge far more reliably.
    TH1D* p_base = (TH1D*)c_eta_p_base.ProjectionX();
    p_base->SetDirectory(nullptr);

    // Fit range capped just below the spectrum peak, so only the well-populated
    // rising part is used.
    float rangemax_p = 30;
    if (p_base->GetBinCenter(p_base->GetMaximumBin()) < rangemax_p) rangemax_p = 0.8 * p_base->GetBinCenter(p_base->GetMaximumBin());
    
    TF1 f_p_base("f_p_base","[0]*([1]+erf((log(x)-[2])/[3]))",0,rangemax_p);
    f_p_base.SetParameter(0,560);
    f_p_base.FixParameter(1,1.0);   // overall offset fixed: only the shape matters
    f_p_base.SetParameter(2,3.50116e+00);
    f_p_base.SetParameter(3,0.60152e+00);

    // "L" = likelihood fit, appropriate for a low-statistics binned spectrum.
    p_base->Fit(&f_p_base, "QRSL", "", 0, rangemax_p);

    // Only the two shape parameters are propagated to the toys; the normalisation
    // is re-derived per eta slice inside fillPredMass.
    double par_p2 = f_p_base.GetParameter(2);
    double par_p3 = f_p_base.GetParameter(3);

    delete p_base;


    // -----------------------------------------------------------------------
    // One pseudo-experiment.
    //
    // Runs in a forked process, so everything it touches must be either copied in
    // or created locally; the only thing that crosses back is the returned tuple.
    // Every histogram is explicitly detached from any directory (gROOT->cd() plus
    // SetDirectory(nullptr)) to avoid ROOT's global directory ownership across the
    // fork boundary.
    //
    // Returns: (predicted mass, predicted mass vs eta, ABCD normalisation,
    //           Ih template used, reweighted 1/p template).
    // -----------------------------------------------------------------------
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
        
        // Setup: shallow copies again, plus value copies of the histograms that
        // will be smeared, so the originals are never modified.
        Region a = A;
        Region b = B;
        Region b_ifIhpSAME = B_ifIhpSAME;
        Region c = C;
        Region bc = BC;

        // Reference (unsmeared) templates and their eta profiles.
        // Note the axis conventions: ProjectionX on ih_eta gives eta,
        // ProjectionY on eta_p gives eta.
        TH2D a_ih_eta_base(*a.ih_eta);
        TH1D* a_eta_base = (TH1D*)a_ih_eta_base.ProjectionX();
        TH2D b_ih_eta_base(*b.ih_eta);
        TH2D b_ifIhpSAME_ih_eta_base(*b_ifIhpSAME.ih_eta);
        TH2D b_ifIhpSAME_eta_p_base(*b_ifIhpSAME.eta_p);
        TH1D* b_ifIhpSAME_eta_base = (TH1D*)b_ifIhpSAME_eta_p_base.ProjectionY();
        TH2D b_eta_p_base(*b.eta_p);
        TH1D* b_eta_base = (TH1D*)b_eta_p_base.ProjectionY();
        TH2D c_eta_p_base(*c.eta_p);

        // Optional Fpixel-correlation corrections, applied before the smearing so
        // that the Poisson fluctuations are drawn around the corrected means.
        if (corrTemplateIh) corrIh(&b_ih_eta_base, etaName);
        if (corrTemplate1oP) corr1oP(&c_eta_p_base, etaName);
        

        // --- 1/p fit function, configured once and refitted per eta slice inside
        // fillPredMass. Two shapes are available: the legacy erf-of-log and the
        // newer symmetric cosh-like parametrisation.
        TH1D* p_base = (TH1D*)c_eta_p_base.ProjectionX();
        float rangemax_p = 30;
        if (p_base->GetBinCenter(p_base->GetMaximumBin()) < rangemax_p) rangemax_p = 0.8 * p_base->GetBinCenter(p_base->GetMaximumBin());
        
        TF1 f_p("f_p", useOld1oPFit ? "[0]*([1]+erf((log(x)-[2])/[3]))": "0.5*(exp([0]*x*x+[1]*x)+exp(-[0]*x*x-[1]*x))-1", 0, rangemax_p);
        if (useOld1oPFit) {
            f_p.SetParLimits(0, 0, 1);
            f_p.FixParameter(1, 1.0);
            f_p.SetParameter(2, par_p2);   // seeded from the global pre-fit
            f_p.SetParameter(3, par_p3);
        }
        else {
            // Tight limits: the new shape is very sensitive and would otherwise
            // wander into unphysical regions on low-statistics slices.
            f_p.SetParLimits(0, 0, 1e-3);
            f_p.SetParLimits(1, -1e-2, 1e-2);
        }

        // --- Ih fit function: a Gaussian anchored on the observed peak, fitted on
        // the falling tail only (from the peak up to 8 MeV/cm).
        TH1D* ih_base = (TH1D*)b_ih_eta_base.ProjectionX();
        float max_ih = ih_base->GetBinCenter(ih_base->GetMaximumBin());

        TF1 f_ihg("f_ihg", "gaus", max_ih, 8);
        f_ihg.SetParameter(0, 0.01*ih_base->Integral());
        f_ihg.SetParameter(1, max_ih);
        f_ihg.SetParameter(2, ih_base->GetStdDev());


        // --- Pseudo-experiment ---------------------------------------------
        // Seed derived from the worker index: reproducible across runs, distinct
        // across toys.
        TRandom3* RNG = new TRandom3(kSeedBase + workerID);
        bc.pred_mass->Reset();
        bc.pred_mass_eta->Reset();

        // Every input template is smeared independently. The eta profiles are
        // smeared separately from the 2D histograms they come from, so the
        // reweighting factors fluctuate as well.
        TH2D* a_ih_eta = poissonHisto(a_ih_eta_base, RNG);
        TH2D* b_ih_eta = poissonHisto(b_ih_eta_base, RNG);
        TH2D* b_ifIhpSAME_ih_eta = poissonHisto(b_ifIhpSAME_ih_eta_base, RNG);
        TH2D* b_eta_p = poissonHisto(b_eta_p_base, RNG);
        TH2D* c_eta_p = poissonHisto(c_eta_p_base, RNG);
        
        TH1D* b_eta = poissonHisto(*b_eta_base, RNG);
        TH1D* b_ifIhpSAME_eta = poissonHisto(*b_ifIhpSAME_eta_base, RNG);
        TH1D* a_eta = poissonHisto(*a_eta_base, RNG);
        
        // Eta reweighting. Two mutually exclusive paths:
        //  - ifIhpSAME: the Ih template is reweighted to the B/A eta profile;
        //  - otherwise: the 1/p template is reweighted to the eta profile of B.
        if(ifIhpSAME) etaReweighingP_X(b_ih_eta, b_ifIhpSAME_eta, a_eta);
        else etaReweighingP_Y(c_eta_p,b_eta);

        // Control template (written out for monitoring), always reweighted B/A.
        etaReweighingP_Y(b_eta_p, b_ifIhpSAME_eta, a_eta);


        // Mass prediction in the BC region
        // The local Region "bc" is pointed at this toy's templates, then the
        // convolution fills its pred_mass / pred_mass_eta.
        bc.eta_p = c_eta_p;
        bc.ih_eta = b_ih_eta;
        bc.fillPredMass(filename, st, calib, f_p, f_ihg, useFit, fitIh, fitP, 
                        useOldIhFit, useOld1oPFit, etaName, saveFits,
                        par_p2, par_p3, workerID);
        
        // ABCD normalisation: yields integrated over the full histograms,
        // under/overflow included.
        double normA = a_ih_eta->Integral(0, a_ih_eta->GetNbinsX()+1, 0, a_ih_eta->GetNbinsY()+1);
        double normB = b_ih_eta->Integral(0, b_ih_eta->GetNbinsX()+1, 0, b_ih_eta->GetNbinsY()+1);
        double normC = c_eta_p->Integral(0, c_eta_p->GetNbinsX()+1, 0, c_eta_p->GetNbinsY()+1);
        double normB_ifIhpSAME = b_ifIhpSAME_ih_eta->Integral(0, b_ifIhpSAME_ih_eta->GetNbinsX()+1, 0, b_ifIhpSAME_ih_eta->GetNbinsY()+1);

        // The expected yield in D under the ABCD assumption of factorisation.
        double normalisationABC = normB * normC / normA;
        if (ifIhpSAME) normalisationABC = normB_ifIhpSAME * normC / normA;


        // TEMPORARY
        // Debug switch: keeps the convolution's own normalisation instead of the
        // ABCD one, which isolates shape effects from normalisation effects.
        //normalisationABC = bc.pred_mass->Integral();
        // TEMPORARY
        

        // The convolution only produces a *shape*; the absolute scale is imposed
        // here. The 1D and 2D predictions are normalised independently, both to
        // the same target, so their projections stay consistent.
        const double itg = bc.pred_mass->Integral();
        const double itg2D = bc.pred_mass_eta->Integral();
        if (itg > 0) bc.pred_mass->Scale(normalisationABC/itg);
        else std::cerr << "toy " << workerID << ": empty prediction" << std::endl;
        if (itg2D > 0) bc.pred_mass_eta->Scale(normalisationABC/itg2D);
        else std::cerr << "toy " << workerID << ": empty prediction" << std::endl;

        // Value copies for the return trip across the process boundary; the
        // pointers themselves cannot cross.
        TH1D out_pred_mass    (*bc.pred_mass);
        TH2D out_pred_mass_eta(*bc.pred_mass_eta);
        TH2D out_ih_eta       (*bc.ih_eta);
        TH2D out_p_eta    (*b_eta_p);
        out_pred_mass.SetDirectory(nullptr);
        out_pred_mass_eta.SetDirectory(nullptr);
        out_ih_eta.SetDirectory(nullptr);
        out_p_eta.SetDirectory(nullptr);

        // End
        // Explicit cleanup of everything allocated in this worker.
        delete a_eta;              delete a_ih_eta;           delete a_eta_base;
        delete b_ih_eta;           delete b_eta_p;            delete b_eta;
        delete b_ifIhpSAME_eta;    delete b_ifIhpSAME_ih_eta; delete b_eta_base;
        delete b_ifIhpSAME_eta_base;
        delete c_eta_p;
        delete p_base;             delete ih_base;
        // bc borrowed these pointers; they have just been deleted, so they are
        // cleared to prevent any later dereference or double free.
        bc.ih_eta = nullptr;
        bc.eta_p  = nullptr;
        delete RNG;

        return {out_pred_mass, out_pred_mass_eta, normalisationABC, out_ih_eta, out_p_eta};
    };

    
    // Bind every constant argument, leaving only the worker index free (_1),
    // which is what TProcessExecutor::Map supplies.
    // Loop on the toys
    auto workItemToRun = std::bind (workItem, _1, filename, B, C, BC, A, ifIhpSAME, B_ifIhpSAME, st, calib,
                                    useFit, useOldIhFit, useOld1oPFit, corrTemplateIh, corrTemplate1oP, etaName, saveFits, fitIh, fitP, 
                                    par_p2, par_p3);
    
    auto vPE = workers.Map(workItemToRun, ROOT::TSeqI(nPE));
    // A crashed worker silently drops its result, so the count is checked
    // explicitly: a partial set of toys would bias the mean and underestimate the
    // spread.
    if (vPE.size() != (size_t)nPE) {
        std::cerr << "ERROR: " << vPE.size() << "/" << nPE
                  << " toys back, some workers crashed" << std::endl;
        return false;
    }


    // ---- Collect the results ------------------------------------------------
    // Get the results
    std::vector<TH1D> histo_pred_mass;
    std::vector<TH2D> histo_pred_mass_eta;
    std::vector<double> normalisations;
    std::vector<TH1D> Ih_eta, oP_eta;

    for (const auto& result: vPE) {
        histo_pred_mass.push_back(std::get<0>(result));     // predicted mass spectrum
        histo_pred_mass_eta.push_back(std::get<1>(result)); // predicted mass vs eta
        normalisations.push_back(std::get<2>(result));      // ABCD normalisation factor
        
        // Control plot: Ih projection of the template actually used (eta on X, so
        // ProjectionY gives Ih). Unique names avoid ROOT name clashes.
        const TH2D& h2 = std::get<3>(result);
        std::unique_ptr<TH1D> proj(h2.ProjectionY(Form("ih_eta_py_%zu", Ih_eta.size())));
        proj->SetDirectory(nullptr);
        Ih_eta.push_back(*proj);
        Ih_eta.back().SetDirectory(nullptr);

        // Control plot: 1/p projection of the reweighted template (eta on Y, so
        // ProjectionX gives 1/p).
        const TH2D& h3 = std::get<4>(result);
        std::unique_ptr<TH1D> proj2(h3.ProjectionX(Form("p_eta_px_%zu", oP_eta.size())));
        proj2->SetDirectory(nullptr);
        oP_eta.push_back(*proj2);
    }

    // Bin-by-bin mean over the toys; the bin errors are the toy-to-toy RMS, i.e.
    // the statistical uncertainty of the prediction.
    TH1D h_ih_eta_mean = MeanOfToys(Ih_eta, ("ih_eta_mean_"+st).c_str());
    TH1D h_oP_eta_mean = MeanOfToys(oP_eta, ("oP_eta_mean_"+st).c_str());
    TH1D h_temp        = MeanOfToys(histo_pred_mass,     ("pred_mass_mean_"+st).c_str());
    TH2D h_temp_eta    = MeanOfToys(histo_pred_mass_eta, ("pred_mass_eta_mean_"+st).c_str());
        
    // bc is a local shallow copy: repointing it here leaves BC untouched. The two
    // stack objects h_temp / h_temp_eta must outlive every use of bc below.
    bc.pred_mass = &h_temp;
    bc.pred_mass_eta = &h_temp_eta;
    double avgNormalisation = std::accumulate(normalisations.begin(), normalisations.end(), 0.0) / normalisations.size();


    // ---- Observed spectrum --------------------------------------------------
    TH1D* d_mass = (TH1D*) d.mass->Clone(("mass_obs_"+st).c_str());
    d_mass->SetDirectory(nullptr);
    
    if(blind) blindMass(d_mass,300);
    // Both spectra get the same overflow treatment so the cumulative ratio is
    // computed on identical ranges.
    overflowLastBin(d_mass);
    overflowLastBin(bc.pred_mass);


    // ---- Output -------------------------------------------------------------
    // Saving histograms
    // Writes mass_obs, mass_predBC and their cumulative ratio mass_predBCR.
    saveHistoRatio(d_mass, bc.pred_mass, ("mass_obs_"+st).c_str(), ("mass_predBC_"+st).c_str(), ("mass_predBCR_"+st).c_str());

    // 2D prediction and its eta profile.
    bc.pred_mass_eta->Write();
    bc.pred_mass_eta->ProjectionY()->Write();

    // Observed mass vs eta in D and C, with their eta profiles, for the
    // closure/validation plots produced by the plotting scripts.
    d.mass_eta->SetName(("mass_eta_D_"+st).c_str());        d.mass_eta->Write();
    d.mass_eta->ProjectionY(("mass_eta_D_"+st+"_py").c_str())->Write();
    c.mass_eta->SetName(("mass_eta_C_"+st).c_str());        c.mass_eta->Write();
    c.mass_eta->ProjectionY(("mass_eta_C_"+st+"_py").c_str())->Write();

    // Raw (unsmeared) input templates, kept so that the inputs of a given run can
    // always be re-inspected from the output file alone.
    A.ih_eta->Write();
    A.eta_p->Write();
    B.ih_eta->Write();
    C.eta_p->Write();
    D.eta_p->Write();
    D.ih_eta->Write();
    
    // Projections of the same, used directly by the control plots.
    D.ih_eta->ProjectionY()->Write();
    B.ih_eta->ProjectionY()->Write();
    D.ih_eta->ProjectionX()->Write();
    B.ih_eta->ProjectionX()->Write();
    D.eta_p->ProjectionX()->Write();
    C.eta_p->ProjectionX()->Write();   
    C.eta_p->ProjectionY()->Write();

    // Mass spectrum of C rescaled to the average ABCD normalisation: shows what
    // the prediction would look like without the convolution, i.e. isolates the
    // effect of the mass convolution itself.
    std::unique_ptr<TH1D> cMassScaled(static_cast<TH1D*>(C.mass->Clone(("mass_C_"+st).c_str())));
    cMassScaled->SetDirectory(nullptr);
    const double itgC = cMassScaled->Integral();
    if (itgC > 0) cMassScaled->Scale(avgNormalisation/itgC);
    cMassScaled->Write();

    // Superposition of every individual toy: a visual check of the spread and of
    // any pathological toy. The rebinning is applied here only, after the means
    // have been computed, so it does not affect the prediction.
    // Draw on one canvas all the histos of the vector histo_pred_mass
    TCanvas* c1  = new TCanvas(("toys_"+st).c_str(), "toys", 800, 800);
    c1->cd();
    for (size_t i=0; i<histo_pred_mass.size(); i++) histo_pred_mass[i].Rebin(10);
    histo_pred_mass[0].Draw("hist");
    for (size_t i=1; i<histo_pred_mass.size(); i++) histo_pred_mass[i].Draw("hist same");
    c1->Write();

    // Distribution of the ABCD closure over the toys: predicted yield divided by
    // the observed one. A distribution centred on 1 means the method closes.
    // fill a new histo with norm/d.mass->integral:
    TH1D* h_norm = new TH1D(("h_norm_"+st).c_str(), "h_norm", 40, 0.9, 1.1);
    h_norm->GetXaxis()->SetTitle("BC/AD");
    h_norm->GetYaxis()->SetTitle("Entries");
    for(size_t i=0; i<normalisations.size(); i++) h_norm->Fill(normalisations[i]/d.mass->Integral());
    h_norm->Write();

    // Template-vs-observation comparisons in the validation region: the toy-mean
    // Ih and 1/p templates against what is actually seen in D.
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