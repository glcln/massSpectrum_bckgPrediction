#pragma once

// =============================================================================
//  Regions.h
// -----------------------------------------------------------------------------
//  Core building blocks of the HSCP data-driven background estimate (ABCD method
//  with mass convolution).
//
//  This header provides:
//    * the dE/dx calibration constants (K, C) used to turn a (p, Ih) pair into a
//      reconstructed mass;
//    * a few helpers acting on the template histograms (unit normalisation,
//      Fpixel-correlation corrections, blinding);
//    * the Region class, which owns the set of histograms describing one ABCD
//      region and implements fillPredMass(), i.e. the mass convolution itself.
//
//  Conventions used throughout the whole package:
//    - "Ih"     : most-probable ionisation energy-loss estimator [MeV/cm]
//    - "1/p"    : inverse momentum, stored *scaled* by kInvPScale (see below)
//    - "eta"    : track pseudorapidity
//    - A / B / C / D are the usual ABCD regions. The prediction is built in the
//      "BC" region by convolving the Ih template of B with the 1/p template of C,
//      then normalised by the ABCD factor N_B * N_C / N_A.
//
//  Axis conventions of the 2D histograms (must match what step1 writes):
//    - eta_p       : X = 1/p (scaled), Y = eta
//    - ih_eta      : X = eta,          Y = Ih
//    - mass_eta    : X = mass,         Y = eta
//    - pred_mass_eta: X = mass,         Y = eta
// =============================================================================

#include <TCanvas.h>
#include <TLegend.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TDirectory.h>
#include <TRatioPlot.h>
#include <THStack.h>
#include <iostream>
#include <map>
#include <string>
#include <memory> 
#include <algorithm>
#include <vector>
#include <cmath>
#include <iterator>
#include <functional>
#include "TF1.h"
#include "TFitResult.h"
#include "Math/Integrator.h"
#include "Math/IntegratorOptions.h"
#include "Math/WrappedTF1.h"

// The X axis of the eta_p histograms does NOT hold the momentum in GeV: it holds
// the inverse momentum multiplied by this scale factor, i.e. x = kInvPScale / p.
// Converting a bin position back to a physical momentum is therefore
//     p [GeV] = kInvPScale / x
// which is exactly what is done in fillPredMass() below.
constexpr float kInvPScale = 10000.f;

// SETUP

// dE/dx calibration pair entering the mass formula m = p * sqrt((Ih - C) / K).
// K is in MeV/cm per (GeV/c)^2-like units, C is the offset in MeV/cm.
struct DeDxCalib { float K, C; };

// Returns the (K, C) calibration constants for a given sample key.
// The key is carried by the Python launcher (dataset dict field "sampleType"),
// so that data and MC of the same year automatically get their own calibration.
// Unknown keys fall back to the 2024 data calibration and emit a warning rather
// than aborting, so a typo degrades the result instead of killing the job.
inline DeDxCalib GetDeDxCalib(const std::string& sample) {
    static const std::map<std::string, DeDxCalib> kCalib = {
        {"data2017", {2.54f,    3.14f}},
        {"data2018", {2.55f,    3.14f}},
        {"data2024", {2.8202f,  2.9784f}},
        {"mc2017",   {2.48f,    3.19f}},
        {"mc2018",   {2.49f,    3.19f}},
        {"mc2024",   {2.83894f, 3.01756f}},
        {"ttbar2024",   {2.83894f, 3.01756f}},
    };
    auto it = kCalib.find(sample);
    if (it == kCalib.end()) {
        std::cerr << "GetDeDxCalib: unknown sample '" << sample << "', will take 2024 param as default ones " << std::endl;
        return {2.8202f, 2.9784f};
    }
    return it->second;
}

// Flat relative systematic uncertainty assigned to the background prediction and
// added in quadrature inside ratioIntegral() (CommonFunctions.h).
// It is deliberately set to 0 while *measuring* the systematics, so that the
// spread between variations is not double counted.
constexpr float systErr_ = 0.; //set to 0 for systematic studies




// Normalise a 1D histogram to unit area, under/overflow bins included.
// A non-positive integral is reported and the histogram left untouched, because
// scaling by 1/0 would silently poison every downstream template.
void scale(TH1D* h) {
    const double itg = h->Integral(0, h->GetNbinsX()+1);
    if (itg <= 0) { std::cerr << "scale: null integral for " << h->GetName() << std::endl; return; }
    h->Scale(1./itg);
}


// ---------------------------------------------------------------------------
// corrIh: removes the residual Ih dependence of the Fpixel selection from the
// B-region template.
//
// The correction is a first-order polynomial in Ih measured beforehand (ratio of
// the Ih spectrum with and without the Fpixel requirement). Its parameters depend
// on the eta range of the region, hence the lookup on the region name.
//
// The template is divided by the correction factor, bin content and bin error
// alike, so the relative statistical uncertainty of each bin is preserved.
// ---------------------------------------------------------------------------
void corrIh(TH2D* ih_eta, const std::string& etaName) {
    constexpr double kIhLow = 2., kIhUp = 10.;   // range over which the pol1 was fitted
    TF1 f_correlation_Ih_Fpix("f_correlation_Ih_Fpix", "pol1", kIhLow, kIhUp);

    // Default = the widest eta range. The tests below must stay ordered from the
    // most specific name to the least specific one, because they use find():
    // "Eta1_2p4" also contains the substring "Eta1".
    float par0 = 1.03878, par1 = -0.0120062;
    if      (etaName.find("Eta2p4")   != std::string::npos) {par0 = 1.03878; par1 = -0.0120062;}
    else if (etaName.find("Eta1_2p4") != std::string::npos) {par0 = 1.12024; par1 = -0.03684;}
    else if (etaName.find("Eta1")     != std::string::npos) {par0 = 1.15082; par1 = -0.0459028;}

    f_correlation_Ih_Fpix.SetParameter(0, par0);
    f_correlation_Ih_Fpix.SetParameter(1, par1);

    const TAxis* ay = ih_eta->GetYaxis();

    // The correction depends on Ih only: it is evaluated once per Ih bin and then
    // applied to every eta bin of that row (under/overflow included).
    for (int bin_ih = 0; bin_ih <= ih_eta->GetNbinsY() + 1; ++bin_ih) {

        // Outside the fitted range the correction is frozen to its edge value
        // (constant extrapolation). A linear extrapolation would diverge quickly
        // and can even change sign, which is why it is deliberately avoided.
        const double x    = std::clamp(ay->GetBinCenter(bin_ih), kIhLow, kIhUp);
        const double corr = f_correlation_Ih_Fpix.Eval(x);

        // A non-positive factor would flip the sign of the template: skip the bin.
        if (corr <= 0.) {
            std::cerr << "corrIh: factor <= 0 (" << corr << ") for Ih=" << x
                      << " -> bin not corrected" << std::endl;
            continue;
        }

        for (int bin_eta = 0; bin_eta <= ih_eta->GetNbinsX() + 1; ++bin_eta) {
            ih_eta->SetBinContent(bin_eta, bin_ih, ih_eta->GetBinContent(bin_eta, bin_ih) / corr);
            ih_eta->SetBinError  (bin_eta, bin_ih, ih_eta->GetBinError  (bin_eta, bin_ih) / corr);
        }
    }
}

// ---------------------------------------------------------------------------
// corr1oP: same idea as corrIh, but for the 1/p template of the C region.
// The correction is a pol1 in the *scaled* 1/p variable (x = kInvPScale / p),
// so the fitted range [0, 200] corresponds to p above ~50 GeV.
// Note the transposed axis convention: here eta is on Y and 1/p on X.
// ---------------------------------------------------------------------------
void corr1oP(TH2D* eta_p, const std::string& etaName) {
    constexpr double k1oPLow = 0., k1oPUp = 200.;
    TF1 f_correlation_1oP_Fpix("f_correlation_1oP_Fpix", "pol1", k1oPLow, k1oPUp);

    // Same ordering caveat as in corrIh: most specific region name first.
    float par0 = 0.871848, par1 = 0.0020839;
    if      (etaName.find("Eta2p4")   != std::string::npos) {par0 = 0.871848; par1 = 0.0020839;}
    else if (etaName.find("Eta1_2p4") != std::string::npos) {par0 = 0.926461; par1 = 0.0019605;}
    else if (etaName.find("Eta1")     != std::string::npos) {par0 = 0.938621; par1 = 0.000665965;}

    f_correlation_1oP_Fpix.SetParameter(0, par0);
    f_correlation_1oP_Fpix.SetParameter(1, par1);

    const TAxis* ax = eta_p->GetXaxis();

    // One evaluation per 1/p bin, applied to all eta bins of that column.
    for (int bin_1oP = 0; bin_1oP <= eta_p->GetNbinsX() + 1; ++bin_1oP) {

        // Constant extrapolation outside the fitted range (see corrIh).
        const double x    = std::clamp(ax->GetBinCenter(bin_1oP), k1oPLow, k1oPUp);
        const double corr = f_correlation_1oP_Fpix.Eval(x);

        if (corr <= 0.) {
            std::cerr << "corr1oP: factor <= 0 (" << corr << ") for 1/p=" << x
                      << " -> bin not corrected" << std::endl;
            continue;
        }

        for (int bin_eta = 0; bin_eta <= eta_p->GetNbinsY() + 1; ++bin_eta) {
            eta_p->SetBinContent(bin_1oP, bin_eta, eta_p->GetBinContent(bin_1oP, bin_eta) / corr);
            eta_p->SetBinError  (bin_1oP, bin_eta, eta_p->GetBinError  (bin_1oP, bin_eta) / corr);
        }
    }
}

// Blinding: empties every mass bin whose *lower edge* is at or above mass_value.
// Using the lower edge means a bin straddling the threshold is removed entirely,
// which is the conservative choice. Under/overflow are covered by the 0..N+1 loop.
void blindMass(TH1D* h_m, float mass_value=300) {
    for(int i=0; i<h_m->GetNbinsX()+2; i++){
        if(h_m->GetBinLowEdge(i)>=mass_value) {
            h_m->SetBinContent(i,0);
            h_m->SetBinError(i,0);
        }
    }
}

// HSCP mass estimator from the Bethe-Bloch-like parametrisation
//     Ih = K * m^2 / p^2 + C   ==>   m = p * sqrt((Ih - C) / K)
// Returns -1 when Ih sits below the offset C (unphysical, no real solution);
// callers must test for that and skip the entry.
inline float GetMass(float p, float ih, float K, float C) {
    if (ih - C < 0) return -1;
    return std::sqrt((ih - C) / K) * p;
}


// ---------------------------------------------------------------------------
// Region
//
// Holds the histogram set describing one ABCD region, plus the mass convolution.
//
// IMPORTANT — ownership: the members are raw pointers and the class relies on the
// implicit copy constructor, so copying a Region is a *shallow* copy: the copy
// shares the histograms with the original. This is used on purpose in
// bckgEstimate() (CommonFunctions.h), where a local copy has its eta_p / ih_eta
// pointers reassigned to per-toy histograms and reset to nullptr afterwards, so
// the original region is never touched. Never let two Regions delete the same
// pointer.
//
// The histograms are either created by initHisto() (step1, owned by the output
// directory) or loaded and cloned by loadHistograms() (step2, detached from any
// directory via SetDirectory(nullptr)).
// ---------------------------------------------------------------------------
class Region{
    public:
        Region();
        ~Region();

        // Convenience constructor: names the region and books its histograms in
        // one go. TDir is templated so that this header does not depend on CMSSW:
        // step1 passes a TFileDirectory, step2 can pass any object exposing
        // make<T>(...).
        template <typename TDir>
        Region(TDir& dir, std::string suffix, int etabins, int ihbins, int pbins, int massbins) {
            suffix_ = std::move(suffix);
            initHisto(dir, etabins, ihbins, pbins, massbins);
        }

        template <typename TDir>
        void initHisto(TDir& dir, int etabins, int ihbins, int pbins, int massbins);

        // Core of the method: convolves the 1/p template (eta_p) with the Ih
        // template (ih_eta), eta bin by eta bin, and fills pred_mass /
        // pred_mass_eta. See the definition below for the full description.
        void fillPredMass(const std::string& filename,
                          const std::string& st,
                          const DeDxCalib& calib,
                          TF1& f_p,
                          TF1& f_ih,
                          const bool useFit,
                          const int& fit_ih_err = 1,
                          const int& fit_p_err = 1,
                          bool useOldIhFit = false,
                          bool useOld1oPFit = false,
                          const std::string& etaName = "",
                          bool saveFits = false,
                          const double par_p2 = 4.70839,
                          const double par_p3 = 1.05005,
                          const UInt_t workerID = 0);

        // Binning definitions, kept as members so that step1 and step2 agree.
        int np;          // number of 1/p bins
        float plow;      // lower edge of the (scaled) 1/p axis
        float pup;       // upper edge of the (scaled) 1/p axis
        int npt;         // pt binning: declared for completeness, unused here
        float ptlow;
        float ptup;
        int nih;         // number of Ih bins
        float ihlow;
        float ihup;
        int neta;        // number of eta bins
        float etalow;
        float etaup;
        int nmass;       // number of mass bins
        float masslow;
        float massup;
        std::string suffix_;   // region tag appended to every histogram name

        TH2D* eta_p         = nullptr;   // X = 1/p (scaled), Y = eta
        TH2D* ih_eta        = nullptr;   // X = eta,          Y = Ih
        TH1D* mass          = nullptr;   // observed mass spectrum
        TH1D* pred_mass     = nullptr;   // predicted mass spectrum (filled here)
        TH2D* mass_eta      = nullptr;   // X = mass, Y = eta (observed)
        TH2D* pred_mass_eta = nullptr;   // X = mass, Y = eta (predicted)
};


Region::Region(){}


// ---------------------------------------------------------------------------
// initHisto: books the histogram set of the region.
//
// Templated on the directory type so that this header stays CMSSW-free: any
// object exposing a make<T>(args...) factory works. The explicit `template`
// keyword in `dir.template make<...>` is mandatory because TDir is a dependent
// type, otherwise `<` is parsed as a comparison.
//
// Sumw2 is forced globally so that every booked histogram carries proper errors
// from the very first fill.
// ---------------------------------------------------------------------------
template <typename TDir>
void Region::initHisto(TDir& dir, int etabins, int ihbins, int pbins, int massbins) {
    TH1::SetDefaultSumw2(kTRUE);
    TH2::SetDefaultSumw2(kTRUE);

    np = pbins;  plow = 0;  pup = 10000;
    nih = ihbins; ihlow = 0; ihup = 20;
    neta = etabins; etalow = -3; etaup = 3;
    nmass = massbins; masslow = 0; massup = 4000;
    const std::string suffix = suffix_;

    // 'dir.template make<...>': the 'template' keyword is required because dir
    // depends on the TDir template parameter.
    eta_p         = dir.template make<TH2D>(("eta_p"        + suffix).c_str(), ";10^{4}/p [GeV^{-1}];#eta",              np,   plow,   pup,     neta, etalow, etaup);
    ih_eta        = dir.template make<TH2D>(("ih_eta"       + suffix).c_str(), ";#eta;I_{h} [MeV/cm]",       neta, etalow, etaup,   nih,  ihlow,  ihup);
    mass          = dir.template make<TH1D>(("mass"         + suffix).c_str(), ";Mass [GeV]",                nmass, masslow, massup);
    pred_mass     = dir.template make<TH1D>(("pred_mass"    + suffix).c_str(), ";Mass [GeV]",                nmass, masslow, massup);
    mass_eta      = dir.template make<TH2D>(("mass_eta"     + suffix).c_str(), ";Mass [GeV];#eta",           nmass, masslow, massup, neta, etalow, etaup);
    pred_mass_eta = dir.template make<TH2D>(("pred_mass_eta"+ suffix).c_str(), ";Mass [GeV];#eta",           nmass, masslow, massup, neta, etalow, etaup);

    // Poisson (Garwood) asymmetric errors on the *observed* spectrum only.
    // Any histogram that is later averaged over toys must keep the default kNormal
    // errors, otherwise the RMS-based toy uncertainty gets contaminated
    // (see MeanOfToys in CommonFunctions.h, which explicitly resets this option).
    mass->SetBinErrorOption(TH1::EBinErrorOpt::kPoisson);
}

Region::~Region(){}



// =============================================================================
// fillPredMass — the mass convolution
// -----------------------------------------------------------------------------
// For each eta bin i:
//   1. project out the 1/p spectrum (from eta_p) and the Ih spectrum (from ih_eta)
//      of that eta slice;
//   2. normalise the 1/p spectrum to unit area (the absolute yield is restored
//      later by the global ABCD normalisation);
//   3. fit the tails of both spectra, where the statistics are too poor to be used
//      bin by bin:
//        - Ih : Gaussian (new) or the legacy shape, above ~the peak;
//        - 1/p: erf-of-log (old) or cosh-like (new), below a fraction of the peak
//                position, i.e. at *high* momentum;
//   4. loop over every (1/p, Ih) bin pair, build the mass from the pair, and fill
//      pred_mass / pred_mass_eta with the product of the two contents (or of the
//      fitted densities where the fit replaces the data).
//
// Parameters:
//   filename, st, etaName: only used to build the debug-fit output file name.
//   calib                : dE/dx constants (K, C) used by GetMass().
//   f_p, f_ih            : the fit functions, pre-configured by the caller and
//                           reused (and refitted) for every eta bin.
//   useFit               : master switch; false = pure bin-by-bin data product.
//   fit_ih_err, fit_p_err: 1 = nominal, 2 = +1 sigma, otherwise -1 sigma. The
//                           shift is computed with TF1::IntegralError using the
//                           MINUIT covariance matrix of that eta bin's fit.
//   useOldIhFit/useOld1oPFit: select the legacy fit shapes and ranges.
//   saveFits             : dump the fitted histograms to DebugFit/ for inspection.
//   par_p2, par_p3       : seed values for the legacy 1/p fit, obtained from a
//                           pre-fit of the full C-region 1/p spectrum. Seeding from
//                           the pre-fit makes the per-toy fits far more stable.
//   workerID             : toy index; enters the debug file name and the log lines
//                           so that a bad fit can be traced back to its toy.
// =============================================================================
void Region::fillPredMass(const std::string& filename,
                          const std::string& st,
                          const DeDxCalib& calib,
                          TF1& f_p,
                          TF1& f_ih,
                          const bool useFit,
                          const int& fit_ih_err,
                          const int& fit_p_err,
                          bool useOldIhFit,
                          bool useOld1oPFit,
                          const std::string& etaName,
                          bool saveFits,
                          const double par_p2,
                          const double par_p3,
                          const UInt_t workerID) {

    // ---- Optional debug output ------------------------------------------------
    // One file per (sample, region, fit flavour, eta range, toy). The workerID in
    // the name is what prevents forked workers from overwriting each other.
    TFile* OutputHisto = nullptr;
    std::string filenameOutputFit = "DebugFit/Fits_" + filename + "_" + st + ((useOldIhFit || useOld1oPFit) ? "_OldFit": "_NewFit") + etaName + "_" + std::to_string(workerID) +  ".root";
    if (saveFits) {
        OutputHisto = new TFile(filenameOutputFit.c_str(), "RECREATE");
        OutputHisto->cd();
    }


    // ---- Setup ---------------------------------------------------------------
    // eta axis reference, taken from ih_eta (where eta is on X). It is used both
    // to drive the eta loop and to print readable eta values in the fit warnings.
    TH1D* eta = (TH1D*) ih_eta->ProjectionX();
    eta->SetDirectory(nullptr);

    const float K = calib.K;
    const float C = calib.C;

    // Per-eta-bin switches: they start from useFit and are turned off locally
    // whenever the corresponding fit fails, so one bad eta slice does not
    // invalidate the whole prediction.
    bool useFitIh = true;
    bool useFitP = true;
    // Tight tolerance for the numerical integration of the fit functions; the
    // default is too loose for the very small tail integrals involved here.
    ROOT::Math::IntegratorOneDimOptions::SetDefaultRelTolerance(1.E-9);


    // ---- Loop over the eta bins (1..N, under/overflow excluded) ---------------
    for(int i=1;i<eta->GetNbinsX()+1;i++) {
        // Reset the per-bin switches to the global setting.
        useFitIh = useFit;
        useFitP = useFit;
        // eta_p  has eta on Y -> ProjectionX(i,i) gives the 1/p spectrum of slice i
        // ih_eta has eta on X -> ProjectionY(i,i) gives the Ih  spectrum of slice i
        // The "e" option propagates the bin errors into the projection.
        // Both regions must share the same eta binning for index i to be consistent.
        std::unique_ptr<TH1D> p (static_cast<TH1D*>(eta_p ->ProjectionX(Form("proj_p_eta%d",  i), i, i, "e")));
        std::unique_ptr<TH1D> ih(static_cast<TH1D*>(ih_eta->ProjectionY(Form("proj_ih_eta%d", i), i, i, "e")));
        p->SetDirectory(nullptr);
        ih->SetDirectory(nullptr);

        // Empty eta slices carry no information: skip them rather than fitting noise.
        if (ih->GetEntries() < 1 || p->GetEntries() < 1) continue;
        if (p->Integral(0, p->GetNbinsX()+1) <= 0) continue;
        // Only the 1/p template is normalised to unit area. The Ih template keeps
        // its yield, so the raw product c_ih * c_p is proportional to the number of
        // events; the absolute scale is fixed afterwards by the ABCD factor.
        scale(p.get());

        // Default fit ranges. endIhFit is in MeV/cm; the 1/p range is in scaled
        // 1/p units and is refined below from the position of the spectrum peak.
        float endIhFit = 6., end1oPFit = 30., start1oPFit = 0;
        

        // ------------------------- Ih fit -------------------------------------
        TFitResultPtr ptr1 = 0;
        float max_ih = ih->GetBinCenter(ih->GetMaximumBin());
        // Legacy fit starts at a fixed 3 MeV/cm; the Gaussian one starts just above
        // the peak so that only the falling tail is fitted.
        float start_fit = (useOldIhFit)? 3: 1.1*max_ih; // 1.1*max_ih; for gauss
        // Index of the last populated bin (the name says "content" but it holds a
        // bin index). Used to detect slices whose spectrum stops before start_fit.
        int lastBinContent = ih->GetNbinsX();
        while (lastBinContent > 1 && ih->GetBinContent(lastBinContent) == 0) --lastBinContent;
        // If the requested start lies beyond the populated range, fall back to the
        // peak position so that the fit range is never empty.
        if(start_fit > ih->GetBinCenter(lastBinContent)) start_fit = max_ih;

        if (useFitIh) {
            // "QRS": quiet, restrict to the given range, and return a TFitResultPtr
            // (needed for the covariance matrix used in the error propagation).
            ptr1 = ih->Fit(&f_ih, "QRS", "", start_fit, endIhFit);
            // Quality criterion: a valid result with a meaningful chi2/ndf.
            bool goodFit = ptr1.Get() && ptr1->Ndf() > 0 && ptr1->Chi2()/ptr1->Ndf() < 6;
            if (!goodFit) {
                // Verbose diagnostics: the toy index and the eta value make it
                // possible to reproduce a failing fit in isolation.
                std::cout << "Bad fit Ih in " << ih->GetName() << " workerID=" << workerID
                        << " eta=" << eta->GetBinCenter(i);
                if (ptr1.Get()) {
                    std::cout << " status=" << ptr1->Status()
                            << " covMatrixStatus=" << ptr1->CovMatrixStatus()
                            << " edm=" << ptr1->Edm()
                            << " chi2/ndf=" << ptr1->Chi2() << "/" << ptr1->Ndf()
                            << " p-value=" << ptr1->Prob();
                } else std::cout << " (fit not performed)";
                std::cout << std::endl;
                if (saveFits) { OutputHisto->cd(); ih->Write(); }
                // Fall back to the raw bin contents for this eta slice.
                useFitIh = false;
            }
            else {                             // Good fit
                if (saveFits) { OutputHisto->cd(); ih->Write(); }
            }
        }

        // Scale factor bringing the fitted function onto the observed yield.
        // intIh is a plain bin-content sum while intFih is a true integral over x:
        // the ratio therefore absorbs the bin width, and since SFih is only ever
        // used to rescale *other* integrals of the same function, the convention is
        // self-consistent (the fit integrated over [3, endIhFit] reproduces the
        // observed number of entries in that window).
        TF1* const f_ih2 = &f_ih;
        double intFih = f_ih2->Integral(3, endIhFit);
        double intIh = ih->Integral(ih->FindBin(3), ih->FindBin(endIhFit));

        double SFih = (intFih > 0)? intIh/intFih: -1;
        if(SFih < 0 && useFitIh) {
            std::cout<<"ERROR > INTEGRAL FIT IH IS <= 0.   ITG = " << intFih << " FOR ETA BIN #" << i << std::endl;
            useFitIh = false;
        }
        

        // Cache the fit parameters and the covariance matrix, needed by
        // TF1::IntegralError for the +/-1 sigma template variations. hasCov guards
        // against a failed MINUIT error computation (empty covariance matrix).
        TF1* f_ih3 = f_ih2;
        bool hasCov = false;
        const double* fit_ih_params = nullptr;
        const double* fit_ih_cov = nullptr;

        if (ptr1.Get() && ptr1->Status() == 0) {
            const TMatrixDSym& cov = ptr1->GetCovarianceMatrix();
            if (cov.GetNrows() > 0) {
                fit_ih_params = ptr1->GetParams();
                fit_ih_cov    = cov.GetMatrixArray();
                hasCov = true;
            }
        }


        // ------------------------- 1/p fit ------------------------------------
        double SFp = 0;
        TFitResultPtr ptr2 = 0;
        bool hasCovP = false;
        const double* fit_p_params = nullptr;
        const double* fit_p_cov    = nullptr;
        int statusFit = 1;   // 1 = no good fit yet, 0 = converged

        // Taking the fit from the Down variation, as it is the one with the thicker bins, thus less statistical fluctuations for the fit to converge
        // Working copy: the fit modifies the histogram's associated function list,
        // and p itself is reused unmodified in the convolution loop below.
        std::unique_ptr<TH1D> p_forfit(static_cast<TH1D*>(p->Clone(Form("forfit_p_eta%d", i))));
        p_forfit->SetDirectory(nullptr);


        // Successive shrinking of the upper fit bound, as a fraction of the peak
        // position. The first fraction that yields a converged fit is kept.
        const float endFracs[5] = {0.9f, 0.8f, 0.7f, 0.6f, 0.5f};
        float peak = p_forfit->GetBinCenter(p_forfit->GetMaximumBin());
        int incrFit_end = 0;   // index of the fraction actually retained
        if (useFitP) {

            // Legacy shape: seed the parameters from the global pre-fit and bound
            // them, which is what makes the per-toy fits reproducible.
            if (useOld1oPFit) {
                f_p.SetParameter(0, p_forfit->GetMaximum());
                f_p.FixParameter(1, 1.0);
                f_p.SetParameter(2, par_p2);
                f_p.SetParameter(3, par_p3);
                f_p.SetParLimits(0, 0, 10*p_forfit->GetMaximum());
                f_p.SetParLimits(2, 0, 10*par_p2);
                f_p.SetParLimits(3, 0, 10*par_p3);
            }

            for (unsigned int incrFit = 0; incrFit < std::size(endFracs); incrFit++) {

                // The new shape uses a single fixed fraction; only the legacy one
                // actually scans the list.
                if (useOld1oPFit)  end1oPFit = endFracs[incrFit] * peak;
                else end1oPFit = 0.6 * peak;

                // Degenerate range (peak in the first bins): nothing to fit.
                if (end1oPFit <= start1oPFit) { statusFit = 1; continue; }

                ptr2 = p_forfit->Fit(&f_p, "QRS", "", start1oPFit, end1oPFit);
                // Convergence: MINUIT status 0 and a small estimated distance to
                // minimum. The chi2 test is relaxed for very few degrees of freedom,
                // where chi2/ndf is not a meaningful quality measure.
                bool converged = ptr2.Get() && ptr2->Status() == 0 && ptr2->Edm() < 1e-2;
                bool chi2ok = ptr2.Get() && (ptr2->Ndf() < 5 || ptr2->Chi2()/ptr2->Ndf() < 5);
                bool goodFitP = converged && chi2ok;


                statusFit = goodFitP ? 0: 1;

                incrFit_end = incrFit;
                if (statusFit == 0) break;  // good fit, we keep it
            }

            if (statusFit == 0) {
                // Same normalisation logic as for Ih, but here the histogram
                // integral uses the "width" option because p was scaled to unit
                // *area*, making it a density that must be compared to the integral of f.
                ROOT::Math::IntegratorOneDim intOneDim_p(f_p, ROOT::Math::IntegrationOneDim::kGAUSS);
                double intFp = intOneDim_p.Integral(start1oPFit, end1oPFit);
                if (intFp <= 0) std::cout << "ERROR > INTEGRAL FIT P IS <= 0.   ITG = " << intFp << std::endl;

                double intP = p_forfit->Integral(p_forfit->FindBin(start1oPFit), p_forfit->FindBin(end1oPFit), "width");
                SFp = (intFp > 0) ? intP / intFp: -1;
                if (SFp < 0) useFitP = false;

                // Covariance matrix for the +/-1 sigma variations of the 1/p template.
                if (ptr2.Get()) {
                    const TMatrixDSym& covP = ptr2->GetCovarianceMatrix();
                    if (covP.GetNrows() > 0) {
                        fit_p_params = ptr2->GetParams();
                        fit_p_cov    = covP.GetMatrixArray();
                        hasCovP      = true;
                    }
                }
            }

        }


        // No fraction produced a usable fit: log it and fall back to the raw
        // 1/p bin contents for this eta slice.
        if (statusFit != 0 && useFit) {
            if (saveFits) { OutputHisto->cd(); p_forfit->Write(); }
            std::cout << "Bad fit 1/p in " << p_forfit->GetName()
                    << " workerID=" << workerID << " eta=" << eta->GetBinCenter(i);
            if (ptr2.Get()) {
                std::cout << " status=" << ptr2->Status()
                        << " covMatrixStatus=" << ptr2->CovMatrixStatus()
                        << " edm=" << ptr2->Edm()
                        << " chi2/ndf=" << ptr2->Chi2() << "/" << ptr2->Ndf()
                        << " p-value=" << ptr2->Prob();
            } else {
                std::cout << " (no fit attempted: empty range)";
            }
            std::cout << std::endl;
            useFitP = false;
        }
        else {                              // Good fit
            if (saveFits) { OutputHisto->cd(); p_forfit->Write(); }
        }

        // Thresholds deciding, bin by bin, whether the data or the fit is used.
        //  - dedx_temp: above this Ih value the observed spectrum is tail-dominated
        //    and the Ih fit takes over.
        //  - mom_temp:  below this *scaled 1/p* value (i.e. above the corresponding
        //    momentum) the 1/p fit takes over. The factor 0.2 is a deliberate
        //    safety margin: the fit is only trusted well inside its fitted range,
        //    never up to its upper bound.
        float dedx_temp = (useOldIhFit)? 3.5: start_fit;
        float mom_temp = 0.2*endFracs[incrFit_end] * peak;

        // ------------- If false: no fit -------------
        // Manual override kept for debugging: uncomment to force the pure
        // bin-by-bin data convolution regardless of the fit quality.

                    //useFitIh = false;
                    //useFitP = false;

        // --------------------------------------------

        // ---- Convolution loop over the (1/p, Ih) bin pairs -------------------
        // Both loops run to N+1, i.e. the overflow bin is included on purpose:
        // dropping it would lose the highest-mass entries. (Bin 0, underflow, is
        // excluded because it corresponds to unphysical values.)
        for(int j=1;j<p->GetNbinsX()+2;j++)
        {
            for(int k=1;k<ih->GetNbinsX()+2;k++)
            {
                // "mom" is the scaled 1/p at the bin lower edge, NOT a momentum:
                // a small mom means a large physical momentum.
                float mom = p->GetBinLowEdge(j);
                float dedx = ih->GetBinLowEdge(k);
                double c_p = p->GetBinContent(j);
                double c_ih = ih->GetBinContent(k);
                float pLowEdge = p->GetBinLowEdge(j);
                float pUpEdge = p->GetBinLowEdge(j+1);
                float dedxLowEdge = ih->GetBinLowEdge(k);
                float dedxUpEdge = ih->GetBinLowEdge(k+1);

                double weight = 0;

                float mom_GeV = 0;
                float mass = -1;
                int bin_mass = 0;

                // When a fit replaces the data, the bin is sub-divided in 5 slices
                // so that the strongly varying mass mapping is sampled finely
                // enough inside a single (coarse) bin.
                float dedx_sampling = (dedxUpEdge-dedxLowEdge)/5.;
                float mom_sampling = (pUpEdge-pLowEdge)/5.;

                
                // --- Case 1: Ih taken from the fit -------------------------------
                // Conditions: the bin is sparsely populated (< 100 entries), it sits
                // in the Ih tail, and the Ih fit succeeded for this eta slice.
                if(c_ih < 100 && dedx > dedx_temp && useFitIh) {
                    for(double divdedx=dedxLowEdge; divdedx<dedxUpEdge; divdedx+=dedx_sampling){
                        
                        // Fitted yield in the sub-slice, brought back to the observed
                        // normalisation by SFih.
                        c_ih = f_ih3->Integral(divdedx,divdedx+dedx_sampling);
                        c_ih *= SFih;
                        if(c_ih==0) continue;
                        
                        // Systematic variation of the Ih template: shift the fitted
                        // yield by +/-1 sigma, propagated from the MINUIT covariance
                        // matrix. The 5e-2 argument is the relative tolerance of the
                        // numerical error integration.
                        if (fit_ih_err != 1 && hasCov) {
                            const double dc_ih = SFih * f_ih3->IntegralError(divdedx, divdedx+dedx_sampling,
                                                                            fit_ih_params, fit_ih_cov, 5e-2);
                            c_ih += (fit_ih_err == 2) ? dc_ih: -dc_ih;
                        }
                        // A downward variation must never produce a negative yield.
                        if (c_ih < 0) c_ih = 0;

                        // --- Case 1a: Ih fit AND 1/p fit (25 sub-samples per bin) --
                        if(mom < mom_temp && mom > 0 && useFitP){
                            for(double divmom=pLowEdge; divmom<pUpEdge; divmom+=mom_sampling){
                                c_p = f_p.Integral(divmom, divmom + mom_sampling);
                                c_p *= SFp;
                                if (c_p == 0) continue;

                                // Same +/-1 sigma treatment for the 1/p template.
                                if (fit_p_err != 1 && hasCovP) {
                                    const double dc_p = SFp * f_p.IntegralError(divmom, divmom + mom_sampling,
                                                                                fit_p_params, fit_p_cov, 5e-2);
                                    c_p += (fit_p_err == 2) ? dc_p: -dc_p;
                                }
                                if (c_p < 0) c_p = 0;
                                
                                // The two variables are assumed uncorrelated inside
                                // an eta slice: the joint weight is the product.
                                weight = c_ih * c_p;
                                
                                // Evaluate both variables at the centre of the
                                // sub-slice, then convert 1/p back to GeV.
                                dedx = divdedx+dedx_sampling/2.;
                                mom_GeV = kInvPScale/(divmom+mom_sampling/2.);
                                mass = GetMass(mom_GeV,dedx,K,C);
                                if (mass < 0) continue;   // Ih below the offset C

                                bin_mass = pred_mass->FindBin(mass);
                                pred_mass->SetBinContent(bin_mass,pred_mass->GetBinContent(bin_mass)+weight);
                                pred_mass_eta->SetBinContent(bin_mass,i,pred_mass_eta->GetBinContent(bin_mass,i)+weight);

                                // NaN guard: a NaN here would silently propagate to
                                // the whole prediction, so it is reported loudly.
                                if( std::isnan(pred_mass->GetBinContent(bin_mass)+weight)) std::cout << "ERROR: BIN CONTENT SET IS NAN ! 1" << std::endl;
                            }
                        }
                        // --- Case 1b: Ih from the fit, 1/p from the data ----------
                        else{
                            c_p = p->GetBinContent(j);
                            weight = c_ih * c_p;
                            dedx = divdedx+dedx_sampling/2.;
                            mom_GeV = kInvPScale/p->GetBinCenter(j);
                            mass = GetMass(mom_GeV,dedx,K,C);
                            if (mass < 0) continue;

                            bin_mass = pred_mass->FindBin(mass);
                            pred_mass->SetBinContent(bin_mass,pred_mass->GetBinContent(bin_mass)+weight);
                            pred_mass_eta->SetBinContent(bin_mass,i,pred_mass_eta->GetBinContent(bin_mass,i)+weight);

                            if( std::isnan(pred_mass->GetBinContent(bin_mass)+weight)) std::cout << "ERROR: BIN CONTENT SET IS NAN ! 2" << std::endl;
                        }
                    }
                }
                else{
                    // --- Case 2: 1/p from the fit, Ih from the data ---------------
                    if(mom < mom_temp && mom > 0 && useFitP){
                        for(double divmom=pLowEdge; divmom<pUpEdge; divmom+=mom_sampling){
                            c_p = f_p.Integral(divmom, divmom + mom_sampling);
                            c_p *= SFp;
                            if (c_p == 0) continue;

                            if (fit_p_err != 1 && hasCovP) {
                                const double dc_p = SFp * f_p.IntegralError(divmom, divmom + mom_sampling,
                                                                            fit_p_params, fit_p_cov, 5e-2);
                                c_p += (fit_p_err == 2) ? dc_p: -dc_p;
                            }
                            if (c_p < 0) c_p = 0;
                            
                            weight = c_ih * c_p;
                            dedx = ih->GetBinCenter(k);
                            mom_GeV = kInvPScale/(divmom+mom_sampling/2.);
                            mass = GetMass(mom_GeV,dedx,K,C);
                            if (mass < 0) continue;

                            bin_mass = pred_mass->FindBin(mass);
                            pred_mass->SetBinContent(bin_mass,pred_mass->GetBinContent(bin_mass)+weight);
                            pred_mass_eta->SetBinContent(bin_mass,i,pred_mass_eta->GetBinContent(bin_mass,i)+weight);

                            if( std::isnan(pred_mass->GetBinContent(bin_mass)+weight)) std::cout << "ERROR: BIN CONTENT SET IS NAN ! 3" << std::endl;
                        }
                    }
                    // --- Case 3: pure data x data (the bulk of the spectrum) ------
                    else{
                        c_p = p->GetBinContent(j);
                        c_ih = ih->GetBinContent(k);
                        weight = c_ih * c_p;
                        
                        dedx = ih->GetBinCenter(k);
                        mom_GeV = kInvPScale/p->GetBinCenter(j);
                        mass = GetMass(mom_GeV,dedx,K,C);
                        if (mass < 0) continue;

                        bin_mass = pred_mass->FindBin(mass);
                        pred_mass->SetBinContent(bin_mass,pred_mass->GetBinContent(bin_mass)+weight);
                        pred_mass_eta->SetBinContent(bin_mass,i,pred_mass_eta->GetBinContent(bin_mass,i)+weight);
                        if(std::isnan(pred_mass->GetBinContent(bin_mass)+weight)) std::cout << "ERROR: BIN CONTENT SET IS NAN ! 4" << std::endl;
                    }
                }
            }
        }
    }
    // The projection was detached from any directory, so it must be freed by hand.
    delete eta;

    if (saveFits && OutputHisto) {
        OutputHisto->Write();
        OutputHisto->Close();
        delete OutputHisto;
    }
}