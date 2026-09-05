#pragma once

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

constexpr float kInvPScale = 10000.f;

// SETUP

struct DeDxCalib { float K, C; };

inline DeDxCalib GetDeDxCalib(const std::string& sample) {
    static const std::map<std::string, DeDxCalib> kCalib = {
        {"data2017", {2.54f,    3.14f}},
        {"data2018", {2.55f,    3.14f}},
        {"data2024", {2.8202f,  2.9784f}},
        {"mc2017",   {2.48f,    3.19f}},
        {"mc2018",   {2.49f,    3.19f}},
        {"mc2024",   {2.83894f, 3.01756f}},
    };
    auto it = kCalib.find(sample);
    if (it == kCalib.end()) {
        std::cerr << "GetDeDxCalib: unknown sample '" << sample << "', will take 2024 param as default ones " << std::endl;
        return {2.8202f, 2.9784f};
    }
    return it->second;
}

//Systematic error due to the background estimate method
constexpr float systErr_ = 0.; //set to 0 for systematic studies




// Scale the 1D-histogram given to the unit 
void scale(TH1D* h) {
    const double itg = h->Integral(0, h->GetNbinsX()+1);
    if (itg <= 0) { std::cerr << "scale: integrale nulle pour " << h->GetName() << std::endl; return; }
    h->Scale(1./itg);
}


void corrIh(TH2D* ih_eta, const std::string& etaName) {
    constexpr double kIhLow = 2., kIhUp = 10.;   // domaine sur lequel le pol1 a ete ajuste
    TF1 f_correlation_Ih_Fpix("f_correlation_Ih_Fpix", "pol1", kIhLow, kIhUp);

    float par0 = 1.03878, par1 = -0.0120062;
    if      (etaName.find("Eta2p4")   != std::string::npos) {par0 = 1.03878; par1 = -0.0120062;}
    else if (etaName.find("Eta1_2p4") != std::string::npos) {par0 = 1.12024; par1 = -0.03684;}
    else if (etaName.find("Eta1")     != std::string::npos) {par0 = 1.15082; par1 = -0.0459028;}

    f_correlation_Ih_Fpix.SetParameter(0, par0);
    f_correlation_Ih_Fpix.SetParameter(1, par1);

    const TAxis* ay = ih_eta->GetYaxis();

    // La correction ne depend que de Ih : une evaluation par bin en Ih
    for (int bin_ih = 0; bin_ih <= ih_eta->GetNbinsY() + 1; ++bin_ih) {

        // Hors du domaine du fit : extrapolation constante (valeur au bord),
        // pas d'extrapolation lineaire qui derive vite.
        const double x    = std::clamp(ay->GetBinCenter(bin_ih), kIhLow, kIhUp);
        const double corr = f_correlation_Ih_Fpix.Eval(x);

        if (corr <= 0.) {
            std::cerr << "corrIh: facteur <= 0 (" << corr << ") a Ih=" << x
                      << " -> bin laisse non corrige" << std::endl;
            continue;
        }

        for (int bin_eta = 0; bin_eta <= ih_eta->GetNbinsX() + 1; ++bin_eta) {
            ih_eta->SetBinContent(bin_eta, bin_ih, ih_eta->GetBinContent(bin_eta, bin_ih) / corr);
            ih_eta->SetBinError  (bin_eta, bin_ih, ih_eta->GetBinError  (bin_eta, bin_ih) / corr);
        }
    }
}

void corr1oP(TH2D* eta_p, const std::string& etaName) {
    constexpr double k1oPLow = 0., k1oPUp = 200.;
    TF1 f_correlation_1oP_Fpix("f_correlation_1oP_Fpix", "pol1", k1oPLow, k1oPUp);

    float par0 = 0.871848, par1 = 0.0020839;
    if      (etaName.find("Eta2p4")   != std::string::npos) {par0 = 0.871848; par1 = 0.0020839;}
    else if (etaName.find("Eta1_2p4") != std::string::npos) {par0 = 0.926461; par1 = 0.0019605;}
    else if (etaName.find("Eta1")     != std::string::npos) {par0 = 0.938621; par1 = 0.000665965;}

    f_correlation_1oP_Fpix.SetParameter(0, par0);
    f_correlation_1oP_Fpix.SetParameter(1, par1);

    const TAxis* ax = eta_p->GetXaxis();

    for (int bin_1oP = 0; bin_1oP <= eta_p->GetNbinsX() + 1; ++bin_1oP) {

        const double x    = std::clamp(ax->GetBinCenter(bin_1oP), k1oPLow, k1oPUp);
        const double corr = f_correlation_1oP_Fpix.Eval(x);

        if (corr <= 0.) {
            std::cerr << "corr1oP: facteur <= 0 (" << corr << ") a 1/p=" << x
                      << " -> bin laisse non corrige" << std::endl;
            continue;
        }

        for (int bin_eta = 0; bin_eta <= eta_p->GetNbinsY() + 1; ++bin_eta) {
            eta_p->SetBinContent(bin_1oP, bin_eta, eta_p->GetBinContent(bin_1oP, bin_eta) / corr);
            eta_p->SetBinError  (bin_1oP, bin_eta, eta_p->GetBinError  (bin_1oP, bin_eta) / corr);
        }
    }
}

void blindMass(TH1D* h_m, float mass_value=300) {
    for(int i=0; i<h_m->GetNbinsX()+2; i++){
        if(h_m->GetBinLowEdge(i)>=mass_value) {
            h_m->SetBinContent(i,0);
            h_m->SetBinError(i,0);
        }
    }
}

inline float GetMass(float p, float ih, float K, float C) {
    if (ih - C < 0) return -1;
    return std::sqrt((ih - C) / K) * p;
}


// class using to definite signal and control regions. 
class Region{
    public:
        Region();
        ~Region();

        template <typename TDir>
        Region(TDir& dir, std::string suffix, int etabins, int ihbins, int pbins, int massbins) {
            suffix_ = std::move(suffix);
            initHisto(dir, etabins, ihbins, pbins, massbins);
        }

        template <typename TDir>
        void initHisto(TDir& dir, int etabins, int ihbins, int pbins, int massbins);

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

        int np;
        float plow;
        float pup;
        int npt;
        float ptlow;
        float ptup;
        int nih;
        float ihlow;
        float ihup;
        int neta;
        float etalow;
        float etaup;
        int nmass;
        float masslow;
        float massup;
        std::string suffix_;

        TH2D* eta_p         = nullptr;
        TH2D* ih_eta        = nullptr;
        TH1D* mass          = nullptr;
        TH1D* pred_mass     = nullptr;
        TH2D* mass_eta      = nullptr;
        TH2D* pred_mass_eta = nullptr;
};


Region::Region(){}


template <typename TDir>
void Region::initHisto(TDir& dir, int etabins, int ihbins, int pbins, int massbins) {
    TH1::SetDefaultSumw2(kTRUE);
    TH2::SetDefaultSumw2(kTRUE);

    np = pbins;  plow = 0;  pup = 10000;
    nih = ihbins; ihlow = 0; ihup = 20;
    neta = etabins; etalow = -3; etaup = 3;
    nmass = massbins; masslow = 0; massup = 4000;
    const std::string suffix = suffix_;

    // 'dir.template make<...>' : le 'template' est obligatoire car dir depend du parametre TDir
    eta_p         = dir.template make<TH2D>(("eta_p"        + suffix).c_str(), ";p [GeV];#eta",              np,   plow,   pup,     neta, etalow, etaup);
    ih_eta        = dir.template make<TH2D>(("ih_eta"       + suffix).c_str(), ";#eta;I_{h} [MeV/cm]",       neta, etalow, etaup,   nih,  ihlow,  ihup);
    mass          = dir.template make<TH1D>(("mass"         + suffix).c_str(), ";Mass [GeV]",                nmass, masslow, massup);
    pred_mass     = dir.template make<TH1D>(("pred_mass"    + suffix).c_str(), ";Mass [GeV]",                nmass, masslow, massup);
    mass_eta      = dir.template make<TH2D>(("mass_eta"     + suffix).c_str(), ";Mass [GeV];#eta",           nmass, masslow, massup, neta, etalow, etaup);
    pred_mass_eta = dir.template make<TH2D>(("pred_mass_eta"+ suffix).c_str(), ";Mass [GeV];#eta",           nmass, masslow, massup, neta, etalow, etaup);

    mass->SetBinErrorOption(TH1::EBinErrorOpt::kPoisson);
}

Region::~Region(){}



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

    // Debug Fit
    TFile* OutputHisto = nullptr;
    std::string filenameOutputFit = "DebugFit/Fits_" + filename + "_" + st + ((useOldIhFit || useOld1oPFit) ? "_OldFit" : "_NewFit") + etaName + "_" + std::to_string(workerID) +  ".root";
    if (saveFits) {
        OutputHisto = new TFile(filenameOutputFit.c_str(), "RECREATE");
        OutputHisto->cd();
    }


    // Setup
    TH1D* eta = (TH1D*) ih_eta->ProjectionX();
    eta->SetDirectory(nullptr);

    const float K = calib.K;
    const float C = calib.C;

    bool useFitIh = true;
    bool useFitP = true;
    ROOT::Math::IntegratorOneDimOptions::SetDefaultRelTolerance(1.E-9);


    // Loop over the eta bins
    for(int i=1;i<eta->GetNbinsX()+1;i++) {
        // Setup
        useFitIh = useFit;
        useFitP = useFit;
        std::unique_ptr<TH1D> p (static_cast<TH1D*>(eta_p ->ProjectionX(Form("proj_p_eta%d",  i), i, i, "e")));
        std::unique_ptr<TH1D> ih(static_cast<TH1D*>(ih_eta->ProjectionY(Form("proj_ih_eta%d", i), i, i, "e")));
        p->SetDirectory(nullptr);
        ih->SetDirectory(nullptr);

        if (ih->GetEntries() < 1 || p->GetEntries() < 1) continue;
        if (p->Integral(0, p->GetNbinsX()+1) <= 0) continue;
        scale(p.get());

        float endIhFit = 6., end1oPFit = 30., start1oPFit = 0;
        

        // Ih fit
        TFitResultPtr ptr1 = 0;
        float max_ih = ih->GetBinCenter(ih->GetMaximumBin());
        float start_fit = (useOldIhFit)? 3 : 1.1*max_ih; // 1.1*max_ih; for gauss
        int lastBinContent = ih->GetNbinsX();
        while (lastBinContent > 1 && ih->GetBinContent(lastBinContent) == 0) --lastBinContent;
        if(start_fit > ih->GetBinCenter(lastBinContent)) start_fit = max_ih;

        if (useFitIh) {
            ptr1 = ih->Fit(&f_ih, "QRS", "", start_fit, endIhFit);
            bool goodFit = ptr1.Get() && ptr1->Ndf() > 0 && ptr1->Chi2()/ptr1->Ndf() < 6;
            if (!goodFit) {
                std::cout << "Bad fit Ih in " << ih->GetName() << " workerID=" << workerID
                        << " eta=" << eta->GetBinCenter(i);
                if (ptr1.Get()) {
                    std::cout << " status=" << ptr1->Status()
                            << " covMatrixStatus=" << ptr1->CovMatrixStatus()
                            << " edm=" << ptr1->Edm()
                            << " chi2/ndf=" << ptr1->Chi2() << "/" << ptr1->Ndf()
                            << " p-value=" << ptr1->Prob();
                } else std::cout << " (fit non effectue)";
                std::cout << std::endl;
                if (saveFits) { OutputHisto->cd(); ih->Write(); }
                useFitIh = false;
            }
            else {                             // Good fit
                if (saveFits) { OutputHisto->cd(); ih->Write(); }
            }
        }

        TF1* const f_ih2 = &f_ih;
        double intFih = f_ih2->Integral(3, endIhFit);
        double intIh = ih->Integral(ih->FindBin(3), ih->FindBin(endIhFit));

        double SFih = (intFih > 0)? intIh/intFih : -1;
        if(SFih < 0 && useFitIh) {
            std::cout<<"ERROR > INTEGRAL FIT IH IS <= 0.   ITG = " << intFih << " FOR ETA BIN #" << i << std::endl;
            useFitIh = false;
        }
        

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


        // 1/p fit
        double SFp = 0;
        TFitResultPtr ptr2 = 0;
        bool hasCovP = false;
        const double* fit_p_params = nullptr;
        const double* fit_p_cov    = nullptr;
        int statusFit = 1;

        // Taking the fit from the Down variation, as it is the one with the thicker bins, thus less statistical fluctuations for the fit to converge
        std::unique_ptr<TH1D> p_forfit(static_cast<TH1D*>(p->Clone(Form("forfit_p_eta%d", i))));
        p_forfit->SetDirectory(nullptr);


        const float endFracs[5] = {0.9f, 0.8f, 0.7f, 0.6f, 0.5f};
        float peak = p_forfit->GetBinCenter(p_forfit->GetMaximumBin());
        int incrFit_end = 0;
        if (useFitP) {

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

                if (useOld1oPFit)  end1oPFit = endFracs[incrFit] * peak;
                else end1oPFit = 0.6 * peak;

                if (end1oPFit <= start1oPFit) { statusFit = 1; continue; }

                ptr2 = p_forfit->Fit(&f_p, "QRS", "", start1oPFit, end1oPFit);
                bool converged = ptr2.Get() && ptr2->Status() == 0 && ptr2->Edm() < 1e-2;
                bool chi2ok = ptr2.Get() && (ptr2->Ndf() < 5 || ptr2->Chi2()/ptr2->Ndf() < 5);
                bool goodFitP = converged && chi2ok;


                statusFit = goodFitP ? 0 : 1;

                incrFit_end = incrFit;
                if (statusFit == 0) break;  // good fit, we keep it
            }

            if (statusFit == 0) {
                ROOT::Math::IntegratorOneDim intOneDim_p(f_p, ROOT::Math::IntegrationOneDim::kGAUSS);
                double intFp = intOneDim_p.Integral(start1oPFit, end1oPFit);
                if (intFp <= 0) std::cout << "ERROR > INTEGRAL FIT P IS <= 0.   ITG = " << intFp << std::endl;

                double intP = p_forfit->Integral(p_forfit->FindBin(start1oPFit), p_forfit->FindBin(end1oPFit), "width");
                SFp = (intFp > 0) ? intP / intFp : -1;
                if (SFp < 0) useFitP = false;

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
                std::cout << " (aucun fit tente : plage vide)";
            }
            std::cout << std::endl;
            useFitP = false;
        }
        else {                              // Good fit
            if (saveFits) { OutputHisto->cd(); p_forfit->Write(); }
        }

        float dedx_temp = (useOldIhFit)? 3.5 : start_fit;
        float mom_temp = 0.2*endFracs[incrFit_end] * peak;

        // ------------- If false: no fit -------------

                    //useFitIh = false;
                    //useFitP = false;

        // --------------------------------------------

        // Loop over the bins in (p,ih)
        for(int j=1;j<p->GetNbinsX()+2;j++)
        {
            for(int k=1;k<ih->GetNbinsX()+2;k++)
            {
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

                float dedx_sampling = (dedxUpEdge-dedxLowEdge)/5.;
                float mom_sampling = (pUpEdge-pLowEdge)/5.;

                
                // use Ih fit
                if(c_ih < 100 && dedx > dedx_temp && useFitIh) {
                    for(double divdedx=dedxLowEdge; divdedx<dedxUpEdge; divdedx+=dedx_sampling){
                        
                        c_ih = f_ih3->Integral(divdedx,divdedx+dedx_sampling);
                        c_ih *= SFih;
                        if(c_ih==0) continue;
                        
                        if (fit_ih_err != 1 && hasCov) {
                            const double dc_ih = SFih * f_ih3->IntegralError(divdedx, divdedx+dedx_sampling,
                                                                            fit_ih_params, fit_ih_cov, 5e-2);
                            c_ih += (fit_ih_err == 2) ? dc_ih : -dc_ih;
                        }
                        if (c_ih < 0) c_ih = 0;

                        // use Ih fit AND 1/p fit
                        if(mom < mom_temp && mom > 0 && useFitP){
                            for(double divmom=pLowEdge; divmom<pUpEdge; divmom+=mom_sampling){
                                c_p = f_p.Integral(divmom, divmom + mom_sampling);
                                c_p *= SFp;
                                if (c_p == 0) continue;

                                if (fit_p_err != 1 && hasCovP) {
                                    const double dc_p = SFp * f_p.IntegralError(divmom, divmom + mom_sampling,
                                                                                fit_p_params, fit_p_cov, 5e-2);
                                    c_p += (fit_p_err == 2) ? dc_p : -dc_p;
                                }
                                if (c_p < 0) c_p = 0;
                                
                                weight = c_ih * c_p;
                                
                                dedx = divdedx+dedx_sampling/2.;
                                mom_GeV = kInvPScale/(divmom+mom_sampling/2.);
                                mass = GetMass(mom_GeV,dedx,K,C);
                                if (mass < 0) continue;

                                bin_mass = pred_mass->FindBin(mass);
                                pred_mass->SetBinContent(bin_mass,pred_mass->GetBinContent(bin_mass)+weight);
                                pred_mass_eta->SetBinContent(bin_mass,i,pred_mass_eta->GetBinContent(bin_mass,i)+weight);

                                if( std::isnan(pred_mass->GetBinContent(bin_mass)+weight)) std::cout << "ERROR : BIN CONTENT SET IS NAN ! 1" << std::endl;
                            }
                        }
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

                            if( std::isnan(pred_mass->GetBinContent(bin_mass)+weight)) std::cout << "ERROR : BIN CONTENT SET IS NAN ! 2" << std::endl;
                        }
                    }
                }
                else{
                    // use 1/p fit
                    if(mom < mom_temp && mom > 0 && useFitP){
                        for(double divmom=pLowEdge; divmom<pUpEdge; divmom+=mom_sampling){
                            c_p = f_p.Integral(divmom, divmom + mom_sampling);
                            c_p *= SFp;
                            if (c_p == 0) continue;

                            if (fit_p_err != 1 && hasCovP) {
                                const double dc_p = SFp * f_p.IntegralError(divmom, divmom + mom_sampling,
                                                                            fit_p_params, fit_p_cov, 5e-2);
                                c_p += (fit_p_err == 2) ? dc_p : -dc_p;
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

                            if( std::isnan(pred_mass->GetBinContent(bin_mass)+weight)) std::cout << "ERROR : BIN CONTENT SET IS NAN ! 3" << std::endl;
                        }
                    }
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
                        if(std::isnan(pred_mass->GetBinContent(bin_mass)+weight)) std::cout << "ERROR : BIN CONTENT SET IS NAN ! 4" << std::endl;
                    }
                }
            }
        }
    }
    delete eta;

    if (saveFits && OutputHisto) {
        OutputHisto->Write();
        OutputHisto->Close();
        delete OutputHisto;
    }
}
