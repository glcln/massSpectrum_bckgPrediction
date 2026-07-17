#include <TCanvas.h>
#include <TLegend.h>
#include "TFile.h"
#include "TH1.h"
#include "TDirectory.h"
#include <TRatioPlot.h>
#include <THStack.h>
#include "CommonTools/UtilAlgos/interface/TFileService.h"

using namespace std::placeholders;
gErrorIgnoreLevel = kFatal;

// SETUP

float K_data2017 = 2.54, C_data2017 = 3.14;
float K_data2018 = 2.55, C_data2018 = 3.14;
float K_data2024 = 2.8202, C_data2024 = 2.9784;

float K_mc2017 = 2.48, C_mc2017 = 3.19;
float K_mc2018 = 2.49, C_mc2018 = 3.19;

//Systematic error due to the background estimate method
float systErr_ = 0.; //set to 0 for systematic studies




// Scale the 1D-histogram given to the unit 
void scale(TH1F* h) {
    h->Scale(1./h->Integral(0,h->GetNbinsX()+1));
}


void corrIh(TH2F* ih_eta) {
    TF1 f_correlation_Ih_Fpix("f_correlation_Ih_Fpix","pol1",2,10);
    f_correlation_Ih_Fpix.SetParameter(0, 1.15515);
    f_correlation_Ih_Fpix.SetParameter(1, -0.0464627);

    for(int bin_eta=0; bin_eta<ih_eta->GetNbinsX(); bin_eta++){
        for(int bin_ih=0; bin_ih<ih_eta->GetNbinsY(); bin_ih++){
            ih_eta->SetBinContent(bin_eta, bin_ih, 
                                  ih_eta->GetBinContent(bin_eta, bin_ih) / f_correlation_Ih_Fpix.Eval(ih_eta->GetYaxis()->GetBinCenter(bin_ih)));
        }
    }
}

void blindMass(TH1F* h_m, float mass_value=300) {
    for(int i=0; i<h_m->GetNbinsX()+2; i++){
        if(h_m->GetBinLowEdge(i)>=mass_value) {
            h_m->SetBinContent(i,0);
            h_m->SetBinError(i,0);
        }
    }
}

float GetMass(float P, float I, float dEdxK, float dEdxC) {
  float& K = dEdxK;
  float& C = dEdxC;

  if (I - C < 0)
    return -1;
  return sqrt((I - C) / K) * P;
}


// class using to definite signal and control regions. 
class Region{
    public:
        Region();
        Region(TFileDirectory &dir, std::string suffix, int& etabins, int& ihbins, int& pbins, int& massbins);
        ~Region();
        void initHisto(TFileDirectory &dir, int etabins, int ihbins, int pbins, int massbins);
        void fill(float& eta, float& p, float& pt, float& pterr, float& ih, float& ias, float& m, float& w);
        void fillPredMass(const std::string& filename,
                          const std::string& st,
                          const std::string& st_sample,
                          TF1& f_p,
                          TF1& f_ih,
                          const bool useFit,
                          const int& fit_ih_err = 1,
                          const int& fit_p_err = 1,
                          float weight_ = -1,
                          bool useOldIhFit = false,
                          bool useOld1oPFit = false,
                          const std::string& etaName = "",
                          bool saveFits = false,
                          const int rebinp = 1,
                          const float MyIhCut = C_data2024,
                          const double par_p2 = 4.70839,
                          const double par_p3 = 1.05005,
                          const UInt_t workerID = 0);
        void write();

        float K_;
        float C_;


        int np;
        float plow;
        float pup;
        int npt;
        float ptlow;
        float ptup;
        int nih;
        float ihlow;
        float ihup;
        int nias;
        float iaslow;
        float iasup;
        int neta;
        float etalow;
        float etaup;
        int nmass;
        float masslow;
        float massup;
        std::vector<double> VectOfBins_P_;
        std::string suffix_;
        TH2F* eta_p;
        TH2F* ih_eta;
        TH2F* ih_p;
        TH2F* ias_p;
        TH2F* ias_pt;
        TH1F* mass;
        TH1F* pred_mass;
        TH2F* mass_eta;
        TH2F* mass_ih;
        TH2F* pred_mass_eta;
        TH1F* pred_mass_fitIh;
        TH1F* pred_mass_fitP;
        TH1F* pred_mass_fitIh_fitP;
        TH1F* pred_mass_noFit;
        TH2F* pt_pterroverpt;
        TH2F* ih_p_cross1D;
        TH2F* ih_p_cross1D_fit;
        TH2F* ih_p_cross1D_corr;
};


Region::Region(){}


Region::Region(TFileDirectory &dir, std::string suffix, int& etabins, int& ihbins, int& pbins, int& massbins) {
    suffix_ = suffix;
    initHisto(dir,etabins,ihbins,pbins,massbins);
} 

Region::~Region(){}


void Region::initHisto(TFileDirectory &dir, int etabins, int ihbins, int pbins, int massbins) {
    TH1::SetDefaultSumw2(kTRUE);
    TH2::SetDefaultSumw2(kTRUE);
    TH3::SetDefaultSumw2(kTRUE);
    np = pbins;
    plow = 0;
    pup = 10000;
    npt = pbins;
    ptlow = 0;
    ptup = 10000; 
    nih = ihbins;
    ihlow = 0;
    ihup = 20;
    nias = ihbins;
    iaslow = 0;
    iasup = 1;
    neta = etabins;
    etalow = -3;
    etaup = 3;
    nmass = massbins;
    masslow = 0;
    massup = 4000;
    std::string suffix = suffix_;

    eta_p = dir.make<TH2F>(("eta_p"+suffix).c_str(),";p [GeV];#eta",np,plow,pup,neta,etalow,etaup); 
    ih_eta = dir.make<TH2F>(("ih_eta"+suffix).c_str(),";#eta;I_{h} [MeV/cm]",neta,etalow,etaup,nih,ihlow,ihup); 
    ih_p = dir.make<TH2F>(("ih_p"+suffix).c_str(),";p [GeV];I_{h} [MeV/cm]",np,plow,pup,nih,ihlow,ihup);
    ih_p_cross1D = dir.make<TH2F>(("ih_p_cross1D"+suffix).c_str(),";p [GeV];I_{h} [MeV/cm]",np,plow,pup,nih,ihlow,ihup);
    ih_p_cross1D_fit = dir.make<TH2F>(("ih_p_cross1D_fit"+suffix).c_str(),";p [GeV];I_{h} [MeV/cm]",np,plow,pup,nih,ihlow,ihup);
    ih_p_cross1D_corr = dir.make<TH2F>(("ih_p_cross1D_corr"+suffix).c_str(),";p [GeV];I_{h} [MeV/cm]",np,plow,pup,nih,ihlow,ihup);
    ias_p = dir.make<TH2F>(("ias_p"+suffix).c_str(),";p [GeV];I_{as}",np,plow,pup,nias,iaslow,iasup); 
    ias_pt = dir.make<TH2F>(("ias_pt"+suffix).c_str(),";pt [GeV];I_{as}",npt,ptlow,ptup,nias,iaslow,iasup);
    mass = dir.make<TH1F>(("mass"+suffix).c_str(),";Mass [GeV]",nmass,masslow,massup); 
    pred_mass = dir.make<TH1F>(("pred_mass"+suffix).c_str(),";Mass [GeV]",nmass,masslow,massup);
    mass_eta = dir.make<TH2F>(("mass_eta"+suffix).c_str(),";Mass [GeV];#eta",nmass,masslow,massup,neta,etalow,etaup);
    mass_ih = dir.make<TH2F>(("mass_ih"+suffix).c_str(),";I_{h} [MeV/cm];Mass [GeV]",nih,ihlow,ihup, nmass,masslow,massup);
    pred_mass_eta = dir.make<TH2F>(("pred_mass_eta"+suffix).c_str(),";Mass [GeV];#eta",nmass,masslow,massup,neta,etalow,etaup); 
    
    mass->SetBinErrorOption(TH1::EBinErrorOpt::kPoisson);
    pred_mass->SetBinErrorOption(TH1::EBinErrorOpt::kPoisson);
    pt_pterroverpt = dir.make<TH2F>(("pt_pterroverpt"+suffix).c_str(),";p_{T} [GeV];#frac{#sigma_{pT}}{p_{T}}",npt,ptlow,ptup,100,0,1);
}


void Region::fill(float& eta, float& p, float& pt, float& pterr, float& ih, float& ias, float& m, float& w) {
   eta_p->Fill(p,eta,w);
   ih_eta->Fill(eta,ih,w);
   ih_p->Fill(p,ih,w);
   ias_p->Fill(p,ias,w);
   ias_pt->Fill(pt,ias,w);
   mass->Fill(m,w);
   mass_eta->Fill(m,eta,w);
   mass_ih->Fill(m,ih,w);
   pt_pterroverpt->Fill(pt,pterr/pt,w);
}


void Region::fillPredMass(const std::string& filename,
                          const std::string& st,
                          const std::string& st_sample,
                          TF1& f_p,
                          TF1& f_ih,
                          const bool useFit,
                          const int& fit_ih_err = 1,
                          const int& fit_p_err = 1,
                          float weight_ = -1,
                          bool useOldIhFit = false,
                          bool useOld1oPFit = false,
                          const std::string& etaName = "",
                          bool saveFits = false,
                          const int rebinp = 1,
                          const float MyIhCut = C_data2024,
                          const double par_p2 = 4.70839,
                          const double par_p3 = 1.05005,
                          const UInt_t workerID = 0) {

    // Debug Fit
    TFile* OutputHisto = nullptr;
    std::string filenameOutputFit = "DebugFit/Fits_" + filename + "_" + st + ((useOldIhFit || useOld1oPFit) ? "_OldFit" : "_NewFit") + etaName + "_" + std::to_string(workerID) + '_' + std::to_string(MyIhCut) +  ".root";
    if (saveFits) {
        OutputHisto = new TFile(filenameOutputFit.c_str(), "RECREATE");
        OutputHisto->cd();
    }


    // Setup
    TH1F* eta = (TH1F*) ih_eta->ProjectionX();

    float K = 2.27, C = 3.16;
    if(st_sample=="data2017"){K = K_data2017; C = C_data2017;}
    else if(st_sample=="data2018"){K = K_data2018; C = C_data2018;}
    else if (st_sample=="data2024"){K = K_data2024; C = C_data2024;}
    else if(st_sample=="mc2017"){K = K_mc2017; C = C_mc2017;}
    else if(st_sample=="mc2018"){K = K_mc2018; C = C_mc2018;}

    bool useFitIh = true;
    bool useFitP = true;
    ROOT::Math::IntegratorOneDimOptions::SetDefaultRelTolerance(1.E-9);


    // Loop over the eta bins
    for(int i=1;i<eta->GetNbinsX()+1;i++)
    {
        // Setup
        useFitIh = useFit;
        useFitP = useFit;
        TH1F* p  = (TH1F*) eta_p->ProjectionX(Form("proj_p_eta%d", i), i, i, "e");
        TH1F* ih = (TH1F*) ih_eta->ProjectionY(Form("proj_ih_eta%d", i), i, i, "e");

        scale(p); //only scale one of the two distributions ih or p --> keep the information of the normalisation 
        if(ih->GetEntries() < 1 || p->GetEntries() < 1) continue;

        float endIhFit = 6., end1oPFit = 30., start1oPFit = 5;
        

        // Ih fit
        TFitResultPtr ptr1 = 0;
        float max_ih = ih->GetBinCenter(ih->GetMaximumBin());
        float start_fit = (useOldIhFit)? 3 : 1.1*max_ih; // 1.1*max_ih; for gauss
        int lastBinContent = ih->GetNbinsX();
        while(ih->GetBinContent(lastBinContent)==0) lastBinContent--;
        if(start_fit > ih->GetBinCenter(lastBinContent)) start_fit = max_ih;

        if (useFitIh) {
            ptr1 = ih->Fit(&f_ih, "QRS", "", start_fit, endIhFit);
            bool goodFit = ptr1.Get() && ptr1->Ndf() > 0 && ptr1->Chi2()/ptr1->Ndf() < 6;
            if (!goodFit) { // Bad fit
                
                cout << "Bad fit Ih in " << ih->GetName() << " workerID=" << std::to_string(workerID)
                          << "    status=" << ptr1->Status()
                          << " covMatrixStatus=" << ptr1->CovMatrixStatus()
                          << " isValid=" << ptr1->IsValid()
                          << " edm=" << ptr1->Edm()
                          << " chi2/ndf=" << ptr1->Chi2() << "/" << ptr1->Ndf()
                          << " eta=" << eta->GetBinCenter(i)
                          << " p-value=" << ptr1->Prob() << endl;
                if (saveFits) { OutputHisto->cd(); ih->Write(); }
                useFitIh = false;
            }
            else {                             // Good fit
                if (saveFits) { OutputHisto->cd(); ih->Write(); }
            }
        }

        TF1* const f_ih2 = &f_ih;
        float intFih = f_ih2->Integral(3, endIhFit);
        float intIh = ih->Integral(ih->FindBin(3), ih->FindBin(endIhFit));

        float SFih = (intFih > 0)? intIh/intFih : -1;
        if(SFih < 0 && useFitIh) {
            cout<<"ERROR > INTEGRAL FIT IH IS <= 0.   ITG = " << intFih << " FOR ETA BIN #" << i << std::endl;
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
        float SFp = 0;
        TF1* f_p3 = &f_p;
        TFitResultPtr ptr2 = 0;
        int statusFit = 1;

        start1oPFit = 0;
        // Taking the fit from the Down variation, as it is the one with the thicker bins, thus less statistical fluctuations for the fit to converge
        TH1F *p_forfit = (TH1F*) p->Clone(Form("forfit_p_eta%d", i));
        p_forfit->SetDirectory(nullptr);

        // for (int b = 1; b <= p_forfit->GetNbinsX(); b++) {
        //     if (p_forfit->GetBinContent(b) > 0) {
        //         start1oPFit = p_forfit->GetBinCenter(b);
        //         break;
        //     }
        // }

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
                TF1* f_p2 = &f_p;
                ROOT::Math::IntegratorOneDim intOneDim_p(*f_p2, ROOT::Math::IntegrationOneDim::kGAUSS);
                float intFp = intOneDim_p.Integral(start1oPFit, end1oPFit);
                if (intFp <= 0) cout << "ERROR > INTEGRAL FIT P IS <= 0.   ITG = " << intFp << std::endl;

                float intP = p_forfit->Integral(p_forfit->FindBin(start1oPFit), p_forfit->FindBin(end1oPFit), "width");
                SFp = (intFp > 0) ? intP / intFp : -1;
                if (SFp < 0) useFitP = false;
                f_p3 = f_p2;
            }
        }


        if (statusFit!=0 && useFit!=0) {      // Bad fit
            if (saveFits) { OutputHisto->cd(); p_forfit->Write(); }
            cout << "Bad fit 1/p in " << p_forfit->GetName() << " workerID=" << std::to_string(workerID)
                 << "    status=" << ptr2->Status()
                 << " covMatrixStatus=" << ptr2->CovMatrixStatus()
                 << " isValid=" << ptr2->IsValid()
                 << " edm=" << ptr2->Edm()
                 << " chi2/ndf=" << ptr2->Chi2() << "/" << ptr2->Ndf()
                 << " eta=" << eta->GetBinCenter(i)
                 << " p-value=" << ptr2->Prob() << endl;
            cout << " param: " << f_p.GetParameter(0) << " " << f_p.GetParameter(1) << " " << f_p.GetParameter(2) << " " << f_p.GetParameter(3) << endl;
            useFitP = false;
        }
        else {                              // Good fit
            if (saveFits) { OutputHisto->cd(); p_forfit->Write(); }
        }

        ROOT::Math::IntegratorOneDim* intOneDimFP3 = nullptr;
        if (useFitP) {
            intOneDimFP3 = new ROOT::Math::IntegratorOneDim(*f_p3, ROOT::Math::IntegrationOneDim::kGAUSS);
        }
        float dedx_temp = (useOldIhFit)? 3.5 : start_fit;
        float mom_temp = 0.2*endFracs[incrFit_end] * peak;

        // ------------- If false: no fit -------------

                    //useFitIh = false;
                    //useFitP = false;

        // --------------------------------------------
        int atLeastOne=0;

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

                float invMom = 0;
                float mass = -1;
                int bin_mass = 0;

                float dedx_sampling = (dedxUpEdge-dedxLowEdge)/5.;
                float mom_sampling = (pUpEdge-pLowEdge)/5.;

                ih_p_cross1D->SetBinContent(j,k,ih_p_cross1D->GetBinContent(j,k)+(p->GetBinContent(j)*ih->GetBinContent(k)));

                // use Ih fit
                if(c_ih < 100 && dedx > dedx_temp && useFitIh) {
                    for(double divdedx=dedxLowEdge; divdedx<dedxUpEdge; divdedx+=dedx_sampling){
                        
                        c_ih = f_ih3->Integral(divdedx,divdedx+dedx_sampling);
                        c_ih *= SFih;
                        if(c_ih==0) continue;
                        
                        if (fit_ih_err == 0 && hasCov) c_ih -= (SFih*f_ih3->IntegralError(divdedx, divdedx+dedx_sampling, fit_ih_params, fit_ih_cov, 5e-2));
                        if (fit_ih_err == 2 && hasCov) c_ih += (SFih*f_ih3->IntegralError(divdedx, divdedx+dedx_sampling, fit_ih_params, fit_ih_cov, 5e-2));
                       
                        // use Ih fit AND 1/p fit
                        if(mom < mom_temp && mom > 0 && useFitP){
                            for(double divmom=pLowEdge; divmom<pUpEdge; divmom+=mom_sampling){
                                c_p = intOneDimFP3->Integral(divmom,divmom+mom_sampling);
                                c_p *= SFp;
                                if(c_p==0) continue;

                                if (fit_p_err == 0) c_p -= (SFp*intOneDimFP3->Error());
                                if (fit_p_err == 2) c_p += (SFp*intOneDimFP3->Error());
                                
                                weight = c_ih * c_p;
                                
                                dedx = divdedx+dedx_sampling/2.;
                                invMom = 10000./(divmom+mom_sampling/2.);
                                mass = GetMass(invMom,dedx,K,C);
                                bin_mass = pred_mass->FindBin(mass);
                                if (dedx >= MyIhCut) pred_mass->SetBinContent(bin_mass,pred_mass->GetBinContent(bin_mass)+weight);
                                pred_mass_eta->SetBinContent(i,bin_mass,pred_mass_eta->GetBinContent(i,bin_mass)+weight);
                                if( std::isnan(pred_mass->GetBinContent(bin_mass)+weight)) cout << "ERROR : BIN CONTENT SET IS NAN ! 1" << std::endl;

                                pred_mass_fitIh_fitP->SetBinContent(bin_mass,pred_mass_fitIh_fitP->GetBinContent(bin_mass)+weight);
                                ih_p_cross1D_fit->SetBinContent(j,k,ih_p_cross1D_fit->GetBinContent(j,k)+weight);
                            }
                        }
                        else{
                            c_p = p->GetBinContent(j);
                            weight = c_ih * c_p;
                            dedx = divdedx+dedx_sampling/2.;
                            invMom = 10000./p->GetBinCenter(j);
                            mass = GetMass(invMom,dedx,K,C);
                            bin_mass = pred_mass->FindBin(mass);
                            if (dedx >= MyIhCut) pred_mass->SetBinContent(bin_mass,pred_mass->GetBinContent(bin_mass)+weight);
                            pred_mass_eta->SetBinContent(i,bin_mass,pred_mass_eta->GetBinContent(i,bin_mass)+weight);
                            if( std::isnan(pred_mass->GetBinContent(bin_mass)+weight)) cout << "ERROR : BIN CONTENT SET IS NAN ! 2" << std::endl;
                            
                            pred_mass_fitIh->SetBinContent(bin_mass,pred_mass_fitIh->GetBinContent(bin_mass)+weight);
                            ih_p_cross1D_fit->SetBinContent(j,k,ih_p_cross1D_fit->GetBinContent(j,k)+weight);
                        }
                    }
                }
                else{
                    // use 1/p fit
                    if(mom < mom_temp && mom > 0 && useFitP){
                        //if (c_ih < 100 && dedx > dedx_temp && atLeastOne==0) { cout << "no Ih but P " << useFitIh << " " << SFih << endl; atLeastOne++;}
                        for(double divmom=pLowEdge; divmom<pUpEdge; divmom+=mom_sampling){
                            c_p = intOneDimFP3->Integral(divmom,divmom+mom_sampling);
                            c_p *= SFp;
                            if(c_p==0) continue;
                            if (fit_p_err == 0) c_p -= (SFp*intOneDimFP3->Error());
                            if (fit_p_err == 2) c_p += (SFp*intOneDimFP3->Error());
                            
                            weight = c_ih * c_p;
                            dedx = ih->GetBinCenter(k);
                            invMom = 10000./(divmom+mom_sampling/2.);
                            mass = GetMass(invMom,dedx,K,C);
                            bin_mass = pred_mass->FindBin(mass);
                            if (dedx >= MyIhCut) pred_mass->SetBinContent(bin_mass,pred_mass->GetBinContent(bin_mass)+weight);
                            pred_mass_eta->SetBinContent(i,bin_mass,pred_mass_eta->GetBinContent(i,bin_mass)+weight);
                            if( std::isnan(pred_mass->GetBinContent(bin_mass)+weight)) cout << "ERROR : BIN CONTENT SET IS NAN ! 3" << std::endl;
                            
                            pred_mass_fitP->SetBinContent(bin_mass,pred_mass_fitP->GetBinContent(bin_mass)+weight);
                            ih_p_cross1D_fit->SetBinContent(j,k,ih_p_cross1D_fit->GetBinContent(j,k)+weight);
                        }
                    }
                    else{
                        c_p = p->GetBinContent(j);
                        c_ih = ih->GetBinContent(k);
                        weight = c_ih * c_p;
                        
                        dedx = ih->GetBinCenter(k);
                        invMom = 10000./p->GetBinCenter(j);
                        mass = GetMass(invMom,dedx,K,C);
                        bin_mass = pred_mass->FindBin(mass);
                        
                        if (dedx >= MyIhCut) pred_mass->SetBinContent(bin_mass,pred_mass->GetBinContent(bin_mass)+weight);
                        pred_mass_eta->SetBinContent(i,bin_mass,pred_mass_eta->GetBinContent(i,bin_mass)+weight);
                        if(std::isnan(pred_mass->GetBinContent(bin_mass)+weight)) cout << "ERROR : BIN CONTENT SET IS NAN ! 4" << std::endl;
                        
                        pred_mass_noFit->SetBinContent(bin_mass,pred_mass_noFit->GetBinContent(bin_mass)+weight);
                        ih_p_cross1D_fit->SetBinContent(j,k,ih_p_cross1D_fit->GetBinContent(j,k)+weight);
                    }
                }
            }
        }
        delete p;
        delete ih;
    }
    delete eta;

    if (saveFits && OutputHisto) {
        OutputHisto->Write();
        OutputHisto->Close();
        delete OutputHisto;
    }
}


void Region::write() {
    eta_p->Write();
    ih_eta->Write();
    ih_p->Write();
    ih_p_cross1D->Write();
    ih_p_cross1D_fit->Write();
    ias_p->Write();
    ias_pt->Write();
    mass->Write();
    pred_mass->Write();
    mass_eta->Write();
    mass_ih->Write();
    pred_mass_eta->Write();
    pred_mass_fitIh->Write();
    pred_mass_fitP->Write();
    pred_mass_fitIh_fitP->Write();
    pred_mass_noFit->Write();
    pt_pterroverpt->Write();
}
