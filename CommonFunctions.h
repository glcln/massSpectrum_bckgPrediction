#include <TCanvas.h>
#include <TLegend.h>
#include "TFile.h"
#include "TH1.h"
#include "TDirectory.h"
#include <TRatioPlot.h>
#include <THStack.h>

#include "Regions.h"

gErrorIgnoreLevel = kFatal;

TH2F* RebinTH2Y_varBins(TH2F* h, int nEta, double* eEta) {
    int nX = h->GetNbinsX();
    const TArrayD* xArr = h->GetXaxis()->GetXbins();
    TH2F* hNew;
    if (xArr->GetSize() > 0)
        hNew = new TH2F(h->GetName(), h->GetTitle(),
                        nX, xArr->GetArray(),
                        nEta, eEta);
    else
        hNew = new TH2F(h->GetName(), h->GetTitle(),
                        nX, h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax(),
                        nEta, eEta);

    for (int ix = 1; ix <= nX; ix++)
        for (int iy = 1; iy <= h->GetNbinsY(); iy++) {
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

TH2F* RebinTH2X_varBins(TH2F* h, int nEta, double* eEta) {
    int nY = h->GetNbinsY();
    const TArrayD* yArr = h->GetYaxis()->GetXbins();
    TH2F* hNew;
    if (yArr->GetSize() > 0)
        hNew = new TH2F(h->GetName(), h->GetTitle(),
                        nEta, eEta,
                        nY, yArr->GetArray());
    else
        hNew = new TH2F(h->GetName(), h->GetTitle(),
                        nEta, eEta,
                        nY, h->GetYaxis()->GetXmin(), h->GetYaxis()->GetXmax());

    for (int ix = 1; ix <= h->GetNbinsX(); ix++)
        for (int iy = 1; iy <= nY; iy++) {
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

TH2F* FoldAbsTH2X(TH2F* h, const std::string& newName) {
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

    TH2F* hf = new TH2F(newName.c_str(), h->GetTitle(),
                        nxPos, xedges.data(), ny, yedges.data());
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

TH2F* FoldAbsTH2Y(TH2F* h, const std::string& newName) {
    int nx = h->GetNbinsX();
    int ny = h->GetNbinsY();
    const TAxis* ay = h->GetYaxis();

    int jzero = ay->FindBin(0.0 + 1e-9);
    int nyPos = ny - jzero + 1;

    std::vector<double> xedges;
    for (int i = 1; i <= nx + 1; ++i) xedges.push_back(h->GetXaxis()->GetBinLowEdge(i));
    std::vector<double> yedges;
    for (int j = jzero; j <= ny + 1; ++j) yedges.push_back(ay->GetBinLowEdge(j));

    TH2F* hf = new TH2F(newName.c_str(), h->GetTitle(),
                        nx, xedges.data(), nyPos, yedges.data());
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


void loadHistograms(Region& r, TFile* f, const std::string& regionName, bool bool_rebin=true, int rebineta=1, int rebinp=1, int rebinih=1, bool TakeAbsEta = false) {

    cout << "loading region " << regionName << endl;

    r.eta_p                             = (TH2F*) f->Get(("eta_1oP_"+regionName).c_str())->Clone();

    r.ih_eta                            = (TH2F*) f->Get(("ih_eta_"+regionName).c_str())->Clone();
    r.ih_p                              = (TH2F*) f->Get(("ih_p_"+regionName).c_str())->Clone();
    r.ih_p_cross1D                      = (TH2F*) r.ih_p->Clone(); r.ih_p_cross1D->Reset(); r.ih_p_cross1D->SetName(("cross1D_"+regionName).c_str());
    r.ih_p_cross1D_fit                  = (TH2F*) r.ih_p->Clone(); r.ih_p_cross1D_fit->Reset(); r.ih_p_cross1D_fit->SetName(("cross1D_fit_"+regionName).c_str());
    r.ih_p_cross1D_corr                 = (TH2F*) r.ih_p->Clone(); r.ih_p_cross1D_corr->Reset(); r.ih_p_cross1D_corr->SetName(("cross1D_corr_"+regionName).c_str());

    r.ias_p                             = (TH2F*) f->Get(("ias_p_"+regionName).c_str())->Clone();
    r.ias_pt                            = (TH2F*) f->Get(("ias_pt_"+regionName).c_str())->Clone();

    r.mass                              = (TH1F*) f->Get(("mass_"+regionName).c_str())->Clone();
    r.mass_eta                          = (TH2F*) f->Get(("mass_eta_"+regionName).c_str())->Clone();

    r.pred_mass                         = (TH1F*) r.mass->Clone(); r.pred_mass->SetName(("pred_mass_"+regionName).c_str()); r.pred_mass->Reset();
    r.pred_mass_eta                     = (TH2F*) r.mass_eta->Clone(); r.pred_mass_eta->SetName(("pred_mass_eta_"+regionName).c_str()); r.pred_mass_eta->Reset();
    r.pred_mass_fitIh                   = (TH1F*) r.pred_mass->Clone(); r.pred_mass_fitIh->SetName(("pred_mass_fitIh_"+regionName).c_str());
    r.pred_mass_fitP                    = (TH1F*) r.pred_mass->Clone(); r.pred_mass_fitP->SetName(("pred_mass_fitP_"+regionName).c_str());
    r.pred_mass_fitIh_fitP              = (TH1F*) r.pred_mass->Clone(); r.pred_mass_fitIh_fitP->SetName(("pred_mass_fitIh_fitP_"+regionName).c_str());
    r.pred_mass_noFit                   = (TH1F*) r.pred_mass->Clone(); r.pred_mass_noFit->SetName(("pred_mass_noFit_"+regionName).c_str());


    std::vector<double> RebinEta_1_2p4_Down = {-2.4, -2.0, -1.6, -1., 1., 1.6, 2.0, 2.4};
    std::vector<double> RebinEta_1_Down = {-2.4, -1., -0.6, -0.2, 0.2, 0.6, 1., 2.4};

    if (bool_rebin) {

        r.ih_p->Rebin2D(rebinp,rebinih);
        r.ias_p->Rebin2D(rebinp,rebinih);
        r.ias_pt->Rebin2D(rebinp,rebinih);

        if (rebineta==8) {
            if (regionName.find("Eta2p4") != std::string::npos) {
                r.eta_p->Rebin2D(rebinp, rebineta);
                r.ih_eta->Rebin2D(rebineta, rebinih);
                r.mass_eta->Rebin2D(1, rebineta);
                r.pred_mass_eta->Rebin2D(rebineta, 1);
            }
            else {
                std::vector<double>& RebinEtaVec = (regionName.find("Eta1_2p4") != std::string::npos)
                                                ? RebinEta_1_2p4_Down
                                                : RebinEta_1_Down;
                int nEta = RebinEtaVec.size() - 1;
                double* eEta  = RebinEtaVec.data();

                // eta_p : X=p (uniforme), Y=eta (variable)
                r.eta_p->RebinX(rebinp);
                TH2F* tmp = RebinTH2Y_varBins(r.eta_p, nEta, eEta);
                delete r.eta_p; r.eta_p = tmp;

                // ih_eta : X=eta (variable), Y=ih (uniforme)
                r.ih_eta->RebinY(rebinih);
                tmp = RebinTH2X_varBins(r.ih_eta, nEta, eEta);
                delete r.ih_eta; r.ih_eta = tmp;

                // mass_eta : X=mass (uniforme), Y=eta (variable)
                r.mass_eta->RebinX(1);
                tmp = RebinTH2Y_varBins(r.mass_eta, nEta, eEta);
                delete r.mass_eta; r.mass_eta = tmp;

                // pred_mass_eta : X=eta (variable), Y=mass (uniforme)
                r.pred_mass_eta->RebinY(1);
                tmp = RebinTH2X_varBins(r.pred_mass_eta, nEta, eEta);
                delete r.pred_mass_eta; r.pred_mass_eta = tmp;
            }
        }
        else {
            r.eta_p->Rebin2D(rebinp,rebineta);
            r.ih_eta->Rebin2D(rebineta,rebinih);
            r.mass_eta->Rebin2D(1,rebineta);
            r.pred_mass_eta->Rebin2D(rebineta,1);
        }
    }

    // ---- Repliement |eta| ----
    // eta_p        : eta sur l'axe Y  -> FoldAbsTH2Y
    // ih_eta       : eta sur l'axe X  -> FoldAbsTH2X
    // mass_eta     : eta sur l'axe Y  -> FoldAbsTH2Y
    // pred_mass_eta: eta sur l'axe X  -> FoldAbsTH2X
    if (TakeAbsEta) {
        TH2F* tmp;

        tmp = FoldAbsTH2Y(r.eta_p, ("eta_1oP_"+regionName+"_absEta").c_str());
        delete r.eta_p; r.eta_p = tmp;

        tmp = FoldAbsTH2X(r.ih_eta, ("ih_eta_"+regionName+"_absEta").c_str());
        delete r.ih_eta; r.ih_eta = tmp;

        tmp = FoldAbsTH2Y(r.mass_eta, ("mass_eta_"+regionName+"_absEta").c_str());
        delete r.mass_eta; r.mass_eta = tmp;

        tmp = FoldAbsTH2X(r.pred_mass_eta, ("pred_mass_eta_"+regionName+"_absEta").c_str());
        delete r.pred_mass_eta; r.pred_mass_eta = tmp;
    }

    return;
}


TH1F* poissonHisto(const TH1F& h,TRandom3* RNG) {
    TH1F* hres = (TH1F*) h.Clone();
    for(int i=0;i<h.GetNbinsX()+1;i++){
        hres->SetBinContent(i,RNG->Poisson(h.GetBinContent(i)));
    }
    return hres;
}


TH2F* poissonHisto(const TH2F& h,TRandom3* RNG) {
    TH2F* hres = (TH2F*) h.Clone();
    for(int i=0;i<h.GetNbinsX()+1;i++){
        for(int j=0;j<h.GetNbinsY()+1;j++){
            hres->SetBinContent(i,j,RNG->Poisson(h.GetBinContent(i,j)));
        }
    }
    return hres;
}


// Function doing the eta reweighing between two 2D-histograms as done in the Hscp background estimate method,
// because of the correlation between variables (momentum & transverse momentum). 
// The first given 2D-histogram is weighted in respect to the 1D-histogram 
void etaReweighingP(TH2F* eta_p_1, const TH1F* eta2_) {
    TH1F* eta1 = (TH1F*) eta_p_1->ProjectionY(); 
    TH1F* eta2 = (TH1F*) eta2_->Clone();
    eta1->Scale(1./eta1->Integral(0,eta1->GetNbinsX()+1));
    eta2->Scale(1./eta2->Integral(0,eta2->GetNbinsX()+1));
    eta2->Divide(eta1);
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
void etaReweighingP(TH2F* ih_eta_C, const TH1F* eta_B_, const TH1F* eta_A_) {
    TH1F* eta_B = (TH1F*) eta_B_->Clone();
    TH1F* eta_A = (TH1F*) eta_A_->Clone();
    eta_B->Scale(1./eta_B->Integral(0,eta_B->GetNbinsX()+1));
    eta_A->Scale(1./eta_A->Integral(0,eta_A->GetNbinsX()+1));
    eta_B->Divide(eta_A);

    for(int i=0; i<ih_eta_C->GetNbinsY()+2; i++)  // ih bins
    {
        for(int j=0; j<ih_eta_C->GetNbinsX()+2; j++)  // eta bins
        {
            ih_eta_C->SetBinContent(j, i, ih_eta_C->GetBinContent(j,i)*eta_B->GetBinContent(j));
            ih_eta_C->SetBinError(j, i, ih_eta_C->GetBinError(j,i)*eta_B->GetBinContent(j));
        }
    }
}

// add the overflow bin to the last one
void overflowLastBin(TH1F* h) {
    h->SetBinContent(h->GetNbinsX(),h->GetBinContent(h->GetNbinsX())+h->GetBinContent(h->GetNbinsX()+1));
    h->SetBinContent(h->GetNbinsX()+1,0);
}


// rebinning histogram according to an array of bins
TH1F* rebinHisto(TH1F* h) {
    double xbins[33]={0.,20.,40.,60.,80.,100.,120.,140.,160.,180.,200.,220.,240.,260.,280.,300.,320.,340.,360.,380.,410.,440.,480.,530.,590.,660.,760.,880.,1030.,1210.,1440.,1730.,2000.};
    std::vector<double> xbins_v;
    for(double i=0.0;i<=1000.0;i+=50) xbins_v.push_back(i);
    std::string newname = h->GetName(); 
    newname += "_rebinned";
    TH1F* hres = (TH1F*) h->Rebin(32,newname.c_str(),xbins);
    return hres;
}


// Function returning the ratio of right integer (from x to infty) for two 1D-histograms
// This function is used in the Hscp data-driven background estimate to test the mass shape prediction
// The argument to use this type of ratio is that we're in case of cut & count experiment 
TH1F* ratioIntegral(TH1F* h1, TH1F* h2) {    
    float SystError = systErr_;
    TH1F* res = (TH1F*) h1->Clone(); res->Reset();
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


TH1F* pull(TH1F* h1, TH1F* h2) {
    float SystError = systErr_;
    TH1F* res = (TH1F*) h2->Clone(); //res->Reset();
    res->Divide(h1);

    return res;
}


void saveHistoRatio(TH1F* h1,TH1F* h2,std::string st1,std::string st2,std::string st3,bool rebin=false) {
    h1->SetName(st1.c_str());
    h2->SetName(st2.c_str());
    if(rebin){
        h1=rebinHisto(h1);
        h2=rebinHisto(h2);
    }
    h1->Write();
    h2->Write();
    TH1F* R = (TH1F*) ratioIntegral(h2,h1)->Clone();
    if(rebin) st3+="_rebinned";
    R->SetName(st3.c_str());
    R->Write();
}


TH1F meanHistoPE(std::vector<TH1F> vPE) {
    TH1F h = TH1F(vPE[0]);
    h.Reset();
    h.SetBinErrorOption(TH1::EBinErrorOpt::kPoisson);
    h.Sumw2();

    for(int i=0;i<h.GetNbinsX()+1;i++)
    {
        float mean = 0, err = 0;

        for(unsigned int pe=0; pe<vPE.size(); pe++) mean += vPE[pe].GetBinContent(i);
        mean /= vPE.size();

        for(unsigned int pe=0; pe<vPE.size(); pe++) err += pow(mean - vPE[pe].GetBinContent(i),2);

        if(vPE.size()>1) err = sqrt(err/(vPE.size()-1));
        else err = sqrt(err);

        h.SetBinContent(i, mean);
        h.SetBinError(i, err);
    }

    return h;
}


TH2F meanHistoPE_2D(std::vector<TH2F> vPE) {
    float SystError = systErr_;
    TH2F h(vPE[0]);  // Copier le premier histogramme
    h.Reset();
    h.SetBinErrorOption(TH1::EBinErrorOpt::kPoisson);

    // Parcours des bins 2D
    for (int i = 0; i <= h.GetNbinsX() + 1; i++) {  // Inclut les underflow et overflow
        for (int j = 0; j <= h.GetNbinsY() + 1; j++)
        {
            float mean = 0, err = 0;

            for (unsigned int pe = 0; pe < vPE.size(); pe++) mean += vPE[pe].GetBinContent(i, j);
            mean /= vPE.size();

            for (unsigned int pe = 0; pe < vPE.size(); pe++) err += pow(mean - vPE[pe].GetBinContent(i, j), 2);

            float fact = 1;
            if (vPE.size() > 1)
            {
                err = sqrt(err / (vPE.size() - 1));
                fact = vPE.size() / (vPE.size() - 1);
            }
            else err = sqrt(err);

            h.SetBinContent(i, j, mean);
            h.SetBinError(i, j, err);
        }
    }

    return h;
}

void bckgEstimate(const std::string& filename,
                  const std::string& st_sample, 
                  const std::string& dirname,
                  const Region& B,
                  const Region& C,
                  const Region& BC,
                  const Region& A,
                  const Region& D,
                  bool ifIhpSAME,
                  const Region& B_ifIhpSAME,
                  const std::string& st,
                  const int& nPE = 200,
                  const bool useFit = true,
                  const bool useOldIhFit = false,
                  const bool useOld1oPFit = false,
                  const bool corrTemplateIh = false,
                  const std::string& etaName = "",
                  const bool saveFits = false,
                  const int& fitIh = 1,
                  const int& fitP = 1,
                  bool blind = false) {

    std::vector<TH1F> vPE_;
    std::vector<TH1F> vPE_corr;
    std::vector<TH2F> vPE_cross1D;
    std::vector<TH2F> vPE_cross1D_corr;  
    ROOT::EnableImplicitMT(true);

    
    Region a = A;
    Region b = B;
    Region c = C;
    Region bc = BC;
    Region d = D;

    TH2F a_ih_eta_base(*a.ih_eta);
    TH2F b_ih_eta_base(*b.ih_eta);
    TH2F c_ih_eta_base(*c.ih_eta);
    TH2F b_eta_p_base(*b.eta_p);
    TH2F c_eta_p_base(*c.eta_p);

    // Pre-fit of 1/p to get the parameters for the next fits in the toys
    TFitResultPtr ptr_pinc = 0;
    TH1F* p_base = (TH1F*)c_eta_p_base.ProjectionX();

    float rangemax_p = 30;
    if (p_base->GetBinCenter(p_base->GetMaximumBin()) < rangemax_p) rangemax_p = 0.8 * p_base->GetBinCenter(p_base->GetMaximumBin());
    
    TF1 f_p_base("f_p_base","[0]*([1]+erf((log(x)-[2])/[3]))",0,rangemax_p);
    f_p_base.SetParameter(0,560);
    f_p_base.FixParameter(1,1.0);
    f_p_base.SetParameter(2,3.50116e+00);
    f_p_base.SetParameter(3,0.60152e+00);

    ptr_pinc = p_base->Fit(&f_p_base,"QRSL","",0,rangemax_p);

    cout << "Did the initial 1/p fit" << endl;
    cout << endl;
    double par_p2 = f_p_base.GetParameter(2);
    double par_p3 = f_p_base.GetParameter(3);


    // Toys lambda function
    auto workItem = [] (UInt_t workerID, const std::string& filename, const std::string& st_sample, const std::string& dirname, const Region& B, const Region& C, 
                        const Region& BC, const Region& A, const Region& D, bool ifIhpSAME, const Region& B_ifIhpSAME, const std::string& st,
                        const int& nPE = 200, const bool useFit = true, const bool useOldIhFit = false, const bool useOld1oPFit = false, const bool corrTemplateIh=false,
                        const std::string& etaName = "", const bool& saveFits = false, const int& fitIh = 1, const int& fitP = 1, bool blind = false, 
                        const double& par_p2 = 4.70839, const double& par_p3 = 1.05005)
                        -> std::tuple<TH1F, TH2F, float>
    {

        //cout<<"workerId: "<<workerID<<endl;
    
        
        // Setup
        Region a = A;
        Region b = B;
        Region b_ifIhpSAME = B_ifIhpSAME;
        Region c = C;
        Region bc = BC;

        TH2F a_ih_eta_base(*a.ih_eta);
        TH1F* a_eta_base = (TH1F*)a_ih_eta_base.ProjectionX();
        TH2F b_ih_eta_base(*b.ih_eta);
        TH2F b_ifIhpSAME_ih_eta_base(*b_ifIhpSAME.ih_eta);
        TH2F b_ifIhpSAME_eta_p_base(*b_ifIhpSAME.eta_p);
        TH1F* b_ifIhpSAME_eta_base = (TH1F*)b_ifIhpSAME_eta_p_base.ProjectionY();
        TH2F c_ih_eta_base(*c.ih_eta);
        TH2F b_eta_p_base(*b.eta_p);
        TH1F* b_eta_base = (TH1F*)b_eta_p_base.ProjectionY();
        TH2F c_eta_p_base(*c.eta_p);

        if(corrTemplateIh) corrIh(&b_ih_eta_base);
        
        
        if (st_sample=="data2017"){bc.K_=K_data2017;bc.C_=C_data2017;}
        else if (st_sample=="data2018"){bc.K_=K_data2018;bc.C_=C_data2018;}
        else if (st_sample=="data2024"){bc.K_=K_data2024;bc.C_=C_data2024;}
        else if (st_sample=="mc2017"){bc.K_=K_mc2017;bc.C_=C_mc2017;}
        else if (st_sample=="mc2018"){bc.K_=K_mc2018;bc.C_=C_mc2018;}

        
        // 1/p fit
        TH1F* p_base = (TH1F*)c_eta_p_base.ProjectionX();
        float rangemax_p = 30;
        if (p_base->GetBinCenter(p_base->GetMaximumBin()) < rangemax_p) rangemax_p = 0.8 * p_base->GetBinCenter(p_base->GetMaximumBin());
        
        TF1 f_p("f_p", useOld1oPFit ? "[0]*([1]+erf((log(x)-[2])/[3]))" : "0.5*(exp([0]*x*x+[1]*x)+exp(-[0]*x*x-[1]*x))-1", 0, rangemax_p);
        if (useOld1oPFit) {
            f_p.SetParLimits(0, 0, 9000);
            f_p.FixParameter(1, 1.0);
            f_p.FixParameter(2, par_p2);
            f_p.FixParameter(3, par_p3);
        }
        else {
            f_p.SetParLimits(0, 0, 1e-3);
            f_p.SetParLimits(1, -1e-2, 1e-2);
        }

        // Ih fit
        TH1F* ih_base = (TH1F*)b_ih_eta_base.ProjectionX();
        float max_ih = ih_base->GetBinCenter(ih_base->GetMaximumBin());

        TF1 f_ihg("f_ihg", "gaus", max_ih, 8);
        f_ihg.SetParameter(0, 0.5*ih_base->Integral());
        f_ihg.SetParameter(1, max_ih);
        f_ihg.SetParameter(2, ih_base->GetStdDev());


        // Toys
        TRandom3* RNG = new TRandom3(workerID);
        bc.pred_mass->Reset();
        bc.pred_mass_eta->Reset();

        TH2F* a_ih_eta = poissonHisto(a_ih_eta_base,RNG);
        TH2F* b_ih_eta = poissonHisto(b_ih_eta_base,RNG);
        TH2F* b_ifIhpSAME_ih_eta = poissonHisto(b_ifIhpSAME_ih_eta_base,RNG);
        TH2F* b_eta_p = poissonHisto(b_eta_p_base,RNG);
        TH2F* c_eta_p = poissonHisto(c_eta_p_base,RNG);
        
        TH1F* b_ih = (TH1F*)b_ih_eta->ProjectionY();
        TH1F* b_eta = poissonHisto(*b_eta_base,RNG);
        TH1F* b_ifIhpSAME_eta = poissonHisto(*b_ifIhpSAME_eta_base,RNG);
        TH1F* a_eta = poissonHisto(*a_eta_base,RNG);
        
        bool bloutaba = false;
        if(ifIhpSAME) etaReweighingP(b_ih_eta, b_ifIhpSAME_eta, a_eta); //bloutaba = true; in C
        else etaReweighingP(c_eta_p,b_eta);


        // Mass prediction in the BC region
        bc.eta_p = c_eta_p;
        bc.ih_eta = b_ih_eta;
        if(st.find("ias") != std::string::npos) bc.fillPredMass(filename, st, st_sample, f_p, f_ihg, useFit, fitIh, fitP, -1, useOldIhFit, useOld1oPFit, etaName, saveFits);
        else bc.fillPredMass(filename, st, st_sample, f_p, f_ihg, useFit, fitIh, fitP, -1, useOldIhFit, useOld1oPFit, etaName, saveFits);

        float normA = a_ih_eta->Integral(0, a_ih_eta->GetNbinsX()+1, 0, a_ih_eta->GetNbinsY()+1);
        float normB = b_ih_eta->Integral(0, b_ih_eta->GetNbinsX()+1, 0, b_ih_eta->GetNbinsY()+1);
        float normC = c_eta_p->Integral(0, c_eta_p->GetNbinsX()+1, 0, c_eta_p->GetNbinsY()+1);

        float normalisationABC = normB*normC/normA;
        if (ifIhpSAME) normalisationABC = b_ifIhpSAME_ih_eta->Integral(0,b_ifIhpSAME_ih_eta->GetNbinsX()+1,0,b_ifIhpSAME_ih_eta->GetNbinsY()+1) * normC / normA;
        
        bc.pred_mass->Scale(normalisationABC/bc.pred_mass->Integral());
        bc.pred_mass_eta->Scale(normalisationABC/bc.pred_mass_eta->Integral());


        // End
        delete a_ih_eta;
        delete b_ih_eta;
        delete b_eta_p;
        delete c_eta_p;
        delete b_eta;
        delete b_ifIhpSAME_ih_eta;

        return {*bc.pred_mass, *bc.pred_mass_eta, normalisationABC};
    };

    
    // Loop on the toys
    ROOT::TProcessExecutor workers(26);
    auto workItemToRun = std::bind (workItem, _1, filename, st_sample, dirname, B, C, BC, A, D, ifIhpSAME, B_ifIhpSAME, st, nPE,
                                    useFit, useOldIhFit, useOldIhFit, corrTemplateIh, etaName, saveFits, fitIh, fitP, blind, par_p2, par_p3);
    auto vPE = workers.Map(workItemToRun, ROOT::TSeqI(nPE));


    // Get the results
    std::vector<TH1F> histo_pred_mass;
    std::vector<TH2F> histo_pred_mass_eta;
    std::vector<float> normalisations;

    for (const auto& result : vPE) {
        histo_pred_mass.push_back(std::get<0>(result));   // Récupère *bc.pred_mass
        histo_pred_mass_eta.push_back(std::get<1>(result)); // Récupère *bc.pred_mass_eta
        normalisations.push_back(std::get<2>(result)); // Récupère normalisationABC
    }


    // Mean histogram of the toys
    TH1F h_temp = meanHistoPE(histo_pred_mass);
    TH2F h_temp_eta = meanHistoPE_2D(histo_pred_mass_eta);
    std::cout<<"Mean histo of all PE entries : " << h_temp.GetEntries() << " , and integral : " << h_temp.Integral() << std::endl;
    std::cout<<"Mean histo of bc.pred_mass BEFORE changing it to NPE : " << bc.pred_mass->GetEntries() << " , and integral : " << bc.pred_mass->Integral() << std::endl;
    if(nPE>1){ bc.pred_mass = &h_temp; bc.pred_mass_eta = &h_temp_eta; }
    float avgNormalisation = std::accumulate(normalisations.begin(), normalisations.end(), 0.0f) / normalisations.size();

    if(blind) blindMass(d.mass,300);
     
    std::cout<<"Before overflowLastBin : bc.pred_mass entries : " << bc.pred_mass->GetEntries() << " and d.mass entries : " << d.mass->GetEntries() << std::endl; 
    std::cout<<"Before overflowLastBin : bc.pred_mass integral : " << bc.pred_mass->Integral() << " and d.mass integral : " << d.mass->Integral() << std::endl; 
    overflowLastBin(d.mass);
    overflowLastBin(bc.pred_mass);


    // Saving histograms
    saveHistoRatio(d.mass,bc.pred_mass,("mass_obs_"+st).c_str(),("mass_predBC_"+st).c_str(),("mass_predBCR_"+st).c_str());

    bc.pred_mass_eta->Write();
    bc.pred_mass_eta->ProjectionY()->Write();
    d.mass_eta->Write();
    d.mass_eta->ProjectionY()->Write();
    c.mass_eta->Write();
    c.mass_eta->ProjectionY()->Write();

    A.ih_eta->Write();
    A.eta_p->Write();
    B.ih_eta->Write();
    C.eta_p->Write();
    D.eta_p->Write();
    D.ih_eta->Write();
    
    b_ih_eta_base.Write();
    b_ih_eta_base.ProjectionY("_ih_py")->Write();
    c_eta_p_base.Write();
    c_eta_p_base.ProjectionX("_p_px")->Write();
    
    D.ih_eta->ProjectionY()->Write();
    B.ih_eta->ProjectionY()->Write();
    D.ih_eta->ProjectionX()->Write();
    B.ih_eta->ProjectionX()->Write();
    D.eta_p->ProjectionX()->Write();
    C.eta_p->ProjectionX()->Write();   
    C.eta_p->ProjectionY()->Write();


    C.mass->Scale(avgNormalisation/C.mass->Integral());
    C.mass->Write();

    // Draw on one canvas all the histos of the vector histo_pred_mass
    TCanvas *c1 = new TCanvas("c1","c1",800,800);
    c1->cd();
    for (int i=0; i<histo_pred_mass.size(); i++) histo_pred_mass[i].Rebin(10);
    histo_pred_mass[0].Draw("hist");
    for (int i=1; i<histo_pred_mass.size(); i++) histo_pred_mass[i].Draw("hist same");
    c1->Write();

    // fill a new histo with norm/d.mass->integral:
    TH1F* h_norm = new TH1F("h_norm", "h_norm", 400, 0.9, 1.1);
    h_norm->GetXaxis()->SetTitle("BC/AD");
    h_norm->GetYaxis()->SetTitle("Entries");
    for(int i=0; i<normalisations.size(); i++) h_norm->Fill(normalisations[i]/d.mass->Integral());
    h_norm->Write();

}