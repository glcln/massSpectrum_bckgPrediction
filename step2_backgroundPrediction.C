    //************************************************************************************
// Usage:
// root -l -q -b step2_backgroundPrediction.C

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "TRandom3.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraphErrors.h>

#include "./Regions.h"

using namespace std;



void step2_backgroundPrediction(){
    ifstream infile;
    infile.open("configFile_readHisto_toLaunch.txt");
    std::string line;
    std::string filename;
    std::string st_sample;
    std::string dirname;
    int nPE, cutIndex;
    int rebineta,rebinih,rebinp,rebinmass;
    bool rebin;
    bool corrTemplateIh, corrTemplateP;
    int fitIh, fitP;
    while(std::getline(infile,line)){
        if(std::strncmp(line.c_str(),"#",1)==0) continue;
        std::cout << line << std::endl;
        std::stringstream ss(line);
        ss >> filename >> st_sample >> dirname >> nPE >> cutIndex >> rebin >> rebineta >> rebinih >> rebinp >> rebinmass >> corrTemplateIh >> corrTemplateP >> fitIh >> fitP;
    }

    std::string commandDir = "mkdir -p "+dirname;
    system(commandDir.c_str());

    std::string outfilename_;
    if(!rebin)outfilename_ = filename+"_cutIndex"+to_string(cutIndex)+"_analysed";
    else outfilename_ = filename +"_cutIndex"+to_string(cutIndex)+"_rebinEta"+to_string(rebineta)+"_rebinIh"+to_string(rebinih)+"_rebinP"+to_string(rebinp)+"_rebinMass"+to_string(rebinmass);

    if(corrTemplateIh) outfilename_ += "_corrTemplateIh";
    if(corrTemplateP) outfilename_ += "_corrTemplateP";
    if(fitIh==0) outfilename_ += "_fitIhDown";
    if(fitIh==2) outfilename_ += "_fitIhUp";
    if(fitP==0) outfilename_ += "_fitPDown";
    if(fitP==2) outfilename_ += "_fitPUp";
    outfilename_ += "_EtaReweighting";

    std::cout << outfilename_ << std::endl;

    bool bool_rebin=rebin;
    
    TFile* ifile = new TFile((filename+".root").c_str());
    TFile* ofile = new TFile((outfilename_+".root").c_str(),"RECREATE");

    std::string Ext = "_ReRunRaph";
    
    // ------------------------------------------------------------------
    //                              If Gstrip
    //
    //                    pT
    //                      |           |      |     |
    //                      |     C     |      |  D  | blind
    //                      |           |      |   VR|
    //                   70 |-----------|------|-----|------
    //                      |           |      |     |
    //                      |     A     |      |  B  |
    //                   55 |___________|______|_____|______
    //                      0         0.018  0.038 0.057     Gstrip
    //                                (50%)  (80%) (90%)
    //
    // ------------------------------------------------------------------
    /*
    Region ra_ias50;

    Region rb_50ias60;
    Region rb_50ias90;
    Region rb_60ias70;
    Region rb_70ias80;
    Region rb_80ias90;

    Region rc_ias50;

    Region rd_50ias60;
    Region rd_50ias90;
    Region rd_60ias70;
    Region rd_70ias80;
    Region rd_80ias90;

    Region rbc_50ias60;
    Region rbc_60ias70;
    Region rbc_70ias80;
    Region rbc_80ias90;
    Region rbc_50ias90;

    // loading histograms used to validate the background estimate method in data --> base on Ias slices
    std::cout << std::endl;
    std::cout << "    Loading..." << std::endl;

    loadHistograms(ra_ias50,ifile,"regionA_ias50"+Ext,bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    
    loadHistograms(rb_50ias60,ifile,"regionB_50ias60"+Ext,bool_rebin,rebineta,rebinp,rebinih,rebinmass);
    loadHistograms(rb_50ias90,ifile,"regionB_50ias90"+Ext,bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rb_60ias70,ifile,"regionB_60ias70"+Ext,bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rb_70ias80,ifile,"regionB_70ias80"+Ext,bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rb_80ias90,ifile,"regionB_80ias90"+Ext,bool_rebin,rebineta,rebinp,rebinih,rebinmass);

    loadHistograms(rc_ias50,ifile,"regionC_ias50"+Ext,bool_rebin,rebineta,rebinp,rebinih,rebinmass); 

    loadHistograms(rd_50ias60,ifile,"regionD_50ias60"+Ext,bool_rebin,rebineta,rebinp,rebinih,rebinmass);
    loadHistograms(rd_50ias90,ifile,"regionD_50ias90"+Ext,bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rd_60ias70,ifile,"regionD_60ias70"+Ext,bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rd_70ias80,ifile,"regionD_70ias80"+Ext,bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rd_80ias90,ifile,"regionD_80ias90"+Ext,bool_rebin,rebineta,rebinp,rebinih,rebinmass);
     
    loadHistograms(rbc_50ias60,ifile,"regionD_50ias60"+Ext,bool_rebin,rebineta,rebinp,rebinih,rebinmass);
    loadHistograms(rbc_50ias90,ifile,"regionD_50ias90"+Ext,bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rbc_60ias70,ifile,"regionD_60ias70"+Ext,bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rbc_70ias80,ifile,"regionD_70ias80"+Ext,bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rbc_80ias90,ifile,"regionD_80ias90"+Ext,bool_rebin,rebineta,rebinp,rebinih,rebinmass);
    
    std::cout << "    Regions loaded" << std::endl;
    std::cout << std::endl;
    std::cout << "    Background estimation..." << std::endl;
    std::cout << std::endl;

    // Estimate the background in different Ias slices, each containing 10% of the statistic 
    bool blind = false;

    //bckgEstimate(st_sample, dirname, rb_50ias60, rc_ias50, rbc_50ias60, ra_ias50, rd_50ias60, "50ias60", nPE, corrTemplateIh, corrTemplateP, fitIh, fitP, blind);
    //bckgEstimate(st_sample, dirname, rb_60ias70, rc_ias50, rbc_60ias70, ra_ias50, rd_60ias70, "60ias70", nPE, corrTemplateIh, corrTemplateP, fitIh, fitP, blind);
    //bckgEstimate(st_sample, dirname, rb_70ias80, rc_ias50, rbc_70ias80, ra_ias50, rd_70ias80, "70ias80", nPE, corrTemplateIh, corrTemplateP, fitIh, fitP, blind);
    bckgEstimate(st_sample, dirname, rb_80ias90, rc_ias50, rbc_80ias90, ra_ias50, rd_80ias90, false, rb_80ias90, "80ias90", nPE, corrTemplateIh, corrTemplateP, fitIh, fitP, blind);
    //bckgEstimate(st_sample, dirname, rb_50ias90, rc_ias50, rbc_50ias90, ra_ias50, rd_50ias90, "50ias90", nPE, corrTemplateIh, corrTemplateP, fitIh, fitP, blind);
    */
    
    // ------------------------------------------------------------------
    //                              If Fpixel
    //
    //                    pT
    //                      |           |          |
    //                      |     C     |     D    | blind   
    //                      |           |        VR|
    //                   70 |-----------|----------|------
    //                      |           |          |
    //                      |     A     |     B    |
    //                   55 |___________|__________|______
    //                     0.3         0.8        0.9      Fpixel
    //
    // ------------------------------------------------------------------
    
    Region ra_3fp8;
    Region rb_8fp9;
    Region rc_3fp8;
    Region rd_8fp9;
    Region rbc_8fp9;
    Region rbc_3fp8;

    std::cout << std::endl;
    std::cout << "    Loading..." << std::endl;

    loadHistograms(ra_3fp8,ifile,"regionA_3fp8"+Ext,bool_rebin,rebineta,rebinp,rebinih,rebinmass);
    loadHistograms(rb_8fp9,ifile,"regionB_8fp9"+Ext,bool_rebin,rebineta,rebinp,rebinih,rebinmass);
    loadHistograms(rc_3fp8,ifile,"regionC_3fp8"+Ext,bool_rebin,rebineta,rebinp,rebinih,rebinmass);
    loadHistograms(rd_8fp9,ifile,"regionD_8fp9"+Ext,bool_rebin,rebineta,rebinp,rebinih,rebinmass);
    loadHistograms(rbc_8fp9,ifile,"regionD_8fp9"+Ext,bool_rebin,rebineta,rebinp,rebinih,rebinmass);
    loadHistograms(rbc_3fp8,ifile,"regionC_3fp8"+Ext,bool_rebin,rebineta,rebinp,rebinih,rebinmass);

    std::cout << "    Regions loaded" << std::endl;
    std::cout << std::endl;
    std::cout << "    Background estimation..." << std::endl;
    std::cout << std::endl;

    // Estimate the background in different Fpixel slices
    bool blind = false;
    bool ifIhpSAME = true; // TRUE: Ih and p templates in C region but B still used for the normalisation
    

    bckgEstimate(st_sample, dirname, rc_3fp8, rc_3fp8, rbc_8fp9, ra_3fp8, rd_8fp9, ifIhpSAME, rb_8fp9, "8fp9", nPE, corrTemplateIh, corrTemplateP, fitIh, fitP, blind);
    
        // Only in C: no Ih reweighting in Region.h to change
    //bckgEstimate(st_sample, dirname, rc_3fp8, rc_3fp8, rbc_3fp8, rc_3fp8, rc_3fp8, ifIhpSAME, rc_3fp8, "3fp8", nPE, corrTemplateIh, corrTemplateP, fitIh, fitP, blind);
    
    
    delete ofile;

    return;
}
