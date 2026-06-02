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

#include "CommonFunctions.h"

using namespace std;
gErrorIgnoreLevel = kFatal;



void step2_backgroundPrediction() {
    ifstream infile;
    infile.open("configFile_readHisto_toLaunch.txt");
    std::string line;
    std::string filename;
    std::string st_sample;
    std::string dirname;
    int nPE;
    int rebineta, rebinih, rebinp;
    bool rebin;
    bool corrTemplateIh;
    int fitIh, fitP;
    bool useFit;
    while(std::getline(infile,line)){
        if(std::strncmp(line.c_str(),"#",1)==0) continue;
        std::cout << line << std::endl;
        std::stringstream ss(line);
        ss >> filename >> st_sample >> dirname >> nPE >> rebin >> rebineta >> rebinih >> rebinp >> corrTemplateIh >> fitIh >> fitP >> useFit;
    }

    std::string commandDir = "mkdir -p "+dirname;
    system(commandDir.c_str());

    std::string outfilename_;
    outfilename_ = filename + "_rebinEta" + to_string(rebineta) + "_rebinIh" + to_string(rebinih) + "_rebinP" + to_string(rebinp);

    if(corrTemplateIh) outfilename_ += "_corrTemplateIh";
    if(fitIh==0) outfilename_ += "_fitIhDown";
    if(fitIh==2) outfilename_ += "_fitIhUp";
    if(fitP==0) outfilename_ += "_fitPDown";
    if(fitP==2) outfilename_ += "_fitPUp";
    outfilename_ += "_EtaReweighting";


    //std::string Ext = "_METanalysis_Eta1";
    //std::string Ext = "_METanalysis_Eta1_2p4";
    std::string Ext = "_METanalysis_Eta2p4";

    std::string etaName = "";
    if (Ext == "_METanalysis_Eta1") etaName = "_Eta1";
    else if (Ext == "_METanalysis_Eta1_2p4") etaName += "_Eta1_2p4";
    else if (Ext == "_METanalysis_Eta2p4") etaName += "_Eta2p4";
    outfilename_ += etaName;


    bool bool_rebin = rebin;
    bool blind = false;
    bool useOldIhFit = false, useOld1oPFit = false;
    bool saveFits = true;
    std::string DataSetName = filename.substr(filename.find_last_of('/') + 1);
    bool TakeAbsEta = true;
    cout << "valeur absolue de eta: " << TakeAbsEta << endl;
    
    //useFit = false; // to be commented !
    cout << "Use fit: " << useFit << endl;
    cout << "useOldIhFit: " << useOldIhFit << ", useOld1oPFit: " << useOld1oPFit << endl;
    cout << "saveFits: " << saveFits << endl;

    if ((useOldIhFit || useOld1oPFit) && useFit) outfilename_ += "_OldFit";
    else if ((!useOldIhFit && !useOld1oPFit) && useFit) outfilename_ += "_NewFit";
    else if (!useFit) outfilename_ += "_NoFit";


    std::cout << "Input file:      " << DataSetName << std::endl;
    std::cout << "Output file:     " << outfilename_ << std::endl;
    TFile* ifile = new TFile((filename+".root").c_str());
    TFile* ofile = new TFile((outfilename_+".root").c_str(),"RECREATE");

    
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

    Region rb_999ias100;
    Region rd_999ias100;
    Region rbc_999ias100;


    // loading histograms used to validate the background estimate method in data --> base on Ias slices
    std::cout << std::endl;
    std::cout << "    Loading..." << std::endl;

    loadHistograms(ra_ias50, ifile, "regionA_ias50" + Ext ,bool_rebin ,rebineta ,rebinp ,rebinih); 
    
    loadHistograms(rb_50ias60, ifile, "regionB_50ias60"+E xt,bool_reb in,rebine ta,rebi np,rebin i);
    loadHistograms(rb_50ias90, ifile, "regionB_50ias90"+E xt,bool_reb in,rebine ta,rebi np,rebin i); 
    loadHistograms(rb_60ias70, ifile, "regionB_60ias70"+E xt,bool_reb in,rebine ta,rebi np,rebin i); 
    loadHistograms(rb_70ias80, ifile, "regionB_70ias80"+E xt,bool_reb in,rebine ta,rebi np,rebin i); 
    loadHistograms(rb_80ias90, ifile, "regionB_80ias90"+E xt,bool_reb in,rebine ta,rebi np,rebin i);

    loadHistograms(rc_ias50, ifile, "regionC_ias50" + Ext ,bool_rebin ,rebineta ,rebinp ,rebinih); 

    loadHistograms(rd_50ias60, ifile, "regionD_50ias60"+E xt,bool_reb in,rebine ta,rebi np,rebin i);
    loadHistograms(rd_50ias90, ifile, "regionD_50ias90"+E xt,bool_reb in,rebine ta,rebi np,rebin i); 
    loadHistograms(rd_60ias70, ifile, "regionD_60ias70"+E xt,bool_reb in,rebine ta,rebi np,rebin i); 
    loadHistograms(rd_70ias80, ifile, "regionD_70ias80"+E xt,bool_reb in,rebine ta,rebi np,rebin i); 
    loadHistograms(rd_80ias90, ifile, "regionD_80ias90"+E xt,bool_reb in,rebine ta,rebi np,rebin i);
     
    loadHistograms(rbc_50ias60, ifile, "regionD_50ias60"+E xt,bool_reb in,rebine ta,rebi np,rebin i);
    loadHistograms(rbc_50ias90, ifile, "regionD_50ias90"+E xt,bool_reb in,rebine ta,rebi np,rebin i); 
    loadHistograms(rbc_60ias70, ifile, "regionD_60ias70"+E xt,bool_reb in,rebine ta,rebi np,rebin i); 
    loadHistograms(rbc_70ias80, ifile, "regionD_70ias80"+E xt,bool_reb in,rebine ta,rebi np,rebin i); 
    loadHistograms(rbc_80ias90, ifile, "regionD_80ias90"+E xt,bool_reb in,rebine ta,rebi np,rebin i);


    loadHistograms(rb_999ias100, ifile, "regionB_999ias100"  + Ext,bool_r ebin,rebi neta,re binp,reb ini); 
    loadHistograms(rd_999ias100, ifile, "regionD_999ias100"  + Ext,bool_r ebin,rebi neta,re binp,reb ini); 
    loadHistograms(rbc_999ias100, ifile, "regionD_999ias100"  + Ext,bool_r ebin,rebi neta,re binp,reb ini);
    
    std::cout << "    Regions loaded" << std::endl;
    std::cout << std::endl;
    std::cout << "    Background estimation..." << std::endl;
    std::cout << std::endl;

    // Estimate the background in different Ias slices, each containing 10% of the statistic 
    bool blind = false;

    //bckgEstimate(st_sample, dirname, rb_50ias60, rc_ias50, rbc_50ias60, ra_ias50, rd_50ias60, "50ias60", nPE, corrTemplateIh, corrTemplateIh, fitIh, fitP, blind);
    //bckgEstimate(st_sample, dirname, rb_60ias70, rc_ias50, rbc_60ias70, ra_ias50, rd_60ias70, "60ias70", nPE, corrTemplateIh, corrTemplateIh, fitIh, fitP, blind);
    //bckgEstimate(st_sample, dirname, rb_70ias80, rc_ias50, rbc_70ias80, ra_ias50, rd_70ias80, "70ias80", nPE, corrTemplateIh, corrTemplateIh, fitIh, fitP, blind);
    //bckgEstimate(st_sample, dirname, rb_80ias90, rc_ias50, rbc_80ias90, ra_ias50, rd_80ias90, false, rb_80ias90, "80ias90", nPE, corrTemplateIh, corrTemplateIh, fitIh, fitP, blind);
    //bckgEstimate(st_sample, dirname, rb_50ias90, rc_ias50, rbc_50ias90, ra_ias50, rd_50ias90, "50ias90", nPE, corrTemplateIh, corrTemplateIh, fitIh, fitP, blind);
    bckgEstimate(st_sample, dirname, rb_999ias100, rc_ias50, rbc_999ias100, ra_ias50, rd_999ias100, false, rb_999ias100, "999ias100", nPE, true, corrTemplateIh, corrTemplateIh, fitIh, fitP, blind);
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
    Region ra_3fp9;
    Region rb_8fp9;
    Region rb_9fp10;
    Region rc_3fp8;
    Region rc_3fp9;
    Region rd_8fp9;
    Region rd_9fp10;
    Region rbc_8fp9;
    Region rbc_9fp10;
    Region rbc_3fp8;

    std::cout << std::endl;
    std::cout << "    Loading..." << std::endl;

    loadHistograms(ra_3fp8, ifile, "regionA_3fp8" + Ext, bool_rebin, rebineta, rebinp, rebinih, TakeAbsEta);
    loadHistograms(ra_3fp9, ifile, "regionA_3fp9" + Ext, bool_rebin, rebineta, rebinp, rebinih, TakeAbsEta);
    loadHistograms(rb_8fp9, ifile, "regionB_8fp9" + Ext, bool_rebin, rebineta, rebinp, rebinih, TakeAbsEta);
    loadHistograms(rb_9fp10, ifile, "regionB_9fp10" + Ext, bool_rebin, rebineta, rebinp, rebinih, TakeAbsEta);
    loadHistograms(rc_3fp8, ifile, "regionC_3fp8" + Ext, bool_rebin, rebineta, rebinp, rebinih, TakeAbsEta);
    loadHistograms(rc_3fp9, ifile, "regionC_3fp9" + Ext, bool_rebin, rebineta, rebinp, rebinih, TakeAbsEta);
    loadHistograms(rd_8fp9, ifile, "regionD_8fp9" + Ext, bool_rebin, rebineta, rebinp, rebinih, TakeAbsEta);
    loadHistograms(rd_9fp10, ifile, "regionD_9fp10" + Ext ,bool_rebin, rebineta, rebinp, rebinih, TakeAbsEta);
    loadHistograms(rbc_8fp9, ifile, "regionD_8fp9" + Ext, bool_rebin, rebineta, rebinp, rebinih, TakeAbsEta);
    loadHistograms(rbc_9fp10, ifile, "regionD_9fp10" + Ext, bool_rebin, rebineta, rebinp, rebinih, TakeAbsEta);
    loadHistograms(rbc_3fp8, ifile, "regionC_3fp8" + Ext, bool_rebin, rebineta, rebinp, rebinih, TakeAbsEta);

    std::cout << "    Regions loaded" << std::endl;
    std::cout << std::endl;
    std::cout << "    Background estimation..." << std::endl;
    std::cout << std::endl;

    // Estimate the background in different Fpixel slices
    bool ifIhpSAME = true; // TRUE: Ih and p templates in C region but B still used for the normalisation
    

        // In 8fp9
    blind = false;
    bckgEstimate(DataSetName, st_sample, dirname, rc_3fp8, rc_3fp8, rbc_8fp9, ra_3fp8, rd_8fp9, ifIhpSAME, rb_8fp9,
                  "8fp9", nPE, useFit, useOldIhFit, useOldIhFit, corrTemplateIh, etaName, saveFits, fitIh, fitP, blind);

        // In 9fp10 with A and C in 3fp9
    blind = true;
    // bckgEstimate(DataSetName, st_sample, dirname, rc_3fp9, rc_3fp9, rbc_9fp10, ra_3fp9, rd_9fp10, ifIhpSAME, rb_9fp10,
    //              "9fp10", nPE, useFit, useOldIhFit, useOldIhFit, corrTemplateIh, etaName, saveFits, fitIh, fitP, blind);
    


        // Only in C: no Ih reweighting in Region.h to change
    //bckgEstimate(st_sample, dirname, rc_3fp8, rc_3fp8, rbc_3fp8, rc_3fp8, rc_3fp8, ifIhpSAME, rc_3fp8, "3fp8", nPE, fitIh, fitP, blind);
    

    delete ofile;
    return;
}
