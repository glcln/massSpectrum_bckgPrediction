#include "CommonFunctions.h"

using namespace std;
gErrorIgnoreLevel = kFatal;



void step2_backgroundPrediction() {
    ifstream infile;
    infile.open("configFile_readHisto_toLaunch.txt");
    std::string line;
    std::string filename;
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
        ss >> filename >> nPE >> rebin >> rebineta >> rebinih >> rebinp >> corrTemplateIh >> fitIh >> fitP >> useFit;
    }

    std::string outfilename_;
    outfilename_ = filename + "_rebinEta" + to_string(rebineta) + "_rebinIh" + to_string(rebinih) + "_rebinP" + to_string(rebinp);

    if(corrTemplateIh) outfilename_ += "_corrTemplateIh";
    if(fitIh==0) outfilename_ += "_fitIhDown";
    if(fitIh==2) outfilename_ += "_fitIhUp";
    if(fitP==0) outfilename_ += "_fitPDown";
    if(fitP==2) outfilename_ += "_fitPUp";
    outfilename_ += "_EtaReweighting";


    std::string st_sample = "data2024";

    std::string Ext = "_METanalysis_TestPUppiMETCut_Eta1";
    //std::string Ext = "_METanalysis_TestPUppiMETCut_Eta1_2p4";
    //std::string Ext = "_METanalysis_TestPUppiMETCut_Eta2p4";
    //std::string Ext = "_METanalysis_TestPUppiMETCut_Eta1p2_2p2";
    //std::string Ext = "_METanalysis_TestPUppiMETCut_Eta1p2_2p4";
    
    //std::string Ext = "_METanalysis_TestPUppiMETCut_Eta1_2p4_IhC";
    // std::string Ext = "_METanalysis_TestPUppiMETCut_Eta1_2p4_Ih3p1";
    //std::string Ext = "_METanalysis_TestPUppiMETCut_Eta1_2p4_Ih3p2";
    //std::string Ext = "_METanalysis_TestPUppiMETCut_Eta1_2p4_Ih3p3";
    // std::string Ext = "_METanalysis_TestPUppiMETCut_Eta1_2p4_Ih3p4";
    // std::string Ext = "_METanalysis_TestPUppiMETCut_Eta1_2p4_Ih3p5";
    //std::string Ext = "_METanalysis_TestPUppiMETCut_Eta1_2p4_Ih3p6";
    // std::string Ext = "_METanalysis_TestPUppiMETCut_Eta1_2p4_Ih3p7";
    // std::string Ext = "_METanalysis_TestPUppiMETCut_Eta1_2p4_Ih3p8";
    // std::string Ext = "_METanalysis_TestPUppiMETCut_Eta1_2p4_Ih3p9";
    //std::string Ext = "_METanalysis_TestPUppiMETCut_Eta1_2p4_Ih4p0";
    
    std::string etaName = "";
    if (Ext == "_METanalysis_TestPUppiMETCut_Eta1") etaName = "_Eta1";
    else if (Ext == "_METanalysis_TestPUppiMETCut_Eta1_2p4") etaName += "_Eta1_2p4";
    else if (Ext == "_METanalysis_TestPUppiMETCut_Eta2p4") etaName += "_Eta2p4";
    else if (Ext == "_METanalysis_TestPUppiMETCut_Eta1p2_2p2") etaName += "_Eta1p2_2p2";
    else if (Ext == "_METanalysis_TestPUppiMETCut_Eta1p2_2p4") etaName += "_Eta1p2_2p4";

    else if (Ext == "_METanalysis_TestPUppiMETCut_Eta1_2p4_IhC") etaName += "_Eta1p2_2p4_IhC";
    else if (Ext == "_METanalysis_TestPUppiMETCut_Eta1_2p4_Ih3p1") etaName += "_Eta1p2_2p4_Ih3p1";
    else if (Ext == "_METanalysis_TestPUppiMETCut_Eta1_2p4_Ih3p2") etaName += "_Eta1p2_2p4_Ih3p2";
    else if (Ext == "_METanalysis_TestPUppiMETCut_Eta1_2p4_Ih3p3") etaName += "_Eta1p2_2p4_Ih3p3";
    else if (Ext == "_METanalysis_TestPUppiMETCut_Eta1_2p4_Ih3p4") etaName += "_Eta1p2_2p4_Ih3p4";
    else if (Ext == "_METanalysis_TestPUppiMETCut_Eta1_2p4_Ih3p5") etaName += "_Eta1p2_2p4_Ih3p5";
    else if (Ext == "_METanalysis_TestPUppiMETCut_Eta1_2p4_Ih3p6") etaName += "_Eta1p2_2p4_Ih3p6";
    else if (Ext == "_METanalysis_TestPUppiMETCut_Eta1_2p4_Ih3p7") etaName += "_Eta1p2_2p4_Ih3p7";
    else if (Ext == "_METanalysis_TestPUppiMETCut_Eta1_2p4_Ih3p8") etaName += "_Eta1p2_2p4_Ih3p8";
    else if (Ext == "_METanalysis_TestPUppiMETCut_Eta1_2p4_Ih3p9") etaName += "_Eta1p2_2p4_Ih3p9";
    else if (Ext == "_METanalysis_TestPUppiMETCut_Eta1_2p4_Ih4p0") etaName += "_Eta1p2_2p4_Ih4p0";
    outfilename_ += etaName;


    bool bool_rebin = rebin;
    bool blind = false;
    bool useOldIhFit = false, useOld1oPFit = true;
    bool saveFits = true;
    std::string DataSetName = filename.substr(filename.find_last_of('/') + 1);
    bool TakeAbsEta = false;
    cout << "valeur absolue de eta: " << TakeAbsEta << endl;
    
    cout << "Use fit: " << useFit << endl;
    cout << "useOldIhFit: " << useOldIhFit << ", useOld1oPFit: " << useOld1oPFit << endl;
    cout << "saveFits: " << saveFits << endl;

    if ((useOldIhFit || useOld1oPFit) && useFit) outfilename_ += "_OldFit";
    else if ((!useOldIhFit && !useOld1oPFit) && useFit) outfilename_ += "_NewFit";
    else if (!useFit) outfilename_ += "_NoFit";

    // ---
    std::vector <float> MyIhCut_values = {2.9784, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6};
    std::vector <std::string> MyIhCut_names  = {"_IhC", "_Ih3p1", "_Ih3p2", "_Ih3p3", "_Ih3p4", "_Ih3p5", "_Ih3p6"};
    float MyIhCut = MyIhCut_values[0];
    outfilename_ += MyIhCut_names[0];

    std::cout << "Input file:      " << DataSetName << std::endl;
    std::cout << "Output file:     " << outfilename_ << std::endl;
    TFile* ifile = new TFile((filename+".root").c_str());
    TFile* ofile = new TFile((outfilename_+".root").c_str(),"RECREATE");

    

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
    bckgEstimate(DataSetName, st_sample, rc_3fp8, rc_3fp8, rbc_8fp9, ra_3fp8, rd_8fp9, ifIhpSAME, rb_8fp9,
                 "8fp9", nPE, useFit, useOldIhFit, useOld1oPFit, corrTemplateIh, etaName, saveFits, fitIh, fitP, rebinp, 
                 MyIhCut, blind);

    // Debug by taking templates in D
    // bckgEstimate(DataSetName, st_sample, rd_8fp9, rd_8fp9, rbc_8fp9, ra_3fp8, rd_8fp9, ifIhpSAME, rb_8fp9,
    //              "8fp9", nPE, useFit, useOldIhFit, useOld1oPFit, corrTemplateIh, etaName, saveFits, fitIh, fitP, rebinp, 
    //              MyIhCut, blind);


        // In 9fp10 with A and C in 3fp9
    // blind = true;
    // bckgEstimate(DataSetName, st_sample, rc_3fp9, rc_3fp9, rbc_9fp10, ra_3fp9, rd_9fp10, ifIhpSAME, rb_9fp10,
    //              "9fp10", nPE, useFit, useOldIhFit, useOld1oPFit, corrTemplateIh, etaName, saveFits, fitIh, fitP, rebinp, 
    //              MyIhCut, blind);
    


        // Only in C: no Ih reweighting in Region.h to change
    //bckgEstimate(st_sample, rc_3fp8, rc_3fp8, rbc_3fp8, rc_3fp8, rc_3fp8, ifIhpSAME, rc_3fp8, "3fp8", nPE, fitIh, fitP, blind);
    

    delete ofile;
    return;
}
