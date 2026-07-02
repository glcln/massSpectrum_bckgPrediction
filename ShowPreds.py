import ROOT
import os

regions = ['8fp9', 
           '8fp9',
           '8fp9']

ifile = ['/safe/ui3_1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/macros/DataMET_2024_V12p31__8fp9_etaAbs_etaRebinPerso_1oPRebinBig_Oldfit__PUppiMETcut/Eta1/JetMET2024_V12p31_rebinEta4_rebinIh4_rebinP4_EtaReweighting_Eta1_OldFit_IhC.root',
         '/safe/ui3_1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/macros/DataMET_2024_V12p31__8fp9_etaAbs_etaRebinPerso_1oPRebinBig_Oldfit__PUppiMETcut/Eta1_2p4/JetMET2024_V12p31_rebinEta4_rebinIh4_rebinP4_EtaReweighting_Eta1_2p4_OldFit_IhC.root',
         '/safe/ui3_1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/macros/DataMET_2024_V12p31__8fp9_etaAbs_etaRebinPerso_1oPRebinBig_Oldfit__PUppiMETcut/Eta2p4/JetMET2024_V12p31_rebinEta4_rebinIh4_rebinP4_EtaReweighting_Eta2p4_OldFit_IhC_bis.root']

odir = ['/safe/ui3_1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/macros/DataMET_2024_V12p31__8fp9_etaAbs_etaRebinPerso_1oPRebinBig_Oldfit__PUppiMETcut/Eta1/Plots_'     + regions[0],
        '/safe/ui3_1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/macros/DataMET_2024_V12p31__8fp9_etaAbs_etaRebinPerso_1oPRebinBig_Oldfit__PUppiMETcut/Eta1_2p4/Plots_' + regions[1],
        '/safe/ui3_1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/macros/DataMET_2024_V12p31__8fp9_etaAbs_etaRebinPerso_1oPRebinBig_Oldfit__PUppiMETcut/Eta2p4/Plots_'   + regions[2]]

labelName = ["METanalysis_TestPUppiMETCut_Eta1",
             "METanalysis_TestPUppiMETCut_Eta1_2p4",
             "METanalysis_TestPUppiMETCut_Eta2p4"]

OnlyNominal = True

for j in range(len(regions)):
    command = "python3 MyMacroMass.py --ifile {} --labelName {} --ofile mass_plot --region {} --odir {} --nom {}".format(
        ifile[j], labelName[j], regions[j], odir[j], OnlyNominal)
    
    print('')
    print('       Running:')
    print("{}\n".format(command))
    print('')
    os.system(command)