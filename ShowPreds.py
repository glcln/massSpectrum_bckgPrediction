import ROOT
import os

regions = ['9fp10', 
           '9fp10',
           '9fp10']

ifile = ['/opt/sbg/cms/safe1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/macros/DataMET_2024_V12p24/Eta1/JetMET2024_V12p24_rebinEta4_rebinIh4_rebinP2_EtaReweighting_Eta1_NewFit.root',
         '/opt/sbg/cms/safe1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/macros/DataMET_2024_V12p24/Eta1_2p4/JetMET2024_V12p24_rebinEta4_rebinIh4_rebinP2_EtaReweighting_Eta1_2p4_NewFit.root',
         '/opt/sbg/cms/safe1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/macros/DataMET_2024_V12p24/Eta2p4/JetMET2024_V12p24_rebinEta4_rebinIh4_rebinP2_EtaReweighting_Eta2p4_NewFit.root']

odir = ['/opt/sbg/cms/safe1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/macros/DataMET_2024_V12p24/Eta1/Plots_NewFit_'     + regions[0],
        '/opt/sbg/cms/safe1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/macros/DataMET_2024_V12p24/Eta1_2p4/Plots_NewFit_' + regions[1],
        '/opt/sbg/cms/safe1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/macros/DataMET_2024_V12p24/Eta2p4/Plots_NewFit_'   + regions[2]]

labelName = ["METanalysis_Eta1",
             "METanalysis_Eta1_2p4",
             "METanalysis_Eta2p4"]

OnlyNominal = True

for j in range(len(regions)):
    command = "python3 MyMacroMass.py --ifile {} --labelName {} --ofile mass_plot --region {} --odir {} --nom {}".format(
        ifile[j], labelName[j], regions[j], odir[j], OnlyNominal)
    
    print('')
    print('       Running:')
    print("{}\n".format(command))
    print('')
    os.system(command)