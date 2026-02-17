import ROOT
import os

regions = ['8fp9', 
           '8fp9',
           '8fp9']

ifile = ['/opt/sbg/cms/safe1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/macros/DataMET_2024_V12p24/JetMET2024_V12p24_cutIndex3_rebinEta4_rebinIh4_rebinP2_rebinMass1_EtaReweighting_Eta1_NewFit.root',
         '/opt/sbg/cms/safe1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/macros/DataMET_2024_V12p24/JetMET2024_V12p24_cutIndex3_rebinEta4_rebinIh4_rebinP2_rebinMass1_EtaReweighting_Eta1_2p4_NewFit.root',
         '/opt/sbg/cms/safe1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/macros/DataMET_2024_V12p24/JetMET2024_V12p24_cutIndex3_rebinEta4_rebinIh4_rebinP2_rebinMass1_EtaReweighting_Eta2p4_NewFit.root']

odir = ['/opt/sbg/cms/safe1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/macros/DataMET_2024_V12p24/Plots_2024_Fpix_Eta1_NewFit_NoC',
        '/opt/sbg/cms/safe1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/macros/DataMET_2024_V12p24/Plots_2024_Fpix_Eta1_2p4_NewFit_NoC',
        '/opt/sbg/cms/safe1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/macros/DataMET_2024_V12p24/Plots_2024_Fpix_Eta2p4_NewFit_NoC']

labelName = ["METanalysis_Eta1",
             "METanalysis_Eta1_2p4",
             "METanalysis_Eta2p4"]


j = 0
for i in regions:
    
    command = "python3 MyMacroMass.py --ifile {} --labelName {} --ofile mass_plot_NoC --region {} --odir {}".format(ifile[j], labelName[j], i, odir[j])
    
    print ('')
    print ('       Running:')
    print("{}\n".format(command))
    print ('')
    os.system(command)

    j += 1
