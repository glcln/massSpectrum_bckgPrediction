import ROOT
import os

regions = '8fp9'
ifile = '/opt/sbg/cms/safe1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/macros/DataMET_2024_V12p24/JetMET2024_V12p24_cutIndex3_rebinEta4_rebinIh4_rebinP2_rebinMass1_EtaReweighting_Eta2p4_NewFit.root'
odir = '/opt/sbg/cms/safe1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/macros/DataMET_2024_V12p24/Plots_2024_Fpix_Eta2p4_NewFit_NoC_Newbinning'
labelName = "METanalysis_Eta2p4"


command = "python3 MyMacroMass.py --ifile {} --labelName {} --ofile mass_plot_NoC --region {} --odir {}".format(ifile, labelName, regions, odir)

print ('')
print ('       Running:')
print("{}\n".format(command))
print ('')
os.system(command)

