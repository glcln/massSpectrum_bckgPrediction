import ROOT
import os

regions = '8fp9'
option = '_etaAbs_chi2cut'

ifile = '/opt/sbg/cms/safe1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/macros/DataMET_2024_V12p24__8fp9' + option + '/Eta2p4/JetMET2024_V12p24_rebinEta4_rebinIh4_rebinP2_EtaReweighting_Eta2p4_NewFit.root'
odir = '/opt/sbg/cms/safe1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/macros/DataMET_2024_V12p24__8fp9' + option + '/Eta2p4/Plots_NewFit_' + regions
labelName = "METanalysis_Eta2p4"


command = "python3 MyMacroMass.py --ifile {} --labelName {} --ofile mass_plot_NoC --region {} --odir {}".format(ifile, labelName, regions, odir)

print ('')
print ('       Running:')
print("{}\n".format(command))
print ('')
os.system(command)

