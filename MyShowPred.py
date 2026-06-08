import ROOT
import os

regions = '8fp9'
option = '_etaAbs_rebinEta' #'_Ihcut'
ofile = 'mass_plot'

ifile = '/safe/ui3_1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/macros/DataMET_2024_V12p24__' + regions + option + '/Eta1_2p4/JetMET2024_V12p24_rebinEta8_rebinIh4_rebinP2_EtaReweighting_Eta1_2p4_NewFit.root'
odir = '/safe/ui3_1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/macros/DataMET_2024_V12p24__' + regions + option + '/Eta1_2p4/Plots_NewFit_' + regions
labelName = "METanalysis_Eta1_2p4"


command = "python3 MyMacroMass.py --ifile {} --labelName {} --ofile {}_NoC --region {} --odir {}".format(ifile, labelName, ofile, regions, odir)

print ('')
print ('       Running:')
print("{}\n".format(command))
print ('')
os.system(command)

