import ROOT
import os

regions = '8fp9'
option = '_etaAbs_Ihcut_v2' #'_Ihcut'

ifile = '/safe/ui3_1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/macros/DataMET_2024_V12p31__8fp9_etaAbs_etaRebinPerso_1oPRebinBig_Oldfit__PUppiMETcut/Eta2p4/JetMET2024_V12p31_rebinEta4_rebinIh4_rebinP4_EtaReweighting_Eta2p4_OldFit_IhC_bis.root'
odir  = '/safe/ui3_1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/macros/DataMET_2024_V12p31__8fp9_etaAbs_etaRebinPerso_1oPRebinBig_Oldfit__PUppiMETcut/Eta2p4/Plots_' + regions
labelName = "METanalysis_TestPUppiMETCut_Eta2p4"


OnlyNominal = True


command = "python3 MyMacroMass.py --ifile {} --labelName {} --ofile mass_plot --region {} --odir {} --nom {}".format(ifile, labelName, regions, odir, OnlyNominal)

print ('')
print ('       Running:')
print("{}\n".format(command))
print ('')
os.system(command)

