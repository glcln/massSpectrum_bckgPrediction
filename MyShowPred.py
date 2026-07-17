import ROOT
import os

regions     = '8fp9'
eta         = 'Eta1'
OnlyNominal = True
isMC        = False
#HistForBkg_MC
#


ifile = '/safe/ui3_1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/macros/DataMET_2024_V12p31__' + regions + '_etaRebinPerso_Oldfit__PUppiMETcut/' + eta + '/JetMET2024_V12p32_rebinEta4_rebinIh4_rebinP2_EtaReweighting_' + eta + '_OldFit_IhC.root'
odir  = '/safe/ui3_1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/macros/DataMET_2024_V12p31__' + regions + '_etaRebinPerso_Oldfit__PUppiMETcut/' + eta + '/Plots_' + regions
labelName = 'METanalysis_TestPUppiMETCut_' + eta



command = "python3 MyMacroMass.py --ifile {} --labelName {} --ofile mass_plot_Cnoabs --region {} --odir {} --nom {} --eta {} --isMC {}".format(ifile, labelName, regions, odir, OnlyNominal, eta, isMC)

print ('')
print ('       Running:')
print("{}\n".format(command))
print ('')
os.system(command)

