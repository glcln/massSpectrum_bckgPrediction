import ROOT
import os

regions     = '9fp10'
eta         = 'Eta2p4'
OnlyNominal = True
isMC        = True


ifile = '/safe/ui3_1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/macros/DataMET_2024_V12p31__' + regions + '_etaRebinPerso_Oldfit__PUppiMETcut/BKG/' + eta + '/HistForBkg_MC_rebinEta4_rebinIh4_rebinP2_EtaReweighting_' + eta + '_OldFit_IhC.root'
odir  = '/safe/ui3_1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/macros/DataMET_2024_V12p31__' + regions + '_etaRebinPerso_Oldfit__PUppiMETcut/BKG/' + eta + '/Plots_' + regions
labelName = 'METanalysis_TestPUppiMETCut_' + eta



command = "python3 MyMacroMass.py --ifile {} --labelName {} --ofile mass_plot_Cnoabs --region {} --odir {} --nom {} --eta {} --isMC {}".format(ifile, labelName, regions, odir, OnlyNominal, eta, isMC)

print ('')
print ('       Running:')
print("{}\n".format(command))
print ('')
os.system(command)