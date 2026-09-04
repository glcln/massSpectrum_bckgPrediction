import ROOT
import os

regions     = '8fp9'
eta         = 'Eta2p4'
version     = '3'
OnlyNominal = True
isMC        = True
option      = 'SigmaPtoverPt_0p5_EoP_0p1'


ifile = '/safe/ui3_1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/macros/MC_2024_V12p35' + '__' + regions + '_' + option + '/' + eta + '/HistForBkg_MC_V' + version + '_rebinEta4_rebinIh4_rebinP2_EtaReweighting_' + eta + '_' + option + '_OldFit_IhC.root'
odir  = '/safe/ui3_1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/macros/MC_2024_V12p35' + '__' + regions + '_' + option + '/' + eta + '/Plots_' + regions
labelName = 'METanalysis_TestPUppiMETCut_' + option + '_' + eta



command = "python3 MyMacroMass.py --ifile {} --labelName {} --ofile mass_plot_Cnoabs --region {} --odir {} --nom {} --eta {} --isMC {}".format(ifile, labelName, regions, odir, OnlyNominal, eta, isMC)

print ('')
print ('       Running:')
print("{}\n".format(command))
print ('')
os.system(command)