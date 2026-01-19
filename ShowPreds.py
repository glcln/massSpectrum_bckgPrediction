import ROOT
import os

regions = ["8fp9"]

ifile = '/opt/sbg/cms/safe1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/macros/DataMET_2024_V12p20/DataMET_2024_test_miniAOD_V12p20_Eta2p4_cutIndex3_rebinEta4_rebinIh4_rebinP2_rebinMass1_EtaReweighting.root'

odir = '/opt/sbg/cms/safe1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/macros/DataMET_2024_V12p20/Plots_Gstrip_Eta1_NewFit_NoC'


for i in regions:
    
    command = "python3 MyMacroMass.py --ifile {} --ofile mass_plot_NoC --region {} --odir {}".format(ifile, i, odir)
    
    print ('')
    print ('       Running:')
    print("{}\n".format(command))
    print ('')
    os.system(command)
