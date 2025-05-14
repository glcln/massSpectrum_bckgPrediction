import ROOT
import os

#regions = ["80ias90"]
#regions = ["3fp8"]
regions = ["8fp9"]

# Fpix VS Gstrip
ifile = '/opt/sbg/cms/ui3_data1/gcoulon/CMSSW_10_6_30/src/HSCPTreeAnalyzer/macros/Fpix_V3p2/MET_2017_2018_massCut_0_pT70_V3p2_Fpix_Eta2p4_cutIndex3_rebinEta4_rebinIh4_rebinP2_rebinMass1_EtaReweighting_NoFit.root'

odir = '/opt/sbg/cms/ui3_data1/gcoulon/CMSSW_10_6_30/src/HSCPTreeAnalyzer/macros/Fpix_V3p2/Plots_Fpix_Eta2p4_NoFit_MassRebin_NoC'


for i in regions:
    
    command = "python2.7 MyMacroMass.py --ifile {} --ofile mass_plot_NoC --region {} --odir {}".format(ifile, i, odir)
    
    print ''
    print '       Running:'
    print("{}\n".format(command))
    print ''
    os.system(command)