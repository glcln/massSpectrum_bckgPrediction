import ROOT
import os

regions = ['9fp10', 
           '9fp10',
           '9fp10',
           '9fp10',
           '9fp10',
           '9fp10',
           '9fp10',
           '9fp10',
           '9fp10',
           '9fp10'
           ]

labelName = '/safe/ui3_1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/macros/DataMET_2024_V12p31__9fp10_etaRebinPerso_Oldfit__PUppiMETcut__IhCuts/'

ifile = [labelName + 'Eta1_2p4_Ih3p1/JetMET2024_V12p31_rebinEta4_rebinIh4_rebinP2_EtaReweighting_Eta1p2_2p4_Ih3p1_OldFit_IhC.root',
         labelName + 'Eta1_2p4_Ih3p2/JetMET2024_V12p31_rebinEta4_rebinIh4_rebinP2_EtaReweighting_Eta1p2_2p4_Ih3p2_OldFit_IhC.root',
         labelName + 'Eta1_2p4_Ih3p3/JetMET2024_V12p31_rebinEta4_rebinIh4_rebinP2_EtaReweighting_Eta1p2_2p4_Ih3p3_OldFit_IhC.root',
         labelName + 'Eta1_2p4_Ih3p4/JetMET2024_V12p31_rebinEta4_rebinIh4_rebinP2_EtaReweighting_Eta1p2_2p4_Ih3p4_OldFit_IhC.root',
         labelName + 'Eta1_2p4_Ih3p5/JetMET2024_V12p31_rebinEta4_rebinIh4_rebinP2_EtaReweighting_Eta1p2_2p4_Ih3p5_OldFit_IhC.root',
         labelName + 'Eta1_2p4_Ih3p6/JetMET2024_V12p31_rebinEta4_rebinIh4_rebinP2_EtaReweighting_Eta1p2_2p4_Ih3p6_OldFit_IhC.root',
         labelName + 'Eta1_2p4_Ih3p7/JetMET2024_V12p31_rebinEta4_rebinIh4_rebinP2_EtaReweighting_Eta1p2_2p4_Ih3p7_OldFit_IhC.root',
         labelName + 'Eta1_2p4_Ih3p8/JetMET2024_V12p31_rebinEta4_rebinIh4_rebinP2_EtaReweighting_Eta1p2_2p4_Ih3p8_OldFit_IhC.root',
         labelName + 'Eta1_2p4_Ih3p9/JetMET2024_V12p31_rebinEta4_rebinIh4_rebinP2_EtaReweighting_Eta1p2_2p4_Ih3p9_OldFit_IhC.root',
         labelName + 'Eta1_2p4_Ih4p0/JetMET2024_V12p31_rebinEta4_rebinIh4_rebinP2_EtaReweighting_Eta1p2_2p4_Ih4p0_OldFit_IhC.root'
        ]

odir = [labelName + 'Eta1_2p4_Ih3p1/',
        labelName + 'Eta1_2p4_Ih3p2/',
        labelName + 'Eta1_2p4_Ih3p3/',
        labelName + 'Eta1_2p4_Ih3p4/',
        labelName + 'Eta1_2p4_Ih3p5/',
        labelName + 'Eta1_2p4_Ih3p6/',
        labelName + 'Eta1_2p4_Ih3p7/',
        labelName + 'Eta1_2p4_Ih3p8/',
        labelName + 'Eta1_2p4_Ih3p9/',
        labelName + 'Eta1_2p4_Ih4p0/'
        ]

labelNameHist = ["METanalysis_TestPUppiMETCut_Eta1_2p4_Ih3p1",
             "METanalysis_TestPUppiMETCut_Eta1_2p4_Ih3p2",
             "METanalysis_TestPUppiMETCut_Eta1_2p4_Ih3p3",
             "METanalysis_TestPUppiMETCut_Eta1_2p4_Ih3p4",
             "METanalysis_TestPUppiMETCut_Eta1_2p4_Ih3p5",
             "METanalysis_TestPUppiMETCut_Eta1_2p4_Ih3p6",
             "METanalysis_TestPUppiMETCut_Eta1_2p4_Ih3p7",
             "METanalysis_TestPUppiMETCut_Eta1_2p4_Ih3p8",
             "METanalysis_TestPUppiMETCut_Eta1_2p4_Ih3p9",
             "METanalysis_TestPUppiMETCut_Eta1_2p4_Ih4p0"
            ]

OnlyNominal = True

for j in range(len(regions)):
    command = "python3 MyMacroMass.py --ifile {} --labelName {} --ofile mass_plot --region {} --odir {} --nom {}".format(
        ifile[j], labelNameHist[j], regions[j], odir[j], OnlyNominal)
    
    print('')
    print('       Running:')
    print("{}\n".format(command))
    print('')
    os.system(command)