import sys, os, time, re
import numpy as np
from optparse import OptionParser
parser = OptionParser(usage="Usage: python %prog codeVersion")
(opt,args) = parser.parse_args()

datasetList = [
    "/safe/ui3_1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/output/JetMET2024_V12/JetMET2024_V12p24",
]

nPE = "200"

#[label, rebinEta, rebinIh, rebinMom, fitIh, fitMom, useFit, corrTemplateIh]
config = [
    ["nominal", "4", "4", "4", "1", "1", "1", "0"],
    #["etaup", "2", "4", "4", "1", "1", "1", "0"],
    #["etadown", "8", "4", "4", "1", "1", "1", "0"],
    #["ihup", "4", "2", "4", "1", "1", "1", "0"],
    #["ihdown", "4", "8", "4", "1", "1", "1", "0"],
    #["momup", "4", "4", "2", "1", "1", "1", "0"],
    #["momdown", "4", "4", "8", "1", "1", "1", "0"],
    #["useFit", "4", "4", "4", "1", "1", "0", "0"],
    #["FitIhUp", "4", "4", "4", "2", "1", "1", "0"],
    #["FitIhDown", "4", "4", "4", "0", "1", "1", "0"],
    #["FitMomUp", "4", "4", "4", "1", "2", "1", "0"],
    #["FitMomDown", "4", "4", "4", "1", "0", "1", "0"],
    #["corrTemplateIh", "4", "4", "4", "1", "1", "1", "1"],
]


i = 0
for dataset in datasetList:
    print('')
    print("Launch on dataset:    " + dataset)
    print('')
    os.system("cp configFile_readHist_template.txt configFile_readHisto_toLaunch_tmp.txt")
    os.system("sed -i 's|sample|" + dataset + "|g' configFile_readHisto_toLaunch_tmp.txt")
    os.system("sed -i 's|nPE|" + nPE + "|g' configFile_readHisto_toLaunch_tmp.txt")
    for conf in config:
        os.system("cp configFile_readHisto_toLaunch_tmp.txt configFile_readHisto_toLaunch.txt")
        os.system("sed -i 's|label|" + conf[0] + "|g' configFile_readHisto_toLaunch.txt")
        os.system("sed -i 's|rebinEta|" + conf[1] + "|g' configFile_readHisto_toLaunch.txt")
        os.system("sed -i 's|rebinIh|" + conf[2] + "|g' configFile_readHisto_toLaunch.txt")
        os.system("sed -i 's|rebinMom|" + conf[3] + "|g' configFile_readHisto_toLaunch.txt")
        os.system("sed -i 's|fitIh|" + conf[4] + "|g' configFile_readHisto_toLaunch.txt")
        os.system("sed -i 's|fitMom|" + conf[5] + "|g' configFile_readHisto_toLaunch.txt")
        os.system("sed -i 's|useFit|" + conf[6] + "|g' configFile_readHisto_toLaunch.txt")
        os.system("sed -i 's|corrTemplateIh|" + conf[7] + "|g' configFile_readHisto_toLaunch.txt")

        os.system("cat configFile_readHisto_toLaunch.txt")
        os.system("time root -l -q -b step2_backgroundPrediction.C")#  &> exit.txt")
        #os.system("root -l -b -q 'exit_displayChi2.C(\"exit.txt\",\"exit_displayChi2_" + conf[0] + ".root\")'")    
    i += 1

