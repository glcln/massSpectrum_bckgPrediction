import sys, os, time, re
import numpy as np
from optparse import OptionParser
parser = OptionParser(usage="Usage: python %prog codeVersion")
(opt,args) = parser.parse_args()

datasetList = [
    "/opt/sbg/cms/safe1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/macros/DataMET_2024_test_miniAOD_V12p20_Eta2p4"
]
tagKC=[
    "data2024"
]
odir=[
    "\/opt\/sbg\/cms\/safe1\/cms\/gcoulon\/CMSSW_15_0_13_patch1\/src\/TupleAnalysis\/macros\/todelete"
]

nPE = "200"

#[label, rebinEta, rebinIh, rebinMom, corrTemplateIh, corrTemplateMom, fitIh, fitMom, useFit]
config=[
    ["nominal", "4", "4", "2", "0", "0", "1", "1", "1"],
    #["etaup", "2", "4", "2", "0", "0", "1", "1", "1"],
    #["etadown", "8", "4", "2", "0", "0", "1", "1", "1"],
    #["ihup", "4", "2", "2", "0", "0", "1", "1", "1"],
    #["ihdown", "4", "8", "2", "0", "0", "1", "1", "1"],
    #["momup", "4", "4", "1", "0", "0", "1", "1", "1"],
    #["momdown", "4", "4", "4", "0", "0", "1", "1", "1"],
    #["useFit", "4", "4", "2", "0", "0", "1", "1", "0"],
    #["corrIh", "4", "4", "2", "1", "0", "1", "1", "1"],
    #["FitIhUp", "4", "4", "2", "0", "0", "2", "1", "1"],
    #["FitIhDown", "4", "4", "2", "0", "0", "0", "1", "1"],
    #["FitMomUp", "4", "4", "2", "0", "0", "1", "2", "1"],
    #["FitMomDown", "4", "4", "2", "0", "0", "1", "0", "1"],
]


i = 0
for dataset in datasetList:
    print('')
    print("Launch on dataset:    " + dataset)
    print('')
    os.system("cp configFile_readHist_template.txt configFile_readHisto_toLaunch_tmp.txt")
    os.system("sed -i 's|sample|" + dataset + "|g' configFile_readHisto_toLaunch_tmp.txt")
    os.system("sed -i 's/tag_KC/" + tagKC[i] + "/g' configFile_readHisto_toLaunch_tmp.txt")
    os.system("sed -i 's/dir/" + odir[i] + "/g' configFile_readHisto_toLaunch_tmp.txt")
    os.system("sed -i 's/nPE/" + nPE + "/g' configFile_readHisto_toLaunch_tmp.txt")
    for conf in config:
        os.system("cp configFile_readHisto_toLaunch_tmp.txt configFile_readHisto_toLaunch.txt")
        os.system("sed -i 's/label/" + conf[0] + "/g' configFile_readHisto_toLaunch.txt")
        os.system("sed -i 's/rebinEta/" + conf[1] + "/g' configFile_readHisto_toLaunch.txt")
        os.system("sed -i 's/rebinIh/" + conf[2] + "/g' configFile_readHisto_toLaunch.txt")
        os.system("sed -i 's/rebinMom/" + conf[3] + "/g' configFile_readHisto_toLaunch.txt")
        os.system("sed -i 's/corrTemplateIh/" + conf[4] + "/g' configFile_readHisto_toLaunch.txt")
        os.system("sed -i 's/corrTemplateMom/" + conf[5] + "/g' configFile_readHisto_toLaunch.txt")
        os.system("sed -i 's/fitIh/" + conf[6] + "/g' configFile_readHisto_toLaunch.txt")
        os.system("sed -i 's/fitMom/" + conf[7] + "/g' configFile_readHisto_toLaunch.txt")
        os.system("sed -i 's/useFit/" + conf[8] + "/g' configFile_readHisto_toLaunch.txt")

        os.system("cat configFile_readHisto_toLaunch.txt")
        os.system("time root -l -q -b step2_backgroundPrediction.C")    
    i += 1

