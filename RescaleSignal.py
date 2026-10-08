import ROOT, sys, os, time, re, numpy
from optparse import OptionParser
from tqdm import tqdm
sys.path.append("/safe/ui3_1/cms/gcoulon/")


# SETUP

version = "20p0" # just the number, like 18p2
SelectionNames = []
intLumi = 108950.0 # 2024 in pb^-1


crossSectionArray = {
"Gluino_Run3_MET_1000_V" + version + ".root" : 4.842E-01,
"Gluino_Run3_MET_1200_V" + version + ".root" : 1.288E-01,
"Gluino_Run3_MET_1400_V" + version + ".root" : 3.883E-02,
"Gluino_Run3_MET_1600_V" + version + ".root" : 1.282E-02,
"Gluino_Run3_MET_1800_V" + version + ".root" : 4.524E-03,
"Gluino_Run3_MET_2000_V" + version + ".root" : 1.684E-03,
"Gluino_Run3_MET_2200_V" + version + ".root" : 6.535E-04,
"Gluino_Run3_MET_2400_V" + version + ".root" : 2.627E-04,
"Gluino_Run3_MET_2600_V" + version + ".root" : 1.089E-04,

"Gluino_Run3_MET_madgraph_1100_V" + version + ".root" : 2.450E-01,
"Gluino_Run3_MET_madgraph_1200_V" + version + ".root" : 1.288E-01,
"Gluino_Run3_MET_madgraph_1300_V" + version + ".root" : 6.978E-02,
"Gluino_Run3_MET_madgraph_1400_V" + version + ".root" : 3.883E-02,
"Gluino_Run3_MET_madgraph_1600_V" + version + ".root" : 1.282E-02,
"Gluino_Run3_MET_madgraph_1800_V" + version + ".root" : 4.524E-03,
"Gluino_Run3_MET_madgraph_2000_V" + version + ".root" : 1.684E-03,
"Gluino_Run3_MET_madgraph_2200_V" + version + ".root" : 6.535E-04,
"Gluino_Run3_MET_madgraph_2400_V" + version + ".root" : 2.627E-04,
"Gluino_Run3_MET_madgraph_2600_V" + version + ".root" : 1.089E-04,

"Stop_Run3_MET_madgraph_700_V" + version + ".root" : 9.905000E-02,
"Stop_Run3_MET_madgraph_800_V" + version + ".root" : 4.193000E-02,
"Stop_Run3_MET_madgraph_900_V" + version + ".root" : 1.903000E-02,
"Stop_Run3_MET_madgraph_1000_V" + version + ".root" : 9.123000E-03,
"Stop_Run3_MET_madgraph_1200_V" + version + ".root" : 2.373000E-03,
"Stop_Run3_MET_madgraph_1400_V" + version + ".root" : 6.976000E-04,
"Stop_Run3_MET_madgraph_1600_V" + version + ".root" : 2.236000E-04,
"Stop_Run3_MET_madgraph_1800_V" + version + ".root" : 7.654000E-05,
"Stop_Run3_MET_madgraph_2000_V" + version + ".root" : 2.752000E-05,
"Stop_Run3_MET_madgraph_2200_V" + version + ".root" : 1.029000E-05,
"Stop_Run3_MET_madgraph_2400_V" + version + ".root" : 3.957000E-06,
"Stop_Run3_MET_madgraph_2600_V" + version + ".root" : 1.556000E-06,

"Stau_Run3_MET_247_V" + version + ".root" : 1.472533E-02,
"Stau_Run3_MET_308_V" + version + ".root" : 6.163604E-03,
"Stau_Run3_MET_432_V" + version + ".root" : 1.474467E-03,
"Stau_Run3_MET_557_V" + version + ".root" : 4.561581E-04,
"Stau_Run3_MET_651_V" + version + ".root" : 2.095727E-04,
"Stau_Run3_MET_745_V" + version + ".root" : 1.036057E-04,
"Stau_Run3_MET_871_V" + version + ".root" : 4.263775E-05,
"Stau_Run3_MET_1029_V" + version + ".root" : 1.526866E-05,
"Stau_Run3_MET_1218_V" + version + ".root" : 4.956323E-06,
"Stau_Run3_MET_1409_V" + version + ".root" : 1.715460E-06,
"Stau_Run3_MET_1599_V" + version + ".root" : 6.262536E-07,
}

# SAMPLES
pathSigGluino = "/safe/ui3_1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/output/Gluino_V19/"
pathSigStop   = "/safe/ui3_1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/output/Stop_V21/"
pathSigStau   = "/safe/ui3_1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/output/Stau_V20/"


'''
BackgroundSamples_HISTO = [
pathSigGluino + "Gluino_Run3_MET_1000_V" + version + ".root",
pathSigGluino + "Gluino_Run3_MET_1200_V" + version + ".root",
pathSigGluino + "Gluino_Run3_MET_1400_V" + version + ".root",
pathSigGluino + "Gluino_Run3_MET_1600_V" + version + ".root",
pathSigGluino + "Gluino_Run3_MET_1800_V" + version + ".root",
pathSigGluino + "Gluino_Run3_MET_2000_V" + version + ".root",
pathSigGluino + "Gluino_Run3_MET_2200_V" + version + ".root",
pathSigGluino + "Gluino_Run3_MET_2400_V" + version + ".root",
pathSigGluino + "Gluino_Run3_MET_2600_V" + version + ".root",

pathSigGluino + "Gluino_Run3_MET_madgraph_1100_V" + version + ".root",
pathSigGluino + "Gluino_Run3_MET_madgraph_1200_V" + version + ".root",
pathSigGluino + "Gluino_Run3_MET_madgraph_1300_V" + version + ".root",
pathSigGluino + "Gluino_Run3_MET_madgraph_1400_V" + version + ".root",
pathSigGluino + "Gluino_Run3_MET_madgraph_1600_V" + version + ".root",
pathSigGluino + "Gluino_Run3_MET_madgraph_1800_V" + version + ".root",
pathSigGluino + "Gluino_Run3_MET_madgraph_2000_V" + version + ".root",
pathSigGluino + "Gluino_Run3_MET_madgraph_2200_V" + version + ".root",
pathSigGluino + "Gluino_Run3_MET_madgraph_2400_V" + version + ".root",
pathSigGluino + "Gluino_Run3_MET_madgraph_2600_V" + version + ".root",

pathSigStop + "Stop_Run3_MET_madgraph_700_V" + version + ".root",
pathSigStop + "Stop_Run3_MET_madgraph_800_V" + version + ".root",
pathSigStop + "Stop_Run3_MET_madgraph_900_V" + version + ".root",
pathSigStop + "Stop_Run3_MET_madgraph_1000_V" + version + ".root",
pathSigStop + "Stop_Run3_MET_madgraph_1200_V" + version + ".root",
pathSigStop + "Stop_Run3_MET_madgraph_1400_V" + version + ".root",
pathSigStop + "Stop_Run3_MET_madgraph_1600_V" + version + ".root",
pathSigStop + "Stop_Run3_MET_madgraph_1800_V" + version + ".root",
pathSigStop + "Stop_Run3_MET_madgraph_2000_V" + version + ".root",
pathSigStop + "Stop_Run3_MET_madgraph_2200_V" + version + ".root",
pathSigStop + "Stop_Run3_MET_madgraph_2400_V" + version + ".root",
pathSigStop + "Stop_Run3_MET_madgraph_2600_V" + version + ".root",

pathSigStau + "Stau_Run3_MET_247_V" + version + ".root",
pathSigStau + "Stau_Run3_MET_308_V" + version + ".root",
pathSigStau + "Stau_Run3_MET_432_V" + version + ".root",
pathSigStau + "Stau_Run3_MET_557_V" + version + ".root",
pathSigStau + "Stau_Run3_MET_651_V" + version + ".root",
pathSigStau + "Stau_Run3_MET_745_V" + version + ".root",
pathSigStau + "Stau_Run3_MET_871_V" + version + ".root",
pathSigStau + "Stau_Run3_MET_1029_V" + version + ".root",
pathSigStau + "Stau_Run3_MET_1218_V" + version + ".root",
pathSigStau + "Stau_Run3_MET_1409_V" + version + ".root",
pathSigStau + "Stau_Run3_MET_1599_V" + version + ".root",
]
'''


BackgroundSamples_HISTO = [
pathSigStau + "Stau_Run3_MET_247_V" + version + ".root",
pathSigStau + "Stau_Run3_MET_308_V" + version + ".root",
pathSigStau + "Stau_Run3_MET_432_V" + version + ".root",
pathSigStau + "Stau_Run3_MET_557_V" + version + ".root",
pathSigStau + "Stau_Run3_MET_651_V" + version + ".root",
pathSigStau + "Stau_Run3_MET_745_V" + version + ".root",
pathSigStau + "Stau_Run3_MET_871_V" + version + ".root",
pathSigStau + "Stau_Run3_MET_1029_V" + version + ".root",
pathSigStau + "Stau_Run3_MET_1218_V" + version + ".root",
pathSigStau + "Stau_Run3_MET_1409_V" + version + ".root",
pathSigStau + "Stau_Run3_MET_1599_V" + version + ".root",
]



# RESCALE BACKGROUND HISTOGRAMS

Bkg_HistoFiles = []
for sample in BackgroundSamples_HISTO:
  if not os.path.exists(sample): 
    print("File not found: {}".format(sample))
    continue
  Bkg_HistoFiles.append(ROOT.TFile.Open(sample,"READ"))


for fileIn in Bkg_HistoFiles:

  base = os.path.basename(fileIn.GetName())
  dirname = os.path.dirname(fileIn.GetName())

  # is there an associated xsec ?
  key = base
  if key not in crossSectionArray:
    print("Skipping (no xsection):", key)
    continue

  # Nbr of pre-trigger events and the weight
  cutflow = fileIn.Get("EventCutflow")
  if not cutflow:
    print("No EventCutflow in", base)
    continue

  nEvetsPreTrig = cutflow.GetBinContent(1)
  if nEvetsPreTrig == 0:
    print("Zero events in", base)
    continue

  weight = intLumi * crossSectionArray[key] / nEvetsPreTrig
  print(f"{base} -> nEvetsPreTrig = {nEvetsPreTrig} -> weight = {weight:.3e}")


  # Dupplicate file with weighted histograms
  weighted_name = base.replace(".root", "_weighted.root")
  weighted_path = os.path.join(dirname, weighted_name)
  fileOut = ROOT.TFile(weighted_path, "RECREATE")
  
  param_w = ROOT.TParameter("double")("weight", float(weight))
  param_L = ROOT.TParameter("double")("lumi", float(intLumi))
  param_w.Write()
  param_L.Write()

  for keyObj in fileIn.GetListOfKeys():
    obj = keyObj.ReadObj()
    if not obj: continue

    name = obj.GetName()

    if obj.InheritsFrom("TH1") or obj.InheritsFrom("TH2"):
      h = obj.Clone()
      h.SetDirectory(fileOut)

      # Reweight only histograms from the selections
      if not any(tag in name for tag in SelectionNames):
        if h.GetSumw2N() == 0:
          h.Sumw2()                          # fSumw2[i] = N_i  ->  err = sqrt(N_i)
        h.SetBinErrorOption(ROOT.TH1.kNormal)
        h.Scale(weight)                      # contenu -> w*N , fSumw2 -> w^2*N

      h.Write()

    else: obj.Write()

  fileOut.Close()
  fileIn.Close()

  print(" -> written:", weighted_path, "\n")
