import ROOT, sys, os, time, re, numpy
from optparse import OptionParser
from tqdm import tqdm
sys.path.append("/safe/ui3_1/cms/gcoulon/")


# SETUP

version = "16p3" # just the number, like 18p2
SelectionNames = [] # "Nm1", "EventCutflow"
intLumi = 108950.0 # 2024 in pb^-1

# cross-section [pb]
crossSectionArray = { 
"QCD2024_mu_pt15to20_V" + version + "1.root" : 2982000.0, # +/- 9502.0
"QCD2024_mu_pt20to30_V" + version + "2.root" : 2679000.0, # +/- 8497.0
"QCD2024_mu_pt30to50_V" + version + "3.root" : 1465000.0, # +/- 4593.0
"QCD2024_mu_pt50to80_V" + version + "4.root" : 409500.0, # +/- 1276.0
"QCD2024_mu_pt80to120_V" + version + "5.root" : 96200.0, # +/- 299.4
"QCD2024_mu_pt120to170_V" + version + "6.root" : 22980.0, #+/- 71.7
"QCD2024_mu_pt170to300_V" + version + "7.root" : 7763.0, # +/- 23.8
"QCD2024_mu_pt300to470_V" + version + "8.root" : 699.1, # +/- 2.13
"QCD2024_mu_pt470to600_V" + version + "9.root" : 68.24, # +/- 0.2042
"QCD2024_mu_pt600to800_V" + version + "10.root" : 21.37, # +/- 0.06407
"QCD2024_mu_pt800to1000_V" + version + "11.root" : 3.913, # +/- 0.01164
"QCD2024_mu_pt1000_V" + version + "12.root" : 1.323, # +/- 0.003969	

"TTbar2024_V" + version + ".root" : 98.1, # +2.4/-3.5 
"TTbarSemiLep2024_V" + version + ".root" : 405, # +2.4/-3.5 

"Wjets2024_1J_pt40to100_V" + version + "1.root" : 4211, # +/- 27.71  # 0.0645822 * 63425.1
"Wjets2024_1J_pt100to200_V" + version + "2.root" : 342.3, # +/- 1.842  # 0.00538203 * 63425.1
"Wjets2024_1J_pt200to400_V" + version + "3.root" : 21.84, # +/- 0.1076	  # 0.000373607 * 63425.1
"Wjets2024_1J_pt400to600_V" + version + "4.root" : 0.6845, # +/- 0.003095	  # 1.26886e-05 * 63425.1
"Wjets2024_1J_pt600_V" + version + "5.root" : 0.07753,	# +/- 0.0003313	  # 1.59701e-06 * 63425.1
"Wjets2024_2J_pt40to100_V" + version + "6.root" : 1581, # +/- 16.18  # 0.0237533 * 63425.1
"Wjets2024_2J_pt100to200_V" + version + "7.root" : 411.1, # +/- 3.841  # 0.00640614 * 63425.1
"Wjets2024_2J_pt200to400_V" + version + "8.root" : 53.59, # +/- 0.4098  # 0.000845196 * 63425.1
"Wjets2024_2J_pt400to600_V" + version + "9.root" : 3.099, # +/- 0.01935  # 4.87396e-05 * 63425.1
"Wjets2024_2J_pt600_V" + version + "10.root" : 0.5259, # +/- 0.002768  # 8.27863e-06 * 63425.1

"WjetMuNu2024_V" + version + ".root" : 20790, # +/- 90.5  #63425.1/3
}


# SAMPLES

pathBKG_QCD = "/safe/ui3_1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/output/QCD2024_V16/"
pathBKG_TTbar = "/safe/ui3_1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/output/TTbar2024_V15/"
pathBKG_Wjets = "/safe/ui3_1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/output/Wjets2024_V14/"

'''
BackgroundSamples_HISTO = [
pathBKG_QCD + "QCD2024_mu_pt15to20_V" + version + "1.root",
pathBKG_QCD + "QCD2024_mu_pt20to30_V" + version + "2.root",
pathBKG_QCD + "QCD2024_mu_pt30to50_V" + version + "3.root",
pathBKG_QCD + "QCD2024_mu_pt50to80_V" + version + "4.root",
pathBKG_QCD + "QCD2024_mu_pt80to120_V" + version + "5.root",
pathBKG_QCD + "QCD2024_mu_pt120to170_V" + version + "6.root",
pathBKG_QCD + "QCD2024_mu_pt170to300_V" + version + "7.root",
pathBKG_QCD + "QCD2024_mu_pt300to470_V" + version + "8.root",
pathBKG_QCD + "QCD2024_mu_pt470to600_V" + version + "9.root",
pathBKG_QCD + "QCD2024_mu_pt600to800_V" + version + "10.root",
pathBKG_QCD + "QCD2024_mu_pt800to1000_V" + version + "11.root",
pathBKG_QCD + "QCD2024_mu_pt1000_V" + version + "12.root",

pathBKG_TTbar + "TTbar2024_V" + version + ".root",
pathBKG_TTbar + "TTbarSemiLep2024_V" + version + ".root",

pathBKG_Wjets + "Wjets2024_1J_pt40to100_V" + version + "1.root",
pathBKG_Wjets + "Wjets2024_1J_pt100to200_V" + version + "2.root",
pathBKG_Wjets + "Wjets2024_1J_pt200to400_V" + version + "3.root",
pathBKG_Wjets + "Wjets2024_1J_pt400to600_V" + version + "4.root",
pathBKG_Wjets + "Wjets2024_1J_pt600_V" + version + "5.root",
pathBKG_Wjets + "Wjets2024_2J_pt40to100_V" + version + "6.root",
pathBKG_Wjets + "Wjets2024_2J_pt100to200_V" + version + "7.root",
pathBKG_Wjets + "Wjets2024_2J_pt200to400_V" + version + "8.root",
pathBKG_Wjets + "Wjets2024_2J_pt400to600_V" + version + "9.root",
pathBKG_Wjets + "Wjets2024_2J_pt600_V" + version + "10.root"

pathBKG_Wjets + "WjetMuNu2024_V" + version + ".root",
]
'''

BackgroundSamples_HISTO = [
pathBKG_QCD + "QCD2024_mu_pt15to20_V" + version + "1.root",
pathBKG_QCD + "QCD2024_mu_pt20to30_V" + version + "2.root",
pathBKG_QCD + "QCD2024_mu_pt30to50_V" + version + "3.root",
pathBKG_QCD + "QCD2024_mu_pt50to80_V" + version + "4.root",
pathBKG_QCD + "QCD2024_mu_pt80to120_V" + version + "5.root",
pathBKG_QCD + "QCD2024_mu_pt120to170_V" + version + "6.root",
pathBKG_QCD + "QCD2024_mu_pt170to300_V" + version + "7.root",
pathBKG_QCD + "QCD2024_mu_pt300to470_V" + version + "8.root",
pathBKG_QCD + "QCD2024_mu_pt470to600_V" + version + "9.root",
pathBKG_QCD + "QCD2024_mu_pt600to800_V" + version + "10.root",
pathBKG_QCD + "QCD2024_mu_pt800to1000_V" + version + "11.root",
pathBKG_QCD + "QCD2024_mu_pt1000_V" + version + "12.root",
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
      if not any(tag in name for tag in SelectionNames): h.Scale(weight)

      h.Write()

    else: obj.Write()

  fileOut.Close()
  fileIn.Close()

  print(" -> written:", weighted_path, "\n")
