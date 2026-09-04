#!/usr/bin/python

import sys, os
sys.argv.append( '-b-' )
import ROOT
import math
import array
import ctypes

from ROOT import TFile, THStack, TCanvas, TLegend, TLatex, TPad, TH1, TH2, TLine, TGraph, TGraphErrors
import tdrstyle
from collections import OrderedDict

HARDCODED_SYST = OrderedDict([
    ("lumi", 1.4),   # %
    ("Fpix", 1.6),   # %
])

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning

tdrstyle.setTDRStyle()

def ratioHisto(h2, h1):
    h3 = h1.Clone()
    h3.Divide(h2)
    return h3

def ratioInt(h1, h2):
    h3 = h2.Clone()
    h3.Reset()
    for i in range(0,h1.GetNbinsX()+1):
        e1 = ctypes.c_double(0.0)
        e2 = ctypes.c_double(0.0)
        a = h1.IntegralAndError(i, h1.GetNbinsX() + 1, e1, "")
        b = h2.IntegralAndError(i, h1.GetNbinsX() + 1, e2, "")
        if b != 0 and a != 0:
            c = math.sqrt((e1.value*e1.value)/(a*a) + (e2.value*e2.value)/(b*b)) * a/b
            h3.SetBinContent(i, a/b)
            h3.SetBinError(i, c)
        else:
            h3.SetBinContent(i, 0)
    return h3

def setColorAndMarker(h1,color,markerstyle):
    h1.SetLineColor(color)
    h1.SetMarkerColor(color)
    h1.SetFillColor(color)
    h1.SetMarkerStyle(markerstyle)
    h1.SetMarkerSize(1.2)
    return h1

def overflowInLastBin(h):
    res = h.Clone()
    res.SetBinContent(h.GetNbinsX(), h.GetBinContent(h.GetNbinsX()) + h.GetBinContent(h.GetNbinsX() + 1))
    res.SetBinContent(h.GetNbinsX()+1, 0)
    return res

def underflowInFirstBin(h):
    res = h.Clone()
    res.SetBinContent(1, h.GetBinContent(0) + h.GetBinContent(1))
    res.SetBinContent(0, 0)
    return res

def lowEdge(h1):
    res=ROOT.TGraph(h1.GetNbinsX()-1)
    for i in range (1,h1.GetNbinsX()+1):
        res.SetPoint(i-1,h1.GetBinLowEdge(i),h1.GetBinContent(i))
    return res

def allSet (h, sizeRebinning, rebinning, st):
    h = h.Rebin(sizeRebinning, st, rebinning)
    h = overflowInLastBin(h)
    h = underflowInFirstBin(h)
    return h

def systMass (nominal, down, up, name, typec, binned, mini=0, minNom=0.0):
    if (binned==0):
        ra1 = ratioInt(nominal, up)
        ra2 = ratioInt(nominal, down)
    elif (binned==1):
        ra1 = ratioHisto(nominal, up)
        ra2 = ratioHisto(nominal, down)

    res = ra1.Clone()
    res.SetName(name)

    for i in range (0, res.GetNbinsX()+1):
        if nominal.GetBinContent(i) <= minNom:
            res.SetBinContent(i, 0)
            continue
        r1 = ra1.GetBinContent(i)
        r2 = ra2.GetBinContent(i)
        s1 = abs(1-r1)
        s2 = abs(1-r2)
        if (typec==0):
            m = max(s1,s2)
            if (mini==1):
                m = min(s1,s2)
        elif (typec==1):
            m = (s1+s2)/2.
        res.SetBinContent(i, 100*m)
    return res

def systFlat(nominal, value, name, minNom=0.0):
    res = nominal.Clone(name)
    res.Reset()
    for i in range (0, res.GetNbinsX()+1):
        if nominal.GetBinContent(i) <= minNom:
            res.SetBinContent(i, 0)
        else:
            res.SetBinContent(i, value)
    return res

def systTotal(list_h):
    res = list_h[0].Clone()
    for i in range (0, res.GetNbinsX()+1):
        systotal = 0
        for h in list_h:
            systotal += h.GetBinContent(i)*h.GetBinContent(i)
        res.SetBinContent(i, math.sqrt(systotal))
    return res



def plotter(predNominal, predPullD, predPullU, legNom, leg1, leg2, outDir, outTitle, sampleName):
    
    c1=TCanvas("c1","c1",800,600)
    t1=TPad("t1","t1", 0.0, 0.25, 0.95, 0.95)

    t1.Draw()
    t1.SetLogy(1)
    #t1.SetGrid()
    t1.SetTopMargin(0.05)
    t1.SetBottomMargin(0.12)
    c1.cd()

    # t2=TPad("t2","t2", 0.0, 0.225, 0.95, 0.375)
    # t2.Draw()
    # t2.SetGridy(1)
    # t2.SetTopMargin(0.05)
    # t2.SetBottomMargin(0.1)

    t3=TPad("t3","t3", 0.0, 0., 0.95, 0.25)
    t3.Draw()
    t3.SetGridy(1)
    t3.SetTopMargin(0.05)
    t3.SetBottomMargin(0.4)

    t1.cd()
    c1.SetLogy(1)

    max_mass = 4000
    min_entries=2e-1
    max_entries=predNominal.GetMaximum()*2

    predNominal = setColorAndMarker(predNominal,1,20)
    predNominal.GetXaxis().SetRangeUser(0,max_mass)
    predNominal.GetYaxis().SetRangeUser(min_entries, max_entries)
    predNominal.SetMinimum(min_entries)
    #predNominal.SetMaximum(max_entries)
    predNominal.SetTitle(";Mass (GeV);Tracks")
    predNominal.GetYaxis().SetTitleSize(0.07)
    predNominal.GetYaxis().SetLabelSize(0.06)
    predNominal.GetYaxis().SetTitleOffset(0.7)
    predNominal.GetXaxis().SetTitle("Mass (GeV)")
    predNominal.GetXaxis().SetTitleOffset(0.9)
    predNominal.GetXaxis().SetLabelSize(0.06)
    predNominal.GetXaxis().SetTitleSize(0.06)
    predNominal.Draw()

    predPullD = setColorAndMarker(predPullD,ROOT.kBlue+1,23)
    predPullD.Draw("same")

    predPullU = setColorAndMarker(predPullU,ROOT.kRed+1,22)
    predPullU.Draw("same")

    
    leg = TLegend(0.2,0.65,0.35,0.92)
    leg.AddEntry(predNominal,legNom,"PE1")
    leg.AddEntry(predPullD,leg1,"PE1")
    leg.AddEntry(predPullU,leg2,"PE1")
    leg.SetBorderSize(0)
    leg.SetFillColor(0)

    
    LineLastBin=TLine(predNominal.GetBinLowEdge(predNominal.FindBin(max_mass)-1),0,predNominal.GetBinLowEdge(predNominal.FindBin(max_mass)-1),max_entries)
    LineLastBin.SetLineStyle(1)
    LineLastBin.SetLineColor(1)
    #LineLastBin.Draw("same")
   
    leg.Draw("same")
    
    c1.cd()
    #t2.cd()
    
    frameR2=ROOT.TH1D("frameR2", "frameR2", 1,0, max_mass)
    frameR2.GetXaxis().SetNdivisions(505)
    frameR2.SetTitle("")
    frameR2.SetStats(0)
    frameR2.GetXaxis().SetTitle("")
    frameR2.GetYaxis().SetTitle("RatioR ")
    frameR2.GetXaxis().SetRangeUser(0,max_mass)
    frameR2.SetMaximum(1.1)
    frameR2.SetMinimum(0.9)
    if(outTitle=="Pile-Up"):
        frameR2.SetMaximum(1.2)
        frameR2.SetMinimum(0.8)
    frameR2.GetYaxis().SetLabelFont(43) #give the font size in pixel (instead of fraction)
    frameR2.GetYaxis().SetLabelSize(23) #font size
    frameR2.GetYaxis().SetTitleFont(43) #give the font size in pixel (instead of fraction)
    frameR2.GetYaxis().SetTitleSize(24) #font size
    frameR2.GetYaxis().SetNdivisions(205)
    frameR2.GetYaxis().SetTitleOffset(2)
    frameR2.GetXaxis().SetNdivisions(510)
    frameR2.GetXaxis().SetLabelFont(43) #give the font size in pixel (instead of fraction)
    frameR2.GetXaxis().SetLabelSize(0) #font size
    frameR2.GetXaxis().SetTitleFont(43) #give the font size in pixel (instead of fraction)
    frameR2.GetXaxis().SetTitleSize(24) #font size
    frameR2.GetXaxis().SetTitleOffset(0)
    #frameR2.Draw("AXIS")

    LineAtOne=TLine(0,1,max_mass,1)
    LineAtOne.SetLineStyle(3)
    LineAtOne.SetLineColor(1)
    #LineAtOne.Draw("same")

    ratioInt1=ratioInt(predNominal,predPullD)
    ratioInt2=ratioInt(predNominal,predPullU)
    ratioInt1.Draw("E0 same")
    ratioInt2.Draw("E0 same")
    
    ratioInt1.Write()
    ratioInt2.Write()

    c1.cd()
    t3.cd()

    frameR3 = frameR2.Clone()
    frameR3.GetXaxis().SetLabelSize(23) #font size
    frameR3.GetYaxis().SetTitleOffset(1.8)
    frameR3.GetYaxis().SetTitle("#frac{Nominal}{var}")
    frameR3.GetXaxis().SetRangeUser(0,max_mass)
    frameR3.Draw("AXIS")
    frameR3.GetXaxis().SetTitle("Mass (GeV)")
    frameR3.GetXaxis().SetTitleOffset(1)

    LineAtOne.Draw("same")
   
    predNominalCl1=predNominal.Clone()
    predNominalCl2=predNominal.Clone()
    predNominalCl1.Divide(predPullD)
    predNominalCl2.Divide(predPullU)

    predNominalCl1 = setColorAndMarker(predNominalCl1,ROOT.kBlue+1,23)
    predNominalCl1.Draw("E0 same")

    predNominalCl2 = setColorAndMarker(predNominalCl2,ROOT.kRed+1,23)
    predNominalCl2.Draw("E0 same")

    c1.cd()
    latex = TLatex(0.15, 0.92, "#scale[1.3]{#it{Private Work (CMS simulation)}}")
    latex.SetNDC()
    latex.SetTextFont(42)
    latex.SetTextSize(0.03)
    latex.Draw()

    mass = sampleName.split("_")[1]
    tex = TLatex(0.69, 0.92, "#scale[1.3]{#bf{m_{#tilde{g}}=" + mass + " GeV}}")
    tex.SetNDC()
    tex.SetTextFont(42)
    tex.SetTextSize(0.04)
    tex.Draw()

    tex2 = TLatex(0.71, 0.36, "#scale[1.3]{#bf{" + outTitle + "}}")
    tex2.SetNDC()
    tex2.SetTextFont(42)
    tex2.SetTextSize(0.04)
    tex2.Draw()

    c1.SaveAs(outDir+"/"+outTitle+".pdf")

def plotSummary(syst_K, syst_C, syst_PU, syst_Fpix, syst_lumi,
                syst_trigger, sys_Jet, sysTot, xtitle, outTitle, outDir, eta):

    syst_K = lowEdge(syst_K)
    syst_C = lowEdge(syst_C)
    syst_PU = lowEdge(syst_PU)
    syst_Fpix = lowEdge(syst_Fpix)
    syst_trigger = lowEdge(syst_trigger)
    sys_Jet = lowEdge(sys_Jet)
    syst_lumi = lowEdge(syst_lumi)

    sysTot = lowEdge(sysTot)

    c2=TCanvas("c2","c2",800,600)
    c2.SetBottomMargin(0.12)
    c2.SetLeftMargin(0.1)

    c2.SetLogy()
    c2.SetGrid()
    syst_K.SetMinimum(1)
    syst_K.SetMaximum(500)
    syst_K.GetXaxis().SetTitle(xtitle)
    syst_K.GetYaxis().SetTitle("Systematic Uncertainty [%]")
    syst_K.GetXaxis().SetNdivisions(510)
    syst_K.GetXaxis().SetLabelFont(43) #give the font size in pixel (instead of fraction)
    syst_K.GetXaxis().SetLabelSize(22) #font size
    syst_K.GetXaxis().SetTitleSize(0.05) #font size
    syst_K.GetYaxis().SetLabelFont(43) #give the font size in pixel (instead of fraction)
    syst_K.GetYaxis().SetLabelSize(22) #font size
    syst_K.GetYaxis().SetTitleSize(0.05) #font size
    syst_K.GetYaxis().SetTitleOffset(1.0)
    syst_K.GetXaxis().SetTitleOffset(1.0)

    # draw latex CMS preliminary
    #latex = TLatex(0.10, 0.96, "#scale[1.3]{#bf{CMS}}#it{Simulation Work in progress}")
    latex = TLatex(0.10, 0.96, "#scale[1.3]{#it{Private work (CMS simulation)}}")
    latex.SetNDC()
    latex.SetTextFont(42)
    latex.SetTextSize(0.04)

    latex3 = TLatex(0.8, 0.88, "#bf{|#eta|<1}")
    if (etaS=='Eta1_2p4'): latex3 = TLatex(0.75, 0.88, "#bf{1#leq|#eta|<2.4}")
    if (etaS=='Eta2p4'): latex3 = TLatex(0.8, 0.88, "#bf{|#eta|<2.4}")
    latex3.SetNDC()
    latex3.SetTextFont(42)
    latex3.SetTextSize(0.07)

    mass = outTitle.split("_")[1]
    tex = TLatex(0.73, 0.96, "#scale[1.3]{#bf{m_{#tilde{g}}=" + mass + " GeV}}")
    tex.SetNDC()
    tex.SetTextFont(42)
    tex.SetTextSize(0.04)

    syst_K.GetXaxis().SetLimits(0, 4000)
    syst_K.SetMinimum(0.1)
    syst_K.SetMaximum(2000)

    syst_K = setColorAndMarker(syst_K, ROOT.kMagenta-9, 21)
    syst_C = setColorAndMarker(syst_C, ROOT.kViolet+1, 22)
    syst_PU = setColorAndMarker(syst_PU, ROOT.kBlue+1, 23)
    syst_Fpix = setColorAndMarker(syst_Fpix, ROOT.kCyan, 29)
    syst_lumi = setColorAndMarker(syst_lumi, ROOT.kGreen+2, 20)
    syst_trigger = setColorAndMarker(syst_trigger, ROOT.kOrange, 39)
    sys_Jet = setColorAndMarker(sys_Jet, ROOT.kOrange+1, 33)

    sysTot = setColorAndMarker(sysTot, ROOT.kRed, 34)

    leg2 = TLegend(0.22,0.75,0.6,0.93)
    leg2.SetNColumns(2)
    leg2.AddEntry(sysTot,"Total","PE1")
    leg2.AddEntry(syst_K,"K","PE1")
    leg2.AddEntry(syst_C,"C","PE1")
    leg2.AddEntry(syst_PU,"PU","PE1")
    leg2.AddEntry(syst_Fpix,"F^{pixel}","PE1")
    leg2.AddEntry(syst_lumi, "Lumi", "PE1")
    leg2.AddEntry(syst_trigger,"Trigger","PE1")
    leg2.AddEntry(sys_Jet,"Jet","PE1")

    syst_K.Draw("AP")
    leg2.Draw("same")
    syst_C.Draw("P")
    syst_PU.Draw("P")
    syst_Fpix.Draw("P")
    syst_trigger.Draw("P")
    sys_Jet.Draw("P")
    syst_lumi.Draw("P")
    sysTot.Draw("P")
    

    latex.Draw()
    latex3.Draw()
    tex.Draw()
    c2.Update()

    c2.SaveAs(outDir + "summary_" + outTitle + "_" + eta + ".pdf")

def filledRange(h):
    """Restreint l'axe aux bins non vides. Renvoie False si l'histo est vide."""
    bins = [i for i in range(1, h.GetNbinsX()+1) if h.GetBinContent(i) > 0]
    if not bins:
        return False
    h.GetXaxis().SetRange(bins[0], bins[-1])
    for i in range(0, h.GetNbinsX()+2):
        h.SetBinError(i, 0)
    return True

def plotTotalAllMasses(dict_sysTot, xtitle, outDir, eta, outTitle="sysTot_allMasses"):

    c3 = TCanvas("c3", "c3", 800, 600)
    c3.SetBottomMargin(0.1)
    c3.SetLeftMargin(0.1)
    c3.SetLogy()
    c3.SetGrid()

    colorMap = {
        "2000": ROOT.kOrange + 8,
        "2400": ROOT.kViolet + 1,
        "2600": ROOT.kGreen - 3,
    }
    markers = [23, 21, 22]

    frame = ROOT.TH1D("frameTot", "", 1, 0, 4000)
    frame.SetStats(0)
    frame.SetMinimum(1)
    frame.SetMaximum(2000)
    frame.GetXaxis().SetTitle(xtitle)
    frame.GetYaxis().SetTitle("Total Systematic Uncertainty [%]")
    frame.GetXaxis().SetNdivisions(510)
    frame.GetXaxis().SetLabelFont(43)
    frame.GetXaxis().SetLabelSize(22)
    frame.GetXaxis().SetTitleSize(0.05)
    frame.GetXaxis().SetTitleOffset(1.0)
    frame.GetYaxis().SetLabelFont(43)
    frame.GetYaxis().SetLabelSize(22)
    frame.GetYaxis().SetTitleSize(0.05)
    frame.GetYaxis().SetTitleOffset(1.0)
    frame.Draw("AXIS")

    leg3 = TLegend(0.3, 0.68, 0.63, 0.93)
    leg3.SetBorderSize(0)

    for i, (mass, h) in enumerate(sorted(dict_sysTot.items(), key=lambda x: int(x[0]))):
        if mass not in colorMap:
            continue
        if not filledRange(h):
            continue
        setColorAndMarker(h, colorMap[mass], markers[i % len(markers)])
        h.SetFillStyle(0)
        h.SetLineWidth(2)
        leg3.AddEntry(h, "m_{#tilde{g}}=" + mass + " GeV", "LP")
        h.Draw("HIST SAME")
        h.Draw("P SAME")

    leg3.Draw("same")

    latex = TLatex(0.10, 0.96, "#scale[1.3]{#it{Private work (CMS simulation)}}")
    latex.SetNDC()
    latex.SetTextFont(42)
    latex.SetTextSize(0.04)
    latex.Draw()

    latex3 = TLatex(0.8, 0.88, "#bf{|#eta|<1}")
    if (eta=='Eta1_2p4'): latex3 = TLatex(0.75, 0.88, "#bf{1#leq|#eta|<2.4}")
    if (eta=='Eta2p4'): latex3 = TLatex(0.8, 0.88, "#bf{|#eta|<2.4}")
    latex3.SetNDC()
    latex3.SetTextFont(42)
    latex3.SetTextSize(0.07)
    latex3.Draw()

    c3.Update()
    c3.SaveAs(outDir + outTitle + "_" + eta + ".pdf")


# ---------------- MAIN ----------------

codeVersion = "19p12"

SignalSamples = [
    "Gluino_Run3_MET_madgraph_1100_V" + codeVersion + ".root",
    "Gluino_Run3_MET_madgraph_1200_V" + codeVersion + ".root",
    "Gluino_Run3_MET_madgraph_1300_V" + codeVersion + ".root",
    "Gluino_Run3_MET_madgraph_1400_V" + codeVersion + ".root",
    "Gluino_Run3_MET_madgraph_1600_V" + codeVersion + ".root",
    "Gluino_Run3_MET_madgraph_1800_V" + codeVersion + ".root",
    "Gluino_Run3_MET_madgraph_2000_V" + codeVersion + ".root",
    "Gluino_Run3_MET_madgraph_2200_V" + codeVersion + ".root",
    "Gluino_Run3_MET_madgraph_2400_V" + codeVersion + ".root",
    "Gluino_Run3_MET_madgraph_2600_V" + codeVersion + ".root",
]

sampleName = {
    "Gluino_Run3_MET_madgraph_1100_V" + codeVersion + ".root": "Gluino_1100",
    "Gluino_Run3_MET_madgraph_1200_V" + codeVersion + ".root": "Gluino_1200",
    "Gluino_Run3_MET_madgraph_1300_V" + codeVersion + ".root": "Gluino_1300",
    "Gluino_Run3_MET_madgraph_1400_V" + codeVersion + ".root": "Gluino_1400",
    "Gluino_Run3_MET_madgraph_1600_V" + codeVersion + ".root": "Gluino_1600",
    "Gluino_Run3_MET_madgraph_1800_V" + codeVersion + ".root": "Gluino_1800",
    "Gluino_Run3_MET_madgraph_2000_V" + codeVersion + ".root": "Gluino_2000",
    "Gluino_Run3_MET_madgraph_2200_V" + codeVersion + ".root": "Gluino_2200",
    "Gluino_Run3_MET_madgraph_2400_V" + codeVersion + ".root": "Gluino_2400",
    "Gluino_Run3_MET_madgraph_2600_V" + codeVersion + ".root": "Gluino_2600",
}


idir   = "/safe/ui3_1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/output/Gluino_V19/"
RegS   = "9fp10"
etaS   = "Eta2p4"
option = ''
hname  = "METanalysis_TestPUppiMETCut_" + option + etaS + "_" + RegS + "_SignalMass"

rebinning = array.array('d',[0.,20.,40.,60.,80.,100.,120.,140.,160.,180.,200.,220.,240.,260.,280.,300.,320.,340.,360.,380.,410.,440.,480.,530.,590.,660.,760.,880.,1030.,1210.,1440.,1730.,2000.,2500.,3200.,4000.])
sizeRebinning = len(rebinning) - 1

sysTot_allMasses = {}

for sample in SignalSamples:

    print(sample)
    
    ifile = ROOT.TFile(idir + sample)
    odir = "systSignal/" + sampleName.get(sample) + "_" + RegS + "_" + option + etaS + "/"
    os.system("mkdir -p " + odir)
    ofile = TFile(odir + "sysTotBinned_signal.root", "RECREATE")

    mass_nominal = ifile.Get(hname + "_nominal")
    mass_PU_up = ifile.Get(hname + "_PUUp")
    mass_PU_down = ifile.Get(hname + "_PUDown")
    mass_Trigger_up = ifile.Get(hname + "_TriggerSFUp")
    mass_Trigger_down = ifile.Get(hname + "_TriggerSFDown")
    mass_Jet_up = ifile.Get(hname + "_JetUp")
    mass_Jet_down = ifile.Get(hname + "_JetDown")
    mass_K_up = ifile.Get(hname + "_KUp")
    mass_K_down = ifile.Get(hname + "_KDown")
    mass_C_up = ifile.Get(hname + "_CUp")
    mass_C_down = ifile.Get(hname + "_CDown")


    mass_nominal = allSet(mass_nominal, sizeRebinning, rebinning, "nominal")
    mass_PU_up = allSet(mass_PU_up, sizeRebinning, rebinning, "PU_up")
    mass_PU_down = allSet(mass_PU_down, sizeRebinning, rebinning, "PU_down")
    mass_Trigger_up = allSet(mass_Trigger_up, sizeRebinning, rebinning, "Trigger_up")
    mass_Trigger_down = allSet(mass_Trigger_down, sizeRebinning, rebinning, "Trigger_down")
    mass_Jet_up = allSet(mass_Jet_up, sizeRebinning, rebinning, "Jet_up")
    mass_Jet_down = allSet(mass_Jet_down, sizeRebinning, rebinning, "Jet_down")
    mass_K_up = allSet(mass_K_up, sizeRebinning, rebinning, "K_up")
    mass_K_down = allSet(mass_K_down, sizeRebinning, rebinning, "K_down")
    mass_C_up = allSet(mass_C_up, sizeRebinning, rebinning, "C_up")
    mass_C_down = allSet(mass_C_down, sizeRebinning, rebinning, "C_down")
    

    # Here: systematics computed as the maximum between abs(nom-up) and abs(nom-down) in each bin
    syst_PU_binned = systMass(mass_nominal, mass_PU_down, mass_PU_up, "", 0, 1)
    syst_Trigger_binned = systMass(mass_nominal, mass_Trigger_down, mass_Trigger_up, "", 0, 1)
    syst_Jet_binned = systMass(mass_nominal, mass_Jet_down, mass_Jet_up, "", 0, 1)
    syst_K_binned = systMass(mass_nominal, mass_K_down, mass_K_up, "", 0, 1)
    syst_C_binned = systMass(mass_nominal, mass_C_down, mass_C_up, "", 0, 1)
    syst_lumi_binned = systFlat(mass_nominal, HARDCODED_SYST["lumi"], "syst_lumi")
    syst_Fpix_binned = systFlat(mass_nominal, HARDCODED_SYST["Fpix"], "syst_Fpix")
    
    
    listOfSyst = [syst_PU_binned, syst_Trigger_binned, syst_Jet_binned,
                  syst_K_binned, syst_C_binned,
                  syst_lumi_binned, syst_Fpix_binned]
    sysTot_binned = systTotal(listOfSyst)
    ofile.cd()
    sysTot_binned.Write()

    mass_point = sampleName.get(sample).split("_")[1]
    if mass_point in ["2000", "2400", "2600"]:
        h_tot = sysTot_binned.Clone("sysTot_" + mass_point)
        h_tot.SetDirectory(0)
        sysTot_allMasses[mass_point] = h_tot

    plotter(mass_nominal, mass_PU_down, mass_PU_up, "Nominal", "Down", "Up", odir, "Pile-Up",sampleName.get(sample))
    plotter(mass_nominal, mass_Trigger_down, mass_Trigger_up, "Nominal", "Down", "Up", odir, "Trigger SF",sampleName.get(sample))
    plotter(mass_nominal, mass_Jet_down, mass_Jet_up, "Nominal", "Down", "Up", odir, "JES",sampleName.get(sample))
    plotter(mass_nominal, mass_K_down, mass_K_up, "Nominal", "Down", "Up", odir, "K_{mass}",sampleName.get(sample))
    plotter(mass_nominal, mass_C_down, mass_C_up, "Nominal", "Down", "Up", odir, "C_{mass}",sampleName.get(sample))
    
    plotSummary(syst_K_binned, syst_C_binned, syst_PU_binned,
                syst_Fpix_binned, syst_lumi_binned,
                syst_Trigger_binned, syst_Jet_binned, sysTot_binned,
                "Mass bin", sampleName.get(sample), odir, etaS)
    
odirTot = "systSignal/"
plotTotalAllMasses(sysTot_allMasses, "Mass (GeV)", odirTot, etaS)