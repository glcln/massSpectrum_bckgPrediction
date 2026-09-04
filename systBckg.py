#!/usr/bin/python

import sys, os
sys.argv.append( '-b-' )
import ROOT
import math
import array
import ctypes

from ROOT import TFile, THStack, TCanvas, TLegend, TLatex, TPad, TH1, TH2, TLine, TGraph, TGraphErrors
import tdrstyle

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning

tdrstyle.setTDRStyle()

def ratioHisto(h2, h1):
    hc = h1.Clone()
    hc.Divide(h2)
    return hc

def ratioInt(h1, h2):
    h3 = h2.Clone()
    h3.Reset()
    for i in range(0, h1.GetNbinsX()+1):
        e1 = ctypes.c_double(0.0)
        e2 = ctypes.c_double(0.0)
        a = h1.IntegralAndError(i, h1.GetNbinsX()+1, e1, "")
        b = h2.IntegralAndError(i, h1.GetNbinsX()+1, e2, "")
        if b != 0 and a != 0:
            c = math.sqrt((e1.value*e1.value)/(a*a) + (e2.value*e2.value)/(b*b)) * a/b
            h3.SetBinContent(i, a/b)
            h3.SetBinError(i, c)
        else:
            h3.SetBinContent(i, 0)
    return h3

def lowEdge(h1):
    res=ROOT.TGraph(h1.GetNbinsX()-1)
    for i in range (1,h1.GetNbinsX()+1):
        res.SetPoint(i-1,h1.GetBinLowEdge(i),h1.GetBinContent(i))
    return res

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

def statErr(h1, name):
    statErr = h1.Clone()
    statErr.SetName(name)
    for i in range (1, statErr.GetNbinsX()):
        if statErr.GetBinContent(i)>0:
            statErr.SetBinContent(i, statErr.GetBinError(i)/statErr.GetBinContent(i))
        else:
            statErr.SetBinContent(i,0)
    return 100*statErr

def statErrRInt(h1, name):
    statErr = h1.Clone()
    statErr.SetName(name)
    e = ctypes.c_double(0.0)
    for i in range(1, statErr.GetNbinsX()):
        c = h1.IntegralAndError(i, h1.GetNbinsX(), e, "")
        if e.value > 0:
            statErr.SetBinContent(i, e.value / c)
        else:
            statErr.SetBinContent(i, 0)
    return 100 * statErr

def systMass(nominal, down, up, name, typec, binned, mini=0, debug=False):
    if (binned==0):
        ra1 = ratioInt(nominal, up)
        ra2 = ratioInt(nominal, down)
    elif (binned==1):
        ra1 = ratioHisto(nominal, up)
        ra2 = ratioHisto(nominal, down)
    res = ra1.Clone()
    res.SetName(name)
    for i in range (0,res.GetNbinsX()+1):
        r1 = ra1.GetBinContent(i)
        r2 = ra2.GetBinContent(i)

        s1 = abs(1-r1)
        s2 = abs(1-r2)
        
        if (typec==0):
            if (debug):
                if (s1 > s2): print(f"bin {i}: up")
                else: print(f"bin {i}: down")
            
            m = max(s1, s2)
            if (mini==1):
                m = min(s1, s2)
            elif (mini==2):
                m /= 2
        elif (typec==1):
            m = (s1 + s2)/2.
        res.SetBinContent(i, 100*m)
    return res

def systMassAll(nominal, down, up, st, mini=0, debug=False):
    '''
    Return the maximum between the up and down variation relative to the nominal, in percentage.
    If mini=1, return the minimum instead of the maximum.
    '''
    if (debug): print ("res")
    res = systMass(nominal, down, up, st, 0, 0, mini, debug)
    if (debug): print ("res_mean")
    res_mean = systMass(nominal, down, up, st, 1, 0, mini, debug)
    if (debug): print ("res_binned")
    res_binned = systMass(nominal, down, up, st, 0, 1, mini, debug)
    #res_binned = systMass(nominal, down, down, st, 0, 1, mini, debug)
    if (debug): print ("res_binned_mean")
    res_binned_mean = systMass(nominal, down, up, st, 1, 1, mini, debug)
    return res, res_mean, res_binned, res_binned_mean

def systTotal(list_h):
    res = list_h[0].Clone()
    for i in range (0,res.GetNbinsX()+1):
        systotal = 0
        for h in list_h:
            systotal += h.GetBinContent(i)*h.GetBinContent(i)
        res.SetBinContent(i,math.sqrt(systotal))
    return res

def allSet(h, sizeRebinning,  rebinning , st):
    h = h.Rebin(sizeRebinning, st, rebinning)
    h = overflowInLastBin(h)
    h = underflowInFirstBin(h)
    h.Scale(1./h.Integral())
    return h

def setColorAndMarker(h1,color,markerstyle):
    h1.SetLineColor(color)
    h1.SetMarkerColor(color)
    h1.SetFillColor(color)
    h1.SetMarkerStyle(markerstyle)
    h1.SetMarkerSize(1.2)
    return h1

def lowEdge(h1):
    res=ROOT.TGraph(h1.GetNbinsX()-1)
    for i in range (1,h1.GetNbinsX()+1):
        res.SetPoint(i-1,h1.GetBinLowEdge(i),h1.GetBinContent(i))
    return res

def plotter(predNominal, predPullD, predPullU, legNom, leg1, leg2, outDir, outTitle):
 
    c1 = TCanvas("c1_" + outTitle, "c1", 800, 800)
    t1 = TPad("t1", "t1", 0.0, 0.40, 0.95, 0.95)
 
    t1.Draw()
    t1.SetLogy(1)
    t1.SetGrid(1)
    t1.SetTopMargin(0.005)
    t1.SetBottomMargin(0.005)
    c1.cd()
 
    t2 = TPad("t2", "t2", 0.0, 0.225, 0.95, 0.375)
    t2.Draw()
    t2.SetGridy(1)
    t2.SetTopMargin(0.05)
    t2.SetBottomMargin(0.1)
 
    t3 = TPad("t3", "t3", 0.0, 0., 0.95, 0.20)
    t3.Draw()
    t3.SetGridy(1)
    t3.SetBottomMargin(0.45)
 
    t1.cd()
 
    max_mass = 4000
    min_entries = 1e-7
    max_entries = predNominal.GetMaximum() * 5
 
    predNominal = setColorAndMarker(predNominal, 1, 20)
    predNominal.GetXaxis().SetRangeUser(0, max_mass)
    predNominal.GetYaxis().SetRangeUser(min_entries, max_entries)
    predNominal.SetMinimum(min_entries)
    predNominal.SetTitle(";Mass (GeV);Normalized tracks")
    predNominal.GetYaxis().SetTitleSize(0.07)
    predNominal.GetYaxis().SetLabelSize(0.06)
    predNominal.GetYaxis().SetTitleOffset(0.9)
    predNominal.Draw()
 
    predPullD = setColorAndMarker(predPullD, 38, 21)
    predPullD.Draw("same")
 
    if predPullU is not None:
        predPullU = setColorAndMarker(predPullU, 46, 21)
        predPullU.Draw("same")
 
    leg = TLegend(0.16, 0.75, 0.35, 0.99)
    leg.AddEntry(predNominal, legNom, "PE1")
    leg.AddEntry(predPullD, leg1, "PE1")
    if predPullU is not None:
        leg.AddEntry(predPullU, leg2, "PE1")
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.Draw("same")
 
    c1.cd()
    t2.cd()
 
    frameR2 = ROOT.TH1D("frameR2_" + outTitle, "frameR2", 1, 0, max_mass)
    frameR2.GetXaxis().SetNdivisions(505)
    frameR2.SetTitle("")
    frameR2.SetStats(0)
    frameR2.GetXaxis().SetTitle("")
    frameR2.GetYaxis().SetTitle("RatioR ")
    frameR2.GetXaxis().SetRangeUser(0, max_mass)
    frameR2.SetMaximum(1.1)
    frameR2.SetMinimum(0.9)
    frameR2.GetYaxis().SetLabelFont(43)
    frameR2.GetYaxis().SetLabelSize(20)
    frameR2.GetYaxis().SetTitleFont(43)
    frameR2.GetYaxis().SetTitleSize(20)
    frameR2.GetYaxis().SetNdivisions(205)
    frameR2.GetYaxis().SetTitleOffset(2)
    frameR2.GetXaxis().SetNdivisions(510)
    frameR2.GetXaxis().SetLabelFont(43)
    frameR2.GetXaxis().SetLabelSize(0)
    frameR2.GetXaxis().SetTitleFont(43)
    frameR2.GetXaxis().SetTitleSize(24)
    frameR2.GetXaxis().SetTitleOffset(0)
    frameR2.Draw("AXIS")
 
    LineAtOne = TLine(0, 1, max_mass, 1)
    LineAtOne.SetLineStyle(3)
    LineAtOne.SetLineColor(1)
    LineAtOne.Draw("same")
 
    ratioInt1 = ratioInt(predNominal, predPullD)
    ratioInt1 = setColorAndMarker(ratioInt1, 38, 21)
    ratioInt1.Draw("E0 same")
    ratioInt2 = None
    if predPullU is not None:
        ratioInt2 = ratioInt(predNominal, predPullU)
        ratioInt2 = setColorAndMarker(ratioInt2, 46, 21)
        ratioInt2.Draw("E0 same")
 
    c1.cd()
    t3.cd()
 
    frameR3 = frameR2.Clone("frameR3_" + outTitle)
    frameR3.GetXaxis().SetLabelSize(20)
    frameR3.GetYaxis().SetTitleOffset(2.1)
    frameR3.GetYaxis().SetTitle("#frac{Nominal}{var}")
    frameR3.GetXaxis().SetRangeUser(0, max_mass)
    frameR3.Draw("AXIS")
    frameR3.GetXaxis().SetTitle("Mass (GeV)")
    frameR3.GetXaxis().SetTitleOffset(1)
 
    LineAtOne.Draw("same")
 
    predNominalCl1 = predNominal.Clone()
    predNominalCl1.Divide(predPullD)
    predNominalCl1 = setColorAndMarker(predNominalCl1, 38, 21)
    predNominalCl1.Draw("E0 same")
 
    predNominalCl2 = None
    if predPullU is not None:
        predNominalCl2 = predNominal.Clone()
        predNominalCl2.Divide(predPullU)
        predNominalCl2 = setColorAndMarker(predNominalCl2, 46, 21)
        predNominalCl2.Draw("E0 same")
 
    c1.cd()
    latex = TLatex(0.15, 0.955, "#scale[1.3]{#it{Private work (CMS data)}}")
    latex.SetNDC()
    latex.SetTextFont(42)
    latex.SetTextSize(0.03)
    latex.Draw()
 
    latex2 = TLatex(0.70, 0.955, "109 fb^{-1} (13.6 TeV)")
    latex2.SetNDC()
    latex2.SetTextFont(42)
    latex2.SetTextSize(0.03)
    latex2.Draw()
 
    c1.SaveAs(outDir + "/" + outTitle + ".pdf")

def plotSummary(syst_stat,
                syst_eta,
                syst_ih,
                syst_p,
                syst_fitDeDx,
                syst_fitP,
                syst_nofit,
                syst_corrIh,
                syst_corr1oP,
                sysTot,
                xtitle,
                outTitle,
                etaname):
    
    syst_stat = lowEdge(syst_stat)
    syst_eta = lowEdge(syst_eta)
    syst_ih = lowEdge(syst_ih)
    syst_p = lowEdge(syst_p)
    syst_fitP = lowEdge(syst_fitP)
    syst_nofit = lowEdge(syst_nofit)
    syst_corrIh = lowEdge(syst_corrIh)
    syst_corr1oP = lowEdge(syst_corr1oP)
    syst_fitDeDx = lowEdge(syst_fitDeDx)
    sysTot = lowEdge(sysTot)

    c2=TCanvas("c2","c2",800,600)
    c2.SetBottomMargin(0.12)
    c2.SetLeftMargin(0.1)

    c2.SetLogy()
    c2.SetGrid()
    syst_stat.SetMinimum(1)
    syst_stat.SetMaximum(500)
    syst_stat.GetXaxis().SetTitle(xtitle)
    syst_stat.GetYaxis().SetTitle("Systematic Uncertainty [%]")
    syst_stat.GetXaxis().SetNdivisions(510)
    syst_stat.GetXaxis().SetLabelFont(43) #give the font size in pixel (instead of fraction)
    syst_stat.GetXaxis().SetLabelSize(24) #font size
    syst_stat.GetXaxis().SetTitleSize(0.05) #font size
    syst_stat.GetYaxis().SetLabelFont(43) #give the font size in pixel (instead of fraction)
    syst_stat.GetYaxis().SetLabelSize(24) #font size
    syst_stat.GetYaxis().SetTitleSize(0.05) #font size
    syst_stat.GetYaxis().SetTitleOffset(1.0)
    syst_stat.GetXaxis().SetTitleOffset(1.0)


    # draw latex CMS preliminary
    #latex = TLatex(0.1, 0.96, "#scale[1.3]{#bf{CMS}}#it{Work in progress}")
    latex = TLatex(0.1, 0.96, "#scale[1.3]{#it{Private work (CMS data)}}")
    latex.SetNDC()
    latex.SetTextFont(42)
    latex.SetTextSize(0.04)

    latex2 = TLatex(0.735, 0.96, "109 fb^{-1} (13.6 TeV)")
    latex2.SetNDC()
    latex2.SetTextFont(42)
    latex2.SetTextSize(0.04)

    latex3 = TLatex(0.8, 0.88, "#bf{|#eta|<1}")
    if (etaname=='Eta1_2p4'): latex3 = TLatex(0.75, 0.88, "#bf{1#leq|#eta|<2.4}")
    if (etaname=='Eta2p4'): latex3 = TLatex(0.8, 0.88, "#bf{|#eta|<2.4}")
    latex3.SetNDC()
    latex3.SetTextFont(42)
    latex3.SetTextSize(0.07)

    syst_stat.GetXaxis().SetRangeUser(0,4000)
    syst_stat.SetMinimum(0.1)
    syst_stat.SetMaximum(2000)


    syst_stat = setColorAndMarker(syst_stat,1,20)
    syst_eta = setColorAndMarker(syst_eta, ROOT.kMagenta-9, 21)
    syst_ih = setColorAndMarker(syst_ih, ROOT.kViolet+1, 22)
    syst_p = setColorAndMarker(syst_p, ROOT.kBlue+1, 23)
    syst_fitP = setColorAndMarker(syst_fitP, ROOT.kCyan, 29)
    syst_nofit = setColorAndMarker(syst_nofit, ROOT.kGreen+2, 20)
    syst_corrIh = setColorAndMarker(syst_corrIh, ROOT.kOrange, 39)
    syst_corr1oP = setColorAndMarker(syst_corr1oP, ROOT.kGreen+2, 30)
    syst_fitDeDx = setColorAndMarker(syst_fitDeDx, ROOT.kOrange+1, 33)

    sysTot = setColorAndMarker(sysTot, ROOT.kRed, 34)

    leg2 = TLegend(0.12,0.7,0.5,0.93)
    leg2.SetNColumns(2)
    leg2.AddEntry(sysTot,"Total","PE1")
    leg2.AddEntry(syst_stat,"Stat.","PE1")
    leg2.AddEntry(syst_eta,"#eta binning","PE1")
    leg2.AddEntry(syst_ih,"I_{h} binning","PE1")
    leg2.AddEntry(syst_p,"p binning","PE1")
    leg2.AddEntry(syst_fitP,"p fit","PE1")
    leg2.AddEntry(syst_fitDeDx,"I_{h} fit","PE1")
    #leg2.AddEntry(syst_nofit,"No fit","PE1")
    leg2.AddEntry(syst_corrIh, "corr template I_{h}", "PE1")
    leg2.AddEntry(syst_corr1oP, "corr template 1/p", "PE1")
   
    syst_stat.Draw("AP")
    leg2.Draw("same")
    syst_eta.Draw("P")
    syst_ih.Draw("P")
    syst_p.Draw("P")
    syst_fitP.Draw("P")
    #syst_nofit.Draw("P")
    syst_corrIh.Draw("P")
    syst_corr1oP.Draw("P")
    syst_fitDeDx.Draw("P")
    sysTot.Draw("P")
    
    latex.Draw()
    latex2.Draw()
    latex3.Draw()
    
    commandMkdir='mkdir -p '+oDir+'pdf '+oDir+'Cfile '+oDir+'rootfile'
    os.system(commandMkdir)

    c2.SaveAs(oDir+"pdf/summary_" + outTitle + "_" + sample + "_" + etaname + ".pdf")
    c2.SaveAs(oDir+"Cfile/summary_"+outTitle+".root")
    c2.SaveAs(oDir+"rootfile/summary_"+outTitle+".C")




# Setup
onlyNominal = False    
version     = "V12p35"
etaname     = "Eta2p4"
year        = "2024"
sample      = "JetMET"
region      = "9fp10"
option      = "_CorrelationAdded_SigmaPtoverPt_0p5_EoP_0p1"
option2     = "_OldFit_IhC"
option3     = ""
option4     = "_SigmaPtoverPt_0p5_EoP_0p1"

directory   = "/safe/ui3_1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/macros/DataMET_2024_" + version + "__" + region + option + "/" + etaname + "/"
oDir        = directory + option3 + "SystCombined/"
ofile       = TFile(oDir + "sysToTBinned_" + year + "_" + region + ".root", "RECREATE")

plotType    = "mass_predBC_"

rebinning = array.array('d',[0.,20.,40.,60.,80.,100.,120.,140.,160.,180.,200.,220.,240.,260.,280.,300.,320.,340.,360.,380.,
                             410.,440.,480.,530.,590.,660.,760.,880.,1030.,1210.,1440.,1730.,2000.,2500.,3200.,4000.])
sizeRebinning = len(rebinning) - 1


# Only stat for the nominal
inputNominal = directory + sample + year + "_" + version + "_rebinEta4_rebinIh4_rebinP2_EtaReweighting_" + etaname + option4 + option2 + ".root" 
ifileNominal = TFile(inputNominal)

predNominal_def = ifileNominal.Get(plotType + region)
predNominal_def = allSet(predNominal_def, sizeRebinning, rebinning, "nominal_def")

syst_stat = statErrRInt(predNominal_def, "Stat")
syst_stat_binned = statErr(predNominal_def, "Stat_binned")

if (onlyNominal):
    print("Only nominal, no systematic variation, exiting.")
    ofile.cd()
    syst_stat.Write()
    syst_stat_binned.Write()
    ofile.Close()
    sys.exit(0)


# Up/Down for the others
inputEtaD = directory + sample + year + "_" + version + "_rebinEta8_rebinIh4_rebinP2_EtaReweighting_" + etaname + option4 + option2 + ".root" 
inputEtaU = directory + sample + year + "_" + version + "_rebinEta2_rebinIh4_rebinP2_EtaReweighting_" + etaname + option4 + option2 + ".root" 
ifileEtaD = TFile(inputEtaD)
ifileEtaU = TFile(inputEtaU)

inputIhD = directory + sample + year + "_" + version + "_rebinEta4_rebinIh8_rebinP2_EtaReweighting_" + etaname + option4 + option2 + ".root" 
inputIhU = directory + sample + year + "_" + version + "_rebinEta4_rebinIh2_rebinP2_EtaReweighting_" + etaname + option4 + option2 + ".root" 
ifileIhD = TFile(inputIhD)
ifileIhU = TFile(inputIhU)

inputPD = directory + sample + year + "_" + version + "_rebinEta4_rebinIh4_rebinP4_EtaReweighting_" + etaname + option4 + option2 + ".root" 
inputPU = directory + sample + year + "_" + version + "_rebinEta4_rebinIh4_rebinP1_EtaReweighting_" + etaname + option4 + option2 + ".root" 
ifilePD = TFile(inputPD)
ifilePU = TFile(inputPU)

inputFitPUp     = directory + sample + year + "_" + version + "_rebinEta4_rebinIh4_rebinP2_fitPUp_EtaReweighting_" + etaname + option4 + option2 + ".root" 
inputFitPDown   = directory + sample + year + "_" + version + "_rebinEta4_rebinIh4_rebinP2_fitPDown_EtaReweighting_" + etaname + option4 + option2 + ".root" 
ifileFitPUp     = TFile(inputFitPUp)
ifileFitPDown   = TFile(inputFitPDown)

inputFitDeDxUp      = directory + sample + year + "_" + version + "_rebinEta4_rebinIh4_rebinP2_fitIhUp_EtaReweighting_" + etaname + option4 + option2 + ".root" 
inputFitDeDxDown    = directory + sample + year + "_" + version + "_rebinEta4_rebinIh4_rebinP2_fitIhDown_EtaReweighting_" + etaname + option4 + option2 + ".root" 
ifileFitDeDxUp      = TFile(inputFitDeDxUp)
ifileFitDeDxDown    = TFile(inputFitDeDxDown)

inputNoFit = directory + sample + year + "_" + version + "_rebinEta4_rebinIh4_rebinP2_EtaReweighting_" + etaname + option4 + "_NoFit" + "_IhC" + ".root" 
ifileNoFit = TFile(inputNoFit)

inputcorrTemplateIh  = directory + sample + year + "_" + version + "_rebinEta4_rebinIh4_rebinP2_corrTemplateIh_EtaReweighting_" + etaname + option4 + option2 + ".root"
ifilecorrTemplateIh  = TFile(inputcorrTemplateIh)
inputcorrTemplate1oP = directory + sample + year + "_" + version + "_rebinEta4_rebinIh4_rebinP2_corrTemplate1oP_EtaReweighting_" + etaname + option4 + option2 + ".root"
ifilecorrTemplate1oP = TFile(inputcorrTemplate1oP)


print(inputNominal)


# Compute the systematics
predEtaD        = ifileEtaD.Get(plotType + region)
predEtaU        = ifileEtaU.Get(plotType + region)

predIhD         = ifileIhD.Get(plotType + region)
predIhU         = ifileIhU.Get(plotType + region)

predPD          = ifilePD.Get(plotType + region)
predPU          = ifilePU.Get(plotType + region)

predFitPUp      = ifileFitPUp.Get(plotType + region)
predFitPDown    = ifileFitPDown.Get(plotType + region)

predFitDeDxUp   = ifileFitDeDxUp.Get(plotType + region)
predFitDeDxDown = ifileFitDeDxDown.Get(plotType + region)

predNoFit = ifileNoFit.Get(plotType + region)

predcorrIh  = ifilecorrTemplateIh.Get(plotType + region)
predcorr1oP = ifilecorrTemplate1oP.Get(plotType + region)



predEtaNominal = allSet(predNominal_def, sizeRebinning, rebinning, "eta_nominal")
predEtaD       = allSet(predEtaD, sizeRebinning, rebinning, "eta_down")
predEtaU       = allSet(predEtaU, sizeRebinning, rebinning, "eta_up")

predIhNominal  = allSet(predNominal_def, sizeRebinning, rebinning, "ih_nominal")
predIhD        = allSet(predIhD, sizeRebinning, rebinning, "ih_down")
predIhU        = allSet(predIhU, sizeRebinning, rebinning, "ih_up")

predPNominal   = allSet(predNominal_def, sizeRebinning, rebinning, "p_nominal")
predPD         = allSet(predPD, sizeRebinning, rebinning, "p_down")
predPU         = allSet(predPU, sizeRebinning, rebinning, "p_up")

predFitPNominal = allSet(predNominal_def, sizeRebinning, rebinning, "fitP_nominal")
predFitPDown    = allSet(predFitPDown, sizeRebinning, rebinning, "fitP_down")
predFitPUp      = allSet(predFitPUp, sizeRebinning, rebinning, "fitP_up")

predFitDeDxNominal  = allSet(predNominal_def, sizeRebinning, rebinning, "fitDeDx_nominal")
predFitDeDxDown     = allSet(predFitDeDxDown, sizeRebinning, rebinning, "fitDeDx_down")
predFitDeDxUp       = allSet(predFitDeDxUp, sizeRebinning, rebinning, "fitDeDx_up")

predNoFitNominal  = allSet(predNominal_def, sizeRebinning, rebinning, "Nofit_nominal")
predNoFitDown     = allSet(predNoFit, sizeRebinning, rebinning, "Nofit_down")
predNoFitUp       = allSet(predNoFit, sizeRebinning, rebinning, "Nofit_up")

predcorrIhNominal  = allSet(predNominal_def, sizeRebinning, rebinning, "corrIh_nominal")
predcorrIhDown     = allSet(predcorrIh, sizeRebinning, rebinning, "corrIh_down")
predcorrIhUp       = allSet(predcorrIh, sizeRebinning, rebinning, "corrIh_up")

predcorr1oPNominal  = allSet(predNominal_def, sizeRebinning, rebinning, "corr1oP_nominal")
predcorr1oPDown     = allSet(predcorr1oP, sizeRebinning, rebinning, "corr1oP_down")
predcorr1oPUp       = allSet(predcorr1oP, sizeRebinning, rebinning, "corr1oP_up")




# Here: systematics computed as the maximum between abs(nom-up) and abs(nom-down) in each bin
(syst_eta, syst_eta_mean, syst_eta_binned, syst_eta_binned_mean) = systMassAll(predEtaNominal, predEtaD, predEtaU, "Eta")
(syst_ih, syst_ih_mean, syst_ih_binned, syst_ih_binned_mean) = systMassAll(predIhNominal, predIhD, predIhU, "Ih")
(syst_p, syst_p_mean, syst_p_binned, syst_p_binned_mean) = systMassAll(predPNominal, predPD, predPU, "P")

(syst_fitP, syst_fitP_mean, syst_fitP_binned, syst_fitP_binned_mean) = systMassAll(predFitPNominal, predFitPDown, predFitPUp, "Fit_p")
(syst_fitDeDx, syst_fitDeDx_mean, syst_fitDeDx_binned, syst_fitDeDx_binned_mean) = systMassAll(predFitDeDxNominal, predFitDeDxDown, predFitDeDxUp, "Fit_dedx_systMassAl")
(syst_nofit, syst_nofit_mean, syst_nofit_binned, syst_nofit_binned_mean) = systMassAll(predNoFitNominal, predNoFitDown, predNoFitUp, "NoFit_systMassAl")
(syst_corrIh, syst_corrIh_mean, syst_corrIh_binned, syst_corrIh_binned_mean) = systMassAll(predcorrIhNominal, predcorrIhDown, predcorrIhUp, "corrIh_systMassAl")
(syst_corr1oP, syst_corr1oP_mean, syst_corr1oP_binned, syst_corr1oP_binned_mean) = systMassAll(predcorr1oPNominal, predcorr1oPDown, predcorr1oPUp, "corr1oP_systMassAl")


odirPlots = oDir + "individualSyst"
os.system("mkdir -p " + odirPlots)
 
plotter(predEtaNominal,     predEtaD,        predEtaU,      "Nominal", "#eta down",  "#eta up",  odirPlots, "plot_eta")
plotter(predIhNominal,      predIhD,         predIhU,       "Nominal", "I_{h} down", "I_{h} up", odirPlots, "plot_ih")
plotter(predPNominal,       predPD,          predPU,        "Nominal", "p down",     "p up",     odirPlots, "plot_p")
plotter(predFitPNominal,    predFitPDown,    predFitPUp,    "Nominal", "Fit p down", "Fit p up", odirPlots, "plot_fitP")
plotter(predFitDeDxNominal, predFitDeDxDown, predFitDeDxUp, "Nominal", "Fit I_{h} down", "Fit I_{h} up", odirPlots, "plot_fitDeDx")
plotter(predNoFitNominal,   predNoFitDown,   None,          "Nominal", "No fit",     "",         odirPlots, "plot_nofit")
plotter(predcorrIhNominal,  predcorrIhDown,  None,          "Nominal", "corr template I_{h}", "", odirPlots, "plot_corrIh")
plotter(predcorr1oPNominal, predcorr1oPDown, predcorr1oPUp, "Nominal", "corr template 1/p down", "corr template 1/p up", odirPlots, "plot_corr1oP")

ofile.cd()
syst_stat.Write()
syst_stat_binned.Write()
syst_eta.Write()
syst_ih.Write()
syst_p.Write()
syst_fitP.Write()
syst_fitDeDx.Write()
syst_nofit.Write()
syst_corrIh.Write()
syst_corr1oP.Write()

#syst_nofit not taken !!
listOfSyst = [syst_stat, syst_eta, syst_ih, syst_p, syst_fitP, syst_fitDeDx, syst_corrIh, syst_corr1oP] #syst_nofit
sysTot = systTotal(listOfSyst)
sysTot.SetName("systTotal")
sysTot.Write()
sysTot.SaveAs(oDir + "sysTot.root")

listOfSyst_binned = [syst_stat_binned,
                     syst_eta_binned,
                     syst_ih_binned, 
                     syst_p_binned, 
                     syst_fitP_binned, 
                     syst_fitDeDx_binned, 
                     #syst_nofit_binned, 
                     syst_corrIh_binned,
                     syst_corr1oP_binned]
sysTot_binned = systTotal(listOfSyst_binned)
sysTot_binned.SetName("systTotalBinned")
sysTot_binned.Write()
sysTot_binned.SaveAs(oDir + "sysTotBinned_" + year + "_" + region + ".root")


ofileAllPred = TFile(oDir+"massShapePred.root","RECREATE")
ofileAllPred.cd()
predNominal_def.Write()
sysTot.Write()
sysTot_binned.Write()

plotSummary(syst_stat_binned,
            syst_eta_binned,
            syst_ih_binned,
            syst_p_binned,
            syst_fitDeDx_binned,
            syst_fitP_binned,
            syst_nofit_binned,
            syst_corrIh_binned,
            syst_corr1oP_binned,
            sysTot_binned,
            "Mass bin",
            "binned_syst",
            etaname)


print("Everything is done, closing the output file.")