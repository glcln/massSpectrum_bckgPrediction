#!/usr/bin/python

import sys, getopt, os
sys.argv.append( '-b-' )
import ROOT
import math
import array
import numpy as np
import ctypes

from ROOT import TFile, THStack, TCanvas, TLegend, TLatex, TPad, TH1, TH2, TLine, TGraph, TGraphErrors
import CMS_lumi, tdrstyle

ROOT.gROOT.SetBatch(True)

tdrstyle.setTDRStyle()

def ratioHisto(h2,h1):
    h3=h1.Clone()
    h3.Divide(h2)
    return h3

def ratioInt(h1, h2):
    h3 = h2.Clone()
    h3.Reset()
    for i in range(0, h1.GetNbinsX()+1):
        e1 = ctypes.c_double(0.0)
        e2 = ctypes.c_double(0.0)
        a = h1.IntegralAndError(i, h1.GetNbinsX(), e1, "")
        b = h2.IntegralAndError(i, h1.GetNbinsX(), e2, "")
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

def systMass(nominal,down,up,name,typec,binned,mini=0):
    if (binned==0):
        ra1=ratioInt(nominal,up)
        ra2=ratioInt(nominal,down)
    elif (binned==1):
        ra1=ratioHisto(nominal,up)
        ra2=ratioHisto(nominal,down)
    res=ra1.Clone()
    res.SetName(name)
    for i in range (0,res.GetNbinsX()+1):
        r1=ra1.GetBinContent(i)
        r2=ra2.GetBinContent(i)

        s1=abs(1-r1)
        s2=abs(1-r2)
        s1=abs(r1-1)
        s2=abs(r2-1)

        if (typec==0):
            m=max(s1,s2)
            if (mini==1):
                m=min(s1,s2)
            elif (mini==2):
                m/=2
        elif (typec==1):
            m=(s1+s2)/2.
        res.SetBinContent(i,100*m)
    return res

def systTotal(list_h):
    res=list_h[0].Clone()
    for i in range (0,res.GetNbinsX()+1):
        systotal=0
        for h in list_h:
            systotal+=h.GetBinContent(i)*h.GetBinContent(i)
        res.SetBinContent(i,math.sqrt(systotal))
    return res

def RebinHisto(h, sizeRebinning,  rebinning , st):
    h=h.Rebin(sizeRebinning, st, rebinning)
    h.Scale(1./h.Integral())
    return h

def systMassAll(nominal,down,up,st,mini=0):
    res = systMass(nominal,down,up,st,0,0,mini)
    res_mean = systMass(nominal,down,up,st,1,0,mini)
    res_binned = systMass(nominal,down,up,st,0,1,mini)
    res_binned_mean = systMass(nominal,down,up,st,1,1,mini)
    return res, res_mean, res_binned, res_binned_mean

def setColorAndMarker(h1,color,markerstyle):
    h1.SetLineColor(color)
    h1.SetMarkerColor(color)
    h1.SetFillColor(color)
    h1.SetMarkerStyle(markerstyle)
    return h1

def lowEdge(h1):
    res=ROOT.TGraph(h1.GetNbinsX()-1)
    for i in range (1,h1.GetNbinsX()+1):
        res.SetPoint(i-1,h1.GetBinLowEdge(i),h1.GetBinContent(i))
    return res

def plotSummary(syst_stat,
                syst_eta,
                syst_ih,
                syst_p,
                syst_fitDeDx,
                syst_fitP,
                syst_nofit,
                syst_corrIh,
                sysTot,
                xtitle,
                outTitle,
                labelRegion,
                label_lowEdge=0):
    syst_stat=lowEdge(syst_stat)
    syst_eta=lowEdge(syst_eta)
    syst_ih=lowEdge(syst_ih)
    syst_p=lowEdge(syst_p)
    syst_fitP=lowEdge(syst_fitP)
    syst_nofit=lowEdge(syst_nofit)
    syst_corrIh=lowEdge(syst_corrIh)
    syst_fitDeDx=lowEdge(syst_fitDeDx)
    sysTot=lowEdge(sysTot)

    c2=TCanvas()

    c2.SetLogy()
    c2.SetGrid()
    syst_stat.SetMinimum(1)
    syst_stat.SetMaximum(500)
    syst_stat.GetXaxis().SetTitle(xtitle)
    syst_stat.GetYaxis().SetTitle("Systematic Uncertainty [%]")
    syst_stat.GetXaxis().SetNdivisions(510)
    syst_stat.GetXaxis().SetLabelFont(43) #give the font size in pixel (instead of fraction)
    syst_stat.GetXaxis().SetLabelSize(16) #font size
    syst_stat.GetXaxis().SetTitleSize(0.04) #font size
    syst_stat.GetYaxis().SetLabelFont(43) #give the font size in pixel (instead of fraction)
    syst_stat.GetYaxis().SetLabelSize(16) #font size
    syst_stat.GetYaxis().SetTitleSize(0.04) #font size


    # draw latex CMS preliminary
    latex = TLatex(0.16, 0.96, "#scale[1.3]{#bf{CMS}}#it{Work in progress}")
    latex.SetNDC()
    latex.SetTextFont(42)
    latex.SetTextSize(0.04)

    latex2 = TLatex(0.7, 0.96, "108.95 fb^{-1} (13.6 TeV)")
    latex2.SetNDC()
    latex2.SetTextFont(42)
    latex2.SetTextSize(0.04)

    syst_stat.GetXaxis().SetRangeUser(0,4000)
    syst_stat.SetMinimum(1)
    syst_stat.SetMaximum(2000)

    syst_stat=setColorAndMarker(syst_stat,1,20)
    syst_eta=setColorAndMarker(syst_eta,30,21)
    syst_ih=setColorAndMarker(syst_ih,38,22)
    syst_p=setColorAndMarker(syst_p,46,23)
    syst_fitP=setColorAndMarker(syst_fitP,47,47)
    syst_nofit=setColorAndMarker(syst_nofit,41,48)
    syst_corrIh=setColorAndMarker(syst_corrIh,13,49)
    syst_fitDeDx=setColorAndMarker(syst_fitDeDx,39,29)
    sysTot=setColorAndMarker(sysTot,28,34)

    leg2=TLegend(0.27,0.7,0.5,0.93)
    leg2.AddEntry(sysTot,"Total","PE1")
    leg2.AddEntry(syst_stat,"Stat.","PE1")
    leg2.AddEntry(syst_eta,"#eta binning","PE1")
    leg2.AddEntry(syst_ih,"I_{h} binning","PE1")
    leg2.AddEntry(syst_p,"p binning","PE1")
    leg2.AddEntry(syst_fitP,"p fit","PE1")
    leg2.AddEntry(syst_fitDeDx,"I_{h} fit","PE1")
    leg2.AddEntry(syst_nofit,"No fit","PE1")
    leg2.AddEntry(syst_corrIh,"corr template I_{h}","PE1")
   
    syst_stat.Draw("AP")
    leg2.Draw("same")
    syst_eta.Draw("P")
    syst_ih.Draw("P")
    syst_p.Draw("P")
    syst_fitP.Draw("P")
    syst_nofit.Draw("P")
    syst_corrIh.Draw("P")
    syst_fitDeDx.Draw("P")
    sysTot.Draw("P")
    
    latex.Draw()
    latex2.Draw()
    
    commandMkdir='mkdir -p '+oDir+'pdf '+oDir+'Cfile '+oDir+'rootfile'
    os.system(commandMkdir)

    c2.SaveAs(oDir+"pdf/summary_"+outTitle+".pdf")
    c2.SaveAs(oDir+"Cfile/summary_"+outTitle+".root")
    c2.SaveAs(oDir+"rootfile/summary_"+outTitle+".C")




# Setup
version     = "V12p24"
etaname     = "Eta2p4"
year        = "2024"
sample      = "JetMET"
region      = "8fp9"
option      = "_etaAbs_chi2cut"
directory   = "/opt/sbg/cms/safe1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/macros/DataMET_2024_V12p24__8fp9" + option + "/" + etaname + "/"
oDir        = directory + "SystCombined/"
ofile       = TFile(oDir + "sysToTBinned_" + year + "_" + region + ".root", "RECREATE")

plotType    = "mass_predBC_"



inputNominal = directory + sample + year + "_" + version + "_rebinEta4_rebinIh4_rebinP2_EtaReweighting_" + etaname + "_NewFit.root" 
ifileNominal = TFile(inputNominal)

inputEtaD = directory + sample + year + "_" + version + "_rebinEta2_rebinIh4_rebinP2_EtaReweighting_" + etaname + "_NewFit.root" 
inputEtaU = directory + sample + year + "_" + version + "_rebinEta8_rebinIh4_rebinP2_EtaReweighting_" + etaname + "_NewFit.root" 
ifileEtaD = TFile(inputEtaD)
ifileEtaU = TFile(inputEtaU)

inputIhD = directory + sample + year + "_" + version + "_rebinEta4_rebinIh2_rebinP2_EtaReweighting_" + etaname + "_NewFit.root" 
inputIhU = directory + sample + year + "_" + version + "_rebinEta4_rebinIh8_rebinP2_EtaReweighting_" + etaname + "_NewFit.root" 
ifileIhD = TFile(inputIhD)
ifileIhU = TFile(inputIhU)

inputPD = directory + sample + year + "_" + version + "_rebinEta4_rebinIh4_rebinP1_EtaReweighting_" + etaname + "_NewFit.root" 
inputPU = directory + sample + year + "_" + version + "_rebinEta4_rebinIh4_rebinP4_EtaReweighting_" + etaname + "_NewFit.root" 
ifilePD = TFile(inputPD)
ifilePU = TFile(inputPU)

inputFitPUp     = directory + sample + year + "_" + version + "_rebinEta4_rebinIh4_rebinP2_fitPUp_EtaReweighting_" + etaname + "_NewFit.root" 
inputFitPDown   = directory + sample + year + "_" + version + "_rebinEta4_rebinIh4_rebinP2_fitPDown_EtaReweighting_" + etaname + "_NewFit.root" 
ifileFitPUp     = TFile(inputFitPUp)
ifileFitPDown   = TFile(inputFitPDown)

inputFitDeDxUp      = directory + sample + year + "_" + version + "_rebinEta4_rebinIh4_rebinP2_fitIhUp_EtaReweighting_" + etaname + "_NewFit.root" 
inputFitDeDxDown    = directory + sample + year + "_" + version + "_rebinEta4_rebinIh4_rebinP2_fitIhDown_EtaReweighting_" + etaname + "_NewFit.root" 
ifileFitDeDxUp      = TFile(inputFitDeDxUp)
ifileFitDeDxDown    = TFile(inputFitDeDxDown)

inputNoFit = directory + sample + year + "_" + version + "_rebinEta4_rebinIh4_rebinP2_EtaReweighting_" + etaname + "_NoFit.root"
ifileNoFit = TFile(inputNoFit)

inputcorrTemplateIh = directory + sample + year + "_" + version + "_rebinEta4_rebinIh4_rebinP2_corrTemplateIh_EtaReweighting_" + etaname + "_NewFit.root"
ifilecorrTemplateIh = TFile(inputcorrTemplateIh)


print(inputNominal)


# Compute the systematics
predNominal_def = ifileNominal.Get(plotType + region)

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

predcorrIh = ifilecorrTemplateIh.Get(plotType + region)



rebinning = array.array('d',[0.,20.,40.,60.,80.,100.,120.,140.,160.,180.,200.,220.,240.,260.,280.,300.,320.,340.,360.,380.,410.,440.,480.,530.,590.,660.,760.,880.,1030.,1210.,1440.,1730.,2000.,2500.,3200.,4000.])
sizeRebinning = len(rebinning)-1

predNominal_def = RebinHisto(predNominal_def, sizeRebinning, rebinning, "nominal_def")

syst_stat = statErrRInt(predNominal_def, "Stat")
syst_stat_binned = statErr(predNominal_def, "Stat")


predEtaNominal = RebinHisto(predNominal_def, sizeRebinning, rebinning, "eta_nominal")
predEtaD       = RebinHisto(predEtaD, sizeRebinning, rebinning, "eta_down")
predEtaU       = RebinHisto(predEtaU, sizeRebinning, rebinning, "eta_up")

predIhNominal  = RebinHisto(predNominal_def, sizeRebinning, rebinning, "ih_nominal")
predIhD        = RebinHisto(predIhD, sizeRebinning, rebinning, "ih_down")
predIhU        = RebinHisto(predIhU, sizeRebinning, rebinning, "ih_up")

predPNominal   = RebinHisto(predNominal_def, sizeRebinning, rebinning, "p_nominal")
predPD         = RebinHisto(predPD, sizeRebinning, rebinning, "p_down")
predPU         = RebinHisto(predPU, sizeRebinning, rebinning, "p_up")

predFitPNominal = RebinHisto(predNominal_def, sizeRebinning, rebinning, "fitP_nominal")
predFitPDown    = RebinHisto(predFitPDown, sizeRebinning, rebinning, "fitP_down")
predFitPUp      = RebinHisto(predFitPUp, sizeRebinning, rebinning, "fitP_up")

predFitDeDxNominal  = RebinHisto(predNominal_def, sizeRebinning, rebinning, "fitDeDx_nominal")
predFitDeDxDown     = RebinHisto(predFitDeDxDown, sizeRebinning, rebinning, "fitDeDx_down")
predFitDeDxUp       = RebinHisto(predFitDeDxUp, sizeRebinning, rebinning, "fitDeDx_up")

predNoFitNominal  = RebinHisto(predNominal_def, sizeRebinning, rebinning, "Nofit_nominal")
predNoFitDown     = RebinHisto(predNoFit, sizeRebinning, rebinning, "Nofit_down")
predNoFitUp       = RebinHisto(predNoFit, sizeRebinning, rebinning, "Nofit_up")

predcorrIhNominal  = RebinHisto(predNominal_def, sizeRebinning, rebinning, "corrIh_nominal")
predcorrIhDown     = RebinHisto(predcorrIh, sizeRebinning, rebinning, "corrIh_down")
predcorrIhUp       = RebinHisto(predcorrIh, sizeRebinning, rebinning, "corrIh_up")



(syst_eta, syst_eta_mean, syst_eta_binned, syst_eta_binned_mean) = systMassAll(predEtaNominal,predEtaD,predEtaU,"Eta")
(syst_ih, syst_ih_mean, syst_ih_binned, syst_ih_binned_mean) = systMassAll(predIhNominal,predIhD,predIhU,"Ih")
(syst_p, syst_p_mean, syst_p_binned, syst_p_binned_mean) = systMassAll(predPNominal,predPD,predPU,"P")

(syst_fitP, syst_fitP_mean, syst_fitP_binned, syst_fitP_binned_mean) = systMassAll(predFitPNominal,predFitPDown,predFitPUp,"Fit_p")
(syst_fitDeDx, syst_fitDeDx_mean, syst_fitDeDx_binned, syst_fitDeDx_binned_mean) = systMassAll(predFitDeDxNominal,predFitDeDxDown,predFitDeDxUp,"Fit_dedx_systMassAl")
(syst_nofit, syst_nofit_mean, syst_nofit_binned, syst_nofit_binned_mean) = systMassAll(predNoFitNominal,predNoFitDown,predNoFitUp,"NoFit_systMassAl")
(syst_corrIh, syst_corrIh_mean, syst_corrIh_binned, syst_corrIh_binned_mean) = systMassAll(predcorrIhNominal,predcorrIhDown,predcorrIhUp,"corrIh_systMassAl")

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


listOfSyst = [syst_stat, syst_eta, syst_ih, syst_p, syst_fitP, syst_fitDeDx, syst_nofit, syst_corrIh]
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
                     syst_nofit_binned, 
                     syst_corrIh_binned]
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
            sysTot_binned,
            "Mass bin",
            "binned_syst",
            region)
