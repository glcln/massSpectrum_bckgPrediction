#!/usr/bin/python

import sys, getopt, os
sys.argv.append( '-b-' )
import ROOT
import math
import array
import numpy as np

from ROOT import TFile, THStack, TCanvas, TLegend, TLatex, TPad, TH1, TH2, TLine, TGraph, TGraphErrors
import CMS_lumi, tdrstyle

ROOT.gROOT.SetBatch(True)

tdrstyle.setTDRStyle()

def ratioHisto(h2,h1):
    h3=h1.Clone()
    h3.Divide(h2)
    return h3

def ratioInt(h1,h2):
    h3=h2.Clone()
    h3.Reset()
    for i in range(0,h1.GetNbinsX()+1):
        e1=ROOT.Double(0.0)
        e2=ROOT.Double(0.0)
        a=h1.IntegralAndError(i,h1.GetNbinsX(),e1,"")
        b=h2.IntegralAndError(i,h1.GetNbinsX(),e2,"")
        if b != 0 and a != 0:
            c=math.sqrt((e1*e1)/(a*a)+(e2*e2)/(b*b))*a/b
            h3.SetBinContent(i,a/b)
            h3.SetBinError(i,c)
        else:
            h3.SetBinContent(i,0)
    return h3

def lowEdge(h1):
    res=ROOT.TGraph(h1.GetNbinsX()-1)
    for i in range (1,h1.GetNbinsX()+1):
        res.SetPoint(i-1,h1.GetBinLowEdge(i),h1.GetBinContent(i))
    return res

def statErr(h1,name):
    statErr=h1.Clone()
    statErr.SetName(name)
    for i in range (1,statErr.GetNbinsX()):
        if statErr.GetBinContent(i)>0:
            statErr.SetBinContent(i,statErr.GetBinError(i)/statErr.GetBinContent(i))
        else:
            statErr.SetBinContent(i,0)
    return 100*statErr

def statErrRInt(h1,name):
    statErr=h1.Clone()
    statErr.SetName(name)
    c=ROOT.Double(0.0)
    e=ROOT.Double(0.0)
    for i in range (1,statErr.GetNbinsX()):
        c=h1.IntegralAndError(i,h1.GetNbinsX(),e,"")
        if e>0:
            statErr.SetBinContent(i,e/c)
        else:
            statErr.SetBinContent(i,0)
    return 100*statErr

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

def setColorAndMarker(h1,color,markerstyle):
    h1.SetLineColor(color)
    h1.SetMarkerColor(color)
    h1.SetFillColor(color)
    h1.SetMarkerStyle(markerstyle)
    return h1

def plotter(predNominal,predPullD,predPullU,legNom,leg1,leg2,outDir,outTitle):
    c1=TCanvas("c1","c1",800,800)
    t1=TPad("t1","t1", 0.0, 0.40, 0.95, 0.95)

    t1.Draw()
    #t1.cd()
    t1.SetLogy(1)
    t1.SetGrid(1)
    t1.SetTopMargin(0.005)
    t1.SetBottomMargin(0.005)
    c1.cd()

    t2=TPad("t2","t2", 0.0, 0.225, 0.95, 0.375)
    t2.Draw()
    #t2.cd()
    t2.SetGridy(1)
    t2.SetTopMargin(0.005)
    t2.SetBottomMargin(0.005)

    t3=TPad("t3","t3", 0.0, 0., 0.95, 0.20)
    t3.Draw()
    #t3.cd()
    t3.SetGridy(1)
    t3.SetTopMargin(0.005)
    t3.SetBottomMargin(0.4)

    '''
    t4=TPad("t4","t4", 0.0, 0.0, 0.95, 0.15)
    t4.Draw()
    #t4.cd()
    t4.SetGridy(1)
    t4.SetTopMargin(0.005)
    t4.SetBottomMargin(0.4)
    '''

    t1.cd()
    c1.SetLogy(1)

    max_mass=2000
    min_entries=predNominal.GetMinimum()/10
    max_entries=predNominal.GetMaximum()*10

    predNominal=setColorAndMarker(predNominal,1,20)
    predNominal.GetXaxis().SetRangeUser(0,max_mass)
    predNominal.GetYaxis().SetRangeUser(min_entries,max_entries)
    predNominal.SetTitle(";Mass (GeV);Tracks")
    predNominal.GetYaxis().SetTitleSize(0.07)
    predNominal.GetYaxis().SetLabelSize(0.05)
    predNominal.Draw("P")

    predPullD=setColorAndMarker(predPullD,38,21)
    predPullD.Draw("P,same")

    predPullU=setColorAndMarker(predPullU,46,21)
    predPullU.Draw("P,same")

    
    leg=TLegend(0.82,0.9,0.4,0.6)
    leg.AddEntry(predNominal,legNom,"PE1")
    leg.AddEntry(predPullD,leg1,"PE1")
    leg.AddEntry(predPullU,leg2,"PE1")

    
    LineLastBin=TLine(predNominal.GetBinLowEdge(predNominal.FindBin(max_mass)-1),0,predNominal.GetBinLowEdge(predNominal.FindBin(max_mass)-1),max_entries)
    LineLastBin.SetLineStyle(1)
    LineLastBin.SetLineColor(1)
    LineLastBin.Draw("same")
   
    t=ROOT.TText(0.95,0.7,"+overflow")
    t.SetNDC(True)
    t.SetTextColor(1)
    t.SetTextFont(43)
    t.SetTextSize(24)
    t.SetTextAngle(90)
    t.Draw("same")

    leg.Draw("same")
    
    c1.cd()
    t2.cd()

    frameR2=ROOT.TH1D("frameR2", "frameR2", 1,0, max_mass)
    frameR2.GetXaxis().SetNdivisions(505)
    frameR2.SetTitle("")
    frameR2.SetStats(0)
    frameR2.GetXaxis().SetTitle("")
    frameR2.GetYaxis().SetTitle("RatioR")
    frameR2.GetXaxis().SetRangeUser(0,max_mass)
    frameR2.SetMaximum(1.5)
    frameR2.SetMinimum(0.5)
    frameR2.GetYaxis().SetLabelFont(43) #give the font size in pixel (instead of fraction)
    frameR2.GetYaxis().SetLabelSize(20) #font size
    frameR2.GetYaxis().SetTitleFont(43) #give the font size in pixel (instead of fraction)
    frameR2.GetYaxis().SetTitleSize(20) #font size
    frameR2.GetYaxis().SetNdivisions(205)
    frameR2.GetYaxis().SetTitleOffset(2)
    frameR2.GetXaxis().SetNdivisions(510)
    frameR2.GetXaxis().SetLabelFont(43) #give the font size in pixel (instead of fraction)
    frameR2.GetXaxis().SetLabelSize(16) #font size
    frameR2.GetXaxis().SetTitleFont(43) #give the font size in pixel (instead of fraction)
    frameR2.GetXaxis().SetTitleSize(24) #font size
    frameR2.GetXaxis().SetTitleOffset(3.75)
    frameR2.Draw("AXIS")

    LineAtOne=TLine(0,1,max_mass,1)
    LineAtOne.SetLineStyle(3)
    LineAtOne.SetLineColor(1)
    LineAtOne.Draw("same")
    
    LineAt1p2=TLine(0,1.2,max_mass,1.2)
    LineAt1p2.SetLineStyle(4)
    LineAt1p2.SetLineColor(1)
    LineAt1p2.Draw("same")
    
    LineAt0p8=TLine(0,0.8,max_mass,0.8)
    LineAt0p8.SetLineStyle(4)
    LineAt0p8.SetLineColor(1)
    LineAt0p8.Draw("same")

    ratioInt1=ratioInt(predPullD,predNominal)
    ratioInt1=setColorAndMarker(ratioInt1,38,21)
    ratioInt2=ratioInt(predPullU,predNominal)
    ratioInt2=setColorAndMarker(ratioInt2,46,21)
    
    ratioInt1.Draw("P E0 same")
    ratioInt2.Draw("P E0 same")

    ofile.cd()
    ratioInt1.Write()
    ratioInt2.Write()

    c1.cd()
    t3.cd()

    frameR3=frameR2.Clone()
    frameR3.GetYaxis().SetTitle("#frac{var}{nominal}")
    frameR3.GetXaxis().SetRangeUser(0,max_mass)
    frameR3.Draw("AXIS")
    frameR3.GetXaxis().SetTitle("Mass (GeV)")
    frameR3.GetXaxis().SetTitleOffset(5)

    LineAtOne.Draw("same")
    LineAt1p2.Draw("same")
    LineAt0p8.Draw("same")
    
    predPullDCl=predPullD.Clone()
    predPullUCl=predPullU.Clone()
    predPullDCl.Divide(predNominal)
    predPullUCl.Divide(predNominal)
    
    predPullDCl=setColorAndMarker(predPullDCl,38,21)
    predPullDCl.Draw("P E0 same")

    predPullUCl=setColorAndMarker(predPullUCl,46,21)
    predPullUCl.Draw("P E0 same")

    cmd='mkdir -p '+outDir
    os.system(cmd)

    c1.SaveAs(outDir+"/"+outTitle+".pdf")
    c1.SaveAs(outDir+"/"+outTitle+".root")
    c1.SaveAs(outDir+"/"+outTitle+".C")

def plotSummary(syst_stat,syst_eta,syst_ih,syst_p,syst_corrDeDx,syst_corrP,syst_fitDeDx,syst_fitP,syst_corrBias,sysTot,xtitle,outTitle,labelRegion,label_lowEdge=0):
    syst_stat=lowEdge(syst_stat)
    syst_eta=lowEdge(syst_eta)
    syst_ih=lowEdge(syst_ih)
    syst_p=lowEdge(syst_p)
    syst_corrDeDx=lowEdge(syst_corrDeDx)
    syst_corrP=lowEdge(syst_corrP)
    syst_fitP=lowEdge(syst_fitP)
    syst_fitDeDx=lowEdge(syst_fitDeDx)
    syst_corrBias=lowEdge(syst_corrBias)
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

    syst_stat.GetXaxis().SetRangeUser(0,2000)

    syst_stat=setColorAndMarker(syst_stat,1,20)
    syst_eta=setColorAndMarker(syst_eta,30,21)
    syst_ih=setColorAndMarker(syst_ih,38,22)
    syst_p=setColorAndMarker(syst_p,46,23)
    syst_corrDeDx=setColorAndMarker(syst_corrDeDx,43,43)
    syst_corrP=setColorAndMarker(syst_corrP,45,45)
    syst_fitP=setColorAndMarker(syst_fitP,47,47)
    syst_fitDeDx=setColorAndMarker(syst_fitDeDx,39,29)
    syst_corrBias=setColorAndMarker(syst_corrBias,34,49)
    sysTot=setColorAndMarker(sysTot,28,34)

    leg2=TLegend(0.17,0.7,0.4,0.93)
    leg2.AddEntry(sysTot,"Total","PE1")
    leg2.AddEntry(syst_stat,"Stat.","PE1")
    leg2.AddEntry(syst_eta,"#eta binning","PE1")
    leg2.AddEntry(syst_ih,"I_{h} binning","PE1")
    leg2.AddEntry(syst_p,"p binning","PE1")
    leg2.AddEntry(syst_corrDeDx,"Template extraction I_{h}","PE1")
    leg2.AddEntry(syst_corrP,"Template extraction p","PE1")
    leg2.AddEntry(syst_fitP,"p fit","PE1")
    leg2.AddEntry(syst_fitDeDx,"I_{h} fit","PE1")
    leg2.AddEntry(syst_corrBias,"Bias correction ("+labelRegion+")","PE1")
    
    syst_stat.Draw("AP")
    leg2.Draw("same")
    syst_eta.Draw("P")
    syst_ih.Draw("P")
    syst_p.Draw("P")
    syst_corrDeDx.Draw("P")
    syst_corrP.Draw("P")
    syst_fitP.Draw("P")
    syst_fitDeDx.Draw("P")
    syst_corrBias.Draw("P")
    sysTot.Draw("P")
    
    commandMkdir='mkdir -p '+oDir+'pdf '+oDir+'Cfile '+oDir+'rootfile'
    os.system(commandMkdir)

    c2.SaveAs(oDir+"pdf/summary_"+outTitle+".pdf")
    c2.SaveAs(oDir+"Cfile/summary_"+outTitle+".root")
    c2.SaveAs(oDir+"rootfile/summary_"+outTitle+".C")

def allSet(h,sizeRebinning,rebinning,st):
    h=h.Rebin(sizeRebinning,st,rebinning)
    h.Scale(1./h.Integral())
    return h

def systMassAll(nominal,down,up,st,mini=0):
    res=systMass(nominal,down,up,st,0,0,mini)
    res_mean=systMass(nominal,down,up,st,1,0,mini)
    res_binned=systMass(nominal,down,up,st,0,1,mini)
    res_binned_mean=systMass(nominal,down,up,st,1,1,mini)
    return res, res_mean, res_binned, res_binned_mean

def BiasCorrection(h1,a_,b_):
    h=h1.Clone()
    for i in range (0,h.GetNbinsX()+1):
        mass = h.GetBinLowEdge(i)
        if(mass<25): continue
        h.SetBinContent(i,h.GetBinContent(i)*(a_*mass+b_))
    return h



# Setup
version = "V3p2"
directory = "/opt/sbg/cms/ui3_data1/gcoulon/CMSSW_10_6_30/src/HSCPTreeAnalyzer/macros/Fpix_"+version+"/"
year = "2017_2018"
region = "8fp9"
plotType = "mass_predBC_"
oDir = "/opt/sbg/cms/ui3_data1/gcoulon/CMSSW_10_6_30/src/HSCPTreeAnalyzer/macros/Fpix_"+version+"/SystCombined/"
ofile = TFile(oDir+"sysToTBinned_"+year+"_"+region+".root","RECREATE")
outTitle = "syst"
labelRegion = region
sample = "MET_2017_2018"




inputNominal = directory + sample+"_massCut_0_pT70_"+version+"_Fpix_Eta2p4_cutIndex3_rebinEta4_rebinIh4_rebinP2_rebinMass1_EtaReweighting.root"
ifileNominal = TFile(inputNominal)

inputEtaD = directory + sample+"_massCut_0_pT70_"+version+"_Fpix_Eta2p4_cutIndex3_rebinEta2_rebinIh4_rebinP2_rebinMass1_EtaReweighting.root"
inputEtaU = directory + sample+"_massCut_0_pT70_"+version+"_Fpix_Eta2p4_cutIndex3_rebinEta8_rebinIh4_rebinP2_rebinMass1_EtaReweighting.root"
ifileEtaD = TFile(inputEtaD)
ifileEtaU = TFile(inputEtaU)

inputIhD = directory + sample+"_massCut_0_pT70_"+version+"_Fpix_Eta2p4_cutIndex3_rebinEta4_rebinIh2_rebinP2_rebinMass1_EtaReweighting.root"
inputIhU = directory + sample+"_massCut_0_pT70_"+version+"_Fpix_Eta2p4_cutIndex3_rebinEta4_rebinIh8_rebinP2_rebinMass1_EtaReweighting.root"
ifileIhD = TFile(inputIhD)
ifileIhU = TFile(inputIhU)

inputPD = directory + sample+"_massCut_0_pT70_"+version+"_Fpix_Eta2p4_cutIndex3_rebinEta4_rebinIh4_rebinP1_rebinMass1_EtaReweighting.root"
inputPU = directory + sample+"_massCut_0_pT70_"+version+"_Fpix_Eta2p4_cutIndex3_rebinEta4_rebinIh4_rebinP4_rebinMass1_EtaReweighting.root"
ifilePD = TFile(inputPD)
ifilePU = TFile(inputPU)

inputCorrDeDxDown = directory + sample+"_massCut_0_pT70_"+version+"_Fpix_Eta2p4_cutIndex3_rebinEta4_rebinIh4_rebinP2_rebinMass1_corrTemplateIh_EtaReweighting.root"
inputCorrDeDxUp = directory + sample+"_massCut_0_pT70_"+version+"_Fpix_Eta2p4_cutIndex3_rebinEta4_rebinIh4_rebinP2_rebinMass1_corrTemplateIh_EtaReweighting.root"
ifileCorrDeDxDown = TFile(inputCorrDeDxDown)
ifileCorrDeDxUp = TFile(inputCorrDeDxUp)

#inputCorrPDown = directory + sample+"_massCut_0_pT70_"+version+"_Fpix_Eta2p4_cutIndex3_rebinEta4_rebinIh4_rebinP2_rebinMass1_corrTemplateP_EtaReweighting.root"
#inputCorrPUp = directory + sample+"_massCut_0_pT70_"+version+"_Fpix_Eta2p4_cutIndex3_rebinEta4_rebinIh4_rebinP2_rebinMass1_corrTemplateP_EtaReweighting.root"
#ifileCorrPDown = TFile(inputCorrPDown)
#ifileCorrPUp = TFile(inputCorrPUp)

inputFitPUp = directory + sample+"_massCut_0_pT70_"+version+"_Fpix_Eta2p4_cutIndex3_rebinEta4_rebinIh4_rebinP2_rebinMass1_fitPUp_EtaReweighting.root"
inputFitPDown = directory + sample+"_massCut_0_pT70_"+version+"_Fpix_Eta2p4_cutIndex3_rebinEta4_rebinIh4_rebinP2_rebinMass1_fitPDown_EtaReweighting.root"
ifileFitPUp = TFile(inputFitPUp)
ifileFitPDown = TFile(inputFitPDown)

inputFitDeDxUp = directory + sample+"_massCut_0_pT70_"+version+"_Fpix_Eta2p4_cutIndex3_rebinEta4_rebinIh4_rebinP2_rebinMass1_fitIhUp_EtaReweighting.root"
inputFitDeDxDown = directory + sample+"_massCut_0_pT70_"+version+"_Fpix_Eta2p4_cutIndex3_rebinEta4_rebinIh4_rebinP2_rebinMass1_fitIhDown_EtaReweighting.root"
ifileFitDeDxUp = TFile(inputFitDeDxUp)
ifileFitDeDxDown = TFile(inputFitDeDxDown)


print(inputNominal)

# Compute the systematics
predNominal_def = ifileNominal.Get(plotType+region)

predEtaD = ifileEtaD.Get(plotType+region)
predEtaU = ifileEtaU.Get(plotType+region)

predIhD = ifileIhD.Get(plotType+region)
predIhU = ifileIhU.Get(plotType+region)

predPD = ifilePD.Get(plotType+region)
predPU = ifilePU.Get(plotType+region)

predCorrDeDxDown = ifileCorrDeDxDown.Get(plotType+region)
predCorrDeDxUp = ifileCorrDeDxUp.Get(plotType+region)

#predCorrPDown = ifileCorrPDown.Get(plotType+region)
#predCorrPUp = ifileCorrPUp.Get(plotType+region)

predFitPUp = ifileFitPUp.Get(plotType+region)
predFitPDown = ifileFitPDown.Get(plotType+region)

predFitDeDxUp = ifileFitDeDxUp.Get(plotType+region)
predFitDeDxDown = ifileFitDeDxDown.Get(plotType+region)


rebinning = array.array('d',[0.,20.,40.,60.,80.,100.,120.,140.,160.,180.,200.,220.,240.,260.,280.,300.,320.,340.,360.,380.,410.,440.,480.,530.,590.,660.,760.,880.,1030.,1210.,1440.,1730.,2000.,2500.,3200.,4000.])
sizeRebinning = len(rebinning)-1

predNominal_def = allSet(predNominal_def, sizeRebinning, rebinning, "nominal_def")

syst_stat = statErrRInt(predNominal_def,"Stat")
syst_stat_binned = statErr(predNominal_def,"Stat")


predEtaNominal = allSet(predNominal_def,sizeRebinning,rebinning,"eta_nominal")
predEtaD = allSet(predEtaD,sizeRebinning,rebinning,"eta_down")
predEtaU = allSet(predEtaU,sizeRebinning,rebinning,"eta_up")

predIhNominal = allSet(predNominal_def,sizeRebinning,rebinning,"ih_nominal")
predIhD = allSet(predIhD,sizeRebinning,rebinning,"ih_down")
predIhU = allSet(predIhU,sizeRebinning,rebinning,"ih_up")

predPNominal = allSet(predNominal_def,sizeRebinning,rebinning,"p_nominal")
predPD = allSet(predPD,sizeRebinning,rebinning,"p_down")
predPU = allSet(predPU,sizeRebinning,rebinning,"p_up")

predCorrDeDxNominal = allSet(predNominal_def,sizeRebinning,rebinning,"corrdedx_nominal")
predCorrDeDxDown = allSet(predCorrDeDxDown,sizeRebinning,rebinning,"corrdedx_down")
predCorrDeDxUp = allSet(predCorrDeDxUp,sizeRebinning,rebinning,"corrdedx_up")

#predCorrPNominal = allSet(predNominal_def,sizeRebinning,rebinning,"corrP_nominal")
#predCorrPDown = allSet(predCorrPDown,sizeRebinning,rebinning,"corrP_down")
#predCorrPUp = allSet(predCorrPUp,sizeRebinning,rebinning,"corrP_up")

predFitPNominal = allSet(predNominal_def,sizeRebinning,rebinning,"fitP_nominal")
predFitPDown = allSet(predFitPDown,sizeRebinning,rebinning,"fitP_down")
predFitPUp = allSet(predFitPUp,sizeRebinning,rebinning,"fitP_up")

predFitDeDxNominal = allSet(predNominal_def,sizeRebinning,rebinning,"fitDeDx_nominal")
predFitDeDxDown = allSet(predFitDeDxDown,sizeRebinning,rebinning,"fitDeDx_down")
predFitDeDxUp = allSet(predFitDeDxUp,sizeRebinning,rebinning,"fitDeDx_up")



(syst_eta, syst_eta_mean, syst_eta_binned, syst_eta_binned_mean) = systMassAll(predEtaNominal,predEtaD,predEtaU,"Eta")
(syst_ih, syst_ih_mean, syst_ih_binned, syst_ih_binned_mean) = systMassAll(predIhNominal,predIhD,predIhU,"Ih")
(syst_p, syst_p_mean, syst_p_binned, syst_p_binned_mean) = systMassAll(predPNominal,predPD,predPU,"P")

(syst_corrDeDx, syst_corrDeDx_mean, syst_corrDeDx_binned, syst_corrDeDx_binned_mean) = systMassAll(predCorrDeDxNominal,predCorrDeDxDown,predCorrDeDxUp,"corrDeDx")
#(syst_corrP, syst_corrP_mean, syst_corrP_binned, syst_corrP_binned_mean) = systMassAll(predCorrPNominal,predCorrPDown,predCorrPUp,"corrP")

(syst_fitP, syst_fitP_mean, syst_fitP_binned, syst_fitP_binned_mean) = systMassAll(predFitPNominal,predFitPDown,predFitPUp,"Fit_p")
(syst_fitDeDx, syst_fitDeDx_mean, syst_fitDeDx_binned, syst_fitDeDx_binned_mean) = systMassAll(predFitDeDxNominal,predFitDeDxDown,predFitDeDxUp,"Fit_dedx_systMassAl")


ofile.cd()
syst_stat.Write()
syst_stat_binned.Write()
syst_eta.Write()
syst_ih.Write()
syst_p.Write()
syst_corrDeDx.Write()
#syst_corrP.Write()
syst_fitDeDx.Write()
syst_fitDeDx_binned.Write()

listOfSyst = [syst_stat,syst_eta,syst_ih,syst_p,syst_corrDeDx,syst_fitDeDx,syst_fitP] #syst_corrP
sysTot = systTotal(listOfSyst)
sysTot.SetName("systTotal")
sysTot.Write()
sysTot.SaveAs(oDir+"sysTot.root")

listOfSyst_binned = [syst_stat_binned,syst_eta_binned,syst_ih_binned,syst_p_binned,syst_corrDeDx_binned,syst_fitDeDx_binned,syst_fitP_binned] #syst_corrP_binned
sysTot_binned = systTotal(listOfSyst_binned)
sysTot_binned.SetName("systTotalBinned")
sysTot_binned.Write()
sysTot_binned.SaveAs(oDir+"sysTotBinned_"+year+"_"+region+".root")


#plotter(predNominal_def,predEtaD,predEtaU,"Nominal (20 bins)","Down (10 bins)","Up (40 bins)",oDir+"syst_Eta","plot_Eta")
#plotter(predNominal_def,predIhD,predIhU,"Nominal (100 bins)","Down (50 bins)","Up (200 bins)",oDir+"syst_Ih","plot_Ih")
#plotter(predNominal_def,predPD,predPU,"Nominal (500 bins)","Down (250 bins)","Up (1000 bins)",oDir+"syst_P","plot_P")
#plotter(predNominal_def,predCorrDeDxDown,predCorrDeDxUp,"Nominal","Down","Up",oDir+"syst_corrDeDx","plot_corrDeDx")
#plotter(predNominal_def,predCorrPDown,predCorrPUp,"Nominal","Down","Up",oDir+"syst_corrP","plot_corrP")

#plotter(predFitPNominal,predFitPDown,predFitPUp,"Nominal","Down","Up",oDir+"syst_FitP","plot_FitP")
#plotter(predFitDeDxNominal,predFitDeDxDown,predFitDeDxUp,"Nominal","Down","Up",oDir+"syst_FitDeDx","plot_FitDeDx")

#xtitle = "Mass bin"
#plotSummary(syst_stat_binned,syst_eta_binned,syst_ih_binned,syst_p_binned,syst_corrDeDx_binned,syst_corrP_binned,syst_fitDeDx_binned,syst_fitP_binned,syst_corrBias_binned,sysTot_binned,xtitle,"binned_"+outTitle,labelRegion)


ofileAllPred = TFile(oDir+"massShapePred.root","RECREATE")
ofileAllPred.cd()
predNominal_def.Write()
sysTot.Write()
sysTot_binned.Write()
