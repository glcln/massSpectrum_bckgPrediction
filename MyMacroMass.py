#!/usr/bin/python

import sys, getopt, os
sys.argv.append(' -b- ')
import ROOT
import math
import array
import numpy as np
sys.path.append("/opt/sbg/cms/ui3_data1/gcoulon")
from USE_DATE import USED_DATE, VERSION

from ROOT import THStack, TCanvas, TLegend, TLatex, TPad, TH1, TH2, TLine
import CMS_lumi, tdrstyle

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptFit(1111)
tdrstyle.setTDRStyle()

a_ = 0
b_ = 1

#year='2017'
year='2018'
#year='2017_2018'
#year='ttbarwjets'
#year='wjets'
#year='ttbar'
region_=""

CMS_lumi.lumi_sqrtS = "13 TeV"
#CMS_lumi.cmsText="CMS"
CMS_lumi.cmsText="Private work"
CMS_lumi.cmsTextFont = 52
CMS_lumi.cmsTextSize   = 0.6
CMS_lumi.cmsTextOffset = -0.05
CMS_lumi.lumi_13TeV = "F_{pixel}"
#CMS_lumi.lumi_13TeV = "2018 - 59.7 fb^{-1}"
#CMS_lumi.lumi_13TeV = "2017 - 41.5 fb^{-1}"
CMS_lumi.lumiTextOffset = 0.1
CMS_lumi.lumiTextSize     = 0.55
CMS_lumi.writeExtraText = False
CMS_lumi.extraText = "Internal"
iPos = 0
if( iPos==0 ): CMS_lumi.relPosX = 0.12
iPeriod = 4


#----------------------------------------------------
#                       Functions
#----------------------------------------------------

def setColorAndMarker(h1,color,markerstyle):
    h1.SetLineColor(color)
    h1.SetMarkerColor(color)
    h1.SetMarkerSize(1.2)
    h1.SetMarkerStyle(markerstyle)
    return h1

def poissonning(h):
    res = h.Clone()
    res.Reset()
    res.Sumw2(0)
    res.SetBinErrorOption(ROOT.TH1.kPoisson)

    for i in range (0,h.GetNbinsX()+1):
        for j in range(0,int(h.GetBinContent(i))):
            res.Fill(h.GetBinCenter(i))
    
    return res

def overflowInLastBin(h, data, mass_max_Display):
    debugprint = False
    if (mass_max_Display > 300):
        if(mass_max_Display < h.GetBinCenter(h.GetNbinsX()) + h.GetBinWidth(h.GetNbinsX())/2):
            if (debugprint): print 'Case where mass_max_Display < edge of histogram: ', mass_max_Display, ' ', h.GetBinCenter(h.GetNbinsX()) + h.GetBinWidth(h.GetNbinsX())/2
            bin_content = 0
            bin_error = 0

            if (debugprint): print ''
            if (debugprint): print '       ', h.GetName()
            if (debugprint): print 'mass_max_Display: ', mass_max_Display, 'histo edge last bin: ', h.GetBinCenter(h.GetNbinsX()) + h.GetBinWidth(h.GetNbinsX())/2
            for i in range (h.FindBin(mass_max_Display)-1, h.GetNbinsX()+1):
                bin_content += h.GetBinContent(i)
                bin_error += h.GetBinError(i)**2

                if (debugprint): print 'Bin #{}, mass {}: Points = {} +/- {}'.format(i,h.GetBinCenter(i),h.GetBinContent(i),h.GetBinError(i))

                h.SetBinContent(i,0)
                h.SetBinError(i,0)
                if (debugprint): print 'Bin #{}, mass {}: Points = {} +/- {}'.format(i,h.GetBinCenter(i),h.GetBinContent(i),h.GetBinError(i))

            h.SetBinContent(h.FindBin(mass_max_Display)-1, bin_content + h.GetBinContent(h.GetNbinsX()+1))
            if (data): h.SetBinError(h.FindBin(mass_max_Display)-1,math.sqrt(bin_error))
            else: h.SetBinError(h.FindBin(mass_max_Display)-1,math.sqrt(bin_error + h.GetBinError(h.GetNbinsX()+1)**2))

        elif (mass_max_Display == h.GetBinCenter(h.GetNbinsX()) + h.GetBinWidth(h.GetNbinsX())/2):
            if (debugprint): print 'Case where mass_max_Display == edge of histogram: ', mass_max_Display, ' ', h.GetBinCenter(h.GetNbinsX()) + h.GetBinWidth(h.GetNbinsX())/2
            h.SetBinContent(h.GetNbinsX(),h.GetBinContent(h.GetNbinsX())+h.GetBinContent(h.GetNbinsX()+1))
            
            if(data): h.SetBinError(h.GetNbinsX(),math.sqrt(h.GetBinContent(h.GetNbinsX())))
            else: h.SetBinError(h.GetNbinsX(),math.sqrt(h.GetBinError(h.GetNbinsX())**2+h.GetBinError(h.GetNbinsX()+1)**2))
            
            h.SetBinContent(h.GetNbinsX()+1,0)
            h.SetBinError(h.GetNbinsX()+1,0)
        
        else: print "Error: mass_max_Display > histo edge last bin"
    else:
        if (debugprint): print 'Case where mass_max_Display = ', mass_max_Display, ' useless to overflowInLastBin in this case'

def underflowInFirstBin(h,data=False):
    h.SetBinContent(1,h.GetBinContent(0)+h.GetBinContent(1))
    h.SetBinContent(0,0)
    if(data): 
        h.SetBinError(1,math.sqrt(h.GetBinContent(1)))
    else:
        h.SetBinError(1,math.sqrt(h.GetBinError(0)**2+h.GetBinError(1)**2))
    h.SetBinError(0,0)

def underflowAndOverflow(h, data, mass_max_Display):
    #underflowInFirstBin(h,data)
    overflowInLastBin(h, data, mass_max_Display)

def binWidth(h1):
    res = h1.Clone()
    for i in range (0,h1.GetNbinsX()+1):
        res.SetBinContent(i,h1.GetBinContent(i)/h1.GetBinWidth(i))
        res.SetBinError(i,h1.GetBinError(i)/h1.GetBinWidth(i))
    return res

def ratioHisto(h1,h2):
    h3=h1.Clone()
    h3.Divide(h2)
    return h3

def ratioIntegral(h1,h2,systErr,upTo=-1):
    h3=h2.Clone()
    h3.Reset()
    if(upTo==-1):
        bornUp=h1.GetNbinsX()+1
    else:
        bornUp=h1.FindBin(upTo)
    for i in range(0,bornUp):
        e1=ROOT.Double(0.0)
        e2=ROOT.Double(0.0)
        if(upTo==-1):
            a=h1.IntegralAndError(i,h1.GetNbinsX()+1,e1,"")
            b=h2.IntegralAndError(i,h1.GetNbinsX()+1,e2,"")
        else:
            a=h1.IntegralAndError(i,bornUp-1,e1,"")
            b=h2.IntegralAndError(i,bornUp-1,e2,"")
        if b != 0 and a != 0:
            c=math.sqrt((e1*e1)/(a*a)+(e2*e2)/(b*b))*a/b
            h3.SetBinContent(i,a/b)
            h3.SetBinError(i,c)
        else:
            h3.SetBinContent(i,0)
    return h3

def pullOfHisto(h2,h1,systErr):
    res=h1.Clone()
    for i in range (1,h1.GetNbinsX()+1):
        Perr=0
        Derr=0
        P=h1.GetBinContent(i)
        D=h2.GetBinContent(i)
      
        Perr=h1.GetBinError(i)
        Derr=h2.GetBinErrorLow(i)
        
        if(P+(Perr*Perr) > 0): res.SetBinContent(i,(D-P)/math.sqrt(P+(Perr*Perr)))
        else: res.SetBinContent(i,0)
        #res.SetBinContent(i,(D-P)/math.sqrt(P+(Perr*Perr)))
        #if (Derr*Derr+Perr*Perr > 0): res.SetBinContent(i,(D-P)/math.sqrt(Derr*Derr+Perr*Perr))
        #else: res.SetBinContent(i,0)
    return res

def addSyst(h,syst):
    res=h.Clone()
    res.Sumw2(0)
    for i in range (0,h.GetNbinsX()+1):
        res.SetBinError(i,math.sqrt(h.GetBinError(i)*h.GetBinError(i)+res.GetBinContent(i)*res.GetBinContent(i)*syst*syst))
    return res

def addHSyst(h,h_syst,hCorrBias):
    res=h.Clone()
    resD=h.Clone()
    resU=h.Clone()
    for i in range (0,h.GetNbinsX()+1):
        syst=h_syst.GetBinContent(i)/100
        j=i
        while (h_syst.GetBinContent(j)/100==0):
            j-=1
            syst=h_syst.GetBinContent(j)/100
        diffCorrBias = abs(hCorrBias.GetBinContent(i)-h.GetBinContent(i))
        errorTotal = math.sqrt(h.GetBinError(i)*h.GetBinError(i)+res.GetBinContent(i)*res.GetBinContent(i)*syst*syst+pow(diffCorrBias,2))
        res.SetBinError(i,errorTotal)
        resD.SetBinContent(i,res.GetBinContent(i)-errorTotal)
        resU.SetBinContent(i,res.GetBinContent(i)+errorTotal)
    return (res,resD,resU)

def addHSystSignal(h,h_syst):
    res=h.Clone()
    for i in range (0,h.GetNbinsX()+1):
        syst=h_syst.GetBinContent(i)/100
        j=i
        while (h_syst.GetBinContent(j)/100==0):
            j-=1
            syst=h_syst.GetBinContent(j)/100
        errorTotal = math.sqrt(h.GetBinError(i)*h.GetBinError(i)+res.GetBinContent(i)*res.GetBinContent(i)*syst*syst)
        res.SetBinError(i,errorTotal)
    return res

def testChi2with1(h,x=-1):
    chi2=0
    ndf=-1
    if(x!=-1):
        upTo=h.FindBin(x)
    else:
        upTo=h.GetNbinsX()
    for i in range (1,upTo):
        chi2+=pow(h.GetBinContent(i)-1,2)
        ndf+=1
    
    return (chi2,ndf,ROOT.TMath.Prob(chi2,ndf))

def BiasCorrection(h):
    for i in range (0,h.GetNbinsX()+1):
        mass = h.GetBinLowEdge(i)
        if(mass<25): continue

        #print(h.GetName(), ' mass: ',  mass, 'value: ', h.GetBinContent(i)*(a_*mass+b_))
        #h.SetBinContent(i,h.GetBinContent(i)*(a_*mass+b_))
        if(a_*mass+b_>0): h.SetBinContent(i,h.GetBinContent(i)*(a_*mass+b_))
        else: h.SetBinContent(i,h.GetBinContent(i))
        #h.SetBinContent(i,h.GetBinContent(i))

def blindAnyUp(h,m):
    for i in range (0,h.GetNbinsX()+1):
        mass = h.GetBinLowEdge(i)
        if(mass>m): 
            h.SetBinContent(i,0)

def blindMassUp(h,m):
    for i in range (0,h.GetNbinsX()+1):
        mass = h.GetBinLowEdge(i)
        if(mass>m): 
            h.SetBinContent(i,0)

def blindMass(h,m):
    for i in range (0,h.GetNbinsX()+1):
        mass = h.GetBinLowEdge(i)
        if(mass<m): 
            h.SetBinContent(i,0)


#----------------------------------------------------
#                       Main
#----------------------------------------------------

def main(argv):
    # -------------- Setup --------------
    outputfile=''
    region=''
    odir=''

    try:
        opts, args = getopt.getopt(argv,"hi:o:r:d",["ifile=","ofile=","region=","odir="])
    except getopt.GetoptError:
        print 'test.py -i <inputfile> -o <outputfile> -r <region> -d <odir>'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'test.py -i <inputfile> -o <outputfile> -r <region> -d <odir>'
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
        elif opt in ("-r", "--region"):
            region = arg
        elif opt in ("-d", "--odir"):
            odir = arg

    os.system('mkdir -p '+odir)
    outputfile = odir+'/'+outputfile

    print ' Input file: ', inputfile
    print 'Output file: ', outputfile
    print '     Region: ', region
    print '' 

    signal = True
    blind = False
    isBinWidth = False
    doRebin = True
    isData = True

    labelRegion=""
    regSignal="VR1"
    
    if(region=="8fp9" and a_==0):
        labelRegion="8fp9"
        regSignal="VR1"
    elif(region=="8fp9" and a_!=0):
        labelRegion="8fp9"  
        regSignal="VR1"

    if(region=="3fp8" and a_==0):
        labelRegion="3fp8"
        regSignal="VR1"
    elif(region=="3fp8" and a_!=0):
        labelRegion="3fp8"  
        regSignal="VR1"
    
    if(region=="50ias60" and a_==0):
        labelRegion="RBF(50-60)"
        regSignal="VR4"
    elif(region=="50ias60" and a_!=0):
        labelRegion="SR(50-60)"  
        regSignal="VR4"

    if(region=="60ias70" and a_==0):
        labelRegion="RBF(60-70)"
        regSignal="VR5"
    elif(region=="60ias70" and a_!=0):
        labelRegion="SR(60-70)"  
        regSignal="VR5"

    if(region=="70ias80" and a_==0):
        labelRegion="RBF(70-80)"
        regSignal="VR6"
    elif(region=="70ias80" and a_!=0):
        labelRegion="SR(70-80)"  
        regSignal="VR6"

    if(region=="80ias90" and a_==0):
        labelRegion="RBF(80-90)"
        regSignal="VR7"
    elif(region=="80ias90" and a_!=0):
        labelRegion="SR(80-90)"  
        regSignal="VR7"

    if(region=="50ias90" and a_==0):
        labelRegion="RBF(50-90)"
        regSignal="VR1"
        signal=True
    elif(region=="50ias90" and a_!=0):
        labelRegion="SR(50-90)" 
        regSignal="VR1"
        signal=True

    if(region=="90ias100" and a_==0):
        labelRegion="RBF(90-100)"
        regSignal="SR1"
        signal=True
    elif(region=="90ias100" and a_!=0):
        labelRegion="SR(90-100)"
        regSignal="SR1"  
        signal=True

    if(region=="99ias100" and a_==0):
        labelRegion="RBF(99-100)"
        regSignal="SR2"
        signal=True
    elif(region=="99ias100" and a_!=0):
        labelRegion="SR(99-100)"   
        regSignal="SR2"
        signal=True
    
    if(region=="999ias100" and a_==0):
        labelRegion="RBF(99.9-100)"
        regSignal="SR3"
        signal=True
    elif(region=="999ias100" and a_!=0):
        labelRegion="Mass approach Signal Region"     
        regSignal="SR3"
        signal=True

    region_ = labelRegion
    signal = False

    ifile = ROOT.TFile(inputfile)
    
    obs = ifile.Get("mass_obs_"+region)
    pred = ifile.Get("mass_predBC_"+region)
    if (region=="80ias90"): C_mass = ifile.Get("mass_regionC_ias50_ReRunRaph")
    elif (region=="8fp9"): C_mass = ifile.Get("mass_regionC_3fp8_ReRunRaph")
    else: C_mass = ifile.Get("mass_obs_"+region)
    pred_noSyst = addSyst(pred,0.0)



    idirSignalGlu = "/opt/sbg/cms/ui3_data1/gcoulon/HSCP_prod/UnblindedProd/V1p1/HSCPgluino_V1p1/"
    idirSignalStau = "/opt/sbg/cms/ui3_data1/gcoulon/HSCP_prod/UnblindedProd/V1p1/HSCPpairStau_V1p1/"

    ifileGl2400 = ROOT.TFile("/opt/sbg/cms/ui3_data1/gcoulon/CMSSW_10_6_30/src/HSCPTreeAnalyzer/macros/Gluino2400_massCut_0_pT70_V2p20_Gstrip_Fpix_Eta2p4_Scale.root")
    ifileGl1600 = ROOT.TFile(idirSignalGlu+"HSCPgluino_M-1600_merged.root")
    ifileGl2000 = ROOT.TFile(idirSignalGlu+"HSCPgluino_M-2000_merged.root")
    ifilePPStau557 = ROOT.TFile(idirSignalStau+"HSCPpairStau_M-557_merged.root")
    ifilePPStau871 = ROOT.TFile(idirSignalStau+"HSCPpairStau_M-871_merged.root")

    ntupleDir = "HSCParticleAnalyzer/BaseName/"

    m_Gl2400 = ifileGl2400.Get("mass_regionD_"+region+"_ReRunRaph")
    m_Gl1600 = ifileGl1600.Get(ntupleDir+"PostS_"+regSignal+"_Mass")
    m_Gl2000 = ifileGl2000.Get(ntupleDir+"PostS_"+regSignal+"_Mass")
    m_ppStau557 = ifilePPStau557.Get(ntupleDir+"PostS_"+regSignal+"_Mass")
    m_ppStau871 = ifilePPStau871.Get(ntupleDir+"PostS_"+regSignal+"_Mass")


    # -------------- Work on histograms --------------


    #rebinning=array.array('d',[0.,20.,40.,60.,80.,100.,120.,140.,160.,180.,200.,220.,240.,260.,280.,300.,320.,340.,360.,380.,410.,440.,480.,530.,590.,660.,760.,880.,1030.,1210.,1440.,1730.,2000.])
    rebinning=array.array('d',[0.,20.,40.,60.,80.,100.,120.,140.,160.,180.,200.,220.,240.,260.,280.,300.,320.,340.,360.,380.,410.,440.,480.,530.,590.,660.,760.,880.,1030.,1210.,1440.,1730.,2000.,2500.,3200.,4000.])

    sizeRebinning=len(rebinning)-1
    
    if(doRebin==True):
        pred=pred.Rebin(sizeRebinning,"pred_new",rebinning)
        pred_noSyst=pred_noSyst.Rebin(sizeRebinning,"pred_noSyst_new",rebinning)
        obs=obs.Rebin(sizeRebinning,"obs_new",rebinning)
        C_mass = C_mass.Rebin(sizeRebinning,"C_mass_new",rebinning)

    normSignal=1
    if(year=="2017_2018"):
        normSignal=1
        #CMS_lumi.lumi_13TeV = "2017+2018 - 101 fb^{-1}"
        CMS_lumi.lumi_13TeV = "101 fb^{-1}"

    if(year=="2018"):
        normSignal=59.7/101
        CMS_lumi.lumi_13TeV = "" #2018 - 59.7 fb^{-1}"
    elif(year=="2017"):
        normSignal=41.5/101
        CMS_lumi.lumi_13TeV = "2017 - 41.5 fb^{-1}"
    elif(year=="wjets"): CMS_lumi.lumi_13TeV = "W+jets - 101 fb^{-1}"
    elif(year=="ttbar"): CMS_lumi.lumi_13TeV = "t#bar{t}+jets - 101 fb^{-1}"
    elif(year=="ttbarwjets"): CMS_lumi.lumi_13TeV = "t#bar{t}/W+jets - 101 fb^{-1}"
    elif(year=="zjets"): CMS_lumi.lumi_13TeV = "Z+jets - 101 fb^{-1}"


    m_Gl2400.Scale(normSignal)
    m_Gl1600.Scale(normSignal)
    m_Gl2000.Scale(normSignal)
    m_ppStau557.Scale(normSignal)
    m_ppStau871.Scale(normSignal)

    if(doRebin==True):
        m_Gl2400=m_Gl2400.Rebin(sizeRebinning,"Gl2400_new",rebinning)
        m_Gl1600=m_Gl1600.Rebin(sizeRebinning,"Gl1600_new",rebinning)
        m_Gl2000=m_Gl2000.Rebin(sizeRebinning,"Gl2000_new",rebinning)
        m_ppStau557=m_ppStau557.Rebin(sizeRebinning,"ppStau557_new",rebinning)
        m_ppStau871=m_ppStau871.Rebin(sizeRebinning,"ppStau871_new",rebinning)

    regionSyst=region
    if(region=="50ias60"): regionSyst="90ias100"
    if(region=="60ias70"): regionSyst="90ias100"
    if(region=="70ias80"): regionSyst="90ias100"
    if(region=="80ias90"): regionSyst="90ias100"
    if(region=="50ias90"): regionSyst="90ias100"
    if("fp" in region): regionSyst="90ias100"
    if(region=="8fp9" and not ("MET" in inputfile)): regionSyst="8fp9"
    yearSyst = year
    if(year=="2017_2018"): yearSyst="2018"


    ifileSyst=ROOT.TFile("systBckg/sysTotBinned_"+yearSyst+"_"+regionSyst+".root")
    histoOfSyst=ifileSyst.Get("systTotalBinned")

    pred_noCorrBias=pred.Clone()
    pred_noBlind=pred.Clone("_prednoBlind")
    obs_noBlind=obs.Clone("_obsnoBlind")

    biasCorr = False
    if(a_==0):
        blindMassUp(obs,300)
        blindMassUp(C_mass,300)
        blindMassUp(pred,300)
    else:
        if(biasCorr==True):
            print "Correction bias: ON"
            BiasCorrection(pred)
            BiasCorrection(pred_noBlind)
        else: print "Correction bias: OFF"
    
    if(a_<0): BiasCorrection(pred_noCorrBias)


    (pred,predD,predU)=addHSyst(pred,histoOfSyst,pred_noCorrBias)
    (pred_noBlind,pred_noBlindU,pred_noBlindD)=addHSyst(pred_noBlind,histoOfSyst,pred_noCorrBias)

    err_obs_m300=ROOT.Double(0)
    err_pred_m300=ROOT.Double(0)
    obs_m300 = obs_noBlind.IntegralAndError(obs_noBlind.FindBin(300),obs_noBlind.GetNbinsX()+1,err_obs_m300)
    pred_m300 = pred_noBlind.IntegralAndError(pred_noBlind.FindBin(300),pred_noBlind.GetNbinsX()+1,err_pred_m300)

    predD.SetName("predD")
    predU.SetName("predU")

    mass_fit=300
    if(a_==0):
        min_mass=0
        max_mass=300
    else:
        min_mass=0
        #max_mass=2000
        max_mass=4000
        if(doRebin==False):    
            max_mass=2500

    underflowAndOverflow(obs,True, max_mass)
    underflowAndOverflow(C_mass,True, max_mass)
    underflowAndOverflow(pred, False, max_mass) 
    underflowAndOverflow(pred_noSyst, False, max_mass) 

    h_syst_gl2400 = ROOT.TFile("systSignal/Gluino_M-2400/SR1_73p0/_sysTot.root").Get("c1_n4").GetPrimitive("")
    h_syst_gl1600 = ROOT.TFile("systSignal/Gluino_M-1600/SR1_73p0/_sysTot.root").Get("c1_n5").GetPrimitive("")
    h_syst_gl2000 = ROOT.TFile("systSignal/Gluino_M-2000/SR1_73p0/_sysTot.root").Get("c1_n9").GetPrimitive("")
    h_syst_ppStau557 = ROOT.TFile("systSignal/ppStau_M-557/SR1_73p0/_sysTot.root").Get("c1_n11").GetPrimitive("")
    h_syst_ppStau871 = ROOT.TFile("systSignal/ppStau_M-871/SR1_73p0/_sysTot.root").Get("c1_n14").GetPrimitive("")

    m_Gl2400 = addHSystSignal(m_Gl2400,h_syst_gl2400)
    m_Gl1600 = addHSystSignal(m_Gl1600,h_syst_gl1600)    
    m_Gl2000 = addHSystSignal(m_Gl2000,h_syst_gl2000)  
    m_ppStau557 = addHSystSignal(m_ppStau557,h_syst_ppStau557)    
    m_ppStau871 = addHSystSignal(m_ppStau871,h_syst_ppStau871)    

    underflowAndOverflow(m_Gl2400, False, max_mass)
    underflowAndOverflow(m_Gl1600, False, max_mass)
    underflowAndOverflow(m_Gl2000, False, max_mass)
    underflowAndOverflow(m_ppStau557, False, max_mass)
    underflowAndOverflow(m_ppStau871, False, max_mass)

    listOfMarkerGluino = [21, 22, 29, 23, 33, 34, 39, 47, 43]
    listOfMarkePPStau = [21, 22, 29, 23, 33, 34, 39, 47, 43, 44]

    m_Gl2400 = setColorAndMarker(m_Gl2400,46,listOfMarkerGluino[4])
    m_Gl1600 = setColorAndMarker(m_Gl1600,28,listOfMarkerGluino[2])
    m_Gl2000 = setColorAndMarker(m_Gl2000,46,listOfMarkerGluino[4])
    m_ppStau557 = setColorAndMarker(m_ppStau557,8,listOfMarkePPStau[3])
    m_ppStau871 = setColorAndMarker(m_ppStau871,31,listOfMarkePPStau[6])



    # Poisson errors for the observed distribution
    obs = poissonning(obs)
    
    KSTEST = obs.KolmogorovTest(pred)
    #print 'Kolmogorov test: ', KSTEST
    
    ratioSimpleH = ratioHisto(obs,pred)
    ratio_massC_obs = ratioHisto(C_mass,pred)

    pull = pullOfHisto(obs,pred,0.)
    pullC = pullOfHisto(C_mass,pred,0.)

    ratioInt = ratioIntegral(obs,pred,0.,-1)
    ratioIntC = ratioIntegral(C_mass,pred,0.,-1)

    if(isBinWidth):
        obs = binWidth(obs)
        pred = binWidth(pred)
        C_mass = binWidth(C_mass)
        pred_noSyst = binWidth(pred_noSyst)
        m_Gl2400 = binWidth(m_Gl2400)
        m_Gl1600 = binWidth(m_Gl1600)
        m_Gl2000 = binWidth(m_Gl2000)
        m_ppStau557 = binWidth(m_ppStau557)
        m_ppStau871 = binWidth(m_ppStau871)

    pred_band=pred.Clone()
    pred_band_noSyst=pred_noSyst.Clone()


    # -------------- Display --------------
       
    c1=TCanvas("c1","c1",700,700)
    t1=TPad("t1","t1", 0.0, 0.45, 0.95, 0.95)
    t1.Draw()
    t1.cd()
    t1.SetLogy(1)
    t1.SetTopMargin(0.003)
    t1.SetBottomMargin(0.005)
    c1.cd()

    t2=TPad("t2","t2", 0.0, 0.3, 0.95, 0.45)
    t2.Draw()
    t2.cd()
    t2.SetGridy(1)
    t2.SetPad(0,0.3,0.95,0.45)
    t2.SetTopMargin(0.1)
    t2.SetBottomMargin(0.02)
    c1.cd()
    
    t3=TPad("t3","t3", 0.0, 0.15, 0.95, 0.3)
    t3.Draw()
    t3.cd()
    t3.SetGridy(1)
    t3.SetPad(0,0.15,0.95,0.3)
    t3.SetTopMargin(0.1)
    t3.SetBottomMargin(0.02)
    c1.cd()

    t4=TPad("t4","t4", 0.0, 0.0, 0.95, 0.15)
    t4.Draw()
    t4.cd()
    t4.SetGridy(1)
    t4.SetPad(0,0.0,0.95,0.15)
    t4.SetTopMargin(0.1)
    t4.SetBottomMargin(0.3)
    t1.cd()

    min_entries=pred.GetBinContent(pred.FindBin(max_mass)-1)/5
    if(doRebin==False):
        min_entries=1e-6
    max_entries=pred.GetMaximum()*100

    min_entries = 1e-4
    max_entries = 5e6

    titleYaxis = "Events / bin"
    if (isBinWidth):
        titleYaxis = "Tracks / bin width"
    

    pred_band.GetXaxis().SetTitle("Mass (GeV)")
    pred_band.GetYaxis().SetTitle(titleYaxis)
    pred_band.GetYaxis().SetLabelFont(43)
    pred_band.GetYaxis().SetLabelSize(20)
    pred_band.GetYaxis().SetTitleFont(43)
    pred_band.GetYaxis().SetTitleSize(20)
    pred_band.GetYaxis().SetTitleOffset(2.)

    pred_band.SetMarkerStyle(22)
    pred_band.SetMarkerColor(5)
    pred_band.SetMarkerSize(1.0)
    pred_band.SetLineColor(5)
    pred_band.SetFillColor(5)
    pred_band.SetFillStyle(1001)
    pred_band.GetXaxis().SetRange(min_mass,max_mass)
    pred_band.GetXaxis().SetRangeUser(min_mass,max_mass)
    pred_band.GetYaxis().SetRangeUser(min_entries,max_entries)
    pred_band.GetXaxis().SetTitle("")
    pred_band.Draw("same E5")
    pred_band.SaveAs(odir+'/predband.root')

    pred_band_noSyst.GetXaxis().SetTitle("Mass (GeV)")
    pred_band_noSyst.GetYaxis().SetTitle(titleYaxis)
    pred_band_noSyst.GetYaxis().SetLabelFont(43)
    pred_band_noSyst.GetYaxis().SetLabelSize(20)
    pred_band_noSyst.GetYaxis().SetTitleFont(43)
    pred_band_noSyst.GetYaxis().SetTitleSize(20)
    pred_band_noSyst.GetYaxis().SetTitleOffset(7)

    pred_band_noSyst.SetMarkerStyle(22)
    pred_band_noSyst.SetMarkerColor(5)
    pred_band_noSyst.SetMarkerSize(0.1)
    pred_band_noSyst.SetLineColor(5)
    pred_band_noSyst.SetFillColor(5)
    pred_band_noSyst.SetFillStyle(1001)
    pred_band_noSyst.GetXaxis().SetRange(min_mass,max_mass)
    pred_band_noSyst.GetXaxis().SetRangeUser(min_mass,max_mass)
    pred_band_noSyst.GetYaxis().SetRangeUser(min_entries,max_entries)
    pred_band_noSyst.GetXaxis().SetTitle("")
    pred_band_noSyst.Draw("same E5")
    pred_band_noSyst.SaveAs(odir+'/predband_nosyst.root')

    pred.SetMarkerStyle(21)
    pred.SetMarkerColor(2)
    pred.SetMarkerSize(1)
    pred.SetLineColor(2)
    pred.SetFillColor(0)
    pred.Draw("same HIST P")
    pred.SaveAs(odir+'/pred.root')

    obs.SetMarkerStyle(20)
    obs.SetMarkerColor(1)
    obs.SetMarkerSize(1.0)
    obs.SetLineColor(1)
    obs.SetFillColor(0)
    obs.GetXaxis().SetRange(min_mass,max_mass)
    obs.GetXaxis().SetRangeUser(min_mass,max_mass)
    C_mass.SetMarkerColor(8)
    C_mass.SetLineColor(8)
    C_mass.SetMarkerStyle(23)
    if (region=="3fp8"): 
        obs.SetMarkerColor(8)
        obs.SetLineColor(8)
        obs.SetMarkerStyle(23)
    if(blind==False):
        obs.Draw("same E1")
        if (region == "8fp9"): m_Gl2400.Draw("same E1")
        if (region == "8fp9" and not ("NoC" in outputfile) and not ("MET" in inputfile)): C_mass.Draw("same E1")

    obs.SaveAs(odir+'/obs.root')
    
    if(signal==True):
        m_Gl1600.Draw("same E1")
        m_Gl2000.Draw("same E1")
        m_ppStau557.Draw("same E1")
        m_ppStau871.Draw("same E1")

    leg=TLegend(0.5,0.6,0.8,0.95)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(43)
    leg.SetTextSize(14)
    leg.SetHeader(labelRegion); 

    pred_leg=pred.Clone()
    pred_leg.SetFillColor(pred_band_noSyst.GetFillColor())
    pred_leg.SetFillStyle(pred_band_noSyst.GetFillStyle())

    #print("INTEGRAL PREDICTION SR3 = {}".format(pred.Integral()))
    if(blind==False and isData==False):
        
        if (region=="3fp8"): leg.AddEntry(obs,"Observed in C","PE1")
        else: leg.AddEntry(obs,"Observed","PE1")
        if (region == "8fp9" and not ("NoC" in outputfile) and not ("MET" in inputfile)): leg.AddEntry(C_mass,"Observed in C","PE1")
    elif(blind==False and isData==True):
        if (region=="3fp8"): leg.AddEntry(obs,"Observed in C","PE1")
        else: leg.AddEntry(obs,"Observed","PE1")
        if (region == "8fp9" and not ("NoC" in outputfile) and not ("MET" in inputfile)): leg.AddEntry(C_mass,"Observed in C","PE1")
    if(a_==0):
        if(biasCorr==True):
            if(isData==True):
                leg.AddEntry(pred_leg,"Data-based pred.","PE")
            else:
                leg.AddEntry(pred_leg,"Pred.","PF")
    else:
        if(isData==True):
            #leg.AddEntry(pred_leg,"Data-based pred. w/ bias corr.","PF")
            entry=leg.AddEntry(pred_leg,"Data-based pred.","PF")
            entry.SetFillColor(5)
            entry.SetFillStyle(1001)
            entry.SetLineColor(5)
            entry.SetLineStyle(1)
            entry.SetLineWidth(1)
            entry.SetMarkerColor(2)
            entry.SetMarkerStyle(21)
            entry.SetMarkerSize(1)
            entry.SetTextFont(43)
 
        else:
            leg.AddEntry(pred_leg,"Pred. w/ bias corr.","PF")

    if (region == "8fp9"): leg.AddEntry(m_Gl2400,"#tilde{g} (M=2400 GeV)","PE1")
    if(signal==True):
        leg.AddEntry(m_Gl1600,"#tilde{g} (M=1600 GeV)","PE1")
        leg.AddEntry(m_Gl2000,"#tilde{g} (M=2000 GeV)","PE1")
        leg.AddEntry(m_ppStau557,"pair. #tilde{#tau} (M=557 GeV)","PE1")
        leg.AddEntry(m_ppStau871,"pair. #tilde{#tau} (M=871 GeV)","PE1")


    LineLastBin=TLine(obs.GetBinLowEdge(obs.FindBin(max_mass)-1),0,obs.GetBinLowEdge(obs.FindBin(max_mass)-1),max_entries)
    LineLastBin.SetLineStyle(3)
    LineLastBin.SetLineColor(1)

    LineFit1=TLine(mass_fit,0,mass_fit,max_entries)
    LineFit1.SetLineStyle(1)
    LineFit1.SetLineColor(1)
    #LineFit1.Draw("same")

    leg.Draw("same")
    
    t=ROOT.TText(0.95,0.7,"+overflow")
    t.SetNDC(True)
    t.SetTextColor(1)
    t.SetTextFont(43)
    t.SetTextSize(24)
    t.SetTextAngle(90)


    LineAtOne=TLine(min_mass,1,max_mass,1)
    LineAtOne.SetLineStyle(3)
    LineAtOne.SetLineColor(1)
    #LineAtOne.Draw("same")



    c1.cd()
    t2.cd()
    
    frameR=ROOT.TH1D("frameR", "frameR", 1,min_mass, max_mass)
    frameR.GetXaxis().SetNdivisions(505)
    frameR.SetTitle("")
    frameR.SetStats(0)
    frameR.GetXaxis().SetTitle("")
    frameR.GetYaxis().SetTitle("RatioR     ")
    frameR.SetMaximum(2.)
    frameR.SetMinimum(0.0)
    frameR.GetYaxis().SetLabelFont(43) #give the font size in pixel (instead of fraction)
    frameR.GetYaxis().SetLabelSize(12) #font size
    frameR.GetYaxis().SetTitleFont(43) #give the font size in pixel (instead of fraction)
    frameR.GetYaxis().SetTitleSize(14) #font size
    frameR.GetYaxis().SetNdivisions(503)
    frameR.GetXaxis().SetNdivisions(505)
    frameR.GetXaxis().SetLabelFont(43) #give the font size in pixel (instead of fraction)
    frameR.GetXaxis().SetLabelSize(20) #font size
    frameR.GetXaxis().SetTitleFont(43) #give the font size in pixel (instead of fraction)
    frameR.GetXaxis().SetTitleSize(20) #font size
    frameR.GetXaxis().SetTitleOffset(3.75)
    frameR.Draw("AXIS")

    ratioInt.SetMarkerStyle(21)
    ratioInt.SetMarkerColor(1)
    ratioInt.SetMarkerSize(0.7)
    ratioInt.SetLineColor(1)
    ratioInt.SetFillColor(0)
    ratioIntC.SetMarkerStyle(23)
    ratioIntC.SetMarkerColor(8)
    ratioIntC.SetLineColor(8)
    if (region=="3fp8"):
        ratioInt.SetMarkerColor(8)
        ratioInt.SetLineColor(8)
        ratioInt.SetMarkerStyle(23)
    if(blind==False):
        ratioInt.Draw("same E0")
        if (region == "8fp9" and not ("NoC" in outputfile) and not ("MET" in inputfile)): ratioIntC.Draw("same E0")
    LineAtOne.Draw("same")


    #ratioInt.Draw("same E0")
    ratioInt.GetXaxis().SetRange(min_mass,max_mass)
    ratioInt.GetXaxis().SetRangeUser(min_mass,max_mass)



    c1.cd()
    t3.cd()
    
    frameR2=ROOT.TH1D("frameR2", "frameR2", 1,min_mass, max_mass)
    frameR2.GetXaxis().SetNdivisions(505)
    frameR2.SetTitle("")
    frameR2.SetStats(0)
    frameR2.GetXaxis().SetTitle("")
    frameR2.GetYaxis().SetTitle("obs / pred")
    frameR2.SetMaximum(2.)
    frameR2.SetMinimum(0.0)
    frameR2.GetYaxis().SetLabelFont(43) #give the font size in pixel (instead of fraction)
    frameR2.GetYaxis().SetLabelSize(12) #font size
    frameR2.GetYaxis().SetTitleFont(43) #give the font size in pixel (instead of fraction)
    frameR2.GetYaxis().SetTitleSize(14) #font size
    frameR2.GetYaxis().SetNdivisions(503)
    frameR2.GetXaxis().SetNdivisions(505)
    frameR2.GetXaxis().SetLabelFont(43) #give the font size in pixel (instead of fraction)
    frameR2.GetXaxis().SetLabelSize(20) #font size
    frameR2.GetXaxis().SetTitleFont(43) #give the font size in pixel (instead of fraction)
    frameR2.GetXaxis().SetTitleSize(20) #font size
    frameR2.GetXaxis().SetTitleOffset(3.75)
    frameR2.Draw("AXIS")

    ratioSimpleH.SetMarkerStyle(21)
    ratioSimpleH.SetMarkerColor(1)
    ratioSimpleH.SetMarkerSize(0.7)
    ratioSimpleH.SetLineColor(1)
    ratioSimpleH.SetFillColor(0)
    ratio_massC_obs.SetMarkerStyle(23)
    ratio_massC_obs.SetMarkerColor(8)
    ratio_massC_obs.SetLineColor(8)
    if (region=="3fp8"):
        ratioSimpleH.SetMarkerColor(8)
        ratioSimpleH.SetLineColor(8)
        ratioSimpleH.SetMarkerStyle(23)

    f=None
    if(blind==False):
        ratioSimpleH.Draw("same E0")
        if (region == "8fp9" and not ("NoC" in outputfile) and not ("MET" in inputfile)): ratio_massC_obs.Draw("same E0")
    #ratioSimpleH.Draw("same E0")
    if(a_==0):
        ratioSimpleH.Fit("pol1","RSQ0","same",50,250)
    ratioSimpleH.GetXaxis().SetRange(min_mass,max_mass)
    ratioSimpleH.GetXaxis().SetRangeUser(min_mass,max_mass)

    LineAtOne.Draw("same")

    LineFit3=TLine(mass_fit,0,mass_fit,2.)
    LineFit3.SetLineStyle(1)
    LineFit3.SetLineColor(1)
    #LineFit3.Draw("same")



    f=None
    if(a_==0):
        ratioSimpleH.Fit("pol1","RSQ0","same",50,250)
    ratioSimpleH.GetXaxis().SetRange(min_mass,max_mass)
    ratioSimpleH.GetXaxis().SetRangeUser(min_mass,max_mass)


    c1.cd()
    t4.cd()

    frameR3=ROOT.TH1D("frameR3", "frameR3", 1,min_mass, max_mass)
    frameR3.GetXaxis().SetNdivisions(505)
    frameR3.SetTitle("")
    frameR3.SetStats(0)
    frameR3.GetXaxis().SetTitle("")
    frameR3.GetXaxis().SetTitle("Mass (GeV)")
    frameR3.GetXaxis().SetTitleOffset(5)
    frameR3.GetYaxis().SetTitle("#frac{Data-pred}{#sigma}")
    frameR3.GetYaxis().SetTickLength(frameR3.GetYaxis().GetTickLength()*2)
    frameR3.SetMaximum(3)
    frameR3.SetMinimum(-3)
    frameR3.GetYaxis().SetLabelFont(43) #give the font size in pixel (instead of fraction)
    frameR3.GetYaxis().SetLabelSize(12) #font size
    frameR3.GetYaxis().SetTitleFont(43) #give the font size in pixel (instead of fraction)
    frameR3.GetYaxis().SetTitleSize(16) #font size
    frameR3.GetYaxis().SetNdivisions(503)
    frameR3.GetXaxis().SetNdivisions(510)
    frameR3.GetXaxis().SetLabelFont(43) #give the font size in pixel (instead of fraction)
    frameR3.GetXaxis().SetLabelSize(14) #font size
    frameR3.GetXaxis().SetTitleFont(43) #give the font size in pixel (instead of fraction)
    frameR3.GetXaxis().SetTitleSize(20) #font size
    frameR3.GetXaxis().SetTitleOffset(3.75)
    frameR3.Draw("AXIS")

    if(blind==False):
        pull.Draw("same HIST")
        if (region == "8fp9" and not ("NoC" in outputfile) and not ("MET" in inputfile)): pullC.Draw("same HIST")

    pull.SetLineColor(1)
    pull.SetFillColor(38)
    if (region=="3fp8"):
        pull.SetFillColorAlpha(8, 0.35)
        pull.SetLineColor(8)
    pullC.SetLineColor(8)
    pullC.SetFillColorAlpha(8, 0.35)
    #blindAnyUp(pull,300)
    t4.RedrawAxis()
    t4.RedrawAxis("G")

    
    LineAtZero=TLine(min_mass,0,max_mass,0)
    LineAtZero.SetLineStyle(1)
    LineAtZero.SetLineColor(1)
    LineAtZero.Draw("same")
    
    LineAt1p0=TLine(min_mass,1.0,max_mass,1.0)
    LineAt1p0.SetLineStyle(4)
    LineAt1p0.SetLineColor(1)
    LineAt1p0.Draw("same")
    
    LineAtMin1p0=TLine(min_mass,-1.0,max_mass,-1.0)
    LineAtMin1p0.SetLineStyle(4)
    LineAtMin1p0.SetLineColor(1)
    LineAtMin1p0.Draw("same")
    
    LineAt2p0=TLine(min_mass,2.0,max_mass,2.0)
    LineAt2p0.SetLineStyle(4)
    LineAt2p0.SetLineColor(1)
    LineAt2p0.Draw("same")
    
    LineAtMin2p0=TLine(min_mass,-2.0,max_mass,-2.0)
    LineAtMin2p0.SetLineStyle(4)
    LineAtMin2p0.SetLineColor(1)
    LineAtMin2p0.Draw("same")

    LineFit4=TLine(mass_fit,-3,mass_fit,3)
    LineFit4.SetLineStyle(1)
    LineFit4.SetLineColor(1)
    #LineFit4.Draw("same")


    CMS_lumi.CMS_lumi(c1, iPeriod, iPos)
    if(a_==0):
        c1.Update()
        c1.SaveAs(outputfile+"_region"+region+"_"+year+"_lowmass_"+VERSION+".root")
        c1.SaveAs(outputfile+"_region"+region+"_"+year+"_lowmass_"+VERSION+".C")
        c1.SaveAs(outputfile+"_region"+region+"_"+year+"_lowmass_"+VERSION+".pdf")
    else:
        c1.Update()
        c1.SaveAs(outputfile+"_region"+region+"_"+year+"_"+VERSION+".pdf")
        c1.SaveAs(outputfile+"_region"+region+"_"+year+"_"+VERSION+".root")
        c1.SaveAs(outputfile+"_region"+region+"_"+year+"_"+VERSION+".C")

    if(a_!=0):
        ratioSimpleH.Fit("pol1","QRS","",75,300)
    f=ratioSimpleH.GetFunction("pol1")

    #Chi2ObsPred = obs.Chi2Test(pred,"UWP")
    #print("Chi2 between prediction and observation  = {}".format(Chi2ObsPred))

    return (f.GetChisquare(),f.GetNDF(),f.GetParameter(1),f.GetParameter(0),f.GetParError(1),f.GetParError(0),labelRegion,region,obs_m300,pred_m300,err_obs_m300,err_pred_m300)




if __name__ == "__main__":
    print '----- First: bias correction'    
    (chi2,ndof,a,b,aerr,berr,region_,reg,obs_m300,pred_m300,err_obs_m300,err_pred_m300) = main(sys.argv[1:])

    print ''
    print '----- Second: full run'
    print 'a = ', a, ' and b = ', b
    a_ = a
    b_ = b
    odir = sys.argv[8]
    (chi2,ndof,a,b,aerr,berr,region_,reg,obs_m300,pred_m300,err_obs_m300,err_pred_m300) = main(sys.argv[1:])
