#!/usr/bin/python

import sys, getopt, os
sys.argv.append(' -b- ')
import ROOT
import math
import array
import numpy as np
import ctypes
sys.path.append("/opt/sbg/cms/safe1/cms/gcoulon")

from ROOT import THStack, TCanvas, TLegend, TLatex, TPad, TH1, TH2, TLine
import CMS_lumi, tdrstyle

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning + 1  # supprime les Info
ROOT.Math.MinimizerOptions.SetDefaultPrintLevel(-1)

tdrstyle.setTDRStyle()


year = '2024'
era = ''

PlotSignal = False


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
            if (debugprint): print ('Case where mass_max_Display < edge of histogram: ', mass_max_Display, ' ', h.GetBinCenter(h.GetNbinsX()) + h.GetBinWidth(h.GetNbinsX())/2)
            bin_content = 0
            bin_error = 0

            if (debugprint): print ('')
            if (debugprint): print ('       ', h.GetName())
            if (debugprint): print ('mass_max_Display: ', mass_max_Display, 'histo edge last bin: ', h.GetBinCenter(h.GetNbinsX()) + h.GetBinWidth(h.GetNbinsX())/2)
            for i in range (h.FindBin(mass_max_Display)-1, h.GetNbinsX()+1):
                bin_content += h.GetBinContent(i)
                bin_error += h.GetBinError(i)**2

                if (debugprint): print ('Bin #{}, mass {}: Points = {} +/- {}'.format(i,h.GetBinCenter(i),h.GetBinContent(i),h.GetBinError(i)))

                h.SetBinContent(i,0)
                h.SetBinError(i,0)
                if (debugprint): print ('Bin #{}, mass {}: Points = {} +/- {}'.format(i,h.GetBinCenter(i),h.GetBinContent(i),h.GetBinError(i)))

            h.SetBinContent(h.FindBin(mass_max_Display)-1, bin_content + h.GetBinContent(h.GetNbinsX()+1))
            if (data): h.SetBinError(h.FindBin(mass_max_Display)-1,math.sqrt(bin_error))
            else: h.SetBinError(h.FindBin(mass_max_Display)-1,math.sqrt(bin_error + h.GetBinError(h.GetNbinsX()+1)**2))

        elif (mass_max_Display == h.GetBinCenter(h.GetNbinsX()) + h.GetBinWidth(h.GetNbinsX())/2):
            if (debugprint): print ('Case where mass_max_Display == edge of histogram: ', mass_max_Display, ' ', h.GetBinCenter(h.GetNbinsX()) + h.GetBinWidth(h.GetNbinsX())/2)
            h.SetBinContent(h.GetNbinsX(),h.GetBinContent(h.GetNbinsX())+h.GetBinContent(h.GetNbinsX()+1))
            
            if(data): h.SetBinError(h.GetNbinsX(),math.sqrt(h.GetBinContent(h.GetNbinsX())))
            else: h.SetBinError(h.GetNbinsX(),math.sqrt(h.GetBinError(h.GetNbinsX())**2+h.GetBinError(h.GetNbinsX()+1)**2))
            
            h.SetBinContent(h.GetNbinsX()+1,0)
            h.SetBinError(h.GetNbinsX()+1,0)
        
        else: print ("Error: mass_max_Display > histo edge last bin")
    else:
        if (debugprint): print ('Case where mass_max_Display = ', mass_max_Display, ' useless to overflowInLastBin in this case')

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
    h3 = h1.Clone()
    h3.Sumw2()
    h2.Sumw2()
    h3.Divide(h2)

    return h3

def ratioIntegral(h1,h2,systErr,upTo=-1):
    h3 = h1.Clone()
    h3.Reset()
    if(upTo==-1):
        bornUp=h1.GetNbinsX()+1
    else:
        bornUp=h1.FindBin(upTo)
    for i in range(0,bornUp):
        e1 = ctypes.c_double(0.0)
        e2 = ctypes.c_double(0.0)

        if upTo == -1:
            a = h1.IntegralAndError(i, h1.GetNbinsX()+1, e1, "")
            b = h2.IntegralAndError(i, h1.GetNbinsX()+1, e2, "")
        else:
            a = h1.IntegralAndError(i, bornUp-1, e1, "")
            b = h2.IntegralAndError(i, bornUp-1, e2, "")

        e1 = e1.value
        e2 = e2.value
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
        
        if (Derr*Derr+Perr*Perr > 0): res.SetBinContent(i,(D-P)/math.sqrt(Derr*Derr+Perr*Perr))
        else: res.SetBinContent(i,0)

    return res

def addSyst(h,syst):
    res = h.Clone()
    res.Sumw2(0)
    for i in range (0,h.GetNbinsX()+1):
        res.SetBinError(i,math.sqrt(h.GetBinError(i)*h.GetBinError(i)+res.GetBinContent(i)*res.GetBinContent(i)*syst*syst))
    return res

def addHSyst(h, h_syst, hCorrBias):
    res = h.Clone()
    resD = h.Clone()
    resU = h.Clone()
    for i in range(0, h.GetNbinsX() + 1):
        syst = h_syst.GetBinContent(i) / 100
        j = i
        while j > 1 and h_syst.GetBinContent(j) == 0:
            j -= 1
            syst = h_syst.GetBinContent(j) / 100
        diffCorrBias = abs(hCorrBias.GetBinContent(i) - h.GetBinContent(i))
        errorTotal = math.sqrt(
            h.GetBinError(i)**2
            + res.GetBinContent(i)**2 * syst**2
            + diffCorrBias**2
        )
        res.SetBinError(i, errorTotal)
        resD.SetBinContent(i, res.GetBinContent(i) - errorTotal)
        resU.SetBinContent(i, res.GetBinContent(i) + errorTotal)
    return (res, resD, resU)

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

def poissonHisto(h,RNG):
    for i in range (0,h.GetNbinsX()):
        h.SetBinContent(i,RNG.Poisson(h.GetBinContent(i)))

def PE_Pred(obs,h,nPE):
    h_chi2=ROOT.TH1F("chi2",";#chi^{2};",100,0,200)
    h_KS=ROOT.TH1F("KS",";Kolmogorov-Smirnov test;",100,0,1e-1)
    RNG=ROOT.TRandom3()
    for i in range(nPE):
        poissonHisto(h,RNG)
        chi2=h.Chi2Test(obs,"CHI2/NDF")
        KS=h.KolmogorovTest(obs,"M")
        h_chi2.Fill(chi2)
        h_KS.Fill(KS)
    return h_chi2, h_KS


#----------------------------------------------------
#                       Main
#----------------------------------------------------

def main(argv):
    # -------------- Setup --------------
    outputfile = ''
    region = ''
    odir = ''
    labelName = ''
    nominalOnly = True

    try:
        opts, args = getopt.getopt(argv,"hi:o:r:d",["ifile=","labelName=","ofile=","region=","odir=", "nom="])
    except getopt.GetoptError:
        print ('test.py -i <inputfile> -e <labelName> -o <outputfile> -r <region> -d <odir> -n <nominalOnly>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('test.py -i <inputfile> -e <labelName> -o <outputfile> -r <region> -d <odir> -n <nominalOnly>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-e", "--labelName"):
            labelName = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
        elif opt in ("-r", "--region"):
            region = arg
        elif opt in ("-d", "--odir"):
            odir = arg
        elif opt in ("-n", "--nom"):
            nominalOnly = arg.lower() in ("true", "1", "yes")

    os.system('mkdir -p ' + odir)
    outputfile = odir + '/' + outputfile

    print (' Input file: ', inputfile)
    print ('Output file: ', outputfile)
    print ('     Region: ', region)

    blind       = False
    isBinWidth  = False
    doRebin     = True
    isData      = True
    labelRegion = region
    signal      = False
    option      = "_etaAbs_chi2cut"


    ifile = ROOT.TFile(inputfile)
    
    obs = ifile.Get("mass_obs_" + region)
    pred = ifile.Get("mass_predBC_" + region)
    
    if (region=="8fp9"): C_mass = ifile.Get("mass_regionC_3fp8_" + labelName)
    else: C_mass = ifile.Get("mass_obs_" + region)
    pred_noSyst = addSyst(pred,0.0)


    ifileGl2000 = ROOT.TFile("/opt/sbg/cms/safe1/cms/gcoulon/CMSSW_10_6_30/src/HSCPTreeAnalyzer/output/Gluino2000_massCut_0_pT70_V5p0_Fpix_Eta2p4_Scale.root")
    m_Gl2000 = ifileGl2000.Get("mass_regionD_"+region+"_METContainingMu")


    # -------------- Work on histograms --------------
    rebinning=array.array('d',[0.,20.,40.,60.,80.,100.,120.,140.,160.,180.,200.,220.,240.,260.,280.,300.,320.,340.,360.,380.,410.,440.,480.,530.,590.,660.,760.,880.,1030.,1210.,1440.,1730.,2000.,2500.,3200.,4000.])


    sizeRebinning = len(rebinning)-1
    
    if(doRebin==True):
        pred        = pred.Rebin(sizeRebinning,"pred_new",rebinning)
        pred_noSyst = pred_noSyst.Rebin(sizeRebinning,"pred_noSyst_new",rebinning)
        obs         = obs.Rebin(sizeRebinning,"obs_new",rebinning)
        C_mass      = C_mass.Rebin(sizeRebinning,"C_mass_new",rebinning)

    normSignal = 1
    if(year=="2024"):
        if (era=="F"): normSignal = 25.40/100
        elif (era=="G"): normSignal = 34.4/100
    
    if (PlotSignal):
        m_Gl2000.Scale(normSignal)

    if(doRebin==True):
        if (PlotSignal):
            m_Gl2000 = m_Gl2000.Rebin(sizeRebinning,"Gl2000_new",rebinning)



    ifileSyst = None
    if (region=="8fp9" and not ("MET" in inputfile)):
        print("faire les syst d'abord !")
    if (("MET" in inputfile) and ("/Eta2p4/" in inputfile)):
        ifileSyst = ROOT.TFile("/opt/sbg/cms/safe1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/macros/DataMET_2024_V12p24__" + region + "/Eta2p4/SystCombined/sysTotBinned_2024_" + region + ".root")
    elif (("MET" in inputfile) and ("/Eta1_2p4/" in inputfile)):
        ifileSyst = ROOT.TFile("/opt/sbg/cms/safe1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/macros/DataMET_2024_V12p24__" + region + "/Eta1_2p4/SystCombined/sysTotBinned_2024_" + region + ".root")
    elif (("MET" in inputfile) and ("/Eta1/" in inputfile)):
        ifileSyst = ROOT.TFile("/opt/sbg/cms/safe1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/macros/DataMET_2024_V12p24__" + region + "/Eta1/SystCombined/sysTotBinned_2024_" + region + ".root")

    if ifileSyst is None:
        raise RuntimeError(f"Aucun fichier syst trouvé pour region='{region}', inputfile='{inputfile}'")

    
    histoOfSyst = ifileSyst.Get("systTotalBinned")

    pred_noCorrBias = pred.Clone()
    pred_noBlind = pred.Clone("_prednoBlind")
    obs_noBlind = obs.Clone("_obsnoBlind")


    if (not nominalOnly):
        (pred, predD, predU) = addHSyst(pred, histoOfSyst, pred_noCorrBias)
        (pred_noBlind, pred_noBlindU, pred_noBlindD) = addHSyst(pred_noBlind, histoOfSyst, pred_noCorrBias)
        print(" syst. file: ", ifileSyst.GetName())
    else:
        print(" /!\ only nominal")

        ifileSystnom = None
        if (("MET" in inputfile) and ("/Eta2p4/" in inputfile)):
            ifileSystnom = ROOT.TFile("/opt/sbg/cms/safe1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/macros/DataMET_2024_V12p24__" + region + option + "/Eta2p4/SystCombined/sysToTBinned_2024_" + region + ".root")
        elif (("MET" in inputfile) and ("/Eta1_2p4/" in inputfile)):
            ifileSystnom = ROOT.TFile("/opt/sbg/cms/safe1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/macros/DataMET_2024_V12p24__" + region + option + "/Eta1_2p4/SystCombined/sysToTBinned_2024_" + region + ".root")
        elif (("MET" in inputfile) and ("/Eta1/" in inputfile)):
            ifileSystnom = ROOT.TFile("/opt/sbg/cms/safe1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/macros/DataMET_2024_V12p24__" + region + option + "/Eta1/SystCombined/sysToTBinned_2024_" + region + ".root")
        print(" syst. file: ", ifileSystnom.GetName())

        histoOfSystnom = ifileSystnom.Get("Stat")
        (pred, predD, predU) = addHSyst(pred, histoOfSystnom, pred_noCorrBias)
        (pred_noBlind, pred_noBlindU, pred_noBlindD) = addHSyst(pred_noBlind, histoOfSystnom, pred_noCorrBias)

    err_obs_m300 = ctypes.c_double(0)
    obs_m300 = obs_noBlind.IntegralAndError(obs_noBlind.FindBin(300), obs_noBlind.GetNbinsX()+1, err_obs_m300)
    err_obs_m300 = err_obs_m300.value

    err_pred_m300 = ctypes.c_double(0)
    pred_m300 = pred_noBlind.IntegralAndError(pred_noBlind.FindBin(300), pred_noBlind.GetNbinsX()+1, err_pred_m300)
    err_pred_m300 = err_pred_m300.value


    mass_fit=4000
    min_mass=0
    max_mass=4000
    if(doRebin==False):    
        max_mass=2500

    underflowAndOverflow(obs,True, max_mass)
    underflowAndOverflow(C_mass,True, max_mass)
    underflowAndOverflow(pred, False, max_mass) 
    underflowAndOverflow(pred_noSyst, False, max_mass) 


    if (PlotSignal):
        underflowAndOverflow(m_Gl2000, False, max_mass)
        m_Gl2000 = setColorAndMarker(m_Gl2000, 28, 22)


    # Poisson errors for the observed distribution
    obs = poissonning(obs)
    



    # Stat. tests
    # h_ch2 = ROOT.TH1F("chi2",";#chi^{2};",100,0,200)
    # h_kolmo = ROOT.TH1F("KS",";Kolmogorov-Smirnov test;",100,0,1e-1)
    # h_ch2, h_kolmo = PE_Pred(obs, pred, 10000)

    # KSTEST = obs.KolmogorovTest(pred)
    # print('Kolmogorov test: ', KSTEST)

    # pull = pullOfHisto(obs,pred,0.)
    # empty = pull.Clone("empty")
    # for i in range(pull.GetNbinsX()):
    #     empty.SetBinContent(i,0)
    #     empty.SetBinError(i,0)

    # Chi2ObsPred = obs.Chi2Test(pred,"UWP")
    # print("Chi2 between prediction and observation  = {}".format(Chi2ObsPred))

    # Chi2OfPull =  pull.Chi2Test(empty, "UU")
    # print("Chi2 of the pull plot = {}".format(Chi2OfPull))
    #
    


    
    ratioSimpleH    = ratioHisto(obs,pred)
    ratio_massC_obs = ratioHisto(C_mass,pred)

    pull  = pullOfHisto(obs,pred,0.)
    pullC = pullOfHisto(C_mass,pred,0.)

    ratioInt  = ratioIntegral(obs,    pred, 0., mass_fit)
    ratioIntC = ratioIntegral(C_mass, pred, 0., mass_fit)

    if(isBinWidth):
        obs = binWidth(obs)
        pred = binWidth(pred)
        C_mass = binWidth(C_mass)
        pred_noSyst = binWidth(pred_noSyst)
        if (PlotSignal):
            m_Gl2000 = binWidth(m_Gl2000)

    pred_band=pred.Clone()
    pred_band_noSyst=pred_noSyst.Clone()


    if blind:
        obs_blind       = obs.Clone("obs_blind")
        C_mass_blind    = C_mass.Clone("C_mass_blind")
        ratioInt_blind  = ratioInt.Clone("ratioInt_blind")
        ratioIntC_blind = ratioIntC.Clone("ratioIntC_blind")
        ratioSimpleH_blind   = ratioSimpleH.Clone("ratioSimpleH_blind")
        ratio_massC_blind    = ratio_massC_obs.Clone("ratio_massC_blind")
        pull_blind  = pull.Clone("pull_blind")
        pullC_blind = pullC.Clone("pullC_blind")
        for h in [obs_blind, C_mass_blind, ratioInt_blind, ratioIntC_blind,
                  ratioSimpleH_blind, ratio_massC_blind, pull_blind, pullC_blind]:
            blindAnyUp(h, mass_fit)
    else:
        obs_blind, C_mass_blind        = obs, C_mass
        ratioInt_blind, ratioIntC_blind = ratioInt, ratioIntC
        ratioSimpleH_blind, ratio_massC_blind = ratioSimpleH, ratio_massC_obs
        pull_blind, pullC_blind        = pull, pullC

    # -------------- Display --------------
       
    c1=TCanvas("c1","c1",700,700)
    t1=TPad("t1","t1", 0.0, 0.45, 0.95, 0.95)
    t1.Draw()
    t1.cd()
    t1.SetLogy(1)
    t1.SetTopMargin(0.003)
    t1.SetBottomMargin(0.005)
    c1.cd()

    t2=TPad("t2","t2", 0.0, 0.32, 0.95, 0.45)
    t2.Draw()
    t2.cd()
    t2.SetGridy(1)
    t2.SetPad(0,0.32,0.95,0.45)
    t2.SetTopMargin(0.1)
    t2.SetBottomMargin(0.02)
    c1.cd()
    
    t3=TPad("t3","t3", 0.0, 0.18, 0.95, 0.32)
    t3.Draw()
    t3.cd()
    t3.SetGridy(1)
    t3.SetPad(0,0.18,0.95,0.32)
    t3.SetTopMargin(0.1)
    t3.SetBottomMargin(0.02)
    c1.cd()

    t4=TPad("t4","t4", 0.0, 0.0, 0.95, 0.18)
    t4.Draw()
    t4.cd()
    t4.SetGridy(1)
    t4.SetPad(0,0.0,0.95,0.18)
    t4.SetTopMargin(0.1)
    t4.SetBottomMargin(0.3)
    t1.cd()

    min_entries=pred.GetBinContent(pred.FindBin(max_mass)-1)/10
    if(doRebin==False):
        min_entries=1e-6
    max_entries=pred.GetMaximum()*100

    #min_entries = 1e-4
    #max_entries = 5e6

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

    obs_blind.SetMarkerStyle(20)
    obs_blind.SetMarkerColor(1)
    obs_blind.SetMarkerSize(1.0)
    obs_blind.SetLineColor(1)
    obs_blind.SetFillColor(0)
    obs_blind.GetXaxis().SetRange(min_mass,max_mass)
    obs_blind.GetXaxis().SetRangeUser(min_mass,max_mass)
    C_mass_blind.SetMarkerColor(8)
    C_mass_blind.SetLineColor(8)
    C_mass_blind.SetMarkerStyle(23)
    if (region=="3fp8"): 
        obs_blind.SetMarkerColor(8)
        obs_blind.SetLineColor(8)
        obs_blind.SetMarkerStyle(23)
    obs_blind.Draw("same E1")
    if (region == "8fp9" and PlotSignal):
        m_Gl2000.Draw("same E1")
    if (region == "8fp9" and not ("NoC" in outputfile) and not ("MET" in inputfile)):
        C_mass_blind.Draw("same E1")

    obs_blind.SaveAs(odir+'/obs_blind.root')
    
    if(signal==True):
        m_Gl2000.Draw("same E1")

    leg=TLegend(0.5,0.75,0.8,0.95)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(43)
    leg.SetTextSize(14)
    if (labelRegion == "8fp9"): leg.SetHeader("CR : " + labelRegion)
    else: leg.SetHeader("Region : " + labelRegion)

    tex1 = ROOT.TLatex(0.85, 0.96, "(13.6 TeV)")
    if (year == "2024"):
        tex1 = ROOT.TLatex(0.66, 0.96, "105.8 fb^{-1} (13.6 TeV)")
        if (era == "F"): tex1 = ROOT.TLatex(0.63, 0.96, "2024F - 25.40 fb^{-1} (13.6 TeV)")
        if (era == "G"): tex1 = ROOT.TLatex(0.63, 0.96, "2024G - 34.4 fb^{-1} (13.6 TeV)")
    tex1.SetNDC()
    tex1.SetTextFont(42)
    tex1.SetLineWidth(2)
    tex1.SetTextSize(0.03)
    c1.cd()
    tex1.Draw()

    tex2 = ROOT.TLatex(0.15, 0.96, "#scale[1.3]{#bf{CMS}}#it{Simulation Work in progress}")
    tex2.SetNDC()
    tex2.SetTextFont(42)
    tex2.SetTextSize(0.03)
    tex2.SetLineWidth(2)
    c1.cd()
    tex2.Draw()


    pred_leg = pred.Clone()
    pred_leg.SetFillColor(pred_band_noSyst.GetFillColor())
    pred_leg.SetFillStyle(pred_band_noSyst.GetFillStyle())


    if (region=="3fp8"): leg.AddEntry(obs_blind, "Observed in C", "PE1")
    else: leg.AddEntry(obs_blind, "Observed", "PE1")
    if (region == "8fp9" and not ("NoC" in outputfile) and not ("MET" in inputfile)):
        leg.AddEntry(C_mass_blind, "Observed in C", "PE1")
    
    if(isData==True):
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
 

    if (region == "8fp9" and PlotSignal):
        leg.AddEntry(m_Gl2000,"#tilde{g} (M=2000 GeV)","PE1")
    if(signal==True):
        leg.AddEntry(m_Gl2000,"#tilde{g} (M=2000 GeV)","PE1")


    LineLastBin=TLine(obs_blind.GetBinLowEdge(obs_blind.FindBin(max_mass)-1),0,obs_blind.GetBinLowEdge(obs_blind.FindBin(max_mass)-1),max_entries)
    LineLastBin.SetLineStyle(3)
    LineLastBin.SetLineColor(1)

    LineFit1=TLine(mass_fit,0,mass_fit,max_entries)
    LineFit1.SetLineStyle(1)
    LineFit1.SetLineColor(1)
    
    t1.cd()
    if (blind): LineFit1.Draw("same")
    leg.Draw("same")
    
    t=ROOT.TText(0.95,0.7,"+overflow")
    t.SetNDC(True)
    t.SetTextColor(1)
    t.SetTextFont(43)
    t.SetTextSize(24)
    t.SetTextAngle(90)



    c1.cd()
    t2.cd()
    
    frameR=ROOT.TH1D("frameR", "frameR", 1,min_mass, max_mass)
    frameR.GetXaxis().SetNdivisions(505)
    frameR.SetTitle("")
    frameR.SetStats(0)
    frameR.GetXaxis().SetTitle("")
    frameR.GetYaxis().SetTitle("RatioR  ")
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

    ratioInt_blind.SetMarkerStyle(21)
    ratioInt_blind.SetMarkerColor(1)
    ratioInt_blind.SetMarkerSize(0.7)
    ratioInt_blind.SetLineColor(1)
    ratioInt_blind.SetFillColor(0)
    ratioIntC_blind.SetMarkerStyle(23)
    ratioIntC_blind.SetMarkerColor(8)
    ratioIntC_blind.SetLineColor(8)
    if (region=="3fp8"):
        ratioInt_blind.SetMarkerColor(8)
        ratioInt_blind.SetLineColor(8)
        ratioInt_blind.SetMarkerStyle(23)
    ratioInt_blind.Draw("same E0")
    if (region == "8fp9" and not ("NoC" in outputfile) and not ("MET" in inputfile)):
        ratioIntC_blind.Draw("same E0")

    LineAtOne=TLine(min_mass,1,max_mass,1)
    LineAtOne.SetLineStyle(3)
    LineAtOne.SetLineColor(1)
    LineAtOne.Draw("same")

    ratioInt_blind.GetXaxis().SetRange(min_mass,max_mass)
    ratioInt_blind.GetXaxis().SetRangeUser(min_mass,max_mass)

    LineFit2 = TLine(mass_fit, 0, mass_fit, 2.)
    LineFit2.SetLineStyle(1)
    LineFit2.SetLineColor(1)
    if (blind): LineFit2.Draw("same")



    c1.cd()
    t3.cd()
    
    frameR2=ROOT.TH1D("frameR2", "frameR2", 1,min_mass, max_mass)
    frameR2.GetXaxis().SetNdivisions(505)
    frameR2.SetTitle("")
    frameR2.SetStats(0)
    frameR2.GetXaxis().SetTitle("")
    frameR2.GetYaxis().SetTitle("obs / pred ")
    frameR2.GetYaxis().SetRangeUser(0.,2.)
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

    ratioSimpleH_blind.Sumw2()
    ratioSimpleH_blind.SetMarkerStyle(21)
    ratioSimpleH_blind.SetMarkerColor(1)
    ratioSimpleH_blind.SetMarkerSize(0.7)
    ratioSimpleH_blind.SetLineColor(1)
    ratioSimpleH_blind.SetFillColor(0)
    ratio_massC_blind.SetMarkerStyle(23)
    ratio_massC_blind.SetMarkerColor(8)
    ratio_massC_blind.SetLineColor(8)
    if (region=="3fp8"):
        ratioSimpleH_blind.SetMarkerColor(8)
        ratioSimpleH_blind.SetLineColor(8)
        ratioSimpleH_blind.SetMarkerStyle(23)


    ratioSimpleH_blind.Draw("same E0")
    if (region == "8fp9" and not ("NoC" in outputfile) and not ("MET" in inputfile)):
        ratio_massC_blind.Draw("same E0")
    #ratioSimpleH_blind.Draw("same E0")
    ratioSimpleH_blind.GetXaxis().SetRange(min_mass,max_mass)
    ratioSimpleH_blind.GetXaxis().SetRangeUser(min_mass,max_mass)

    LineAtOne.Draw("same")

    LineFit3=TLine(mass_fit,0,mass_fit,2.)
    LineFit3.SetLineStyle(1)
    LineFit3.SetLineColor(1)
    if (blind): LineFit3.Draw("same")


    ratioSimpleH_blind.GetXaxis().SetRange(min_mass,max_mass)
    ratioSimpleH_blind.GetXaxis().SetRangeUser(min_mass,max_mass)


    c1.cd()
    t4.cd()

    frameR3=ROOT.TH1D("frameR3", "frameR3", 1,min_mass, max_mass)
    frameR3.GetXaxis().SetNdivisions(505)
    frameR3.SetTitle("")
    frameR3.SetStats(0)
    frameR3.GetXaxis().SetTitle("")
    frameR3.GetXaxis().SetTitle("Mass (GeV)")
    frameR3.GetYaxis().SetTitleOffset(1.3)
    frameR3.GetYaxis().SetTitle("#frac{Data-pred}{#sigma} ")
    frameR3.GetYaxis().SetTickLength(frameR3.GetYaxis().GetTickLength()*2)
    frameR3.SetMaximum(3)
    frameR3.SetMinimum(-3)
    frameR3.GetYaxis().SetLabelFont(43) #give the font size in pixel (instead of fraction)
    frameR3.GetYaxis().SetLabelSize(12) #font size
    frameR3.GetYaxis().SetTitleFont(43) #give the font size in pixel (instead of fraction)
    frameR3.GetYaxis().SetTitleSize(14) #font size
    frameR3.GetYaxis().SetNdivisions(503)
    frameR3.GetXaxis().SetNdivisions(510)
    frameR3.GetXaxis().SetLabelFont(43) #give the font size in pixel (instead of fraction)
    frameR3.GetXaxis().SetLabelSize(14) #font size
    frameR3.GetXaxis().SetTitleFont(43) #give the font size in pixel (instead of fraction)
    frameR3.GetXaxis().SetTitleSize(15) #font size
    frameR3.GetXaxis().SetTitleOffset(.85)
    frameR3.Draw("AXIS")

    pull_blind.Draw("same HIST")
    if (region == "8fp9" and not ("NoC" in outputfile) and not ("MET" in inputfile)):
        pullC_blind.Draw("same HIST")

    pull_blind.SetLineColor(1)
    pull_blind.SetFillColor(38)
    if (region=="3fp8"):
        pull_blind.SetFillColorAlpha(8, 0.35)
        pull_blind.SetLineColor(8)
    pullC_blind.SetLineColor(8)
    pullC_blind.SetFillColorAlpha(8, 0.35)
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
    if (blind): LineFit4.Draw("same")


    
    c1.Update()
    c1.SaveAs(outputfile + "_region" + region + "_" + year + "_" + ("onlyNominal" if nominalOnly else "") + ".pdf")
    c1.SaveAs(outputfile + "_region" + region + "_" + year + "_" + ("onlyNominal" if nominalOnly else "") + ".root")
    c1.SaveAs(outputfile + "_region" + region + "_" + year + "_" + ("onlyNominal" if nominalOnly else "") + ".C")

    print("   Saved in: {}".format(odir))
    print('')

    #Chi2ObsPred = obs.Chi2Test(pred,"UWP")
    #print("Chi2 between prediction and observation  = {}".format(Chi2ObsPred))

    return (region, obs_m300, pred_m300, err_obs_m300, err_pred_m300)




if __name__ == "__main__":

    odir = sys.argv[8]
    (reg, obs_m300, pred_m300, err_obs_m300, err_pred_m300) = main(sys.argv[1:])
