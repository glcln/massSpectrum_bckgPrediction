#!/usr/bin/python

# =============================================================================
#  Final mass-spectrum plot: observation vs data-driven prediction.
# -----------------------------------------------------------------------------
#  Reads one step2 output file and draws a four-pad figure:
#    t1 : the spectra (log y)  - prediction with its uncertainty band, observed
#         points, optional MC stack, optional gluino signal overlays
#    t2 : cumulative ratio (only when doYouWantRratio is on)
#    t3 : bin-by-bin ratio obs / pred
#    t4 : pull, (obs - pred) / sigma
#
#  The uncertainty band is built either from the combined-systematics file
#  produced by systBckg.py (--systfile) or, failing that, from the statistical
#  uncertainty of the prediction alone (--nom true).
#
#  In the search region (9fp10) everything above mass_fit = 300 GeV is blinded
#  in the *drawn* histograms only; the unblinded copies are kept so the yields
#  above 300 GeV can still be returned to the caller.
#
#  Invoked by ShowPlots.py, one call per eta range.
# =============================================================================

import sys, getopt, os
import ROOT
import math
import array
import ctypes
sys.path.append("/safe/ui3_1/cms/gcoulon")

from ROOT import THStack, TCanvas, TLegend, TPad, TH1, TLine
import tdrstyle

ROOT.gROOT.SetBatch(True)                   # no X11: the script only writes files
ROOT.gErrorIgnoreLevel = ROOT.kWarning + 1  # drops the Info messages
ROOT.Math.MinimizerOptions.SetDefaultPrintLevel(-1)

tdrstyle.setTDRStyle()


#----------------------------------------------------
#                       Functions
#----------------------------------------------------

# Command-line booleans arrive as strings ("True", "1", ...) since they are
# passed through a shell by ShowPlots.py.
def asBool(s):
    return str(s).strip().lower() in ("true", "1", "yes", "y", "on")

# Mutates and returns its argument, so callers that need the original style
# preserved must pass a clone.
def setColorAndMarker(h1,color,markerstyle):
    h1.SetLineColor(color)
    h1.SetMarkerColor(color)
    h1.SetMarkerSize(1.2)
    h1.SetMarkerStyle(markerstyle)
    return h1

# Rebuilds the observed spectrum by re-filling it entry by entry, so that ROOT
# computes proper Poisson (Garwood) asymmetric error bars instead of sqrt(N).
# That matters in the tail, where a bin with 1 or 2 events has a very asymmetric
# interval.
#
# Note that int(GetBinContent(i)) truncates: this only makes sense on an
# unweighted, integer-content histogram, i.e. on data.
def poissonning(h):
    res = h.Clone()
    res.Reset()
    res.Sumw2(0)                                  # drop the weight array...
    res.SetBinErrorOption(ROOT.TH1.kPoisson)      # ...so the Poisson option applies

    for i in range (0,h.GetNbinsX()+1):
        for j in range(0,int(h.GetBinContent(i))):
            res.Fill(h.GetBinCenter(i))
    
    return res

# Folds everything above mass_max_Display into the last displayed bin, so that
# no entry silently disappears off the right edge of the plot.
#
# Three cases, depending on where the display limit sits relative to the
# histogram's own last edge:
#   - display limit strictly inside the histogram -> sum the bins beyond it,
#     including the overflow, into the last displayed bin and clear the rest;
#   - display limit exactly on the last edge -> the usual overflow fold;
#   - display limit beyond the histogram -> that is a configuration error.
#
# `data` selects the error convention: for data the error of the merged bin is
# recomputed from the summed content (Poisson), for the prediction the errors
# are added in quadrature.
def overflowInLastBin(h, data, mass_max_Display):
    debugprint = False
    # Below 300 GeV nothing is ever folded: the blinding threshold already cuts
    # the spectrum there.
    if (mass_max_Display > 300):
        if(mass_max_Display < h.GetBinCenter(h.GetNbinsX()) + h.GetBinWidth(h.GetNbinsX())/2):
            if (debugprint): print ('Case where mass_max_Display < edge of histogram: ', mass_max_Display, ' ', h.GetBinCenter(h.GetNbinsX()) + h.GetBinWidth(h.GetNbinsX())/2)
            bin_content = 0
            bin_error = 0

            if (debugprint): print ('')
            if (debugprint): print ('       ', h.GetName())
            if (debugprint): print ('mass_max_Display: ', mass_max_Display, 'histo edge last bin: ', h.GetBinCenter(h.GetNbinsX()) + h.GetBinWidth(h.GetNbinsX())/2)
            # Accumulate then clear every bin from the display limit onwards.
            # The -1 targets the bin *containing* the limit, which becomes the
            # last visible one.
            for i in range (h.FindBin(mass_max_Display)-1, h.GetNbinsX()+1):
                bin_content += h.GetBinContent(i)
                bin_error += h.GetBinError(i)**2

                if (debugprint): print ('Bin #{}, mass {}: Points = {} +/- {}'.format(i,h.GetBinCenter(i),h.GetBinContent(i),h.GetBinError(i)))

                h.SetBinContent(i,0)
                h.SetBinError(i,0)
                if (debugprint): print ('Bin #{}, mass {}: Points = {} +/- {}'.format(i,h.GetBinCenter(i),h.GetBinContent(i),h.GetBinError(i)))

            # The overflow is added on top of the accumulated content. Its error
            # only enters for the prediction: on data the error is left as the
            # quadratic sum of the visible bins.
            h.SetBinContent(h.FindBin(mass_max_Display)-1, bin_content + h.GetBinContent(h.GetNbinsX()+1))
            if (data): h.SetBinError(h.FindBin(mass_max_Display)-1,math.sqrt(bin_error))
            else: h.SetBinError(h.FindBin(mass_max_Display)-1,math.sqrt(bin_error + h.GetBinError(h.GetNbinsX()+1)**2))

        elif (mass_max_Display == h.GetBinCenter(h.GetNbinsX()) + h.GetBinWidth(h.GetNbinsX())/2):
            if (debugprint): print ('Case where mass_max_Display == edge of histogram: ', mass_max_Display, ' ', h.GetBinCenter(h.GetNbinsX()) + h.GetBinWidth(h.GetNbinsX())/2)
            # Plain overflow fold into the last bin.
            h.SetBinContent(h.GetNbinsX(),h.GetBinContent(h.GetNbinsX())+h.GetBinContent(h.GetNbinsX()+1))
            
            # Data: Poisson error from the merged content.
            # Prediction: errors added in quadrature.
            if(data): h.SetBinError(h.GetNbinsX(),math.sqrt(h.GetBinContent(h.GetNbinsX())))
            else: h.SetBinError(h.GetNbinsX(),math.sqrt(h.GetBinError(h.GetNbinsX())**2+h.GetBinError(h.GetNbinsX()+1)**2))
            
            h.SetBinContent(h.GetNbinsX()+1,0)
            h.SetBinError(h.GetNbinsX()+1,0)
        
        else: print ("Error: mass_max_Display > histo edge last bin")
    else:
        if (debugprint): print ('Case where mass_max_Display = ', mass_max_Display, ' useless to overflowInLastBin in this case')

# Single entry point for the under/overflow handling. The underflow fold is
# currently disabled: the first mass bin starts at 0, so there is nothing below.
def underflowAndOverflow(h, data, mass_max_Display):
    #underflowInFirstBin(h,data)
    overflowInLastBin(h, data, mass_max_Display)

# Divides every bin by its width, turning yields into a density. Needed only
# when isBinWidth is on, since the mass binning is strongly non-uniform and the
# raw yields then exaggerate the wide bins.
def binWidth(h1):
    res = h1.Clone()
    for i in range (0,h1.GetNbinsX()+1):
        res.SetBinContent(i,h1.GetBinContent(i)/h1.GetBinWidth(i))
        res.SetBinError(i,h1.GetBinError(i)/h1.GetBinWidth(i))
    return res

# Bin-by-bin ratio. Sumw2 is forced on both so ROOT propagates the errors
# instead of falling back to sqrt(N) on the result.
def ratioHisto(h1,h2):
    h3 = h1.Clone()
    h3.Sumw2()
    h2.Sumw2()
    h3.Divide(h2)

    return h3

# Cumulative ratio: bin i holds the ratio of the yields integrated from bin i
# upwards. This is the cut-and-count form, i.e. exactly what a mass threshold
# selects, as opposed to the differential ratio of ratioHisto.
#
# upTo = -1 integrates to the overflow; otherwise the integration stops at the
# bin containing upTo, which is how the blinded region is excluded.
#
# The error assumes h1 and h2 independent; they are not (both derive from the
# same data), so this over-estimates it, i.e. it is conservative.
def ratioIntegral(h1,h2,upTo=-1):
    h3 = h1.Clone()
    h3.Reset()
    if(upTo==-1):
        bornUp=h1.GetNbinsX()+1
    else:
        bornUp=h1.FindBin(upTo)
    for i in range(0,bornUp):
        # ctypes doubles: IntegralAndError returns its error through a C++
        # reference, which a plain Python float cannot bind to.
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
        # Bins where either integral vanishes are left empty: the ratio there is
        # 0/0 and would otherwise be drawn as a spurious point.
        if b != 0 and a != 0:
            c=math.sqrt((e1*e1)/(a*a)+(e2*e2)/(b*b))*a/b
            h3.SetBinContent(i,a/b)
            h3.SetBinError(i,c)
        else:
            h3.SetBinContent(i,0)
    return h3

# Pull: (observed - predicted) / sigma, with sigma the quadratic sum of the two
# uncertainties. GetBinErrorLow is used on the observation because it carries
# asymmetric Poisson errors after poissonning(); taking the lower error is the
# conservative choice for a positive excess.
def pullOfHisto(h2,h1):
    res=h1.Clone()
    for i in range (1,h1.GetNbinsX()+2):
        Perr=0
        Derr=0
        P=h1.GetBinContent(i)     # prediction
        D=h2.GetBinContent(i)     # observation
      
        Perr=h1.GetBinError(i)
        Derr=h2.GetBinErrorLow(i)
        
        # A bin with no uncertainty at all carries no information: set 0 rather
        # than divide by zero.
        if (Derr*Derr+Perr*Perr > 0): res.SetBinContent(i,(D-P)/math.sqrt(Derr*Derr+Perr*Perr))
        else: res.SetBinContent(i,0)

    return res

# Adds a flat relative systematic in quadrature with the existing bin errors.
# Sumw2(0) then SetBinError is the ROOT idiom to rewrite the error array from
# scratch. Called with syst = 0 to produce the "no systematics" band.
def addSyst(h,syst):
    res = h.Clone()
    res.Sumw2(0)
    for i in range (0,h.GetNbinsX()+1):
        res.SetBinError(i,math.sqrt(h.GetBinError(i)*h.GetBinError(i)+res.GetBinContent(i)*res.GetBinContent(i)*syst*syst))
    return res

# Builds the total uncertainty band of the prediction from three contributions:
#   - the statistical error already carried by h (toy-to-toy RMS from step2);
#   - the relative systematic read bin by bin from h_syst, which is in percent;
#   - a bias term, the absolute difference between h and hCorrBias.
#
# The while loop propagates the last non-zero systematic downwards: in the tail
# the systematics histogram has empty bins, and leaving them at zero would make
# the band collapse exactly where it should be widest.
#
# Returns the histogram with its new errors, plus the down and up envelopes as
# separate histograms.
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

# Blinding: zeroes the content of every bin whose lower edge is above m.
# The errors are deliberately left untouched, so a blinded point disappears from
# the plot without leaving a stray error bar at zero.
def blindAnyUp(h,m):
    for i in range (0,h.GetNbinsX()+1):
        mass = h.GetBinLowEdge(i)
        if(mass>m): 
            h.SetBinContent(i,0)

# Folds the overflow into the last bin, contents only. Used inside allSet, where
# the histogram is about to be normalised and the errors are irrelevant.
def MyoverflowInLastBin(h):
    res = h.Clone()
    res.SetBinContent(h.GetNbinsX(), h.GetBinContent(h.GetNbinsX()) + h.GetBinContent(h.GetNbinsX() + 1))
    res.SetBinContent(h.GetNbinsX()+1, 0)
    return res

# Same for the underflow.
def MyunderflowInFirstBin(h):
    res = h.Clone()
    res.SetBinContent(1, h.GetBinContent(0) + h.GetBinContent(1))
    res.SetBinContent(0, 0)
    return res

# Relative statistical uncertainty, in percent, bin by bin.
# Empty bins are set to 0 rather than left as a division by zero.
def statErr(h1, name):
    statErr = h1.Clone()
    statErr.SetName(name)
    for i in range (1, statErr.GetNbinsX()):
        if statErr.GetBinContent(i)>0:
            statErr.SetBinContent(i, statErr.GetBinError(i)/statErr.GetBinContent(i))
        else:
            statErr.SetBinContent(i,0)
    return 100*statErr

# Standard preparation of a spectrum: rebin onto the analysis binning, fold the
# under/overflow in, normalise to unit area. Used only for the fallback
# systematics below, where a shape is enough.
def allSet(h, sizeRebinning,  rebinning , st):
    h = h.Rebin(sizeRebinning, st, rebinning)
    h = MyoverflowInLastBin(h)
    h = MyunderflowInFirstBin(h)
    if h.Integral() > 0: h.Scale(1./h.Integral())
    return h

def getNominalSyst(ifile, plotType, region, sizeRebinning, rebinning):
    """Relative statistical uncertainty (%) of the nominal prediction, bin by bin.
    Stands in for the systematics when the combined-systematics file is absent."""
    pred = ifile.Get(plotType + region)
    if not pred:
        raise RuntimeError("'{}{}' missing from the file".format(plotType, region))

    pred = allSet(pred, sizeRebinning, rebinning, "nominal_def")
    syst_stat_binned = statErr(pred, "Stat_binned")
    syst_stat_binned.SetDirectory(0)
    return syst_stat_binned   # do NOT close ifile: pred/obs depend on it


#----------------------------------------------------
#                       Main
#----------------------------------------------------

def main(argv):
    # -------------- Setup --------------
    outputfile  = ''
    inputfile   = ''
    cuts        = ''
    systfile    = ''
    region      = ''
    odir        = ''
    eta         = ''
    nominalOnly = True
    isMC        = False

    # Settings overridable from the launcher
    isTTbar = False       # MC: plot the dileptonic ttbar only
    Vsignal = "19p12"     # version of the gluino samples
    year    = "2024"
    era     = ""          # "" (whole year) | "F" | "G"

    try:
        opts, args = getopt.getopt(argv, "", ["ifile=", "cuts=", "ofile=",
                                              "region=", "odir=", "nom=",
                                              "eta=", "isMC=", "systfile=",
                                              "vsignal=", "isTTbar=",
                                              "year=", "era="])
    except getopt.GetoptError:
        print("PlottingMacro.py --ifile <f> --cuts <c> --ofile <o> --region <r> "
              "--odir <d> --nom <bool> --eta <e> --isMC <bool> --systfile <f> "
              "[--vsignal <v>] [--isTTbar <bool>] [--year <y>] [--era <e>]")
        sys.exit(2)

    for o, arg in opts:
        if   o == "--ifile":    inputfile   = arg
        elif o == "--cuts":     cuts        = arg
        elif o == "--ofile":    outputfile  = arg
        elif o == "--region":   region      = arg
        elif o == "--odir":     odir        = arg
        elif o == "--nom":      nominalOnly = asBool(arg)
        elif o == "--eta":      eta         = arg
        elif o == "--isMC":     isMC        = asBool(arg)
        elif o == "--systfile": systfile    = arg
        elif o == "--vsignal":  Vsignal     = arg
        elif o == "--isTTbar":  isTTbar     = asBool(arg)
        elif o == "--year":     year        = arg
        elif o == "--era":      era         = arg

    os.system('mkdir -p ' + odir)
    outputfile = odir + '/' + outputfile

    # Plot-level switches, fixed here rather than exposed on the command line.
    isBinWidth  = False        # divide the yields by the bin width
    doRebin     = True         # rebin onto the analysis binning
    PlotSignal  = True         # overlay the gluino mass points
    doYouWantRratio = False    # show the cumulative-ratio pad t2
    blind       = (region == "9fp10")   # the search region is blinded above 300 GeV
    # Base name of the step1 histograms; must reproduce step1's naming exactly.
    mcTag       = "METanalysis_TestPUppiMETCut" + cuts + "_" + eta

    print (' Input file: ', inputfile)
    print ('Output file: ', outputfile)
    print ('     Region: ', region, '   eta:', eta, '   blind:', blind)
    print ('     Signal: V' + Vsignal, '   year:', year, ('era ' + era) if era else '')


    ifile = ROOT.TFile(inputfile)
    if (not ifile) or ifile.IsZombie():
        print("Error: can't open " + inputfile)
        sys.exit(1)

    # The two histograms written by bckgEstimate() for this region.
    pred   = ifile.Get("mass_predBC_" + region)
    obs    = ifile.Get("mass_obs_"    + region)
    if (not pred) or (not obs):
        # Almost always means step2 was run for the other region.
        print("Error: 'mass_predBC_{0}' or 'mass_obs_{0}' missing from {1}".format(
              region, inputfile))
        print("  -> did you run step2 with runVR/runSR for this region?")
        sys.exit(1)
    # Detached so they survive any file being closed later.
    for h in (pred, obs):
        h.SetDirectory(0)
    

    # -------------- MC mode --------------
    # In MC the "observation" is replaced by the sum of the simulated processes,
    # which lets the ABCD closure be tested where the truth is known.
    if isMC:
        ifileWjet  = ROOT.TFile.Open("/safe/ui3_1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/output/Wjets2024_V14/WjetMuNu2024_V14p14_weighted.root")
        ifileTTbar = ROOT.TFile.Open("/safe/ui3_1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/output/TTbar2024_V15/TTbar2024_V15p10_weighted.root")
        ifileTTbarSemiLep = ROOT.TFile.Open("/safe/ui3_1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/output/TTbar2024_V15/TTbarSemiLep2024_V22p1_weighted.root")
        ifileQCD   = ROOT.TFile.Open("/safe/ui3_1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/output/QCD2024_V16/QCD2024_mu_V16p3_weighted.root")

        # The search region reads region D, the validation region reads region C.
        prefix = "mass_regionD_" if region == "9fp10" else "mass_regionC_"
        obsQCD          = ifileQCD.Get(prefix + region + "_" + mcTag)
        obsWjet         = ifileWjet.Get(prefix + region + "_" + mcTag)
        obsTTbar        = ifileTTbar.Get(prefix + region + "_" + mcTag)
        obsTTbarSemiLep = ifileTTbarSemiLep.Get(prefix + region + "_" + mcTag)

        if (not obsQCD) or (not obsWjet) or (not obsTTbar) or (not obsTTbarSemiLep):
            print("Error: one of the MC histograms is None. Check the input files and histogram names.")
            sys.exit(1)

        # Semi-transparent fills with a black outline, so the stack stays
        # readable where the components overlap.
        obsQCD.SetFillColorAlpha(ROOT.kGreen-4, 0.5)
        obsWjet.SetFillColorAlpha(ROOT.kBlue-7, 0.5)
        obsTTbar.SetFillColorAlpha(ROOT.kRed, 0.5)
        obsTTbarSemiLep.SetFillColorAlpha(ROOT.kRed+2, 0.5)
        for h in (obsQCD, obsWjet, obsTTbar, obsTTbarSemiLep):
            h.SetLineColor(ROOT.kBlack)
            h.SetLineWidth(1)

        # isTTbar restricts the "observation" to the dileptonic ttbar alone,
        # which is the cleanest sample to test the method on.
        obs = (obsTTbar if isTTbar else obsQCD).Clone("obs_MC")
        obs.SetDirectory(0)
        if not isTTbar:
            obs.Add(obsWjet)
            obs.Add(obsTTbar)
            obs.Add(obsTTbarSemiLep)

    # Prediction with its statistical error only, drawn as the inner band.
    pred_noSyst = addSyst(pred,0.0)


    # -------------- Signal overlays --------------
    # Three gluino mass points, read from the weighted step1 files.
    ifileGl2000 = ROOT.TFile("/safe/ui3_1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/output/Gluino_V19/Gluino_Run3_MET_madgraph_2000_V" + Vsignal + "_weighted.root")
    m_Gl2000 = ifileGl2000.Get(mcTag + "_" + region + "_SignalMass_nominal")
    ifileGl2400 = ROOT.TFile("/safe/ui3_1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/output/Gluino_V19/Gluino_Run3_MET_madgraph_2400_V" + Vsignal + "_weighted.root")
    #ifileGl2400 = ROOT.TFile("/safe/ui3_1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/output/Stau_V20/Stau_Run3_MET_871_V20p0_weighted.root")
    m_Gl2400 = ifileGl2400.Get(mcTag + "_" + region + "_SignalMass_nominal")
    ifileGl2600 = ROOT.TFile("/safe/ui3_1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/output/Gluino_V19/Gluino_Run3_MET_madgraph_2600_V" + Vsignal + "_weighted.root")
    #ifileGl2600 = ROOT.TFile("/safe/ui3_1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/output/Stop_V21/Stop_Run3_MET_madgraph_1800_V21p0_weighted.root")
    m_Gl2600 = ifileGl2600.Get(mcTag + "_" + region + "_SignalMass_nominal")
    
    if (not m_Gl2400) or (not m_Gl2000) or (not m_Gl2600):
        print("Error: one of the signal histograms is None. Check the input files and histogram names.")
        sys.exit(1)



    # -------------- Work on histograms --------------
    # Analysis mass binning: fine at low mass, growing towards the tail so every
    # bin keeps a usable population. Must stay identical to the xbins array of
    # rebinHisto() in CommonFunctions.h and to REBINNING in the systematics
    # scripts, otherwise the uncertainty band would not line up with the points.
    rebinning = array.array('d',[0.,20.,40.,60.,80.,100.,120.,140.,160.,180.,200.,220.,240.,260.,
                                 280.,300.,320.,340.,360.,380.,410.,440.,480.,530.,590.,660.,760.,
                                 880.,1030.,1210.,1440.,1730.,2000.,2500.,3200.,4000.])


    sizeRebinning = len(rebinning)-1
    
    if(doRebin==True):
        pred        = pred.Rebin(sizeRebinning,"pred_new",rebinning)
        pred_noSyst = pred_noSyst.Rebin(sizeRebinning,"pred_noSyst_new",rebinning)
        obs         = obs.Rebin(sizeRebinning,"obs_new",rebinning)

    if isMC and doRebin:
        obsQCD   = obsQCD.Rebin(sizeRebinning,  "obsQCD_new",   rebinning)
        obsWjet  = obsWjet.Rebin(sizeRebinning, "obsWjet_new",  rebinning)
        obsTTbar = obsTTbar.Rebin(sizeRebinning,"obsTTbar_new", rebinning)
        obsTTbarSemiLep = obsTTbarSemiLep.Rebin(sizeRebinning,"obsTTbarSemiLep_new", rebinning)


    # Integrated luminosity (fb-1) per era; the signal samples are generated
    # for 109 fb-1.
    LUMI = {"": 109.0, "F": 25.40, "G": 34.4}
    if year != "2024":
        raise RuntimeError("Luminosity not defined for year " + year)
    lumi = LUMI[era]
    # Rescale the signal from the luminosity it was generated with to the one
    # actually being plotted.
    normSignal = lumi / 109.0
    
    if (PlotSignal):
        m_Gl2400.Scale(normSignal)
        m_Gl2000.Scale(normSignal)
        m_Gl2600.Scale(normSignal)

    if(doRebin==True):
        if (PlotSignal):
            m_Gl2400 = m_Gl2400.Rebin(sizeRebinning,"Gl1400_new",rebinning)
            m_Gl2000 = m_Gl2000.Rebin(sizeRebinning,"Gl2000_new",rebinning)
            m_Gl2600 = m_Gl2600.Rebin(sizeRebinning,"Gl2600_new",rebinning)


    # Snapshots taken before the uncertainty band and the blinding are applied.
    # pred_noCorrBias feeds the bias term of addHSyst; the _noBlind copies keep
    # the yields above the blinding threshold for the numbers returned at the end.
    pred_noCorrBias = pred.Clone()
    pred_noBlind = pred.Clone("_prednoBlind")
    obs_noBlind = obs.Clone("_obsnoBlind")


    # -------------- Uncertainty band --------------
    if not nominalOnly:
        # Full systematics, read from the file produced by systBckg.py.
        if not systfile:
            raise RuntimeError("--syst requested but --systfile is empty")
        ifileSyst = ROOT.TFile(systfile)
        if ifileSyst.IsZombie():
            raise RuntimeError("Unreadable systematics file: " + systfile)
        print(" syst. file: " + systfile)

        histoOfSyst = ifileSyst.Get("systTotalBinned")
        if not histoOfSyst:
            raise RuntimeError("'systTotalBinned' missing from " + systfile)

        (pred, predD, predU) = addHSyst(pred, histoOfSyst, pred_noCorrBias)
        (pred_noBlind, pred_noBlindD, pred_noBlindU) = addHSyst(
            pred_noBlind, histoOfSyst, pred_noCorrBias)
    else:
        # Fallback: the statistical uncertainty of the prediction stands in for
        # the systematics, so the plot can be made before the systematics run.
        print(" /!\\ only nominal")
        histoOfSystnom = getNominalSyst(ifile, "mass_predBC_", region,
                                        sizeRebinning, rebinning)
        (pred, predD, predU) = addHSyst(pred, histoOfSystnom, pred_noCorrBias)
        (pred_noBlind, pred_noBlindD, pred_noBlindU) = addHSyst(
            pred_noBlind, histoOfSystnom, pred_noCorrBias)
        
    # Yields above 300 GeV, computed on the unblinded copies and returned to the
    # caller. This is the cut-and-count number the analysis ultimately quotes.
    err_obs_m300 = ctypes.c_double(0)
    obs_m300 = obs_noBlind.IntegralAndError(obs_noBlind.FindBin(300), obs_noBlind.GetNbinsX()+1, err_obs_m300)
    err_obs_m300 = err_obs_m300.value

    err_pred_m300 = ctypes.c_double(0)
    pred_m300 = pred_noBlind.IntegralAndError(pred_noBlind.FindBin(300), pred_noBlind.GetNbinsX()+1, err_pred_m300)
    err_pred_m300 = err_pred_m300.value


    # Display range. mass_fit is both the blinding threshold and the position of
    # the vertical line drawn on every pad.
    mass_fit=300
    min_mass=0
    max_mass=4000
    if(doRebin==False):    
        max_mass=2500

    # Fold everything above the display limit into the last visible bin.
    # The `data` flag selects the error convention (see overflowInLastBin).
    underflowAndOverflow(obs,True, max_mass)
    underflowAndOverflow(pred, False, max_mass) 
    underflowAndOverflow(pred_noSyst, False, max_mass) 

    if isMC:
        for h in (obsQCD, obsWjet, obsTTbar, obsTTbarSemiLep):
            underflowAndOverflow(h, False, max_mass)


    if (PlotSignal):
        underflowAndOverflow(m_Gl2400, False, max_mass)
        m_Gl2400 = setColorAndMarker(m_Gl2400, ROOT.kViolet+1, 21)
        m_Gl2400.SetLineWidth(2)

        underflowAndOverflow(m_Gl2000, False, max_mass)
        m_Gl2000 = setColorAndMarker(m_Gl2000, ROOT.kOrange+8, 22)
        m_Gl2000.SetLineWidth(2)

        underflowAndOverflow(m_Gl2600, False, max_mass)
        m_Gl2600 = setColorAndMarker(m_Gl2600, ROOT.kGreen-3, 23)
        m_Gl2600.SetLineWidth(2)


    # Poisson errors for the observed distribution
    # Only on real data: the MC "observation" is weighted, so re-filling it entry
    # by entry would destroy the weights.
    if (not isMC): obs = poissonning(obs)


    # MC stack, built bottom-up. The order here is the order of the layers, not
    # the order of the legend.
    if isMC:
        stackMC = THStack("stackMC", "")
        stackMC.Add(obsTTbar)   # stacking order: bottom to top
        if not isTTbar:
            stackMC.Add(obsWjet)
            stackMC.Add(obsTTbarSemiLep)
            stackMC.Add(obsQCD)
    
    
    # The three comparisons drawn in the lower pads.
    ratioSimpleH    = ratioHisto(obs,pred)
    pull  = pullOfHisto(obs,pred)
    ratioInt  = ratioIntegral(obs,    pred, max_mass)

    # In the search region the cumulative ratio stops at the blinding threshold,
    # so no information leaks in from above it.
    if (blind):
        ratioInt  = ratioIntegral(obs,    pred, 300)

    if(isBinWidth):
        obs = binWidth(obs)
        pred = binWidth(pred)
        pred_noSyst = binWidth(pred_noSyst)
        if (PlotSignal):
            m_Gl2400 = binWidth(m_Gl2400)
            m_Gl2000 = binWidth(m_Gl2000)
            m_Gl2600 = binWidth(m_Gl2600)

    # Copies drawn as filled error bands; the originals stay as markers.
    pred_band=pred.Clone()
    pred_band_noSyst=pred_noSyst.Clone()


    # Blinding applies to the drawn copies only. Outside the search region the
    # names are simply aliased to the originals, so the drawing code below does
    # not need to know which case it is in.
    if blind:
        obs_blind       = obs.Clone("obs_blind")
        ratioInt_blind  = ratioInt.Clone("ratioInt_blind")
        ratioSimpleH_blind   = ratioSimpleH.Clone("ratioSimpleH_blind")
        pull_blind  = pull.Clone("pull_blind")
        for h in [obs_blind, ratioInt_blind, ratioSimpleH_blind, pull_blind]:
            blindAnyUp(h, mass_fit)
    else:
        obs_blind          = obs
        ratioInt_blind     = ratioInt
        ratioSimpleH_blind = ratioSimpleH
        pull_blind         = pull

    # -------------- Display --------------
    # Four pads sharing the x axis. t3 and t4 grow upwards when the cumulative
    # ratio pad t2 is not shown, so the figure stays full either way.
       
    c1=TCanvas("c1","c1",700,700)
    t1=TPad("t1","t1", 0.0, 0.45, 0.95, 0.95)     # spectra
    t1.Draw()
    t1.cd()
    t1.SetLogy(1)
    t1.SetTopMargin(0.003)
    t1.SetBottomMargin(0.04)
    c1.cd()

    t2=TPad("t2","t2", 0.0, 0.32, 0.95, 0.45)     # cumulative ratio
    t2.Draw()
    t2.cd()
    t2.SetGridy(1)
    t2.SetTopMargin(0.1)
    t2.SetBottomMargin(0.06)
    c1.cd()
    
    t3=TPad("t3","t3", 0.0, 0.18, 0.95, 0.32) if doYouWantRratio else TPad("t3","t3", 0.0, 0.27, 0.95, 0.45)   # obs / pred
    t3.Draw()
    t3.cd()
    t3.SetGridy(1)
    t3.SetTopMargin(0.1)
    t3.SetBottomMargin(0.07)
    c1.cd()

    t4=TPad("t4","t4", 0.0, 0.0, 0.95, 0.18) if doYouWantRratio else TPad("t4","t4", 0.0, 0.0, 0.95, 0.27)     # pull
    t4.Draw()
    t4.cd()
    t4.SetGridy(1)
    t4.SetTopMargin(0.1)
    t4.SetBottomMargin(0.32)


    t1.cd()

    # Fixed log-scale range, so plots of different eta ranges can be compared
    # side by side.
    min_entries = 1e-3
    max_entries = 2e5

    titleYaxis = "Events / bin"
    if (isBinWidth):
        titleYaxis = "Tracks / bin width"
    

    # Outer band: total uncertainty of the prediction. Drawn first so everything
    # else sits on top of it. Font code 43 = absolute pixel sizes, so the labels
    # keep the same size across the differently sized pads.
    pred_band.GetXaxis().SetTitle("")
    pred_band.GetYaxis().SetTitle(titleYaxis)
    pred_band.GetYaxis().SetLabelFont(43)
    pred_band.GetYaxis().SetLabelSize(22)
    pred_band.GetYaxis().SetTitleFont(43)
    pred_band.GetYaxis().SetTitleSize(22)
    pred_band.GetYaxis().SetTitleOffset(1.8)

    pred_band.SetMarkerStyle(22)
    pred_band.SetMarkerColor(5)
    pred_band.SetMarkerSize(1.0)
    pred_band.SetLineColor(5)
    pred_band.SetFillColorAlpha(5,0.8)
    pred_band.SetFillStyle(1001)
    pred_band.GetXaxis().SetRange(min_mass,max_mass)
    pred_band.GetXaxis().SetRangeUser(min_mass,max_mass)
    pred_band.GetYaxis().SetRangeUser(min_entries,max_entries)
    pred_band.GetXaxis().SetTitle("")
    pred_band.GetXaxis().SetLabelSize(0)
    pred_band.Draw("same E5")      # E5 = filled error band


    if isMC:
        stackMC.Draw("same hist")


    # Inner band: statistical uncertainty only. Same colour, drawn on top, so the
    # visible difference between the two bands is the systematic contribution.
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
    pred_band_noSyst.SetFillColorAlpha(5,0.8)
    pred_band_noSyst.SetFillStyle(1001)
    pred_band_noSyst.GetXaxis().SetRange(min_mass,max_mass)
    pred_band_noSyst.GetXaxis().SetRangeUser(min_mass,max_mass)
    pred_band_noSyst.GetYaxis().SetRangeUser(min_entries,max_entries)
    pred_band_noSyst.GetXaxis().SetTitle("")
    pred_band_noSyst.Draw("same E5")

    # Central value of the prediction, as red markers over the bands.
    pred.SetMarkerStyle(21)
    pred.SetMarkerColor(2)
    pred.SetMarkerSize(1)
    pred.SetLineColor(2)
    pred.SetFillColor(0)
    pred.Draw("same HIST P")

    # Observation: black points with error bars, blinded above 300 GeV in the
    # search region. Region 3fp8 is drawn in green to mark that it is a
    # cross-check region rather than the nominal one.
    obs_blind.SetMarkerStyle(20)
    obs_blind.SetMarkerColor(1)
    obs_blind.SetMarkerSize(1.0)
    obs_blind.SetLineColor(1)
    obs_blind.SetFillColor(0)
    obs_blind.GetXaxis().SetRange(min_mass,max_mass)
    obs_blind.GetXaxis().SetRangeUser(min_mass,max_mass)
    if (region=="3fp8"): 
        obs_blind.SetMarkerColor(8)
        obs_blind.SetLineColor(8)
        obs_blind.SetMarkerStyle(23)
    obs_blind.Draw("same E1")
    if (PlotSignal):
        m_Gl2400.Draw("same hist")
        m_Gl2000.Draw("same hist")
        m_Gl2600.Draw("same hist")
    

    # The legend starts lower when the signal curves add three more entries.
    leg=TLegend(0.65,0.5 if PlotSignal else 0.6,1.1,1)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(43)
    leg.SetTextSize(16)

    # Standard CMS-style header: luminosity on the right, provenance on the left.
    lumiText = "{:.4g} fb^{{-1}} (13.6 TeV)".format(lumi)
    if era: lumiText = year + era + " - " + lumiText
    tex1 = ROOT.TLatex(0.63 if era else 0.68, 0.96, lumiText)
    tex1.SetNDC()
    tex1.SetTextFont(42)
    tex1.SetLineWidth(2)
    tex1.SetTextSize(0.03)
    c1.cd()
    tex1.Draw()

    tex2 = ROOT.TLatex(0.15, 0.96, "#it{Private work (CMS data)}")
    if (isMC): tex2 = ROOT.TLatex(0.15, 0.96, "#it{Private work (CMS simulation)}")
    tex2.SetNDC()
    tex2.SetTextFont(42)
    tex2.SetTextSize(0.03)
    tex2.SetLineWidth(2)
    c1.cd()
    tex2.Draw()

    # Which Fpixel slice this plot covers.
    tex3 = ROOT.TLatex(0.18, 0.92, "#bf{Validation Region: 0.8< F_{pixel}#leq0.9}")
    if (region == "9fp10"): tex3 = ROOT.TLatex(0.18, 0.92, "#bf{Signal Region: 0.9< F_{pixel}#leq1.0}")
    tex3.SetNDC()
    tex3.SetTextFont(42)
    tex3.SetTextSize(0.03)
    c1.cd()
    tex3.SetLineWidth(2)
    tex3.Draw()

    # Eta range label; defaults to |eta|<1 when the range is not recognised.
    tex4 = ROOT.TLatex(0.18, 0.88, "#bf{|#eta|<1}")
    if (eta == 'Eta1_2p4'): tex4 = ROOT.TLatex(0.18, 0.88, "#bf{1#leq|#eta|<2.4}")
    if (eta == 'Eta2p4'): tex4 = ROOT.TLatex(0.18, 0.88, "#bf{|#eta|<2.4}")
    tex4.SetNDC()
    tex4.SetTextFont(42)
    tex4.SetTextSize(0.03)
    tex4.SetLineWidth(2)
    c1.cd()
    tex4.Draw()


    # Dedicated legend entry combining the marker of `pred` with the fill of the
    # band, so one entry describes both.
    pred_leg = pred.Clone()
    pred_leg.SetFillColor(pred_band_noSyst.GetFillColor())
    pred_leg.SetFillStyle(pred_band_noSyst.GetFillStyle())


    if (region=="3fp8"): leg.AddEntry(obs_blind, "Observed in C", "PE1")
    else: leg.AddEntry(obs_blind, "Observed", "PE1")
    
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

    if isMC:
        if not isTTbar:
            leg.AddEntry(obsWjet,  "W+jets",         "f")
            leg.AddEntry(obsQCD,   "QCD",            "f")
            leg.AddEntry(obsTTbarSemiLep, "t#bar{t}#rightarrow l#nu",       "f")
        leg.AddEntry(obsTTbar, "t#bar{t}#rightarrow 2l2#nu",       "f")
 

    if (PlotSignal):
        leg.AddEntry(m_Gl2000,"#tilde{g} (M=2000 GeV)","l")
        leg.AddEntry(m_Gl2400,"#tilde{g} (M=2400 GeV)","l")
        leg.AddEntry(m_Gl2600,"#tilde{g} (M=2600 GeV)","l")
        #leg.AddEntry(m_Gl2400,"#tilde{#tau}_{L,R} (M=871 GeV)","l")
        #leg.AddEntry(m_Gl2600,"#tilde{t} (M=1800 GeV)","l")

    

    # Vertical line at the blinding threshold, repeated on every pad.
    LineFit1=TLine(mass_fit, min_entries, mass_fit, max_entries)
    LineFit1.SetLineStyle(1)
    LineFit1.SetLineColor(1)
    
    t1.cd()
    t1.RedrawAxis()      # the filled bands would otherwise cover the frame
    if (blind): LineFit1.Draw("same")
    leg.Draw("same")
    
    # Reference line at 1, shared by both ratio pads.
    LineAtOne=TLine(min_mass,1,max_mass,1)
    LineAtOne.SetLineStyle(3)
    LineAtOne.SetLineColor(1)

    # -------------- t2: cumulative ratio (optional) --------------
    if (doYouWantRratio):
        c1.cd()
        t2.cd()
        
        # Empty frame fixing the axes; the ratio is drawn on top with "same".
        frameR=ROOT.TH1D("frameR", "frameR", 1,min_mass, max_mass)
        frameR.GetXaxis().SetNdivisions(505)
        frameR.SetTitle("")
        frameR.SetStats(0)
        frameR.GetXaxis().SetTitle("")
        frameR.GetYaxis().SetTitle("RatioR ")
        frameR.SetMaximum(2.)
        frameR.SetMinimum(0.0)
        frameR.GetYaxis().SetLabelFont(43) #give the font size in pixel (instead of fraction)
        frameR.GetYaxis().SetLabelSize(14) #font size
        frameR.GetYaxis().SetTitleFont(43) #give the font size in pixel (instead of fraction)
        frameR.GetYaxis().SetTitleSize(18) #font size
        frameR.GetYaxis().SetNdivisions(503)
        frameR.GetXaxis().SetLabelSize(0) #font size
        frameR.GetXaxis().SetTitleOffset(3.75)
        frameR.Draw("AXIS")

        ratioInt_blind.SetMarkerStyle(21)
        ratioInt_blind.SetMarkerColor(1)
        ratioInt_blind.SetMarkerSize(0.7)
        ratioInt_blind.SetLineColor(1)
        ratioInt_blind.SetFillColor(0)
        if (region=="3fp8"):
            ratioInt_blind.SetMarkerColor(8)
            ratioInt_blind.SetLineColor(8)
            ratioInt_blind.SetMarkerStyle(23)
        ratioInt_blind.Draw("same E0")

        LineAtOne.Draw("same")

        ratioInt_blind.GetXaxis().SetRange(min_mass,max_mass)
        ratioInt_blind.GetXaxis().SetRangeUser(min_mass,max_mass)

        LineFit2 = TLine(mass_fit, 0, mass_fit, 2.)
        LineFit2.SetLineStyle(1)
        LineFit2.SetLineColor(1)
        if (blind): LineFit2.Draw("same")



    # -------------- t3: bin-by-bin ratio obs / pred --------------
    c1.cd()
    t3.cd()
    
    frameR2=ROOT.TH1D("frameR2", "frameR2", 1,min_mass, max_mass)
    frameR2.GetXaxis().SetNdivisions(505)
    frameR2.SetTitle("")
    frameR2.SetStats(0)
    frameR2.GetXaxis().SetTitle("")
    frameR2.GetYaxis().SetTitle("obs / pred")
    frameR2.GetYaxis().SetRangeUser(0.,2.)
    frameR2.GetYaxis().SetLabelFont(43) #give the font size in pixel (instead of fraction)
    frameR2.GetYaxis().SetLabelSize(22) #font size
    frameR2.GetYaxis().SetTitleFont(43) #give the font size in pixel (instead of fraction)
    frameR2.GetYaxis().SetTitleSize(22) #font size
    frameR2.GetYaxis().SetNdivisions(503)
    frameR2.GetXaxis().SetLabelSize(0) #font size
    frameR2.GetXaxis().SetTitleOffset(3.75)
    frameR2.GetYaxis().SetTitleOffset(1.4)
    frameR2.Draw("AXIS")

    ratioSimpleH_blind.Sumw2()
    ratioSimpleH_blind.SetMarkerStyle(21)
    ratioSimpleH_blind.SetMarkerColor(1)
    ratioSimpleH_blind.SetMarkerSize(0.7)
    ratioSimpleH_blind.SetLineColor(1)
    ratioSimpleH_blind.SetFillColor(0)
    if (region=="3fp8"):
        ratioSimpleH_blind.SetMarkerColor(8)
        ratioSimpleH_blind.SetLineColor(8)
        ratioSimpleH_blind.SetMarkerStyle(23)


    ratioSimpleH_blind.Draw("same E0")
    ratioSimpleH_blind.GetXaxis().SetRange(min_mass,max_mass)
    ratioSimpleH_blind.GetXaxis().SetRangeUser(min_mass,max_mass)

    LineAtOne.Draw("same")

    LineFit3=TLine(mass_fit,0,mass_fit,2.)
    LineFit3.SetLineStyle(1)
    LineFit3.SetLineColor(1)
    if (blind): LineFit3.Draw("same")


    ratioSimpleH_blind.GetXaxis().SetRange(min_mass,max_mass)
    ratioSimpleH_blind.GetXaxis().SetRangeUser(min_mass,max_mass)


    # -------------- t4: pull --------------
    # Bottom pad, the only one carrying the x axis title and labels.
    c1.cd()
    t4.cd()

    frameR3=ROOT.TH1D("frameR3", "frameR3", 1,min_mass, max_mass)
    frameR3.GetXaxis().SetNdivisions(505)
    frameR3.SetTitle("")
    frameR3.SetStats(0)
    frameR3.GetXaxis().SetTitle("Mass (GeV)")
    frameR3.GetYaxis().SetTitleOffset(1.4)
    frameR3.GetYaxis().SetTitle("#frac{M_{obs}-M_{pred}}{#sigma}") if doYouWantRratio else frameR3.GetYaxis().SetTitle("#frac{M_{obs}-M_{pred}}{#sigma} ")
    frameR3.GetYaxis().SetTickLength(frameR3.GetYaxis().GetTickLength()*2)
    frameR3.SetMaximum(3)
    frameR3.SetMinimum(-3)
    frameR3.GetYaxis().SetLabelFont(43) #give the font size in pixel (instead of fraction)
    frameR3.GetYaxis().SetLabelSize(22) #font size
    frameR3.GetYaxis().SetTitleFont(43) #give the font size in pixel (instead of fraction)
    frameR3.GetYaxis().SetTitleSize(22) #font size
    frameR3.GetYaxis().SetNdivisions(503)
    frameR3.GetXaxis().SetNdivisions(510)
    frameR3.GetXaxis().SetLabelFont(43) #give the font size in pixel (instead of fraction)
    frameR3.GetXaxis().SetLabelSize(22) #font size
    frameR3.GetXaxis().SetTitleFont(43) #give the font size in pixel (instead of fraction)
    frameR3.GetXaxis().SetTitleSize(22) #font size
    frameR3.GetXaxis().SetTitleOffset(1)
    frameR3.Draw("AXIS")

    pull_blind.Draw("same HIST")

    pull_blind.SetLineColor(1)
    pull_blind.SetFillColor(38)
    if (region=="3fp8"):
        pull_blind.SetFillColorAlpha(8, 0.35)
        pull_blind.SetLineColor(8)
    t4.RedrawAxis()        # the filled pull histogram covers the frame
    t4.RedrawAxis("G")     # and the grid

    
    # Guide lines at 0, +/-1 and +/-2 sigma, so the size of a deviation can be
    # read off directly.
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


    
    # Three formats: pdf to look at, .root and .C so the figure can be reopened
    # and restyled without rerunning the chain. The file name encodes every
    # switch, so two configurations cannot overwrite each other.
    c1.Update()
    c1.SaveAs(outputfile + "_region" + region + "_" + year + "_" + ("onlyNominal" if nominalOnly else "") + ("_wRatioR" if doYouWantRratio else "") + ("_MC" if isMC else "") + "_" + eta + ".pdf")
    c1.SaveAs(outputfile + "_region" + region + "_" + year + "_" + ("onlyNominal" if nominalOnly else "") + ("_wRatioR" if doYouWantRratio else "") + ("_MC" if isMC else "") + "_" + eta + ".root")
    c1.SaveAs(outputfile + "_region" + region + "_" + year + "_" + ("onlyNominal" if nominalOnly else "") + ("_wRatioR" if doYouWantRratio else "") + ("_MC" if isMC else "") + "_" + eta + ".C")

    print("   Saved in: {}".format(odir))
    print('')

    #Chi2ObsPred = obs.Chi2Test(pred,"UWP")
    #print("Chi2 between prediction and observation  = {}".format(Chi2ObsPred))

    # Returned for a caller that would aggregate the above-300 GeV yields across
    # eta ranges; currently unpacked and dropped by the __main__ block below.
    return (region, obs_m300, pred_m300, err_obs_m300, err_pred_m300)




if __name__ == "__main__":

    (reg, obs_m300, pred_m300, err_obs_m300, err_pred_m300) = main(sys.argv[1:])