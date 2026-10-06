#!/usr/bin/env python3
"""Systematics of the signal samples (gluinos).

Reads the step1 mass histograms (one per variation) and writes, for each mass
point, the bin-by-bin quadratic sum:

    <odir>/Gluino_<mass>_<region><cuts>_<eta>/sysTotBinned_signal.root
                                                    -> key 'systTotalBinned'

    python3 systSignal.py
    python3 systSignal.py --etas Eta1,Eta2p4 --region 9fp10
    python3 systSignal.py --masses 2000,2400,2600
    python3 systSignal.py --no-plots
"""

# =============================================================================
#  Step 3b — systematic envelope of the gluino signal samples.
# -----------------------------------------------------------------------------
#  Counterpart of systBckg.py, on the signal side. The structure is deliberately
#  parallel (same helpers, same output key), but three things differ:
#
#    1. Input layout. The background variations live in separate files, one per
#       label; here all the variations of a mass point sit in the *same* step1
#       file, distinguished by a suffix on the histogram name.
#    2. No normalisation. The spectra keep their absolute yields, since the flat
#       systematics (luminosity, Fpixel efficiency) are multiplicative on yields.
#    3. No statistical term. Only the systematic sources are summed here.
#
#  Sources: two-sided experimental variations (dE/dx K and C, pile-up, trigger
#  scale factors, jet energy scale) plus flat normalisation uncertainties.
# =============================================================================

import os, sys, math, array, ctypes
from optparse import OptionParser
from collections import OrderedDict

import ROOT
from ROOT import TFile, TCanvas, TLegend, TLatex, TPad, TLine
import tdrstyle

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning
# Detach histograms from their TFile so they survive ifile.Close().
ROOT.TH1.AddDirectory(False)
tdrstyle.setTDRStyle()


# ==================================================================
#   Settings
# ==================================================================
IDIR    = "/safe/ui3_1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/output/Gluino_V19/"
VSIGNAL = "19p12"
# Weighted files carry the cross-section and luminosity weights; the raw ones
# hold plain event counts. Must match what the plotting macros use.
WEIGHTED = True          # True -> files '..._weighted.root', like PlottingMacro.py
REGION  = "9fp10"
ETAS    = "Eta2p4"
CUTS    = ""             # same as ShowPlots.py: "" | "_SigmaPtoverPt_0p5_EoP_0p1"
ODIR    = "systSignal"

# Gluino mass points available in the production, in GeV.
MASSES  = [1100, 1200, 1300, 1400, 1600, 1800, 2000, 2200, 2400, 2600]

# Mass points repeated on the all-masses summary plot.
# An OrderedDict, not a plain dict, so the legend order is reproducible; the
# keys also act as the filter deciding which points reach that plot.
HIGHLIGHT = OrderedDict([(2000, ROOT.kOrange + 8),
                         (2400, ROOT.kViolet + 1),
                         (2600, ROOT.kGreen - 3)])

# Prefix of the step1 histograms: identical to PlottingMacro.py's 'mcTag'.
# If step1 changes its naming, this and histoName() are what must follow.
HNAME_PREFIX = "METanalysis_TestPUppiMETCut"

# Two-sided systematics: suffix of the step1 histograms.
# 'key' names the output histogram, 'title' is the label drawn on the individual
# plots, 'legend' the shorter one used in the summary legend.
SYSTEMATICS = [
    dict(key="K",       up="KUp",         down="KDown",
         legend="K",          title="K_{mass}",   color=ROOT.kMagenta - 9, marker=21),
    dict(key="C",       up="CUp",         down="CDown",
         legend="C",          title="C_{mass}",   color=ROOT.kViolet + 1,  marker=22),
    dict(key="PU",      up="PUUp",        down="PUDown",
         legend="PU",         title="Pile-Up",    color=ROOT.kBlue + 1,    marker=23),
    dict(key="Trigger", up="TriggerSFUp", down="TriggerSFDown",
         legend="Trigger",    title="Trigger SF", color=ROOT.kOrange,      marker=39),
    dict(key="Jet",     up="JetUp",       down="JetDown",
         legend="Jet",        title="JES",        color=ROOT.kOrange + 1,  marker=33),
]

# Flat systematics, in %.
# Constant across the whole spectrum, applied wherever the nominal is populated:
# integrated-luminosity uncertainty and Fpixel selection efficiency.
FLAT_SYSTEMATICS = [
    dict(key="lumi", value=1.4, legend="Lumi",      color=ROOT.kGreen + 2, marker=20),
    dict(key="Fpix", value=1.6, legend="F^{pixel}", color=ROOT.kCyan,      marker=29),
]

# Suffix of the reference histogram inside each step1 file.
NOMINAL_SUFFIX = "nominal"

# Same analysis binning as systBckg.py and as rebinHisto() in CommonFunctions.h.
# The three must stay identical for signal and background to be combined.
REBINNING = array.array('d', [0., 20., 40., 60., 80., 100., 120., 140., 160., 180.,
                              200., 220., 240., 260., 280., 300., 320., 340., 360., 380.,
                              410., 440., 480., 530., 590., 660., 760., 880., 1030., 1210.,
                              1440., 1730., 2000., 2500., 3200., 4000.])
SIZE_REBINNING = len(REBINNING) - 1

MAX_MASS = 4000


# ==================================================================
#   Histogram helpers
# ==================================================================
# Same helper set as systBckg.py; kept duplicated so each script stands alone.
_UID = [0]


def uniq(stem):
    # Unique names avoid ROOT's "Replacing existing TH1", which silently
    # destroys an object that may still be referenced.
    _UID[0] += 1
    return "{}_{}".format(stem, _UID[0])


def overflowInLastBin(h):
    # Folds the overflow into the last visible bin, errors in quadrature, and
    # empties bin N+1 so later 0..N+1 loops see zero there.
    res = h.Clone(uniq(h.GetName() + "_ovf"))
    n = h.GetNbinsX()
    res.SetBinContent(n, h.GetBinContent(n) + h.GetBinContent(n + 1))
    res.SetBinError(n, math.sqrt(h.GetBinError(n)**2 + h.GetBinError(n + 1)**2))
    res.SetBinContent(n + 1, 0)
    res.SetBinError(n + 1, 0)
    return res


def underflowInFirstBin(h):
    # Same for the underflow.
    res = h.Clone(uniq(h.GetName() + "_udf"))
    res.SetBinContent(1, h.GetBinContent(0) + h.GetBinContent(1))
    res.SetBinError(1, math.sqrt(h.GetBinError(0)**2 + h.GetBinError(1)**2))
    res.SetBinContent(0, 0)
    res.SetBinError(0, 0)
    return res


def allSet(h, name):
    """Variable Rebin + under/overflow folding. No normalisation: the
    absolute yields are used for flat-systematic corrections (lumi)."""
    # This is the one real difference with systBckg.allSet: no Scale(1/integral).
    res = h.Rebin(SIZE_REBINNING, uniq(name + "_reb"), REBINNING)
    res = overflowInLastBin(res)
    res = underflowInFirstBin(res)
    res.SetName(name)
    res.SetDirectory(0)
    return res


def ratioInt(num, den):
    # Ratio of the cumulative yields [i, overflow], bin by bin. Errors assume
    # num and den independent, which over-estimates them here since both come
    # from the same sample.
    res = num.Clone(uniq("ratioInt"))
    res.SetDirectory(0)
    res.Reset()
    last = num.GetNbinsX() + 1
    for i in range(0, num.GetNbinsX() + 1):
        # ctypes doubles: IntegralAndError returns its error through a C++
        # reference, which a Python float cannot bind to.
        e1 = ctypes.c_double(0.0)
        e2 = ctypes.c_double(0.0)
        a = num.IntegralAndError(i, last, e1, "")
        b = den.IntegralAndError(i, last, e2, "")
        if a != 0 and b != 0:
            res.SetBinContent(i, a / b)
            res.SetBinError(i, math.sqrt((e1.value / a)**2 + (e2.value / b)**2) * a / b)
        else:
            res.SetBinContent(i, 0)
    return res


def ratioBinned(num, den):
    # Differential ratio; this is the form actually used for the signal, whose
    # limits are set from the binned shape.
    res = num.Clone(uniq("ratioBin"))
    res.SetDirectory(0)
    res.Divide(den)
    return res


def systMass(nominal, down, up, name, binned=1, typec=0, mini=0, minNom=0.0):
    """Relative deviation (%) between the variations and the nominal value,
    bin by bin (binned=1) or based on the cumulative integral (binned=0).
    Bins where the nominal is empty are set to 0."""
    # Signal spectra are narrow peaks, so most bins are empty; minNom is the
    # threshold below which a bin is considered unpopulated. Without it, the
    # ratio 0/0 would be reported as a 100 % deviation over the whole tail.
    # Unlike the background version, `binned` defaults to 1 and there is no
    # RATIO_VAR_OVER_NOM switch: the convention is always var/nom.
    var = up if up is not None else down
    if binned:
        ra1 = ratioBinned(var,  nominal)
        ra2 = ratioBinned(down, nominal)
    else:
        ra1 = ratioInt(var,  nominal)
        ra2 = ratioInt(down, nominal)

    last = nominal.GetNbinsX() + 1
    res  = ra1.Clone(name)
    res.SetDirectory(0)
    res.Reset()
    for i in range(0, res.GetNbinsX() + 2):
        ref = nominal.GetBinContent(i) if binned else nominal.Integral(i, last)
        if ref <= minNom:
            res.SetBinContent(i, 0.0)
            continue
        s1 = abs(1 - ra1.GetBinContent(i))
        s2 = abs(1 - ra2.GetBinContent(i))
        if typec == 1:
            m = 0.5 * (s1 + s2)
        else:
            # Default: the envelope, i.e. the larger of the two deviations.
            m = min(s1, s2) if mini == 1 else max(s1, s2)
        res.SetBinContent(i, 100.0 * m)
        res.SetBinError(i, 0.0)
    return res


def systFlat(nominal, value, name, minNom=0.0):
    # Builds a constant systematic shaped like the nominal: `value` percent in
    # every populated bin, 0 elsewhere. Zeroing the empty bins keeps the flat
    # terms from inflating the total where there is no signal at all.
    res = nominal.Clone(name)
    res.SetDirectory(0)
    res.Reset()
    for i in range(0, res.GetNbinsX() + 2):
        res.SetBinContent(i, value if nominal.GetBinContent(i) > minNom else 0.0)
        res.SetBinError(i, 0.0)
    return res


def systTotal(list_h, name):
    # Quadratic sum, assuming the sources are uncorrelated.
    res = list_h[0].Clone(name)
    res.SetDirectory(0)
    res.Reset()
    for i in range(0, res.GetNbinsX() + 2):
        tot = sum(h.GetBinContent(i)**2 for h in list_h)
        res.SetBinContent(i, math.sqrt(tot))
        res.SetBinError(i, 0.0)
    return res


def binCenters(h):
    # Un marker au centre de chaque bin non vide, sans barre.
    g = ROOT.TGraph()
    for i in range(1, h.GetNbinsX() + 1):
        c = h.GetBinContent(i)
        if c <= 0:
            continue
        g.SetPoint(g.GetN(), h.GetBinCenter(i), c)
    return g


def setColorAndMarker(obj, color, markerstyle):
    # Mutates and returns its argument, so callers that need the original style
    # preserved must pass a clone.
    obj.SetLineColor(color)
    obj.SetMarkerColor(color)
    obj.SetFillColor(color)
    obj.SetMarkerStyle(markerstyle)
    obj.SetMarkerSize(1.2)
    return obj


def etaLatex(eta):
    # Eta-range label; x position tuned per string so the longer ones still fit.
    labels = {"Eta1":       (0.80, "#bf{|#eta|<1}"),
              "Eta1_2p4":   (0.75, "#bf{1#leq|#eta|<2.4}"),
              "Eta2p4":     (0.80, "#bf{|#eta|<2.4}"),
              "Eta1p2_2p2": (0.72, "#bf{1.2#leq|#eta|<2.2}"),
              "Eta1p2_2p4": (0.72, "#bf{1.2#leq|#eta|<2.4}")}
    x, txt = labels.get(eta, (0.80, "#bf{" + eta + "}"))
    lat = TLatex(x, 0.88, txt)
    lat.SetNDC()
    lat.SetTextFont(42)
    lat.SetTextSize(0.07)
    return lat


def safeName(s):
    """File name without LaTeX characters."""
    # The plot titles are TLatex strings; stripping braces, '#', slashes and
    # spaces makes them usable as file names.
    for ch in "{}#\\/ ":
        s = s.replace(ch, "")
    return s


# ==================================================================
#   Plots
# ==================================================================
def plotter(nominal, down, up, legNom, legDown, legUp, outDir, title, mass):
    # Two pads: the three spectra on a log scale, and the bin-by-bin ratio.
    # Unlike the background version there is no cumulative-ratio pad and `up` is
    # assumed present, since every signal systematic here is two-sided.
    c1 = TCanvas(uniq("c1"), "c1", 800, 600)
    t1 = TPad(uniq("t1"), "t1", 0.0, 0.25, 0.95, 0.95)
    t1.Draw()
    t1.SetLogy(1)
    t1.SetTopMargin(0.05)
    t1.SetBottomMargin(0.12)
    c1.cd()

    t3 = TPad(uniq("t3"), "t3", 0.0, 0., 0.95, 0.25)
    t3.Draw()
    t3.SetGridy(1)
    t3.SetTopMargin(0.05)
    t3.SetBottomMargin(0.4)

    t1.cd()
    # Floor at 0.2 event: below that the weighted spectra are pure noise.
    min_entries = 2e-1
    max_entries = nominal.GetMaximum() * 2

    nominal = setColorAndMarker(nominal, 1, 20)
    nominal.GetXaxis().SetRangeUser(0, MAX_MASS)
    nominal.GetYaxis().SetRangeUser(min_entries, max_entries)
    nominal.SetMinimum(min_entries)
    nominal.SetTitle(";Mass (GeV);Tracks")
    nominal.GetYaxis().SetTitleSize(0.07)
    nominal.GetYaxis().SetLabelSize(0.06)
    nominal.GetYaxis().SetTitleOffset(0.7)
    nominal.GetXaxis().SetTitle("Mass (GeV)")
    nominal.GetXaxis().SetTitleOffset(0.9)
    nominal.GetXaxis().SetLabelSize(0.06)
    nominal.GetXaxis().SetTitleSize(0.06)
    nominal.Draw()

    down = setColorAndMarker(down, ROOT.kBlue + 1, 23)
    down.Draw("same")
    up = setColorAndMarker(up, ROOT.kRed + 1, 22)
    up.Draw("same")

    leg = TLegend(0.2, 0.65, 0.35, 0.92)
    leg.AddEntry(nominal, legNom, "PE1")
    leg.AddEntry(down, legDown, "PE1")
    leg.AddEntry(up, legUp, "PE1")
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.Draw("same")

    c1.cd()
    t3.cd()

    # Empty frame fixing the ratio axes; the ratios are drawn on top with "same".
    frameR = ROOT.TH1D(uniq("frameR"), "frameR", 1, 0, MAX_MASS)
    frameR.SetDirectory(0)
    frameR.SetTitle("")
    frameR.SetStats(0)
    frameR.GetYaxis().SetTitle("#frac{var}{Nominal}")
    # Pile-up moves the spectrum more than the other sources, so its ratio pad
    # gets a wider window rather than clipped points.
    frameR.SetMaximum(1.2 if title == "Pile-Up" else 1.1)
    frameR.SetMinimum(0.8 if title == "Pile-Up" else 0.9)
    # Font code 43 = absolute pixel sizes, so labels do not shrink with the pad.
    frameR.GetYaxis().SetLabelFont(43)
    frameR.GetYaxis().SetLabelSize(23)
    frameR.GetYaxis().SetTitleFont(43)
    frameR.GetYaxis().SetTitleSize(24)
    frameR.GetYaxis().SetNdivisions(205)
    frameR.GetYaxis().SetTitleOffset(1.8)
    frameR.GetXaxis().SetNdivisions(510)
    frameR.GetXaxis().SetLabelFont(43)
    frameR.GetXaxis().SetLabelSize(23)
    frameR.GetXaxis().SetTitleFont(43)
    frameR.GetXaxis().SetTitleSize(24)
    frameR.GetXaxis().SetTitle("Mass (GeV)")
    frameR.GetXaxis().SetTitleOffset(1)
    frameR.Draw("AXIS")
    frameR.Draw("SAME AXIG")

    lineAtOne = TLine(0, 1, MAX_MASS, 1)
    lineAtOne.SetLineStyle(3)
    lineAtOne.SetLineColor(1)
    lineAtOne.Draw("same")

    # Local variables keep these alive until SaveAs; a bare expression would be
    # garbage-collected and the pad would come out empty.
    rDown = setColorAndMarker(ratioBinned(down, nominal), ROOT.kBlue + 1, 23)
    rDown.Draw("E0 same")
    rUp = setColorAndMarker(ratioBinned(up, nominal), ROOT.kRed + 1, 22)
    rUp.Draw("E0 same")

    c1.cd()
    # "CMS simulation" here, against "CMS data" in systBckg.py.
    latex = TLatex(0.15, 0.92, "#scale[1.3]{#it{Private work (CMS simulation)}}")
    latex.SetNDC()
    latex.SetTextFont(42)
    latex.SetTextSize(0.03)
    latex.Draw()

    tex = TLatex(0.69, 0.92, "#scale[1.3]{#bf{m_{#tilde{g}}=" + str(mass) + " GeV}}")
    tex.SetNDC()
    tex.SetTextFont(42)
    tex.SetTextSize(0.04)
    tex.Draw()

    tex2 = TLatex(0.71, 0.36, "#scale[1.3]{#bf{" + title + "}}")
    tex2.SetNDC()
    tex2.SetTextFont(42)
    tex2.SetTextSize(0.04)
    tex2.Draw()

    c1.SaveAs(os.path.join(outDir, safeName(title) + ".pdf"))


def plotSummary(entries, total, xtitle, outDir, eta, mass):
    """entries: list of (histo, legend, color, marker)."""
    # All sources of one mass point on a single log plot, in percent, total in red.
    c2 = TCanvas(uniq("c2"), "c2", 800, 600)
    c2.SetBottomMargin(0.12)
    c2.SetLeftMargin(0.1)
    c2.SetLogy()
    c2.SetGrid()

    frame = ROOT.TH1D(uniq("frameSummary"), "", 1, 0, MAX_MASS)
    frame.SetDirectory(0)
    frame.SetStats(0)
    frame.SetMinimum(0.1)
    frame.SetMaximum(2000)
    frame.GetXaxis().SetTitle(xtitle)
    frame.GetYaxis().SetTitle("Systematic Uncertainty [%]")
    frame.GetXaxis().SetNdivisions(510)
    frame.GetXaxis().SetLabelFont(43)
    frame.GetXaxis().SetLabelSize(24)
    frame.GetXaxis().SetTitleSize(0.05)
    frame.GetXaxis().SetTitleOffset(1.0)
    frame.GetYaxis().SetLabelFont(43)
    frame.GetYaxis().SetLabelSize(24)
    frame.GetYaxis().SetTitleSize(0.05)
    frame.GetYaxis().SetTitleOffset(1.0)
    frame.Draw("AXIS")
    frame.Draw("SAME AXIG")

    leg = TLegend(0.22, 0.75, 0.6, 0.93)
    leg.SetNColumns(2)
    leg.SetBorderSize(1)

    # `graphs` prevents the TGraphs from being collected before SaveAs.
    graphs = []
    gTot = setColorAndMarker(binCenters(total), ROOT.kRed, 34)
    leg.AddEntry(gTot, "Total", "P")
    for (h, legend, color, marker) in entries:
        g = setColorAndMarker(binCenters(h), color, marker)
        leg.AddEntry(g, legend, "P")
        g.Draw("P")
        graphs.append(g)
    gTot.Draw("P")                 # last, so the total sits on top
    leg.Draw("same")

    latex = TLatex(0.10, 0.96, "#scale[1.3]{#it{Private work (CMS simulation)}}")
    latex.SetNDC()
    latex.SetTextFont(42)
    latex.SetTextSize(0.04)
    latex.Draw()

    tex = TLatex(0.73, 0.96, "#scale[1.3]{#bf{m_{#tilde{g}}=" + str(mass) + " GeV}}")
    tex.SetNDC()
    tex.SetTextFont(42)
    tex.SetTextSize(0.04)
    tex.Draw()

    latex3 = etaLatex(eta)
    latex3.Draw()

    c2.Update()
    c2.SaveAs(os.path.join(outDir, "summary_Gluino_{}_{}.pdf".format(mass, eta)))


def filledRange(h):
    """Restricts the axis to non-empty bins. Returns False if the histogram is empty."""
    # A signal spectrum only populates a narrow mass window; without this the
    # plot would be dominated by a long flat line at zero. The error bars are
    # also cleared, since these histograms hold an envelope, not a measurement.
    bins = [i for i in range(1, h.GetNbinsX() + 1) if h.GetBinContent(i) > 0]
    if not bins:
        return False
    h.GetXaxis().SetRange(bins[0], bins[-1])
    for i in range(0, h.GetNbinsX() + 2):
        h.SetBinError(i, 0)
    return True


def plotTotalAllMasses(dict_sysTot, xtitle, outDir, eta, outTitle="sysTot_allMasses"):
    # Total uncertainty of the highlighted mass points, overlaid, to show how it
    # evolves with the gluino mass. dict_sysTot only ever contains keys present
    # in HIGHLIGHT, which is what makes the HIGHLIGHT[mass] lookup below safe.
    if not dict_sysTot:
        print("No mass point selected -> no summary plot")
        return

    c3 = TCanvas(uniq("c3"), "c3", 800, 600)
    c3.SetBottomMargin(0.1)
    c3.SetLeftMargin(0.1)
    c3.SetLogy()
    c3.SetGrid()

    # Cycled by enumeration index, so two overlapping curves stay separable.
    markers = [23, 21, 22]

    frame = ROOT.TH1D(uniq("frameTot"), "", 1, 0, MAX_MASS)
    frame.SetDirectory(0)
    frame.SetStats(0)
    # Floor at 1 % here rather than 0.1 %: the total is never smaller than the
    # flat terms, so a decade of empty axis would be wasted.
    frame.SetMinimum(1)
    frame.SetMaximum(2000)
    frame.GetXaxis().SetTitle(xtitle)
    frame.GetYaxis().SetTitle("Total Systematic Uncertainty [%]")
    frame.GetXaxis().SetNdivisions(510)
    frame.GetXaxis().SetLabelFont(43)
    frame.GetXaxis().SetLabelSize(25)
    frame.GetXaxis().SetTitleSize(0.05)
    frame.GetXaxis().SetTitleOffset(1.0)
    frame.GetYaxis().SetLabelFont(43)
    frame.GetYaxis().SetLabelSize(25)
    frame.GetYaxis().SetTitleSize(0.05)
    frame.GetYaxis().SetTitleOffset(1.0)
    frame.Draw("AXIS")
    frame.Draw("SAME AXIG")

    leg = TLegend(0.3, 0.68, 0.63, 0.93)
    leg.SetBorderSize(0)
    #leg.SetFillStyle(0)

    # Sorted so the legend goes from the lightest to the heaviest mass point.
    for i, mass in enumerate(sorted(dict_sysTot)):
        h = dict_sysTot[mass]
        # Empty spectra are skipped rather than drawn as a flat zero.
        if not filledRange(h):
            continue
        setColorAndMarker(h, HIGHLIGHT[mass], markers[i % len(markers)])
        h.SetFillStyle(0)          # setColorAndMarker set a fill colour
        h.SetLineWidth(2)
        leg.AddEntry(h, "m_{#tilde{g}}=" + str(mass) + " GeV", "LP")
        # Line then markers: the curve is readable and the bins stay identifiable.
        h.Draw("HIST SAME")
        h.Draw("P SAME")

    leg.Draw("same")

    latex = TLatex(0.10, 0.96, "#scale[1.3]{#it{Private work (CMS simulation)}}")
    latex.SetNDC()
    latex.SetTextFont(42)
    latex.SetTextSize(0.04)
    latex.Draw()

    latex3 = etaLatex(eta)
    latex3.Draw()

    c3.Update()
    c3.SaveAs(os.path.join(outDir, "{}_{}.pdf".format(outTitle, eta)))


# ==================================================================
#   Access to the step1 histograms
# ==================================================================
def signalPath(mass, vsignal, weighted):
    # One file per gluino mass point; the "_weighted" flavour carries the
    # cross-section and luminosity weights.
    return os.path.join(IDIR, "Gluino_Run3_MET_madgraph_{}_V{}{}.root".format(
                              mass, vsignal, "_weighted" if weighted else ""))


def histoName(cuts, eta, region):
    """Same as 'mcTag' of PlottingMacro.py."""
    # Full base name; getVariation() then appends "_" + <variation suffix>.
    # It has to reproduce step1's naming exactly, cuts included.
    return "{}{}_{}_{}_SignalMass".format(HNAME_PREFIX, cuts, eta, region)


def getVariation(ifile, hname, suffix, name):
    # Returns None (and says which histogram is missing) instead of raising, so
    # one absent variation only drops that systematic.
    h = ifile.Get(hname + "_" + suffix)
    if not h:
        print("  '{}_{}' missing in {}".format(hname, suffix, ifile.GetName()))
        return None
    return allSet(h, name)


# ==================================================================
#   Processing of one mass point
# ==================================================================
def processMass(mass, eta, opts):
    # Returns the total systematic histogram, or None if the point was skipped.
    path = signalPath(mass, opts.vsignal, opts.weighted)
    if not os.path.isfile(path):
        print("  missing file: {}".format(path))
        return None

    ifile = TFile.Open(path)
    if (not ifile) or ifile.IsZombie():
        print("  unreadable file: {}".format(path))
        return None

    hname   = histoName(opts.cuts, eta, opts.region)
    nominal = getVariation(ifile, hname, NOMINAL_SUFFIX, "nominal")
    if nominal is None:
        # A missing nominal almost always means the region/eta combination was
        # not produced by step1, hence the explicit hint.
        print("  -> Have you set the {} region and the {} range for step1 correctly?".format(
              opts.region, eta))
        ifile.Close()
        return None

    # One directory per (mass, region, cuts, eta): nothing can overwrite
    # anything else even when several configurations are run in a row.
    outDir = os.path.join(opts.odir, "Gluino_{}_{}{}_{}".format(
                          mass, opts.region, opts.cuts, eta))
    if not os.path.isdir(outDir):
        os.makedirs(outDir)

    # written: saved to the output file
    # inTotal: entering the quadratic sum (here: everything computed)
    # summary: (histogram, legend, colour, marker) for the summary plot
    written, inTotal, summary = [nominal.Clone("mass_nominal")], [], []

    # --- two-sided experimental systematics ---
    for syst in SYSTEMATICS:
        hUp   = getVariation(ifile, hname, syst["up"],   syst["key"] + "_up")
        hDown = getVariation(ifile, hname, syst["down"], syst["key"] + "_down")
        # Both sides are required: keeping a single one would halve the envelope
        # without any visible sign.
        if hUp is None or hDown is None:
            print("  incomplete systematic '{}' -> ignored".format(syst["key"]))
            continue

        h = systMass(nominal, hDown, hUp, "syst_" + syst["key"], binned=1)
        written.append(h)
        inTotal.append(h)
        summary.append((h, syst["legend"], syst["color"], syst["marker"]))

        if not opts.noPlots:
            # nominal is cloned because plotter() restyles what it receives.
            plotter(nominal.Clone(uniq("nomForPlot")), hDown, hUp,
                    "Nominal", "Down", "Up", outDir, syst["title"], mass)

    # --- flat normalisation systematics ---
    # No input histogram needed: they are shaped from the nominal itself.
    for flat in FLAT_SYSTEMATICS:
        h = systFlat(nominal, flat["value"], "syst_" + flat["key"])
        written.append(h)
        inTotal.append(h)
        summary.append((h, flat["legend"], flat["color"], flat["marker"]))

    # Cannot happen in practice since the flat terms always apply, but the guard
    # protects systTotal() from an empty list.
    if not inTotal:
        print("  no systematic computed -> mass {} ignored".format(mass))
        ifile.Close()
        return None

    # 'systTotalBinned' is the key read downstream; same name as on the
    # background side so both can be consumed identically.
    sysTotBin = systTotal(inTotal, "systTotalBinned")
    written.append(sysTotBin)

    outPath = os.path.join(outDir, "sysTotBinned_signal.root")
    ofile = TFile.Open(outPath, "RECREATE")
    if (not ofile) or ofile.IsZombie():
        print("  can't write {}".format(outPath))
        ifile.Close()
        return None
    ofile.cd()
    for h in written:
        h.Write()
    ofile.Close()
    ifile.Close()

    print("  m = {:>4} GeV: max total = {:6.1f} %  -> {}".format(
          mass,
          max(sysTotBin.GetBinContent(i) for i in range(1, sysTotBin.GetNbinsX() + 1)),
          outPath))

    if not opts.noPlots:
        plotSummary(summary, sysTotBin, "Mass (GeV)", outDir, eta, mass)

    return sysTotBin


# ==================================================================
#   Main
# ==================================================================
def main():
    parser = OptionParser(usage="Usage: python3 %prog [options]")
    parser.add_option("--etas", dest="etas", default=ETAS,
                      help="eta range, ex. Eta1,Eta1_2p4,Eta2p4")
    parser.add_option("--region", dest="region", default=REGION,
                      help="8fp9 | 9fp10")
    parser.add_option("--cuts", dest="cuts", default=CUTS,
                      help="step1 selection, ex. _SigmaPtoverPt_0p5_EoP_0p1")
    parser.add_option("--masses", dest="masses", default=None,
                      help="mass points, separated by commas, ex. 1000,1500,2000. Default: {}".format(MASSES))
    parser.add_option("--vsignal", dest="vsignal", default=VSIGNAL,
                      help="signal version")
    parser.add_option("--odir", dest="odir", default=ODIR,
                      help="output directory for the .root and .pdf files")
    # store_false with default=WEIGHTED: --raw flips the module-level setting.
    parser.add_option("--raw", dest="weighted", action="store_false", default=WEIGHTED,
                      help="read the unweighted step1 files (no cross-section nor luminosity applied)")
    parser.add_option("--no-plots", dest="noPlots", action="store_true", default=False,
                      help="only .root output")
    (opts, args) = parser.parse_args()

    # The cut string is a fragment of the histogram name and must start with '_',
    # so both "--cuts _EoP_0p1" and "--cuts EoP_0p1" are accepted.
    if opts.cuts and not opts.cuts.startswith("_"):
        opts.cuts = "_" + opts.cuts

    masses = MASSES
    if opts.masses:
        try:
            masses = [int(m) for m in opts.masses.split(",") if m.strip()]
        except ValueError:
            print("--masses expects integers separated by commas, e.g. 1000,1500,2000")
            sys.exit(1)

    etas = [e.strip() for e in opts.etas.split(",") if e.strip()]
    if not etas:
        print("--etas is empty -> nothing to do")
        sys.exit(1)

    # Echoing the resolved histogram name is the quickest way to diagnose a
    # naming mismatch with step1.
    print("Signal    : V{}{}".format(opts.vsignal, "" if opts.weighted else " (without weights)"))
    print("Histos    : {}_<variation>".format(histoName(opts.cuts, etas[0], opts.region)))
    print("Output    : {}".format(opts.odir))

    total = 0
    for eta in etas:
        print("\n=========== {} / {} ===========".format(eta, opts.region))
        # Collected per eta range, so the recap plot compares mass points at
        # fixed eta.
        sysTot_allMasses = OrderedDict()
        for mass in masses:
            h = processMass(mass, eta, opts)
            if h is None:
                continue
            total += 1
            # Cloned because filledRange() later restyles and truncates it.
            if mass in HIGHLIGHT:
                sysTot_allMasses[mass] = h.Clone(uniq("sysTot_{}".format(mass)))

        if not opts.noPlots:
            if not os.path.isdir(opts.odir):
                os.makedirs(opts.odir)
            plotTotalAllMasses(sysTot_allMasses, "Mass (GeV)", opts.odir, eta)

    print("\n{} mass point(s) done.".format(total))
    # Non-zero exit when nothing at all was produced.
    if total == 0:
        sys.exit(1)


if __name__ == "__main__":
    main()