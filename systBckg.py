#!/usr/bin/env python3
"""Combination of the background prediction systematics.

Reads one .root per systematic label (written by BkgPrediction.C, filed by
LaunchBkgPred.py) and produces the file expected by ShowPlots.py:

    <indir>/<systdir>/sysTotBinned_<eta>_<region>.root   -> key 'systTotalBinned'

    python3 systBckg.py
    python3 systBckg.py --etas Eta1,Eta1_2p4,Eta2p4
    python3 systBckg.py --region 9fp10
    python3 systBckg.py --nominal-only
    python3 systBckg.py --only Eta,Ih          # subset of systematics
"""

# =============================================================================
#  Step 3a — systematic envelope of the background prediction.
# -----------------------------------------------------------------------------
#  Input: one .root per systematic label, produced by BkgPrediction.C and filed
#          by LaunchBkgPred.py. Each file holds the predicted mass spectrum under
#          one variation, as "mass_predBC_<region>".
#  Output: a single .root per eta range holding, bin by bin, the relative
#          deviation of every variation with respect to the nominal, plus their
#          quadratic sum under the key 'systTotalBinned' that ShowPlots.py reads.
#
#  The comparison is done on *shapes*: every spectrum is rebinned onto the
#  analysis binning, its under/overflow folded in, and normalised to unit area.
#  The overall normalisation is not a systematic here, it is fixed by the ABCD
#  factor and cancels in the ratio.
# =============================================================================

import os, re, sys, math, array, ctypes
from optparse import OptionParser

import ROOT
from ROOT import TFile, TCanvas, TLegend, TLatex, TPad, TLine
import tdrstyle

ROOT.gROOT.SetBatch(True)              # no X11: the script only writes files
ROOT.gErrorIgnoreLevel = ROOT.kWarning
# Histograms are not owned by any TFile, so they stay valid after f.Close().
# Without this, every histogram read below would be destroyed with its file.
ROOT.TH1.AddDirectory(False)
tdrstyle.setTDRStyle()


# ==================================================================
#   Settings: MUST match those of ShowPlots.py
# ==================================================================
# These reproduce, by hand, the path convention built by LaunchBkgPred.py.
BASE       = "/safe/ui3_1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/macros"
#DATASET   = "/safe/ui3_1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/output/JetMET2024_V12/JetMET2024_V12p35"
DATASET    = "/safe/ui3_1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/output/HistForBkg_MC_V3"
SAMPLETYPE = "mc2024"     # data2017|data2018|data2024|mc2017|mc2018|mc2024
# SUFFIX must reproduce the launcher's final labelDir, i.e. including the cut
# prefix the launcher prepends when both eopCut and sigmaPtCut are set.
SUFFIX     = "v2"
# CUTS is the step1 selection as it appears inside the *file name*, not in the
# directory name; it is inserted between the eta range and the systematic label.
CUTS       = ""           # step1 cuts: "" | "_SigmaPtoverPt_0p5_EoP_0p1" | ...
REGION     = "8fp9"       # 8fp9 (VR) | 9fp10 (SR)
SYSTDIR    = "SystCombined"
YEAR       = "2024"
ERA        = ""           # "" (all year) | "F" | "G"

# Integrated luminosity (fb-1) per era, for the plot label only.
LUMI = {"": 109.0, "F": 25.40, "G": 34.4}

# Prefix of the prediction histogram written by BkgPrediction.C
# ("mass_predBC_<region>"). The region suffix is appended at read time.
PLOTTYPE = "mass_predBC_"

# Convention for the variation / nominal ratio.
#   True : |var/nom - 1|, for the cumulative form AND the bin-by-bin form
#   False: reproduces the old code, which used nom/var for the cumulative form
#           and var/nom bin by bin (the two only agree to first order)
# Kept as a switch so an old result can be reproduced exactly; True is correct.
RATIO_VAR_OVER_NOM = True

# Analysis mass binning. Must stay identical to the xbins array of rebinHisto()
# in CommonFunctions.h, otherwise the systematics and the prediction they apply
# to would live on different binnings.
REBINNING = array.array('d', [0., 20., 40., 60., 80., 100., 120., 140., 160., 180.,
                              200., 220., 240., 260., 280., 300., 320., 340., 360., 380.,
                              410., 440., 480., 530., 590., 660., 760., 880., 1030., 1210.,
                              1440., 1730., 2000., 2500., 3200., 4000.])
SIZE_REBINNING = len(REBINNING) - 1

MAX_MASS = 4000


# ==================================================================
#   Systematics table
# ==================================================================
# 'down' / 'up' are the labels from LaunchBkgPred.py, they directly give the
# file name. up = None -> one-sided variation.
# inTotal = False -> computed and plotted, but left out of the quadratic sum.
#
# One-sided entries are symmetrised inside systMass(): the single variation is
# used as both sides, so max(|up|, |down|) reduces to that one deviation.
#
# NoFit is deliberately excluded from the total: it is a cross-check of how much
# the tail fits contribute at all, not an independent source of uncertainty.
SYSTEMATICS = [
    dict(key="Eta",      down="binEtaDown",      up="binEtaUp",
         legend="#eta binning",       legDown="#eta down",  legUp="#eta up",
         color=ROOT.kMagenta - 9, marker=21, inTotal=True),
    dict(key="Ih",       down="binIhDown",       up="binIhUp",
         legend="I_{h} binning",      legDown="I_{h} down", legUp="I_{h} up",
         color=ROOT.kViolet + 1,  marker=22, inTotal=True),
    dict(key="P",        down="binMomDown",      up="binMomUp",
         legend="p binning",          legDown="p down",     legUp="p up",
         color=ROOT.kBlue + 1,    marker=23, inTotal=True),
    dict(key="FitIh",    down="fitIhDown",       up="fitIhUp",
         legend="I_{h} fit",          legDown="Fit I_{h} down", legUp="Fit I_{h} up",
         color=ROOT.kOrange + 1,  marker=33, inTotal=True),
    dict(key="FitP",     down="fitMomDown",      up="fitMomUp",
         legend="p fit",              legDown="Fit p down", legUp="Fit p up",
         color=ROOT.kCyan,        marker=29, inTotal=True),
    dict(key="NoFit",    down="noFit",           up=None,
         legend="No fit",             legDown="No fit",     legUp="",
         color=ROOT.kGreen + 2,   marker=20, inTotal=False),
    dict(key="CorrIh",   down="corrTemplateIh",  up=None,
         legend="corr template I_{h}", legDown="corr template I_{h}", legUp="",
         color=ROOT.kOrange,      marker=39, inTotal=True),
    dict(key="Corr1oP",  down="corrTemplate1oP", up=None,
         legend="corr template 1/p",  legDown="corr template 1/p",   legUp="",
         color=ROOT.kGreen + 2,   marker=30, inTotal=True),
]

# Label of the reference run; must match the first entry of the launcher's list.
NOMINAL_LABEL = "nominal"
# Style of the statistical uncertainty, which is drawn alongside the systematics
# and, unlike NoFit, does enter the total.
STAT_STYLE    = dict(legend="Stat.", color=ROOT.kBlack, marker=20)


# ==================================================================
#   Histogram helpers
# ==================================================================
# Global counter feeding uniq(); a list is used so the closure can mutate it
# without a `global` statement.
_UID = [0]


def uniq(stem):
    """Unique name: avoid 'Replacing existing TH1' in ROOT."""
    # ROOT keeps a global name registry: reusing a name silently replaces the
    # previous object, which here would corrupt histograms still in use.
    _UID[0] += 1
    return "{}_{}".format(stem, _UID[0])


def overflowInLastBin(h):
    # Folds the overflow into the last visible bin and empties it, so that no
    # entry is lost when the spectrum is later integrated over 1..N.
    # Errors are added in quadrature.
    res = h.Clone(uniq(h.GetName() + "_ovf"))
    n = h.GetNbinsX()
    res.SetBinContent(n, h.GetBinContent(n) + h.GetBinContent(n + 1))
    res.SetBinError(n, math.sqrt(h.GetBinError(n)**2 + h.GetBinError(n + 1)**2))
    res.SetBinContent(n + 1, 0)
    res.SetBinError(n + 1, 0)
    return res


def underflowInFirstBin(h):
    # Same for the underflow. Zeroing bins 0 and N+1 matters downstream: several
    # loops below run over 0..N+1 and would otherwise pick up raw yields.
    res = h.Clone(uniq(h.GetName() + "_udf"))
    res.SetBinContent(1, h.GetBinContent(0) + h.GetBinContent(1))
    res.SetBinError(1, math.sqrt(h.GetBinError(0)**2 + h.GetBinError(1)**2))
    res.SetBinContent(0, 0)
    res.SetBinError(0, 0)
    return res


def allSet(h, name, normalise=True):
    """Rebin, fold under/overflow, normalisation to 1."""
    # Standard preparation applied to every spectrum read from disk, so that the
    # nominal and all the variations are strictly comparable.
    # The unit normalisation is what turns the comparison into a pure shape
    # comparison: the ABCD normalisation is common to all variations and would
    # otherwise leak into the systematic.
    res = h.Rebin(SIZE_REBINNING, uniq(name + "_reb"), REBINNING)
    res = overflowInLastBin(res)
    res = underflowInFirstBin(res)
    res.SetName(name)
    res.SetDirectory(0)
    if normalise:
        integral = res.Integral()
        if integral > 0:
            res.Scale(1.0 / integral)
        else:
            print("  /!\\ empty histogram: '{}' not normalised".format(name))
    return res


def ratioInt(num, den):
    """num/den of cumulated integrals [i, overflow], bin by bin, with errors."""
    # Cut-and-count form: bin i holds the ratio of the yields above the lower
    # edge of bin i, which is the quantity a mass threshold actually selects.
    # The error assumes num and den are independent; they are not (both derive
    # from the same data), so this is an over-estimate, i.e. conservative.
    res = num.Clone(uniq("ratioInt"))
    res.SetDirectory(0)
    res.Reset()
    last = num.GetNbinsX() + 1
    for i in range(0, num.GetNbinsX() + 1):
        e1 = ctypes.c_double(0.0)
        e2 = ctypes.c_double(0.0)
        # ctypes doubles are required: IntegralAndError writes the error through
        # a C++ reference, which a plain Python float cannot provide.
        a = num.IntegralAndError(i, last, e1, "")
        b = den.IntegralAndError(i, last, e2, "")
        if a != 0 and b != 0:
            res.SetBinContent(i, a / b)
            res.SetBinError(i, math.sqrt((e1.value / a)**2 + (e2.value / b)**2) * a / b)
        else:
            res.SetBinContent(i, 0)
    return res


def ratioBinned(num, den):
    """num/den bin by bin."""
    # Differential form. TH1::Divide propagates the errors assuming independence.
    res = num.Clone(uniq("ratioBin"))
    res.SetDirectory(0)
    res.Divide(den)
    return res


def statErr(h, name):
    """Statistical relative error (%) bin by bin."""
    # The input error bars are the toy-to-toy RMS produced by MeanOfToys, so this
    # is the statistical uncertainty of the prediction, not of the observation.
    res = h.Clone(name)
    res.SetDirectory(0)
    for i in range(0, h.GetNbinsX() + 2):
        c = h.GetBinContent(i)
        res.SetBinContent(i, 100.0 * h.GetBinError(i) / c if c > 0 else 0.0)
        res.SetBinError(i, 0.0)
    return res


def statErrRInt(h, name):
    """Statistical relative error (%) on the cumulated integral [i, overflow]."""
    # Cumulative counterpart of statErr, matching the cut-and-count use case.
    res = h.Clone(name)
    res.SetDirectory(0)
    last = h.GetNbinsX() + 1
    for i in range(0, h.GetNbinsX() + 2):
        e = ctypes.c_double(0.0)
        c = h.IntegralAndError(i, last, e, "")
        res.SetBinContent(i, 100.0 * e.value / c if c > 0 else 0.0)
        res.SetBinError(i, 0.0)
    return res


def systMass(nominal, down, up, name, binned, typec=0, mini=0):
    """Relative difference (%) between nominal and the variations.

    binned=0: on the cumulated integral   binned=1: bin by bin
    typec=0 : max(|up|,|down|)            typec=1 : mean of the two
    mini=1  : min instead of max          mini=2  : max/2
    up=None : variation on one side

    Bins where the nominal value is empty are set to 0: the ratio there is 0/0,
    which would result in a purely artificial deviation of 100%.
    """
    # One-sided case: `var` falls back to `down`, so ra1 and ra2 are identical
    # and the max/min/mean combinations all collapse to that single deviation.
    var = up if up is not None else down
    if binned:
        ra1 = ratioBinned(var,  nominal)
        ra2 = ratioBinned(down, nominal)
    elif RATIO_VAR_OVER_NOM:
        ra1 = ratioInt(var,  nominal)
        ra2 = ratioInt(down, nominal)
    else:
        ra1 = ratioInt(nominal, var)
        ra2 = ratioInt(nominal, down)

    last = nominal.GetNbinsX() + 1
    res  = ra1.Clone(name)
    res.SetDirectory(0)
    res.Reset()
    # 0..N+1: under/overflow are covered, and are zero thanks to allSet.
    for i in range(0, res.GetNbinsX() + 2):
        # The reference is taken on the same footing as the ratio (cumulative or
        # differential), so the emptiness test matches the quantity being used.
        ref = nominal.GetBinContent(i) if binned else nominal.Integral(i, last)
        if ref <= 0:
            res.SetBinContent(i, 0.0)
            continue
        s1 = abs(1 - ra1.GetBinContent(i))
        s2 = abs(1 - ra2.GetBinContent(i))
        if typec == 1:
            m = 0.5 * (s1 + s2)
        else:
            m = min(s1, s2) if mini == 1 else max(s1, s2)
            if mini == 2:
                m /= 2.0
        # Stored as a percentage; the error bar is meaningless on an envelope
        # and is zeroed so it cannot be drawn or propagated by mistake.
        res.SetBinContent(i, 100.0 * m)
        res.SetBinError(i, 0.0)
    return res


def systTotal(list_h, name):
    """Sum in quadrature, bin by bin."""
    # Assumes the sources are uncorrelated. The binning of every entry must be
    # identical, which allSet guarantees.
    res = list_h[0].Clone(name)
    res.SetDirectory(0)
    res.Reset()
    for i in range(0, res.GetNbinsX() + 2):
        tot = sum(h.GetBinContent(i)**2 for h in list_h)
        res.SetBinContent(i, math.sqrt(tot))
        res.SetBinError(i, 0.0)
    return res


def lowEdge(h):
    """TGraph (bin low edge, content): avoids the histogram staircase."""
    # On a strongly non-uniform binning the staircase of a TH1 is misleading;
    # one marker per bin, placed at its lower edge, reads much better.
    g = ROOT.TGraph(h.GetNbinsX())
    for i in range(1, h.GetNbinsX() + 1):
        g.SetPoint(i - 1, h.GetBinLowEdge(i), h.GetBinContent(i))
    return g


def setColorAndMarker(obj, color, markerstyle):
    # Returns the object so it can be used inline; note that it *mutates* its
    # argument, hence the .Clone() calls at the call sites that need to keep the
    # original style intact.
    obj.SetLineColor(color)
    obj.SetMarkerColor(color)
    obj.SetFillColor(color)
    obj.SetMarkerStyle(markerstyle)
    obj.SetMarkerSize(1.2)
    return obj


def etaLatex(eta):
    # Plot label for the eta range. The x position is tuned per label so that
    # the longer strings still fit inside the pad. Unknown ranges fall back to
    # the raw name rather than failing.
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


# ==================================================================
#   Plots
# ==================================================================
def plotter(nominal, down, up, legNom, legDown, legUp, outDir, outTitle, lumiLabel):
    """Nominal + variations, with ratios."""
    # Three stacked pads: the spectra (log y), the cumulative ratio, and the
    # bin-by-bin ratio. Showing both ratios makes it obvious when a variation
    # only moves the tail without changing the integral.
    c1 = TCanvas(uniq("c1"), "c1", 800, 800)
    t1 = TPad(uniq("t1"), "t1", 0.0, 0.40, 0.95, 0.95)
    t1.Draw()
    t1.SetLogy(1)
    t1.SetGrid(1)
    t1.SetTopMargin(0.005)
    t1.SetBottomMargin(0.005)
    c1.cd()

    t2 = TPad(uniq("t2"), "t2", 0.0, 0.225, 0.95, 0.375)
    t2.Draw()
    t2.SetGridy(1)
    t2.SetTopMargin(0.05)
    t2.SetBottomMargin(0.1)

    t3 = TPad(uniq("t3"), "t3", 0.0, 0., 0.95, 0.20)
    t3.Draw()
    t3.SetGridy(1)
    t3.SetBottomMargin(0.45)

    t1.cd()
    # Fixed floor and a headroom factor above the nominal maximum, so that the
    # log-scale range is comparable from one systematic to the next.
    min_entries = 1e-7
    max_entries = nominal.GetMaximum() * 5

    nominal = setColorAndMarker(nominal, 1, 20)
    nominal.GetXaxis().SetRangeUser(0, MAX_MASS)
    nominal.GetYaxis().SetRangeUser(min_entries, max_entries)
    nominal.SetMinimum(min_entries)
    nominal.SetTitle(";Mass (GeV);Normalized tracks")
    nominal.GetYaxis().SetTitleSize(0.07)
    nominal.GetYaxis().SetLabelSize(0.06)
    nominal.GetYaxis().SetTitleOffset(0.9)
    nominal.Draw()

    down = setColorAndMarker(down, 38, 21)
    down.Draw("same")
    # One-sided systematics have no `up` histogram; every use of it is guarded.
    if up is not None:
        up = setColorAndMarker(up, 46, 21)
        up.Draw("same")

    leg = TLegend(0.16, 0.75, 0.35, 0.99)
    leg.AddEntry(nominal, legNom, "PE1")
    leg.AddEntry(down, legDown, "PE1")
    if up is not None:
        leg.AddEntry(up, legUp, "PE1")
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.Draw("same")

    c1.cd()
    t2.cd()

    # Empty frame defining the axes of the ratio pad; the ratios themselves are
    # then drawn on top with "same". Drawing a frame rather than the first ratio
    # keeps the y range identical across all systematics.
    frameR2 = ROOT.TH1D(uniq("frameR2"), "frameR2", 1, 0, MAX_MASS)
    frameR2.SetDirectory(0)
    frameR2.GetXaxis().SetNdivisions(505)
    frameR2.SetTitle("")
    frameR2.SetStats(0)
    frameR2.GetXaxis().SetTitle("")
    frameR2.GetYaxis().SetTitle("RatioR ")
    frameR2.GetXaxis().SetRangeUser(0, MAX_MASS)
    frameR2.SetMaximum(1.1)
    frameR2.SetMinimum(0.9)
    # Font code 43 = absolute sizes in pixels, so the labels keep the same size
    # in the small pads as in the big one.
    frameR2.GetYaxis().SetLabelFont(43)
    frameR2.GetYaxis().SetLabelSize(20)
    frameR2.GetYaxis().SetTitleFont(43)
    frameR2.GetYaxis().SetTitleSize(20)
    frameR2.GetYaxis().SetNdivisions(205)
    frameR2.GetYaxis().SetTitleOffset(2)
    frameR2.GetXaxis().SetNdivisions(510)
    frameR2.GetXaxis().SetLabelFont(43)
    frameR2.GetXaxis().SetLabelSize(0)     # x labels only on the bottom pad
    frameR2.GetXaxis().SetTitleFont(43)
    frameR2.GetXaxis().SetTitleSize(24)
    frameR2.GetXaxis().SetTitleOffset(0)
    frameR2.Draw("AXIS")
    frameR2.Draw("SAME AXIG")

    lineAtOne = TLine(0, 1, MAX_MASS, 1)
    lineAtOne.SetLineStyle(3)
    lineAtOne.SetLineColor(1)
    lineAtOne.Draw("same")

    # `keep` holds a Python reference to every drawn object: without it the
    # garbage collector destroys them before SaveAs and the pads come out empty.
    keep = []
    r1 = setColorAndMarker(ratioInt(down, nominal), 38, 21)
    r1.Draw("E0 same")
    keep.append(r1)
    if up is not None:
        r2 = setColorAndMarker(ratioInt(up, nominal), 46, 21)
        r2.Draw("E0 same")
        keep.append(r2)

    c1.cd()
    t3.cd()

    # The bottom pad reuses the same frame, with the x labels and title enabled.
    frameR3 = frameR2.Clone(uniq("frameR3"))
    frameR3.SetDirectory(0)
    frameR3.GetXaxis().SetLabelSize(20)
    frameR3.GetYaxis().SetTitleOffset(2.1)
    frameR3.GetYaxis().SetTitle("#frac{var}{Nominal}")
    frameR3.GetXaxis().SetRangeUser(0, MAX_MASS)
    frameR3.Draw("AXIS")
    frameR3.Draw("SAME AXIG")
    frameR3.GetXaxis().SetTitle("Mass (GeV)")
    frameR3.GetXaxis().SetTitleOffset(1)
    lineAtOne.Draw("same")

    b1 = setColorAndMarker(ratioBinned(down, nominal), 38, 21)
    b1.Draw("E0 same")
    keep.append(b1)
    if up is not None:
        b2 = setColorAndMarker(ratioBinned(up, nominal), 46, 21)
        b2.Draw("E0 same")
        keep.append(b2)

    c1.cd()
    latex = TLatex(0.15, 0.955, "#scale[1.3]{#it{Private work (CMS data)}}")
    latex.SetNDC()
    latex.SetTextFont(42)
    latex.SetTextSize(0.03)
    latex.Draw()

    latex2 = TLatex(0.70, 0.955, lumiLabel)
    latex2.SetNDC()
    latex2.SetTextFont(42)
    latex2.SetTextSize(0.03)
    latex2.Draw()

    c1.SaveAs(os.path.join(outDir, outTitle + ".pdf"))


def plotSummary(entries, total, xtitle, outDir, outTitle, eta, sampleTag, lumiLabel):
    """entries: list of (histo, legend, color, marker)."""
    # All systematics on one log-scale plot, in percent, with the total in red.
    # Note that `entries` contains every computed source, including those flagged
    # inTotal=False: they are shown for information but are not in `total`.
    c2 = TCanvas(uniq("c2"), "c2", 800, 600)
    c2.SetBottomMargin(0.12)
    c2.SetLeftMargin(0.1)
    c2.SetLogy()
    c2.SetGrid()

    frame = ROOT.TH1D(uniq("frameSummary"), "", 1, 0, MAX_MASS)
    frame.SetDirectory(0)
    frame.SetStats(0)
    # 0.1 % to 2000 %: wide enough that a pathological variation stays on-plot.
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

    leg = TLegend(0.12, 0.7, 0.5, 0.93)
    leg.SetNColumns(2)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)

    # `graphs` keeps the TGraphs alive until SaveAs, same reason as `keep` above.
    graphs = []
    gTot = setColorAndMarker(lowEdge(total), ROOT.kRed, 34)
    leg.AddEntry(gTot, "Total", "PE1")
    for (h, legend, color, marker) in entries:
        g = setColorAndMarker(lowEdge(h), color, marker)
        leg.AddEntry(g, legend, "PE1")
        g.Draw("P")
        graphs.append(g)
    gTot.Draw("P")                 # drawn last so it sits on top
    leg.Draw("same")

    latex = TLatex(0.1, 0.96, "#scale[1.3]{#it{Private work (CMS data)}}")
    latex.SetNDC()
    latex.SetTextFont(42)
    latex.SetTextSize(0.04)
    latex.Draw()

    latex2 = TLatex(0.735, 0.96, lumiLabel)
    latex2.SetNDC()
    latex2.SetTextFont(42)
    latex2.SetTextSize(0.04)
    latex2.Draw()

    latex3 = etaLatex(eta)
    latex3.Draw()

    for sub in ("pdf", "root", "Cfile"):
        d = os.path.join(outDir, sub)
        if not os.path.isdir(d):
            os.makedirs(d)

    # Three formats: pdf to look at, .root and .C so the plot can be reopened
    # and restyled without rerunning the whole chain.
    stem = "summary_{}_{}_{}".format(outTitle, sampleTag, eta)
    c2.SaveAs(os.path.join(outDir, "pdf",   stem + ".pdf"))
    c2.SaveAs(os.path.join(outDir, "root",  stem + ".root"))
    c2.SaveAs(os.path.join(outDir, "Cfile", stem + ".C"))


# ==================================================================
#   Paths (same rules as ShowPlots.py)
# ==================================================================
def getVersion(stem):
    # "..._V12p35" -> "12p35". Slightly looser than the launcher's regex: the
    # version may be followed by an underscore instead of ending the name.
    m = re.search(r"_V([0-9A-Za-z]+?)(?:_|$)", stem)
    if not m:
        print("Can't extract version of '{}' ('_V' expected)".format(stem))
        sys.exit(1)
    return m.group(1)


def predictionPath(indir, eta, stem, cuts, label):
    # Mirrors the output name built by BkgPrediction.C:
    #     <dataset stem> + "_" + <eta range> + <cuts> + "_" + <label>
    # inside the per-eta sub-directory created by the launcher.
    return "{}/{}/{}_{}{}_{}.root".format(indir, eta, stem, eta, cuts, label)


def getPrediction(path, region, name):
    """Open, take mass_predBC_<region>, rebin, normalise, close."""
    # Every failure mode is reported and returns None rather than raising, so a
    # single missing variation does not abort the whole eta range.
    if not os.path.isfile(path):
        print("  missing file: {}".format(path))
        return None
    f = TFile.Open(path)
    if (not f) or f.IsZombie():
        print("  unreadable file: {}".format(path))
        return None
    h = f.Get(PLOTTYPE + region)
    if not h:
        print("  '{}{}' missing in {}".format(PLOTTYPE, region, path))
        f.Close()
        return None
    # allSet detaches the result, which is why closing the file here is safe.
    res = allSet(h, name)
    f.Close()
    return res


# ==================================================================
#   Processing of one eta range
# ==================================================================
def processEta(eta, indir, stem, opts, lumiLabel):
    # Returns True on success. One output file per (eta, region).
    print("\n=========== {} / {} ===========".format(eta, opts.region))

    # Without the nominal there is nothing to compare against: bail out early.
    nominal = getPrediction(predictionPath(indir, eta, stem, opts.cuts, NOMINAL_LABEL),
                            opts.region, "nominal_def")
    if nominal is None:
        print("  can't find nominal -> {} ignored".format(eta))
        return False

    outDir = os.path.join(indir, opts.systdir)
    if not os.path.isdir(outDir):
        os.makedirs(outDir)
    # This exact name and location is what ShowPlots.py expects.
    outPath = os.path.join(outDir, "sysTotBinned_{}_{}.root".format(eta, opts.region))

    # --- statistical uncertainty of the nominal ---
    # It enters the total on the same footing as the systematics, so the key
    # 'systTotalBinned' is really a *total* uncertainty, not a purely systematic
    # one. Keep that in mind when combining it downstream.
    syst_stat        = statErrRInt(nominal, "Stat")
    syst_stat_binned = statErr(nominal, "Stat_binned")

    # written   : everything saved to the output file
    # inTotal   : cumulative-form sources entering the quadratic sum
    # inTotalBin: bin-by-bin counterparts
    # summary   : (histogram, legend, colour, marker) tuples for the plot
    written    = [nominal.Clone("mass_predBC_nominal"), syst_stat, syst_stat_binned]
    inTotal    = [syst_stat]
    inTotalBin = [syst_stat_binned]
    summary    = [(syst_stat_binned, STAT_STYLE["legend"],
                   STAT_STYLE["color"], STAT_STYLE["marker"])]

    if not opts.nominalOnly:
        plotDir = os.path.join(outDir, "individualSyst", eta)
        if not os.path.isdir(plotDir):
            os.makedirs(plotDir)

        for syst in SYSTEMATICS:
            # --only restricts the loop to a subset of keys.
            if opts.only and syst["key"] not in opts.only:
                continue

            hDown = getPrediction(predictionPath(indir, eta, stem, opts.cuts, syst["down"]),
                                  opts.region, syst["key"] + "_down")
            hUp = None
            if syst["up"] is not None:
                hUp = getPrediction(predictionPath(indir, eta, stem, opts.cuts, syst["up"]),
                                    opts.region, syst["key"] + "_up")

            # A two-sided systematic missing one side is skipped entirely: using
            # only the surviving side would silently halve the envelope.
            if hDown is None or (syst["up"] is not None and hUp is None):
                print("  systematic '{}' missing -> ignored".format(syst["key"]))
                continue

            # Both forms are computed and stored: the cumulative one feeds the
            # cut-and-count interpretation, the binned one the shape fit.
            h_int = systMass(nominal, hDown, hUp, syst["key"], binned=0)
            h_bin = systMass(nominal, hDown, hUp, syst["key"] + "_binned", binned=1)
            written += [h_int, h_bin]
            if syst["inTotal"]:
                inTotal.append(h_int)
                inTotalBin.append(h_bin)
            # Plotted regardless of inTotal, so excluded sources stay visible.
            summary.append((h_bin, syst["legend"], syst["color"], syst["marker"]))

            print("  {:<10s} max = {:6.1f} %  (bin by bin)".format(
                  syst["key"], max(h_bin.GetBinContent(i)
                                   for i in range(1, h_bin.GetNbinsX() + 1))))

            if not opts.noPlots:
                # The nominal is cloned because plotter() restyles what it gets.
                plotter(nominal.Clone(uniq("nomForPlot")), hDown, hUp,
                        "Nominal", syst["legDown"], syst["legUp"],
                        plotDir, "plot_" + syst["key"], lumiLabel)

    # 'systTotalBinned' is the key ShowPlots.py reads; do not rename it.
    sysTot     = systTotal(inTotal,    "systTotal")
    sysTotBin  = systTotal(inTotalBin, "systTotalBinned")
    written   += [sysTot, sysTotBin]

    ofile = TFile.Open(outPath, "RECREATE")
    if (not ofile) or ofile.IsZombie():
        print("  can't write {}".format(outPath))
        return False
    ofile.cd()
    for h in written:
        h.Write()
    ofile.Close()
    print("  -> {}".format(outPath))

    if not opts.noPlots:
        plotSummary(summary, sysTotBin, "Mass (GeV)", outDir,
                    "binned_syst", eta, stem, lumiLabel)
    return True


# ==================================================================
#   Main
# ==================================================================
def main():
    # Every module-level constant can be overridden from the command line, so a
    # one-off run never requires editing the file.
    parser = OptionParser(usage="Usage: python3 %prog [options]")
    parser.add_option("--etas", dest="etas", default="Eta2p4",
                      help="eta range, separated using commas")
    parser.add_option("--region", dest="region", default=REGION,
                      help="8fp9 | 9fp10")
    parser.add_option("--suffix", dest="suffix", default=SUFFIX,
                      help="override SUFFIX")
    parser.add_option("--cuts", dest="cuts", default=None,
                      help="override CUTS")
    parser.add_option("--indir", dest="indir", default=None,
                      help="override the working directory")
    parser.add_option("--systdir", dest="systdir", default=SYSTDIR,
                      help="output directory")
    parser.add_option("--only", dest="only", default=None,
                      help="systematic keys, separated using commas")
    parser.add_option("--nominal-only", dest="nominalOnly", action="store_true",
                      default=False, help="only nominal computed")
    parser.add_option("--no-plots", dest="noPlots", action="store_true",
                      default=False, help="only the .root in output")
    parser.add_option("--era", dest="era", default=ERA, help="'' | F | G")
    (opts, args) = parser.parse_args()

    # None (flag absent) means "use CUTS"; an explicit "" means "no cut", which
    # is why the default is None rather than CUTS itself.
    opts.cuts = opts.cuts if opts.cuts is not None else CUTS
    # The cut string is a file-name fragment and must start with '_'; adding it
    # here makes both "--cuts _EoP_0p1" and "--cuts EoP_0p1" work.
    if opts.cuts and not opts.cuts.startswith("_"):
        opts.cuts = "_" + opts.cuts

    # Validate --only against the table before doing any work, and print the
    # valid keys so the user does not have to open the file.
    if opts.only:
        keep    = [s.strip() for s in opts.only.split(",") if s.strip()]
        known   = [s["key"] for s in SYSTEMATICS]
        unknown = [k for k in keep if k not in known]
        if unknown:
            print("Unknown key(s): " + ", ".join(unknown))
            print("Available: " + ", ".join(known))
            sys.exit(1)
        opts.only = keep

    if opts.era not in LUMI:
        print("Unknown era: '{}' (expected: {})".format(opts.era, ", ".join(repr(e) for e in LUMI)))
        sys.exit(1)
    lumiLabel = "{:g} fb^{{-1}} (13.6 TeV)".format(LUMI[opts.era])

    # stem is the dataset base name, used both to find the version and as the
    # prefix of every prediction file.
    stem = os.path.basename(DATASET)
    if opts.indir:
        # --indir bypasses the whole convention; the version is then unknown and
        # only used for display.
        indir, version = opts.indir, "?"
    else:
        version = getVersion(stem)
        indir   = "{}/{}_V{}__{}_{}".format(BASE, SAMPLETYPE, version,
                                            opts.region, opts.suffix)

    print("Repository: {}".format(indir))
    print("Prefix    : {}   version {}   region {}".format(stem, version, opts.region))
    if opts.nominalOnly:
        print("/!\\ --nominal-only: no systematics will be computed")

    etas = [e.strip() for e in opts.etas.split(",") if e.strip()]
    if not etas:
        print("--etas is empty -> nothing to do")
        sys.exit(1)
    # Each eta range is independent: a failure on one does not stop the others.
    done = [eta for eta in etas if processEta(eta, indir, stem, opts, lumiLabel)]

    print("\n{}/{} eta range(s) done.".format(len(done), len(etas)))
    # Non-zero exit only when nothing at all was produced.
    if not done:
        sys.exit(1)


if __name__ == "__main__":
    main()