"""Runs the background estimate.

    python3 LaunchBkgPred.py              -> only nominal
    python3 LaunchBkgPred.py --all        -> all systematics
    python3 LaunchBkgPred.py --only etaup,ihdown
"""

# =============================================================================
#  Driver for the background estimate.
# -----------------------------------------------------------------------------
#  For every (dataset, systematic variation) pair it:
#    1. merges three layers of settings into a flat dict
#       (common session settings -> nominal binning/fit values -> the variation);
#    2. writes that dict as a "key = value" config file;
#    3. runs BkgPrediction.C on it through ROOT;
#    4. collects the freshly written .root files and moves them to a structured
#       output directory.
#
#  Design rule: one config file = one run = one output file, with the systematic
#  label baked into the output name. Nothing in the C++ ever has to be edited to
#  change a variation.
# =============================================================================

import os, sys, glob, time, re, shutil
from optparse import OptionParser

# --- Command line ---------------------------------------------------------
# Default (no flag) runs the nominal only, which is the fast path used while
# debugging; --all is the full systematic sweep.
parser = OptionParser(usage="Usage: python %prog [options]")
parser.add_option("--all",  dest="runAll", action="store_true", default=False,
                  help="all systematics (default: nominal)")
parser.add_option("--only", dest="only", default=None,
                  help="labels, separated by commas")
(opt, args) = parser.parse_args()

# --- Inputs ---------------------------------------------------------------
# (path without the .root extension, sample type key for the dE/dx calibration).
# The second element is what selects (K, C) in GetDeDxCalib, so it must track the
# first: reading MC with the data calibration silently biases the mass spectrum.
datasetList = [
    #("/safe/ui3_1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/output/JetMET2024_V12/JetMET2024_V12p35", "data2024"),
    #("/safe/ui3_1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/output/HistForBkg_MC_V3", "mc2024"),
    ("/safe/ui3_1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/output/TTbar2024_V15/TTbar2024_V15p10", "ttbar2024"),
]

# Root of the directory tree the outputs are filed into, plus the pieces used to
# build the sub-directory names.
outputDir      = "/safe/ui3_1/cms/gcoulon/CMSSW_15_0_13_patch1/src/outputHist"
etaRangeName   = "Eta2p4"
sampleTypeName = "ttbar2024"
labelDir       = "v2"

# Settings shared by every launch of this session.
common = dict(
    # Number of pseudo-experiments (toys). Any integer > 0. 200 = nominal.
    # Run time scales linearly with it: use 20 for a quick test.
    nPE             = 200,

    # Apply the rebinning. 0|1. At 0 the step1 histograms are taken as they are
    # (rebinEta / rebinIh / rebinMom are then ignored).
    rebin           = 1,

    # dE/dx calibration (K, C) used for M = sqrt((Ih-C)/K) * p.
    # data2017 | data2018 | data2024 | mc2017 | mc2018 | mc2024 | ttbar2024
    # Must match the dataset: MC read with the data calibration gives a wrong
    # mass spectrum, with no warning whatsoever.
    sampleType      = sampleTypeName,

    # Eta range; must exist on the step1 side.
    # Eta1 | Eta1_2p4 | Eta2p4 | Eta1p2_2p2 | Eta1p2_2p4
    etaRange        = etaRangeName,

    # step1 E/p cut. "" (no cut) | "0p1"
    eopCut          = "",

    # step1 sigma(pT)/pT cut. "" (no cut) | "0p5"
    sigmaPtCut      = "",

    # Shape of the Ih template fit. 0|1
    #   0 = Gaussian fitted from 1.1 * max
    #   1 = legacy fit, fixed start at Ih = 3 MeV/cm
    useOldIhFit     = 0,

    # Shape of the 1/p template fit. 0|1
    #   0 = 0.5*(exp(ax^2+bx) + exp(-ax^2-bx)) - 1, range [0, 0.6*peak]
    #   1 = [0]*([1]+erf((log(x)-[2])/[3])), range widened by successive tries
    # Careful: this no longer appears in the output name, so two runs differing
    # only by this setting overwrite each other.
    useOld1oPFit    = 1,

    # Write the fits into DebugFit/. 0|1
    # At 1: nPE .root files PER systematic (200 x 14 = 2800 with --all).
    saveFits        = 1,

    # Fold the spectrum onto |eta| before the estimate. 0|1
    takeAbsEta      = 0,

    # Validation region 0.8 < Fpixel <= 0.9, unblinded. 0|1
    runVR           = 1,

    # Search region 0.9 < Fpixel <= 1.0, blinded above 300 GeV. 0|1
    runSR           = 0,
    # runVR and runSR write into the SAME file (_8fp9 / _9fp10 suffixes on the
    # histogram names). Set both to 1 to get them together: in two separate
    # passes the second overwrites the first (RECREATE).

    # Number of processes forked for the toys. <= number of available cores.
    nWorkers        = 25,
)

# When both cuts are active, record them in the destination directory name so
# that two selections cannot land in the same folder.
if (common["eopCut"] != "" and common["sigmaPtCut"] != ""): labelDir = "SigmaPtoverPt_" + common["sigmaPtCut"] + "_EoP_" + common["eopCut"] + "_" + labelDir

# One entry per systematic: only the keys that change are listed.
# The first one is the nominal, it is the one that runs without --all.
#
# The label is not cosmetic: it is what the macro appends to the output file
# name, and an empty or duplicated label makes runs overwrite each other.
# Up/Down pairs:
#   binEta/binIh/binMom: binning granularity of the templates
#   fitIh/fitMom       : +/-1 sigma on the fitted template, from the MINUIT
#                         covariance matrix (2 = up, 0 = down, 1 = nominal)
#   noFit              : pure bin-by-bin convolution, no tail fit at all
#   corrTemplate*      : Fpixel correlation correction applied to a template
config = [
    dict(label="nominal"),
    dict(label="binEtaUp",        rebinEta=2),
    dict(label="binEtaDown",      rebinEta=8),
    dict(label="binIhUp",         rebinIh=2),
    dict(label="binIhDown",       rebinIh=8),
    dict(label="binMomUp",        rebinMom=1),
    dict(label="binMomDown",      rebinMom=4),
    dict(label="fitIhUp",         fitIh=2),
    dict(label="fitIhDown",       fitIh=0),
    dict(label="fitMomUp",        fitMom=2),
    dict(label="fitMomDown",      fitMom=0),
    dict(label="noFit",           useFit=0),
    dict(label="corrTemplateIh",  corrTemplateIh=1),
    dict(label="corrTemplate1oP", corrTemplate1oP=1),
]

# Nominal binning and fit values: an entry in 'config' only overrides what it
# mentions, everything else falls back here.
#
# Merge order below is common -> nominal -> variation, so a key present in both
# 'common' and 'nominal' would be won by 'nominal'. Keep the two sets disjoint.
nominal = dict(rebinEta=4, rebinIh=4, rebinMom=2,
               fitIh=1, fitMom=1, useFit=1,
               corrTemplateIh=0, corrTemplate1oP=0)


# --- selection of the configs to run ---
# --only wins over --all. Unknown labels abort immediately rather than silently
# running nothing, and the list of valid labels is printed to help.
if opt.only:
    keep    = [s.strip() for s in opt.only.split(",")]
    unknown = [k for k in keep if k not in [c["label"] for c in config]]
    if unknown:
        print("Unknown label(s): " + ", ".join(unknown))
        print("Available: " + ", ".join(c["label"] for c in config))
        sys.exit(1)
    toRun = [c for c in config if c["label"] in keep]
elif opt.runAll:
    toRun = config
else:
    toRun = config[:1]

# Version extracted from the dataset: "..._V12p35" -> "12p35"
# Used to build the destination directory, so that outputs from two step1
# productions never mix. A dataset name without the "_V" pattern is a hard error.
def getVersion(dataset):
    m = re.search(r"_V([0-9A-Za-z]+?)$", os.path.basename(dataset))
    if not m:
        print("Can't extract the version from '{}' ('_V' pattern expected)".format(dataset))
        sys.exit(1)
    return m.group(1)

# Region covered by the run, deduced from runVR / runSR.
# Note this is an elif chain: with both flags at 1 the directory is labelled
# "8fp9" even though the file also contains the SR histograms.
if common["runVR"]:                     regionName = "8fp9"
elif common["runSR"]:                     regionName = "9fp10"
else:
    print("Neither runVR nor runSR -> nothing to do")
    sys.exit(1)


# Writes the flat settings dict in the "key = value" format read by Config.
# Keys are sorted so that two runs with the same settings produce byte-identical
# config files, which makes them easy to diff.
def write_config(path, settings):
    with open(path, "w") as f:
        for key in sorted(settings):
            f.write("{:<16} = {}\n".format(key, settings[key]))


# Labels whose run failed; reported at the end and used as the exit status.
failed = []

for dataset, sampleTypeName in datasetList:
    print("\nLaunch on dataset:    " + dataset + "\n")
    for conf in toRun:
        # Three-layer merge: session settings, then nominal values, then the
        # variation's own overrides.
        settings = dict(common)
        settings.update(nominal)
        settings.update(conf)
        settings["sample"] = dataset
        settings["sampleType"] = sampleTypeName

        label    = settings["label"]
        cfgPath  = "configFile_{}.txt".format(label)

        write_config(cfgPath, settings)
        print("---- {} ----".format(label))
        os.system("cat " + cfgPath)

        # The macro writes its output next to the input file, because the
        # "sample" key is a full path and the output name is built from it.
        outdir = os.path.dirname(dataset) or "."

        # t0 is both the timer origin and the cut-off used to identify the files
        # this run produced.
        t0 = time.time()
        rc = os.system("root -l -q -b 'BkgPrediction.C+(\"{}\")'".format(cfgPath))
        elapsed = time.time() - t0

        # Files written since the run started. Note that ROOT creates the output
        # file at RECREATE time, i.e. before the macro can fail, so a non-empty
        # 'fresh' list does not by itself prove the estimate completed.
        fresh = [p for p in glob.glob(os.path.join(outdir, "*.root"))
                 if os.path.getmtime(p) >= t0]

        print("  {}: {:.0f} s, {} file(s) created".format(label, elapsed, len(fresh)))
        if rc != 0 or not fresh:
            # The config file is kept so the run can be replayed identically.
            print("FAILED on '{}' (code {}): no files written in {}"
                  .format(label, os.WEXITSTATUS(rc), outdir))
            failed.append(label)
            continue

        os.remove(cfgPath)

        # Move files to outputDir
        # Destination layout: <outputDir>/<sampleType>_V<version>__<region>_<labelDir>/<etaRange>/
        # Everything that distinguishes one production from another is in the
        # path, so outputs never silently overwrite each other.
        version = getVersion(dataset)
        destdir = os.path.join(outputDir,
                               "{}_V{}__{}_{}".format(sampleTypeName, version,
                                                      regionName, labelDir),
                               etaRangeName)
        if not os.path.isdir(destdir):
            os.makedirs(destdir)

        for f in fresh:
            shutil.move(f, os.path.join(destdir, os.path.basename(f)))
        print("  -> {} files moved to {}".format(len(fresh), destdir))

# Cleanup of the ACLiC artefacts, once at the very end.
# Doing it here rather than per run means the macro is compiled once and reused
# by every subsequent launch.
for f in glob.glob("BkgPrediction_C*"):
    os.remove(f)

# Non-zero exit status so that a wrapper script can detect a partial sweep.
if failed:
    print("\n{} failed run(s): {}".format(len(failed), ", ".join(failed)))
    sys.exit(1)