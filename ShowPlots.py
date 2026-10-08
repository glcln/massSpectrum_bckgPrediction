#!/usr/bin/python
"""Runs the mass spectrum plotting (MyMacroMass.py).

    python3 ShowPlots.py
    python3 ShowPlots.py --etas Eta1,Eta1_2p4,Eta2p4
    python3 ShowPlots.py --label binEtaUp
    python3 ShowPlots.py --syst
"""

# =============================================================================
#  Step 3c - driver of the final mass-spectrum plots.
# -----------------------------------------------------------------------------
#  This script does no physics at all: it only resolves paths and shells out to
#  the plotting macro, once per eta range.
#
#  For each eta range it:
#    1. rebuilds the path of the step2 output file for the requested systematic
#       label (the same convention the launcher used when filing it);
#    2. optionally locates the combined-systematics file written by systBckg.py;
#    3. builds the command line and runs the plotting macro.
#
#  The path convention is duplicated by hand here, in systBckg.py and in the
#  launcher. Nothing enforces that the three agree, so a mismatch shows up as a
#  "file not found" listing at the end rather than as a wrong plot.
# =============================================================================

import os, re, sys
from optparse import OptionParser

# Root of the working directories.
BASE = "/safe/ui3_1/cms/gcoulon/CMSSW_15_0_13_patch1/src/outputHist"

# --- to be set by hand ----------------------------------------------------
# Same path as in the step2 launcher: it gives the prefix of the .root files
# produced and the version (whatever follows "_V").
DATASET    = "/safe/ui3_1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/output/JetMET2024_V12/JetMET2024_V12p35"
#DATASET    = "/safe/ui3_1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/output/HistForBkg_MC_V3"
#DATASET     = "/safe/ui3_1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/output/TTbar2024_V15/TTbar2024_V15p10"
SAMPLETYPE = "data2024"   # data2017|data2018|data2024|mc2017|mc2018|mc2024|ttbar2024
SUFFIX     = "v2"         # free-form suffix of the working directory
CUTS       = ""#_SigmaPtoverPt_0p5_EoP_0p1"           # step1 selection: "" | "_SigmaPtoverPt_0p5_EoP_0p1" | ...
VSIGNAL    = "19p12"   # version of the gluino samples
REGION     = "9fp10"
YEAR       = "2024"
ERA        = ""        # "" | "F" | "G"
ISTTBAR    = False     # MC: plot the dileptonic ttbar only
# --------------------------------------------------------------------------

# Every constant above can be overridden from the command line, so a one-off
# plot never requires editing this file.
parser = OptionParser(usage="Usage: python %prog [options]")
parser.add_option("--etas",    dest="etas",    default="Eta2p4",
                  help="eta ranges, separated by commas")
parser.add_option("--label",   dest="label",   default="nominal",
                  help="systematic to plot (nominal, binEtaUp, ...)")
parser.add_option("--suffix",  dest="suffix",  default=None,
                  help="override SUFFIX")
parser.add_option("--cuts",    dest="cuts",    default=None,
                  help="override CUTS")
parser.add_option("--syst",    dest="syst",    action="store_true", default=False,
                  help="include the systematics (default: nominal only)")
parser.add_option("--systdir", dest="systdir", default="SystCombined",
                  help="sub-directory holding the systematics files")
parser.add_option("--ofile",   dest="ofile",   default="mass_plot_Cnoabs",
                  help="prefix of the output files")
parser.add_option("--dry-run", dest="dryrun",  action="store_true", default=False,
                  help="print the commands without running them")
parser.add_option("--vsignal", dest="vsignal", default=VSIGNAL)
parser.add_option("--year",    dest="year",    default=YEAR)
parser.add_option("--era",     dest="era",     default=ERA,
                  help="'' (whole year) | F | G")
parser.add_option("--ttbar",   dest="ttbar",   action="store_true", default=ISTTBAR,
                  help="MC: plot the dileptonic ttbar only")
(opt, args) = parser.parse_args()

# None (flag absent) means "use the module constant"; an explicit "" means "no
# cut", which is why the defaults are None rather than the constants themselves.
suffix = opt.suffix if opt.suffix is not None else SUFFIX
cuts   = opt.cuts   if opt.cuts   is not None else CUTS
# The cut string is a path fragment and must start with '_'; adding it here
# makes both "--cuts _EoP_0p1" and "--cuts EoP_0p1" work.
if cuts and not cuts.startswith("_"):
    cuts = "_" + cuts

# stem = prefix of the .root files written by step2 (= basename of the dataset)
stem = os.path.basename(DATASET)

# "..._V12p35" -> "12p35". The version goes into the working-directory name so
# that outputs of two step1 productions never mix. Same regex as systBckg.py:
# the version may end the name or be followed by an underscore.
m = re.search(r"_V([0-9A-Za-z]+?)(?:_|$)", stem)
if not m:
    print("Can't extract the version from '{}' ('_V' pattern expected)".format(stem))
    sys.exit(1)
version = m.group(1)

# The sample type is enough to know whether we are plotting MC.
isMC  = SAMPLETYPE.startswith("mc") or SAMPLETYPE.startswith("ttbar")
# Working directory, matching what the launcher built when it filed the outputs:
# <sampleType>_V<version>__<region><cuts>_<suffix>
indir = "{}/{}_V{}__{}{}_{}".format(BASE, SAMPLETYPE, version, REGION, cuts, suffix)

print("Directory : {}".format(indir))
print("Prefix    : {}   version {}   {}".format(
      stem, version, "MC" if isMC else "data"))


etas    = [e.strip() for e in opt.etas.split(",")    if e.strip()]

# missing  : files that could not be found, listed once at the end
# launched : number of plotting jobs actually started
missing, launched = [], 0

for eta in etas:
    # Mirrors the output name built by BkgPrediction.C:
    #     <dataset stem> + "_" + <eta range> + <cuts> + "_" + <label>
    # inside the per-eta sub-directory created by the launcher.
    ifile = "{}/{}/{}_{}{}_{}.root".format(indir, eta, stem, eta, cuts, opt.label)
    # A missing file is collected rather than fatal, so one absent eta range does
    # not stop the others.
    if not os.path.isfile(ifile):
        missing.append(ifile)
        continue

    region = REGION                       # a single region per working directory
    odir   = "{}/{}/Plots_{}".format(indir, eta, region)

    # Combined systematics written by systBckg.py. Same name and location it
    # produced them under; if it is missing we skip rather than silently fall
    # back to the nominal-only plot, which would look identical but mean
    # something different.
    systfile = ""
    if opt.syst:
        systfile = "{}/{}/sysTotBinned_{}_{}.root".format(
                    indir, opt.systdir, eta, region)
        if not os.path.isfile(systfile):
            print("Missing systematics: {} -> region skipped".format(systfile))
            continue

    # --nom is the plotting macro's "nominal only" switch, so it is the negation
    # of --syst here. --systfile is only appended when a file was actually found.
    command = ("python3 PlottingMacro.py"
                " --ifile {} --cuts '{}' --ofile {} --region {}"
                " --odir {} --nom {} --eta {} --isMC {}"
                " --vsignal {} --isTTbar {} --year {}").format(
                    ifile, cuts, opt.ofile, region, odir,
                    not opt.syst, eta, isMC,
                    opt.vsignal, opt.ttbar, opt.year)
    if opt.era:      command += " --era {}".format(opt.era)
    if systfile:     command += " --systfile {}".format(systfile)

    # The command is always printed, so a failing run can be replayed by hand.
    print("\n       Running:\n{}\n".format(command))
    if not opt.dryrun:
        rc = os.system(command)
        # A plotting failure aborts the whole loop: it usually means a histogram
        # name changed, which would affect every eta range identically.
        if rc != 0:
            print("FAILED on eta={} region={} (code {})".format(
                    eta, region, os.WEXITSTATUS(rc)))
            sys.exit(1)
    launched += 1

# Reported at the end rather than inline, so the list is easy to read.
if missing:
    print("\n{} file(s) not found:".format(len(missing)))
    for f in missing:
        print("  " + f)

# Non-zero exit when nothing at all was plotted, so a wrapper script can react.
if launched == 0:
    print("\nNothing was plotted.")
    sys.exit(1)