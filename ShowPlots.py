#!/usr/bin/env python3
"""Lance le plotting du spectre en masse (MyMacroMass.py).

    python ShowPlots.py
    python ShowPlots.py --etas Eta1,Eta1_2p4,Eta2p4
    python ShowPlots.py --label binEtaUp
    python ShowPlots.py --syst
"""

import os, re, sys
from optparse import OptionParser

# Racine des repertoires de travail.
BASE = "/safe/ui3_1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/macros"

# --- a regler a la main ---------------------------------------------------
# Meme chemin que dans le launcher de step2 : il donne le prefixe des .root
# produits et la version (ce qui suit "_V").
#DATASET    = "/safe/ui3_1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/output/JetMET2024_V12/JetMET2024_V12p35"
DATASET    = "/safe/ui3_1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/output/HistForBkg_MC_V3"
SAMPLETYPE = "mc2024"   # data2017|data2018|data2024|mc2017|mc2018|mc2024
SUFFIX     = "v2"         # suffixe libre du repertoire de travail
CUTS       = ""           # selection du step1 : "" | "_SigmaPtoverPt_0p5_EoP_0p1" | ...
VSIGNAL    = "19p12"   # version des echantillons gluino
REGION     = "8fp9"
YEAR       = "2024"
ERA        = ""        # "" | "F" | "G"
ISTTBAR    = False     # MC : ne tracer que le ttbar dileptonique
# --------------------------------------------------------------------------

parser = OptionParser(usage="Usage: python %prog [options]")
parser.add_option("--etas",    dest="etas",    default="Eta2p4",
                  help="plages en eta separees par des virgules")
parser.add_option("--label",   dest="label",   default="nominal",
                  help="systematique a tracer (nominal, binEtaUp, ...)")
parser.add_option("--suffix",  dest="suffix",  default=None,
                  help="surcharge SUFFIX")
parser.add_option("--cuts",    dest="cuts",    default=None,
                  help="surcharge CUTS")
parser.add_option("--syst",    dest="syst",    action="store_true", default=False,
                  help="inclure les systematiques (defaut: nominal seul)")
parser.add_option("--systdir", dest="systdir", default="SystCombined",
                  help="sous-repertoire des fichiers de systematiques")
parser.add_option("--ofile",   dest="ofile",   default="mass_plot_Cnoabs",
                  help="prefixe des fichiers de sortie")
parser.add_option("--dry-run", dest="dryrun",  action="store_true", default=False,
                  help="affiche les commandes sans les executer")
parser.add_option("--vsignal", dest="vsignal", default=VSIGNAL)
parser.add_option("--year",    dest="year",    default=YEAR)
parser.add_option("--era",     dest="era",     default=ERA,
                  help="'' (toute l'annee) | F | G")
parser.add_option("--ttbar",   dest="ttbar",   action="store_true", default=ISTTBAR,
                  help="MC : ne tracer que le ttbar dileptonique")
(opt, args) = parser.parse_args()

suffix = opt.suffix if opt.suffix is not None else SUFFIX
cuts   = opt.cuts   if opt.cuts   is not None else CUTS
if cuts and not cuts.startswith("_"):
    cuts = "_" + cuts

# stem = prefixe des .root ecrits par step2 (= basename du dataset)
stem = os.path.basename(DATASET)

m = re.search(r"_V([0-9A-Za-z]+?)(?:_|$)", stem)
if not m:
    print("Impossible d'extraire la version de '{}' (motif '_V' attendu)".format(stem))
    sys.exit(1)
version = m.group(1)

# Le type d'echantillon suffit a savoir si on trace du MC.
isMC  = SAMPLETYPE.startswith("mc")
indir = "{}/{}_V{}__{}_{}".format(BASE, SAMPLETYPE, version, REGION, suffix)

print("Repertoire : {}".format(indir))
print("Prefixe    : {}   version {}   {}".format(
      stem, version, "MC" if isMC else "data"))


etas    = [e.strip() for e in opt.etas.split(",")    if e.strip()]

missing, launched = [], 0

for eta in etas:
    ifile = "{}/{}/{}_{}{}_{}.root".format(indir, eta, stem, eta, cuts, opt.label)
    if not os.path.isfile(ifile):
        missing.append(ifile)
        continue

    region = REGION                       # une seule region par repertoire de travail
    odir   = "{}/{}/Plots_{}".format(indir, eta, region)

    systfile = ""
    if opt.syst:
        systfile = "{}/{}/sysTotBinned_{}_{}.root".format(
                    indir, opt.systdir, eta, region)
        if not os.path.isfile(systfile):
            print("Systematiques absentes : {} -> region ignoree".format(systfile))
            continue

    command = ("python3 PlottingMacro.py"
                " --ifile {} --cuts '{}' --ofile {} --region {}"
                " --odir {} --nom {} --eta {} --isMC {}"
                " --vsignal {} --isTTbar {} --year {}").format(
                    ifile, cuts, opt.ofile, region, odir,
                    not opt.syst, eta, isMC,
                    opt.vsignal, opt.ttbar, opt.year)
    if opt.era:      command += " --era {}".format(opt.era)
    if systfile:     command += " --systfile {}".format(systfile)

    print("\n       Running:\n{}\n".format(command))
    if not opt.dryrun:
        rc = os.system(command)
        if rc != 0:
            print("ECHEC sur eta={} region={} (code {})".format(
                    eta, region, os.WEXITSTATUS(rc)))
            sys.exit(1)
    launched += 1

if missing:
    print("\n{} fichier(s) introuvable(s) :".format(len(missing)))
    for f in missing:
        print("  " + f)

if launched == 0:
    print("\nRien n'a ete trace.")
    sys.exit(1)