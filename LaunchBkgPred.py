"""Lance l'estimation de bruit de fond.

    python LaunchBkgPred.py              -> nominale seulement
    python LaunchBkgPred.py --all        -> toutes les systematiques
    python LaunchBkgPred.py --only etaup,ihdown
"""

import os, sys, glob, time, re, shutil
from optparse import OptionParser

parser = OptionParser(usage="Usage: python %prog [options]")
parser.add_option("--all",  dest="runAll", action="store_true", default=False,
                  help="toutes les systematiques (defaut: nominale seulement)")
parser.add_option("--only", dest="only", default=None,
                  help="labels precis, separes par des virgules")
(opt, args) = parser.parse_args()

datasetList = [
    #("/safe/ui3_1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/output/JetMET2024_V12/JetMET2024_V12p35", "data2024"),
    ("/safe/ui3_1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/output/HistForBkg_MC_V3", "mc2024"),
]

outputDir      = "/safe/ui3_1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/macros"
etaRangeName   = "Eta2p4"
sampleTypeName = "mc2024"
labelDir       = "v2"

# Reglages communs a tous les lancements de cette session.
common = dict(
    # Nombre de pseudo-experiences (toys). Tout entier > 0. 200 = nominal.
    # Le temps de calcul est proportionnel : 20 pour un test rapide.
    nPE             = 200,

    # Appliquer le rebinning. 0|1. A 0, les histos du step1 sont pris tels quels
    # (rebinEta / rebinIh / rebinMom sont alors ignores).
    rebin           = 1,

    # Calibration dE/dx (K, C) utilisee pour M = sqrt((Ih-C)/K) * p.
    # data2017 | data2018 | data2024 | mc2017 | mc2018 | mc2024
    # Doit correspondre au dataset : un MC lu avec la calibration data
    # donne un spectre en masse faux, sans aucun avertissement.
    sampleType      = sampleTypeName,

    # Plage en eta, doit exister cote step1.
    # Eta1 | Eta1_2p4 | Eta2p4 | Eta1p2_2p2 | Eta1p2_2p4
    etaRange        = etaRangeName,

    # Coupure E/p du step1. "" (pas de coupure) | "0p1"
    eopCut          = "",

    # Coupure sigma(pT)/pT du step1. "" (pas de coupure) | "0p5"
    sigmaPtCut      = "",

    # Forme du fit du template Ih. 0|1
    #   0 = gaussienne ajustee a partir de 1.1 * max
    #   1 = ancien fit, depart fixe a Ih = 3 MeV/cm
    useOldIhFit     = 0,

    # Forme du fit du template 1/p. 0|1
    #   0 = 0.5*(exp(ax^2+bx) + exp(-ax^2-bx)) - 1, plage [0, 0.6*pic]
    #   1 = [0]*([1]+erf((log(x)-[2])/[3])), plage elargie par essais successifs
    # Attention : n'apparait plus dans le nom de sortie, deux runs qui ne
    # different que par ce reglage s'ecrasent mutuellement.
    useOld1oPFit    = 1,

    # Ecrire les fits dans DebugFit/. 0|1
    # A 1 : nPE fichiers .root PAR systematique (200 x 14 = 2800 sur --all).
    saveFits        = 1,

    # Replier le spectre sur |eta| avant l'estimation. 0|1
    takeAbsEta      = 0,

    # Region de validation 0.8 < Fpixel <= 0.9, non blindee. 0|1
    runVR           = 1,

    # Region de recherche 0.9 < Fpixel <= 1.0, blindee au-dessus de 300 GeV. 0|1
    runSR           = 0,
    # runVR et runSR ecrivent dans le MEME fichier (suffixes _8fp9 / _9fp10 sur
    # les histos). Mets les deux a 1 pour les avoir ensemble : en deux passes,
    # la seconde ecrase la premiere (RECREATE).

    # Nombre de process forkes pour les toys. <= nombre de coeurs disponibles.
    nWorkers        = 25,
)

# Une entree par systematique : seules les cles qui changent sont listees.
# La premiere est la nominale, c'est elle qui tourne sans --all.
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

# Valeurs nominales du binning et des fits : une entree de 'config' ne surcharge
# que ce qu'elle mentionne, le reste retombe ici.
nominal = dict(rebinEta=4, rebinIh=4, rebinMom=2,
               fitIh=1, fitMom=1, useFit=1,
               corrTemplateIh=0, corrTemplate1oP=0)


# --- selection des configs a lancer ---
if opt.only:
    keep    = [s.strip() for s in opt.only.split(",")]
    unknown = [k for k in keep if k not in [c["label"] for c in config]]
    if unknown:
        print("Label(s) inconnu(s) : " + ", ".join(unknown))
        print("Disponibles : " + ", ".join(c["label"] for c in config))
        sys.exit(1)
    toRun = [c for c in config if c["label"] in keep]
elif opt.runAll:
    toRun = config
else:
    toRun = config[:1]

# Version extraite du dataset : "..._V12p35" -> "12p35"
def getVersion(dataset):
    m = re.search(r"_V([0-9A-Za-z]+?)$", os.path.basename(dataset))
    if not m:
        print("Impossible d'extraire la version de '{}' (motif '_V' attendu)".format(dataset))
        sys.exit(1)
    return m.group(1)

# Region couverte par le run, deduite de runVR / runSR.
if common["runVR"]:                     regionName = "8fp9"
elif common["runSR"]:                     regionName = "9fp10"
else:
    print("Ni runVR ni runSR -> rien a lancer")
    sys.exit(1)


def write_config(path, settings):
    with open(path, "w") as f:
        for key in sorted(settings):
            f.write("{:<16} = {}\n".format(key, settings[key]))


failed = []

for dataset, sampleTypeName in datasetList:
    print("\nLaunch on dataset:    " + dataset + "\n")
    for conf in toRun:
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

        outdir = os.path.dirname(dataset) or "."

        t0 = time.time()
        rc = os.system("root -l -q -b 'BkgPrediction.C+(\"{}\")'".format(cfgPath))
        elapsed = time.time() - t0

        fresh = [p for p in glob.glob(os.path.join(outdir, "*.root"))
                 if os.path.getmtime(p) >= t0]

        print("  {} : {:.0f} s, {} fichier(s) produit(s)".format(label, elapsed, len(fresh)))
        if rc != 0 or not fresh:
            # On garde le fichier de config pour pouvoir rejouer a l'identique.
            print("ECHEC sur '{}' (code {}) : aucun fichier ecrit dans {}"
                  .format(label, os.WEXITSTATUS(rc), outdir))
            failed.append(label)
            continue

        os.remove(cfgPath)

        # Moove files in outputDir
        version = getVersion(dataset)
        destdir = os.path.join(outputDir,
                               "{}_V{}__{}_{}".format(sampleTypeName, version,
                                                      regionName, labelDir),
                               etaRangeName)
        if not os.path.isdir(destdir):
            os.makedirs(destdir)

        for f in fresh:
            shutil.move(f, os.path.join(destdir, os.path.basename(f)))
        print("  -> {} fichier(s) deplace(s) dans {}".format(len(fresh), destdir))

# Nettoyage des artefacts ACLiC, une seule fois a la fin.
for f in glob.glob("BkgPrediction_C*"):
    os.remove(f)

if failed:
    print("\n{} run(s) en echec : {}".format(len(failed), ", ".join(failed)))
    sys.exit(1)