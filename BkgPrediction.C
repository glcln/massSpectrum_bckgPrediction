#include <fstream>
#include <iostream>
#include <map>
#include <string>

#include "TFile.h"
#include "TError.h"
#include "TSystem.h"

#include "CommonFunctions.h"

// ---------- mini parseur "cle = valeur" ----------
// Format non positionnel : ajouter un parametre ne peut plus decaler les autres.
class Config {
  public:
    bool load(const std::string& path) {
        std::ifstream in(path);
        if (!in.is_open()) { std::cerr << "Config: impossible d'ouvrir " << path << std::endl; return false; }
        std::string line;
        while (std::getline(in, line)) {
            const std::size_t hash = line.find('#');
            if (hash != std::string::npos) line.erase(hash);   // commentaire de fin de ligne
            const std::size_t eq = line.find('=');
            if (eq == std::string::npos) continue;             // ligne vide / sans '='
            kv_[trim(line.substr(0, eq))] = trim(line.substr(eq + 1));
        }
        return true;
    }

    std::string str(const std::string& k, const std::string& d = "") const {
        auto it = kv_.find(k);
        if (it == kv_.end()) { warn(k, d); return d; }
        return it->second;
    }
    int getInt(const std::string& k, int d) const {
        auto it = kv_.find(k);
        if (it == kv_.end() || it->second.empty()) { warn(k, std::to_string(d)); return d; }
        return std::stoi(it->second);
    }
    bool getBool(const std::string& k, bool d) const {
        auto it = kv_.find(k);
        if (it == kv_.end() || it->second.empty()) { warn(k, d ? "1" : "0"); return d; }
        const std::string& v = it->second;
        return (v == "1" || v == "true" || v == "True" || v == "yes");
    }

    void dump() const {
        std::cout << "---- configuration ----" << std::endl;
        for (const auto& p : kv_) std::cout << "  " << p.first << " = " << p.second << std::endl;
        std::cout << "-----------------------" << std::endl;
    }

  private:
    static std::string trim(std::string s) {
        const char* ws = " \t\r\n";
        const std::size_t a = s.find_first_not_of(ws);
        if (a == std::string::npos) return "";
        return s.substr(a, s.find_last_not_of(ws) - a + 1);
    }
    static void warn(const std::string& k, const std::string& d) {
        std::cerr << "Config: cle '" << k << "' absente -> defaut '" << d << "'" << std::endl;
    }
    std::map<std::string, std::string> kv_;
};


void BkgPrediction(const char* configPath = "configFile_readHisto_toLaunch.txt") {

    gErrorIgnoreLevel = kFatal;   // dans le corps : une affectation a portee globale ne compile pas

    Config cfg;
    if (!cfg.load(configPath)) return;
    cfg.dump();

    // ---- entrees, binning, systematiques (pilotes par la boucle python) ----
    const std::string filename        = cfg.str("sample");
    const int         nPE             = cfg.getInt("nPE", 200);
    const bool        bool_rebin      = cfg.getBool("rebin", true);
    const int         rebineta        = cfg.getInt("rebinEta", 4);
    const int         rebinih         = cfg.getInt("rebinIh",  4);
    const int         rebinp          = cfg.getInt("rebinMom", 2);
    const int         fitIh           = cfg.getInt("fitIh",  1);
    const int         fitP            = cfg.getInt("fitMom", 1);
    const bool        useFit          = cfg.getBool("useFit", true);
    const bool        corrTemplateIh  = cfg.getBool("corrTemplateIh",  false);
    const bool        corrTemplate1oP = cfg.getBool("corrTemplate1oP", false);

    // ---- ce qui etait commente/decommente a la main ----
    const std::string st_sample    = cfg.str("sampleType", "data2024");
    const std::string etaRange     = cfg.str("etaRange",   "Eta1");
    const std::string eopCut       = cfg.str("eopCut",     "");
    const std::string sigPtCut     = cfg.str("sigmaPtCut", "");
    const std::string ihLabel      = cfg.str("ihLabel",    "");
    const bool        useOldIhFit  = cfg.getBool("useOldIhFit",  false);
    const bool        useOld1oPFit = cfg.getBool("useOld1oPFit", true);
    const bool        saveFits     = cfg.getBool("saveFits",     false);
    const bool        TakeAbsEta   = cfg.getBool("takeAbsEta",   false);
    const bool        runVR        = cfg.getBool("runVR", false);
    const bool        runSR        = cfg.getBool("runSR", true);
    const unsigned    nWorkers     = static_cast<unsigned>(cfg.getInt("nWorkers", 25));

    if (filename.empty()) { std::cerr << "Config: cle 'sample' vide -> abandon" << std::endl; return; }
    if (!runVR && !runSR) { std::cerr << "Config: ni runVR ni runSR -> rien a faire" << std::endl; return; }

    const DeDxCalib calib = GetDeDxCalib(st_sample);
    std::cout << "dE/dx calibration (" << st_sample << "): K=" << calib.K << " C=" << calib.C << std::endl;

    // Ext     : suffixe des histos dans le fichier d'entree (ordre du step1 : coupures puis eta)
    // etaName : etiquette interne (ordre inverse : eta puis coupures)
    std::string Ext = "_METanalysis_TestPUppiMETCut";
    if (!sigPtCut.empty()) Ext += "_SigmaPtoverPt_" + sigPtCut;
    if (!eopCut.empty())   Ext += "_EoP_" + eopCut;
    Ext += "_" + etaRange;
    if (!ihLabel.empty())  Ext += "_" + ihLabel;

    std::string etaName = "_" + etaRange;
    if (!sigPtCut.empty()) etaName += "_SigmaPtoverPt_" + sigPtCut;
    if (!eopCut.empty())   etaName += "_EoP_" + eopCut;
    if (!ihLabel.empty())  etaName += "_" + ihLabel;

    // Nom de sortie : <dataset>_<etaName>_<label>
    // Le label vient du launcher, il identifie la systematique a lui seul.
    const std::string label = cfg.str("label", "");
    if (label.empty()) {
        std::cerr << "Config: cle 'label' vide -> abandon (risque d'ecraser un autre run)" << std::endl;
        return;
    }
    const std::string outfilename_ = filename + etaName + "_" + label;

    const std::string DataSetName = filename.substr(filename.find_last_of('/') + 1);
    std::cout << "Input file:      " << DataSetName << std::endl;
    std::cout << "Output file:     " << outfilename_ << std::endl;
    std::cout << "Ext:             " << Ext     << std::endl;
    std::cout << "etaName:         " << etaName << std::endl;
    std::cout << "abs(eta):        " << TakeAbsEta << std::endl;

    if (saveFits && gSystem->AccessPathName("DebugFit")) gSystem->mkdir("DebugFit", kTRUE);

    TFile* ifile = TFile::Open((filename + ".root").c_str(), "READ");
    if (!ifile || ifile->IsZombie()) {
        std::cerr << "Impossible d'ouvrir " << filename << ".root" << std::endl;
        return;
    }
    TFile* ofile = TFile::Open((outfilename_ + ".root").c_str(), "RECREATE");
    if (!ofile || ofile->IsZombie()) {
        std::cerr << "Impossible de creer " << outfilename_ << ".root" << std::endl;
        ifile->Close();
        return;
    }

    // ------------------------------------------------------------------
    //                              If Fpixel
    //
    //                    pT
    //                      |           |          |
    //                      |     C     |     D    | blind
    //                      |           |        VR|
    //                   70 |-----------|----------|------
    //                      |           |          |
    //                      |     A     |     B    |
    //                   55 |___________|__________|______
    //                     0.3         0.8        0.9      Fpixel
    //
    // ------------------------------------------------------------------

    const bool ifIhpSAME = true;   // TRUE: templates Ih et p pris en C, B garde la normalisation

    BckgOptions opt;
    opt.nPE             = nPE;
    opt.useFit          = useFit;
    opt.useOldIhFit     = useOldIhFit;
    opt.useOld1oPFit    = useOld1oPFit;
    opt.corrTemplateIh  = corrTemplateIh;
    opt.corrTemplate1oP = corrTemplate1oP;
    opt.etaName         = etaName;
    opt.saveFits        = saveFits;
    opt.fitIh           = fitIh;
    opt.fitP            = fitP;
    opt.nWorkers        = nWorkers;


    bool done = false;

    if (runVR) {
        std::cout << "\n    Loading validation region..." << std::endl;
        Region ra_3fp8, rb_8fp9, rc_3fp8, rd_8fp9, rbc_8fp9;
        bool ok = true;
        ok &= loadHistograms(ra_3fp8,  ifile, "regionA_3fp8"  + Ext, bool_rebin, rebineta, rebinp, rebinih, TakeAbsEta);
        ok &= loadHistograms(rb_8fp9,  ifile, "regionB_8fp9"  + Ext, bool_rebin, rebineta, rebinp, rebinih, TakeAbsEta);
        ok &= loadHistograms(rc_3fp8,  ifile, "regionC_3fp8"  + Ext, bool_rebin, rebineta, rebinp, rebinih, TakeAbsEta);
        ok &= loadHistograms(rd_8fp9,  ifile, "regionD_8fp9"  + Ext, bool_rebin, rebineta, rebinp, rebinih, TakeAbsEta);
        ok &= loadHistograms(rbc_8fp9, ifile, "regionD_8fp9"  + Ext, bool_rebin, rebineta, rebinp, rebinih, TakeAbsEta);

        if (!ok) std::cerr << "VR: chargement incomplet -> region ignoree" << std::endl;
        else {
            ofile->cd();
            opt.blind = false;
            std::cout << "    Background estimation, VR 8fp9..." << std::endl;
            done = bckgEstimate(DataSetName, calib, rc_3fp8, rc_3fp8, rbc_8fp9, ra_3fp8, rd_8fp9,
                         ifIhpSAME, rb_8fp9, "8fp9", opt);
        }
    }

    if (runSR) {
        std::cout << "\n    Loading search region..." << std::endl;
        Region ra_3fp9, rb_9fp10, rc_3fp9, rd_9fp10, rbc_9fp10;
        bool ok = true;
        ok &= loadHistograms(ra_3fp9,   ifile, "regionA_3fp9"  + Ext, bool_rebin, rebineta, rebinp, rebinih, TakeAbsEta);
        ok &= loadHistograms(rb_9fp10,  ifile, "regionB_9fp10" + Ext, bool_rebin, rebineta, rebinp, rebinih, TakeAbsEta);
        ok &= loadHistograms(rc_3fp9,   ifile, "regionC_3fp9"  + Ext, bool_rebin, rebineta, rebinp, rebinih, TakeAbsEta);
        ok &= loadHistograms(rd_9fp10,  ifile, "regionD_9fp10" + Ext, bool_rebin, rebineta, rebinp, rebinih, TakeAbsEta);
        ok &= loadHistograms(rbc_9fp10, ifile, "regionD_9fp10" + Ext, bool_rebin, rebineta, rebinp, rebinih, TakeAbsEta);

        if (!ok) std::cerr << "SR: chargement incomplet -> region ignoree" << std::endl;
        else {
            ofile->cd();
            opt.blind = true;
            std::cout << "    Background estimation, SR 9fp10..." << std::endl;
            done = bckgEstimate(DataSetName, calib, rc_3fp9, rc_3fp9, rbc_9fp10, ra_3fp9, rd_9fp10,
                         ifIhpSAME, rb_9fp10, "9fp10", opt);
        }
    }

    ofile->Close();
    ifile->Close();
    delete ofile;
    delete ifile;
    if (!done) { std::cerr << "Echec de l'estimation -> fichier de sortie vide" << std::endl; return; }
    std::cout << "\nDone: " << outfilename_ << ".root" << std::endl;

    return;
}
