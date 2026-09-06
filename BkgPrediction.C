// =============================================================================
//  BkgPrediction.C  -  step 2 of the HSCP background estimate
// -----------------------------------------------------------------------------
//  Entry point of the analysis chain on the step1 output. It:
//    1. reads a key = value configuration file (one file per systematic
//       variation, written by the Python launcher);
//    2. resolves the dE/dx calibration and builds the input/output names;
//    3. loads the ABCD regions from the step1 ROOT file;
//    4. calls bckgEstimate() for the validation region and/or the search region;
//    5. writes everything into a single output file.
//
//  Invoked as:  root -l -q -b 'BkgPrediction.C+("configFile_<label>.txt")'
//
//  The macro is deliberately "one config file = one run = one output file": the
//  systematic variations are driven entirely from outside, which keeps this file
//  free of the manual comment/uncomment editing it used to require.
//
//  On success the last line printed is "Done: <output>.root". That sentinel is
//  what the launcher greps for: checking the existence or the mtime of the output
//  file is not reliable, because TFile::Open(..., "RECREATE") creates the file on
//  disk before the macro has had any chance to fail.
// =============================================================================

#include <fstream>
#include <iostream>
#include <map>
#include <string>

#include "TFile.h"
#include "TError.h"
#include "TSystem.h"

#include "CommonFunctions.h"

// ---------- minimal "key = value" parser ----------
// Non-positional format: adding a parameter can no longer shift the others.
//
// This replaces the previous positional parsing, where the Nth line of the file
// had to match the Nth variable in the macro. Any mismatch there was silent and
// produced a plausible-looking but wrong result.
//
// Accepted syntax:
//    key = value        (whitespace around key and value is stripped)
//    # comment          (everything after '#' is dropped, end-of-line comments
//                        included)
//    <blank line>       (ignored, as is any line without '=')
//
// A missing key is not fatal: the default is used and a warning is printed, so a
// config file written by an older launcher still runs.
class Config {
  public:
    // Reads the whole file into the key/value map. Returns false only if the file
    // cannot be opened; malformed lines are skipped silently.
    bool load(const std::string& path) {
        std::ifstream in(path);
        if (!in.is_open()) { std::cerr << "Config: can't open " << path << std::endl; return false; }
        std::string line;
        while (std::getline(in, line)) {
            const std::size_t hash = line.find('#');
            if (hash != std::string::npos) line.erase(hash);   // end-of-line comment
            const std::size_t eq = line.find('=');
            if (eq == std::string::npos) continue;             // blank line / no '='
            kv_[trim(line.substr(0, eq))] = trim(line.substr(eq + 1));
        }
        return true;
    }

    // String accessor. Note that, unlike getInt/getBool, an explicitly empty
    // value ("key =") is returned as an empty string rather than replaced by the
    // default: that is how the optional cut labels below are switched off.
    std::string str(const std::string& k, const std::string& d = "") const {
        auto it = kv_.find(k);
        if (it == kv_.end()) { warn(k, d); return d; }
        return it->second;
    }
    // Integer accessor. Missing *or empty* falls back to the default.
    int getInt(const std::string& k, int d) const {
        auto it = kv_.find(k);
        if (it == kv_.end() || it->second.empty()) { warn(k, std::to_string(d)); return d; }
        return std::stoi(it->second);
    }
    // Boolean accessor. Anything that is not one of the accepted "true" spellings
    // counts as false, so a typo turns an option off rather than on.
    bool getBool(const std::string& k, bool d) const {
        auto it = kv_.find(k);
        if (it == kv_.end() || it->second.empty()) { warn(k, d ? "1": "0"); return d; }
        const std::string& v = it->second;
        return (v == "1" || v == "true" || v == "True" || v == "yes");
    }

    // Echoes the parsed configuration into the log. Since the launcher deletes
    // the config file after a successful run, this print is the only record of
    // what a given output file was produced with.
    void dump() const {
        std::cout << "---- configuration ----" << std::endl;
        for (const auto& p: kv_) std::cout << "  " << p.first << " = " << p.second << std::endl;
        std::cout << "-----------------------" << std::endl;
    }

  private:
    // Strips leading and trailing whitespace, including '\r' so that a file
    // written on Windows does not produce keys with a trailing carriage return.
    static std::string trim(std::string s) {
        const char* ws = " \t\r\n";
        const std::size_t a = s.find_first_not_of(ws);
        if (a == std::string::npos) return "";
        return s.substr(a, s.find_last_not_of(ws) - a + 1);
    }
    static void warn(const std::string& k, const std::string& d) {
        std::cerr << "Config: key '" << k << "' missing -> default '" << d << "'" << std::endl;
    }
    std::map<std::string, std::string> kv_;
};


// =============================================================================
// BkgPrediction - main entry point.
//
// configPath: path to the key = value file describing this single run. The
// default value is only a convenience for interactive use; the launcher always
// passes an explicit path.
// =============================================================================
void BkgPrediction(const char* configPath = "configFile_readHisto_toLaunch.txt") {

    // Silences everything below "fatal". The per-eta-bin fits legitimately fail
    // on sparse slices and are already reported by fillPredMass; ROOT's own
    // warnings would drown that in noise.
    gErrorIgnoreLevel = kFatal;   // inside the body: a global-scope assignment does not compile

    Config cfg;
    if (!cfg.load(configPath)) return;
    cfg.dump();

    // ---- inputs, binning, systematics (driven by the python loop) ----
    // These are exactly the knobs the launcher varies to build the systematic
    // envelope: binning granularities, fit variations, template corrections.
    const std::string filename        = cfg.str("sample");
    const int         nPE             = cfg.getInt("nPE", 200);
    const bool        bool_rebin      = cfg.getBool("rebin", true);
    const int         rebineta        = cfg.getInt("rebinEta", 4);
    const int         rebinih         = cfg.getInt("rebinIh",  4);
    const int         rebinp          = cfg.getInt("rebinMom", 2);
    const int         fitIh           = cfg.getInt("fitIh",  1);   // 1 = nominal, 2 = +1s, 0 = -1s
    const int         fitP            = cfg.getInt("fitMom", 1);
    const bool        useFit          = cfg.getBool("useFit", true);
    const bool        corrTemplateIh  = cfg.getBool("corrTemplateIh",  false);
    const bool        corrTemplate1oP = cfg.getBool("corrTemplate1oP", false);

    // ---- what used to be commented/uncommented by hand ----
    // Everything below was previously edited directly in the source before each
    // run, which made a run impossible to reproduce from its output alone.
    const std::string st_sample    = cfg.str("sampleType", "data2024");
    const std::string etaRange     = cfg.str("etaRange",   "Eta1");
    const std::string eopCut       = cfg.str("eopCut",     "");
    const std::string sigPtCut     = cfg.str("sigmaPtCut", "");
    const bool        useOldIhFit  = cfg.getBool("useOldIhFit",  false);
    const bool        useOld1oPFit = cfg.getBool("useOld1oPFit", true);
    const bool        saveFits     = cfg.getBool("saveFits",     false);
    const bool        TakeAbsEta   = cfg.getBool("takeAbsEta",   false);
    const bool        runVR        = cfg.getBool("runVR", false);
    const bool        runSR        = cfg.getBool("runSR", true);
    const unsigned    nWorkers     = static_cast<unsigned>(cfg.getInt("nWorkers", 25));

    // Fail fast on a configuration that cannot produce anything useful, rather
    // than crashing later inside the ROOT I/O.
    if (filename.empty()) { std::cerr << "Config: empty key 'sample' -> aborting" << std::endl; return; }
    if (!runVR && !runSR) { std::cerr << "Config: neither runVR nor runSR -> nothing to do" << std::endl; return; }

    // dE/dx calibration selected from the sample type carried by the config, so
    // that MC is never read with the data constants (which would shift the whole
    // mass spectrum without any warning).
    const DeDxCalib calib = GetDeDxCalib(st_sample);
    std::cout << "dE/dx calibration (" << st_sample << "): K=" << calib.K << " C=" << calib.C << std::endl;

    // Ext    : suffix of the histograms in the input file (step1 ordering: cuts then eta)
    // etaName: internal label (reverse ordering: eta then cuts)
    //
    // The two orderings are genuinely different and must not be merged: Ext has
    // to reproduce byte for byte the names written by step1, while etaName is
    // what the correction functions (corrIh / corr1oP) and the output file name
    // are built from. Any optional cut that is left empty simply drops out of
    // both strings.
    //
    // The leading "_METanalysis_TestPUppiMETCut" is the step1 selection tag; it
    // must be updated whenever step1 changes its naming.
    std::string Ext = "_METanalysis_TestPUppiMETCut";
    if (!sigPtCut.empty()) Ext += "_SigmaPtoverPt_" + sigPtCut;
    if (!eopCut.empty())   Ext += "_EoP_" + eopCut;
    Ext += "_" + etaRange;

    std::string etaName = "_" + etaRange;
    if (!sigPtCut.empty()) etaName += "_SigmaPtoverPt_" + sigPtCut;
    if (!eopCut.empty())   etaName += "_EoP_" + eopCut;

    // Output name: <dataset>_<etaName>_<label>
    // The label comes from the launcher, it identifies the systematic on its own.
    //
    // This is why an empty label is fatal: two systematic variations sharing a
    // name would silently overwrite each other (the output is opened RECREATE),
    // and the loss would only show up much later as a missing variation.
    const std::string label = cfg.str("label", "");
    if (label.empty()) {
        std::cerr << "Config: empty key 'label' -> aborting (could overwrite files)" << std::endl;
        return;
    }
    // Note that `filename` is a full path, so `outfilename_` is one too: the
    // output lands next to the input, not in the current directory.
    const std::string outfilename_ = filename + etaName + "_" + label;

    // Basename only, used to name the debug-fit files (a full path there would
    // produce an unusable file name).
    const std::string DataSetName = filename.substr(filename.find_last_of('/') + 1);
    std::cout << "Input file:      " << DataSetName << std::endl;
    std::cout << "Output file:     " << outfilename_ << std::endl;
    std::cout << "Ext:             " << Ext     << std::endl;
    std::cout << "etaName:         " << etaName << std::endl;
    std::cout << "abs(eta):        " << TakeAbsEta << std::endl;

    // gSystem->AccessPathName returns TRUE when the path does NOT exist, hence
    // the apparently inverted test: create DebugFit/ only when it is missing.
    if (saveFits && gSystem->AccessPathName("DebugFit")) gSystem->mkdir("DebugFit", kTRUE);

    TFile* ifile = TFile::Open((filename + ".root").c_str(), "READ");
    if (!ifile || ifile->IsZombie()) {
        std::cerr << "Can't open " << filename << ".root" << std::endl;
        return;
    }
    TFile* ofile = TFile::Open((outfilename_ + ".root").c_str(), "RECREATE");
    if (!ofile || ofile->IsZombie()) {
        std::cerr << "Can't create " << outfilename_ << ".root" << std::endl;
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
    // ABCD plane: Fpixel on the horizontal axis, pT on the vertical one.
    // A and B sit in the low-pT band (55-70 GeV), C and D above 70 GeV.
    // The signal, if any, populates D only; the ABCD assumption is that the two
    // variables factorise, so N_D = N_B * N_C / N_A.
    //
    // The region names in the input file encode the Fpixel window, e.g.
    // "regionA_3fp8" = region A with 0.3 < Fpixel <= 0.8. The validation region
    // uses the 0.8-0.9 slice (unblinded), the search region the 0.9-1.0 slice.

    // TRUE: Ih and p templates taken in C, B keeps the normalisation
    //
    // Concretely this changes which region feeds what:
    //   - both the Ih and the 1/p templates come from C (see the two identical
    //     arguments in the bckgEstimate calls below);
    //   - the separate region passed as B_ifIhpSAME is used only for the ABCD
    //     normalisation factor and for the eta reweighting.
    // This is the configuration where the Ih and p cuts coincide, so B and C
    // would otherwise supply the same template twice.
    const bool ifIhpSAME = true;

    // All the knobs bundled once; only `blind` differs between VR and SR below.
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


    // Success flag reported at the very end. NOTE: it holds the result of the
    // *last* estimate that actually ran, so when both regions are requested the
    // SR outcome overwrites the VR one.
    bool done = false;

    // ---- Validation region: 0.8 < Fpixel <= 0.9, not blinded ----------------
    if (runVR) {
        std::cout << "\n    Loading validation region..." << std::endl;
        Region ra_3fp8, rb_8fp9, rc_3fp8, rd_8fp9, rbc_8fp9;
        // Every region is loaded with the same binning options so that the
        // templates remain bin-compatible with each other.
        // rbc_8fp9 is loaded from region D on purpose: the "BC" region only needs
        // the right binning for its (empty) prediction histograms.
        bool ok = true;
        ok &= loadHistograms(ra_3fp8,  ifile, "regionA_3fp8"  + Ext, bool_rebin, rebineta, rebinp, rebinih, TakeAbsEta);
        ok &= loadHistograms(rb_8fp9,  ifile, "regionB_8fp9"  + Ext, bool_rebin, rebineta, rebinp, rebinih, TakeAbsEta);
        ok &= loadHistograms(rc_3fp8,  ifile, "regionC_3fp8"  + Ext, bool_rebin, rebineta, rebinp, rebinih, TakeAbsEta);
        ok &= loadHistograms(rd_8fp9,  ifile, "regionD_8fp9"  + Ext, bool_rebin, rebineta, rebinp, rebinih, TakeAbsEta);
        ok &= loadHistograms(rbc_8fp9, ifile, "regionD_8fp9"  + Ext, bool_rebin, rebineta, rebinp, rebinih, TakeAbsEta);

        if (!ok) std::cerr << "VR: incomplete load -> region skipped" << std::endl;
        else {
            ofile->cd();          // every Write() inside bckgEstimate targets this file
            opt.blind = false;    // the VR is the closure test, it must stay visible
            std::cout << "    Background estimation, VR 8fp9..." << std::endl;
            // Argument order is (filename, calib, B, C, BC, A, D, ifIhpSAME,
            // B_ifIhpSAME, st, opt): rc_3fp8 is passed twice because ifIhpSAME
            // makes both templates come from C, while rb_8fp9 supplies the
            // normalisation. "8fp9" is the suffix appended to every output
            // histogram name.
            done = bckgEstimate(DataSetName, calib, rc_3fp8, rc_3fp8, rbc_8fp9, ra_3fp8, rd_8fp9,
                         ifIhpSAME, rb_8fp9, "8fp9", opt);
        }
    }

    // ---- Search region: 0.9 < Fpixel <= 1.0, blinded above 300 GeV ----------
    if (runSR) {
        std::cout << "\n    Loading search region..." << std::endl;
        Region ra_3fp9, rb_9fp10, rc_3fp9, rd_9fp10, rbc_9fp10;
        bool ok = true;
        ok &= loadHistograms(ra_3fp9,   ifile, "regionA_3fp9"  + Ext, bool_rebin, rebineta, rebinp, rebinih, TakeAbsEta);
        ok &= loadHistograms(rb_9fp10,  ifile, "regionB_9fp10" + Ext, bool_rebin, rebineta, rebinp, rebinih, TakeAbsEta);
        ok &= loadHistograms(rc_3fp9,   ifile, "regionC_3fp9"  + Ext, bool_rebin, rebineta, rebinp, rebinih, TakeAbsEta);
        ok &= loadHistograms(rd_9fp10,  ifile, "regionD_9fp10" + Ext, bool_rebin, rebineta, rebinp, rebinih, TakeAbsEta);
        ok &= loadHistograms(rbc_9fp10, ifile, "regionD_9fp10" + Ext, bool_rebin, rebineta, rebinp, rebinih, TakeAbsEta);

        if (!ok) std::cerr << "SR: incomplete load -> region skipped" << std::endl;
        else {
            ofile->cd();
            opt.blind = true;     // observed spectrum zeroed above 300 GeV
            std::cout << "    Background estimation, SR 9fp10..." << std::endl;
            done = bckgEstimate(DataSetName, calib, rc_3fp9, rc_3fp9, rbc_9fp10, ra_3fp9, rd_9fp10,
                         ifIhpSAME, rb_9fp10, "9fp10", opt);
        }
    }

    // Both regions write into the same output file, distinguished by the "_8fp9"
    // and "_9fp10" suffixes on the histogram names.
    ofile->Close();
    ifile->Close();
    delete ofile;
    delete ifile;
    if (!done) { std::cerr << "Can't perform the prediction -> empty output file" << std::endl; return; }
    // Sentinel line parsed by the launcher to decide whether the run succeeded.
    // It must remain the last thing printed on a successful run.
    std::cout << "\nDone: " << outfilename_ << ".root" << std::endl;

    return;
}