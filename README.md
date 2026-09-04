# massSpectrum_bckgPrediction

Data-driven background estimate of the **HSCP mass spectrum** for the Run-3 (2024) Heavy Stable
Charged Particle search, together with the plotting macros and the computation of the background
and signal systematic uncertainties.

The repository takes as input the histograms produced upstream by `TupleAnalysis`, and produces:

1. the predicted mass spectrum in the validation and signal regions (ROOT files),
2. publication-style mass-spectrum plots (observed vs. predicted, ratio, pull),
3. the relative systematic uncertainties on the background prediction and on the signal.

---

## 1. Method in a nutshell

The prediction relies on an **ABCD method** in the two-dimensional plane
(`Fpixel`, `pT`), which are uncorrelated for SM background but strongly correlated for signal.

```
   pT
    |            |             |
    |     C      |      D'     |     D        D  = signal region (blind)
    |    (CR)    |     (VR)    |  (blind)     D' = validation region
 70 |------------|-------------|--------      C  = control region: templates
    |            |             |                   are extracted here
    |     A      |      B'     |     B
 55 |____________|_____________|________
   0.3          0.8           0.9      1.0    Fpixel
```

**Normalisation.** The yield in the target region is obtained from the three other regions:

```
yield_D = yield_B * yield_C / yield_A
```

**Shape.** The mass spectrum is built by convolving two templates extracted in the control region
C: the momentum template `(1/p, η)` and the ionisation template `(η, Ih)`. Because `p`, `pT` and `η`
are correlated, the convolution is performed in bins of `η`. For every bin triplet `(i, j, k)` =
`(η, 1/p, Ih)` the mass and its weight are

```
M_ijk = p_ji * sqrt( (Ih_ik - C_mass) / K_mass )
w_ijk = ( n_ji / Σ_m n_mi ) * h_ik
```

**η reweighting.** The `η` distribution in C does not match the one in the target region, so the
templates in C are reweighted by `η_B / η_A` before the convolution.

**Template tails.** To suppress statistical fluctuations, the tail of the `Ih` template is replaced
by a Gaussian fit and the low tail of the `1/p` template by an error function (CDF of a Gaussian).

**Statistical uncertainty.** The whole prediction is repeated for `nPE` pseudo-experiments (default
200). In each toy, every template bin — and the control-region yield used for the normalisation — is
resampled from a Poisson distribution. The mean over the toys is the central value, the RMS is the
statistical uncertainty band.

The calibration constants `K_mass` and `C_mass` used in the mass formula are defined at the top of
`Regions.h` (`K_data2024 = 2.8202`, `C_data2024 = 2.9784`, `K_mc2024 = 2.83894`,
`C_mc2024 = 3.01756`).

---

## 2. Repository content

| File | Role |
| --- | --- |
| `LaunchBkgOn1Syst.py` | Launcher — runs **one** configuration (usually the nominal one). |
| `LaunchBkgOnAllSyst.py` | Launcher — runs **all** configurations (nominal + every systematic variation) sequentially. |
| `configFile_readHist_template.txt` | Template of the configuration line, filled in by the launchers via `sed`. |
| `step2_backgroundPrediction.C` | Main ROOT macro: reads the configuration, loads the regions, runs the estimate. |
| `Regions.h` | `Region` class (histogram container), mass formula, template correction functions, and the core `fillPredMass` convolution. |
| `CommonFunctions.h` | Helpers (rebinning, η reweighting, Poisson toys, mean/RMS over toys) and the `bckgEstimate` driver that loops over the pseudo-experiments. |
| `MyShowPred.py` | Wrapper that calls `MyMacroMass.py` on **data**. |
| `MyShowPred_MC.py` | Same, for **MC** samples (stacked W+jets / tt̄ / QCD). |
| `MyMacroMass.py` | Mass-spectrum plotter: observed vs. prediction, uncertainty band, signal overlays, ratio and pull panels. |
| `systBckg.py` | Combines the systematic variations of the **background** prediction. |
| `systSignal.py` | Computes the systematic uncertainties on the **signal** mass spectrum. |
| `tdrstyle.py`, `CMS_lumi.py` | Standard CMS plotting style helpers. |

---

## 3. Setup

The code is run inside a CMSSW area (`Regions.h` includes `CommonTools/UtilAlgos/interface/TFileService.h`,
so a `cmsenv` is required):

```bash
export SCRAM_ARCH=<arch matching the release>
cmsrel CMSSW_15_0_13_patch1
cd CMSSW_15_0_13_patch1/src/
cmsenv
```

Cloning requires an SSH key associated with your GitHub account
(see [connecting-to-github-with-ssh-key](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent)):

```bash
git clone -b master git@github.com:dapparu/massSpectrum_bckgPrediction.git massSpectrum_bckgPrediction
cd massSpectrum_bckgPrediction
mkdir -p DebugFit       # required if saveFits = true
```

---

## 4. Overall workflow

```
   TupleAnalysis output (.root with per-region histograms)
                        |
                        v
  [1]  LaunchBkgOn1Syst.py  /  LaunchBkgOnAllSyst.py
            -> configFile_readHisto_toLaunch.txt
            -> root -l -q -b step2_backgroundPrediction.C
                        |
                        v
       one .root file per configuration (nominal + systematic variations)
                        |
            +-----------+------------+
            |                        |
            v                        v
  [2] MyShowPred(.py/_MC.py)   [3] systBckg.py    (background systematics)
       -> MyMacroMass.py           systSignal.py  (signal systematics)
       -> mass spectrum plots      -> sysTotBinned_*.root + summary plots
```

---

## 5. Step 1 — Running the background prediction

### Single configuration

```bash
python LaunchBkgOn1Syst.py
```

Everything is steered from two files.

**In `LaunchBkgOn1Syst.py`:**

- `datasetList` — list of input samples. Give the path to the ROOT file **without** the `.root`
  extension (the macro appends it).
- `nPE` — number of pseudo-experiments used for the statistical uncertainty (default `200`).
- `config` — the configurations to run. Uncomment the line(s) you want; each entry is

  ```
  [label, rebinEta, rebinIh, rebinMom, fitIh, fitMom, useFit, corrTemplateIh, corrTemplate1oP]
  ```

  | Field | Meaning |
  | --- | --- |
  | `rebinEta`, `rebinIh`, `rebinMom` | Rebinning factors of the `η`, `Ih` and `1/p` templates. Nominal: `4, 4, 2`. A **smaller** factor means finer bins (`up`), a **larger** factor coarser bins (`down`). |
  | `fitIh`, `fitMom` | `1` = nominal fit, `0` = fitted value shifted **down** by its uncertainty, `2` = shifted **up**. |
  | `useFit` | `1` = use the tail fits, `0` = use the raw templates (no fit). |
  | `corrTemplateIh`, `corrTemplate1oP` | `1` = apply the linear correction for the residual `Fpixel` dependence of the template. |

The launcher copies `configFile_readHist_template.txt`, substitutes each placeholder with `sed`, and
writes `configFile_readHisto_toLaunch.txt`, which the ROOT macro then reads. The configuration line is

```
sample  nPE  rebin  rebinEta  rebinIh  rebinMom  fitIh  fitMom  useFit  corrTemplateIh  corrTemplate1oP
```

**In `step2_backgroundPrediction.C`:**

- `st_sample` — `"data2024"` or `"mc2024"`; selects the `K_mass` / `C_mass` constants.
- `Ext` — **one of the `Ext` lines must be uncommented.** It selects the η category and the
  pre-selection variant, and is appended to the histogram names read from the input file
  (`_Eta1`, `_Eta1_2p4`, `_Eta2p4`, plus optional `EoP` / `SigmaPtoverPt` variants). `etaName` is
  derived from it and propagated to the output file name and to the template-correction parameters.
- The `bckgEstimate(...)` call at the bottom selects the region configuration. Two are available:
  - **validation region** (`8fp9`): templates from `C(0.3–0.8)`, normalisation from `A(0.3–0.8)`,
    `B(0.8–0.9)`, prediction compared to `D'(0.8–0.9)`, `blind = false`;
  - **signal region** (`9fp10`): templates from `C(0.3–0.9)`, normalisation from `A(0.3–0.9)`,
    `B(0.9–1.0)`, prediction compared to `D(0.9–1.0)`, `blind = true` (data above 300 GeV is masked).
- `saveFits` — writes every template fit into `DebugFit/` for inspection.
- `MyIhCut` — optional additional `Ih` threshold; index `0` (`_IhC`) is the nominal one.

The toys are run in parallel with `ROOT::TProcessExecutor`; the number of workers is set in
`CommonFunctions.h` (`ROOT::TProcessExecutor workers(26);`).

### All configurations at once

```bash
python LaunchBkgOnAllSyst.py
```

This script embeds the same code as a string and re-runs it once per `config` line, uncommenting them
one at a time. It takes several hours — launch it inside a `screen`/`tmux` session.

### Output

One ROOT file per configuration, named from the full input path plus the configuration:

```
<input_path>_rebinEta4_rebinIh4_rebinP2_EtaReweighting_<etaName>_OldFit_IhC.root
```

with `_corrTemplateIh`, `_corrTemplate1oP`, `_fitIhUp/Down`, `_fitPUp/Down` inserted for the
corresponding variations, and `_NoFit` instead of `_OldFit` when `useFit = 0`. Since the name is
built from the input path, the file is created **next to the input sample**; the plotting and
systematics scripts expect these files to be collected under
`.../macros/DataMET_2024_<version>__<region><option>/<etaName>/`.

The file contains, among others:

| Histogram | Content |
| --- | --- |
| `mass_predBC_<region>` | Predicted mass spectrum (mean of the toys, RMS as error). |
| `mass_obs_<region>` | Observed mass spectrum (blinded above 300 GeV in the SR). |
| `mass_predBCR_<region>` | Right-integral ratio observed/predicted. |
| `h_norm` | Distribution of the ABCD normalisation over the toys. |
| `ih_eta_mean`, `oP_eta_mean`, `ih_VR`, `eta_VR` | Mean toy templates and their target-region counterparts, used for the correlation studies. |

---

## 6. Step 2 — Plotting the mass spectra

```bash
python MyShowPred.py        # data
python MyShowPred_MC.py     # MC (stacked W+jets, tt->lv, tt->2l2v, QCD)
```

Both are thin wrappers: set `regions`, `eta`, `version`, `option`, `OnlyNominal` and `isMC` at the
top, and they build the input path and call

```bash
python3 MyMacroMass.py --ifile <input> --labelName <label> --ofile <name> \
                       --region <region> --odir <outdir> --nom <bool> \
                       --eta <etaName> --isMC <bool>
```

| Option | Meaning |
| --- | --- |
| `--ifile` | Input ROOT file produced at step 1. |
| `--labelName` | Histogram-name suffix used in the input file (`METanalysis_..._<eta>`). |
| `--ofile` | Output file label. |
| `--region` | `8fp9` (VR) or `9fp10` (SR). |
| `--odir` | Output directory (created if missing). |
| `--nom` | `True` = statistical band only; `False` = add the total systematic band. |
| `--eta` | `Eta1`, `Eta1_2p4` or `Eta2p4`; drives the label drawn on the plot. |
| `--isMC` | `True` = stack the individual MC processes as the "observed" points. |

The output is a `.pdf` / `.root` / `.C` canvas with the mass spectrum, the data-driven prediction and
its uncertainty band, the signal overlays, and two lower panels: the bin-by-bin ratio
`obs/pred` and the pull `(N_obs − N_pred)/sqrt(N_obs + σ_pred²)`.

> **Note.** `MyMacroMass.py` contains hard-coded paths to the signal samples (`Gluino_V19`), to the MC
> samples used in the stack, and to the systematics file `sysTotBinned_2024_<region>.root`. Update
> them when the sample versions change.

---

## 7. Step 3 — Systematic uncertainties

### Background (`systBckg.py`)

```bash
python systBckg.py
```

Configure the block at the bottom of the file (`version`, `etaname`, `year`, `sample`, `region`,
`option*`, `directory`). The script opens the nominal file **and every variation file** produced at
step 1, rebins them onto the analysis mass binning, normalises them, and computes for each source

```
syst = 100 * max( |1 - varUp/nominal| , |1 - varDown/nominal| )   [%]
```

both bin-by-bin (`*_binned`) and on the right integral. The considered sources are:

| Source | Variations used |
| --- | --- |
| Statistical | RMS of the toys of the nominal prediction |
| `η` binning | `rebinEta = 2 / 8` |
| `Ih` binning | `rebinIh = 2 / 8` |
| `1/p` binning | `rebinP = 1 / 4` |
| `Ih` tail fit | `fitIhUp` / `fitIhDown` |
| `1/p` tail fit | `fitPUp` / `fitPDown` |
| `Ih` template correlation | `corrTemplateIh` (single variation, symmetrised) |
| `1/p` template correlation | `corrTemplate1oP` (single variation, symmetrised) |

The `NoFit` variation is computed but **not** included in the total (see `listOfSyst`).
The total is the quadratic sum of the individual sources.

Outputs, written under `<directory>/SystCombined/`:

- `sysToTBinned_<year>_<region>.root` — all individual sources + `systTotal` / `systTotalBinned`;
- `massShapePred.root` — nominal prediction together with the totals;
- `individualSyst/` — one nominal-vs-up-vs-down plot per source;
- `pdf/`, `Cfile/`, `rootfile/` — the summary plot of all sources vs. mass bin.

If `onlyNominal = True`, only the statistical uncertainty is computed and the script exits early —
useful before all the variation files have been produced.

### Signal (`systSignal.py`)

```bash
python systSignal.py
```

Configure `codeVersion`, `SignalSamples`, `idir`, `RegS`, `etaS` and `option`. Unlike the background
case, the variations are **not** re-derived here: they are read directly from the signal ROOT files as

```
METanalysis_TestPUppiMETCut_<option><eta>_<region>_SignalMass_<variation>
```

with `<variation>` in `nominal`, `PUUp/PUDown`, `TriggerSFUp/TriggerSFDown`, `JetUp/JetDown`,
`KUp/KDown`, `CUp/CDown`. Two normalisation-only uncertainties are hard-coded in `HARDCODED_SYST`
at the top of the file:

- luminosity: **1.4 %**
- `Fpixel` (pixel calibration): **1.6 %**

For each gluino mass point the script writes `systSignal/<sample>_<region>_<eta>/sysTotBinned_signal.root`,
one up/down comparison plot per source, and a summary plot of all sources. A final plot,
`systSignal/sysTot_allMasses_<eta>.pdf`, compares the total uncertainty for
`m_gluino = 2000, 2400, 2600 GeV`.

---

## 8. Expected input histograms

`loadHistograms` (in `CommonFunctions.h`) reads, for each region name
`region{A,B,C,D}_<FpixelRange><Ext>`:

| Name | Type | Axes |
| --- | --- | --- |
| `eta_1oP_<region>` | `TH2F` | X = `10^4/p` [GeV⁻¹], Y = `η` |
| `ih_eta_<region>` | `TH2F` | X = `η`, Y = `Ih` [MeV/cm] |
| `ih_p_<region>` | `TH2F` | X = `p`, Y = `Ih` |
| `ias_p_<region>`, `ias_pt_<region>` | `TH2F` | X = `p` / `pT`, Y = `Ias` |
| `mass_<region>` | `TH1F` | mass |
| `mass_eta_<region>` | `TH2F` | X = mass, Y = `η` |
| `mass_ih_<region>` | `TH2F` | X = `Ih`, Y = mass |

with `<FpixelRange>` among `3fp8` (0.3 < Fpixel ≤ 0.8), `3fp9`, `8fp9` and `9fp10`.

The `η` rebinning uses variable-width bin edges defined in `loadHistograms`, selected from the region
name (`Eta1`, `Eta1_2p4`, `Eta2p4`, `Eta1p2_2p2`, `Eta1p2_2p4`). Note that the `1/p` rebinning factor
is remapped internally: `rebinMom = 1 → 4`, `2 → 6`, `4 → 8`.

---

## 9. Practical notes

- The number of parallel workers (26) and the `nPE` value drive the runtime and the memory usage;
  reduce them when running on a shared machine.
- `saveFits = true` writes one ROOT file per worker in `DebugFit/`; the directory must exist
  beforehand and the files add up quickly.
- Bad fits are not fatal: the code prints a diagnostic (`Bad fit Ih ...` / `Bad fit 1/p ...`) and
  falls back to the raw template for that `η` slice.
- `TakeAbsEta` in `step2_backgroundPrediction.C` folds the templates onto `|η|`; it is off by default.
- `systErr_` in `Regions.h` is set to `0` so that the ratio plots contain the statistical uncertainty
  only during the systematic studies.