# massSpectrum_bckgPrediction

Data-driven background estimate of the **HSCP mass spectrum** for the Run-3 (2024) Heavy Stable
Charged Particle search, together with the systematic uncertainties (background and signal) and the
final mass-spectrum plots.

The repository takes as input the per-region histograms produced upstream by `TupleAnalysis`
(step 1) and produces:

1. the predicted mass spectrum in the validation and/or search region — one ROOT file per
   systematic variation;
2. the combined relative systematic uncertainties, on the background prediction and on the gluino
   signal samples;
3. publication-style mass-spectrum plots (observation vs. prediction, uncertainty band, signal
   overlays, ratio and pull panels).

Two helper scripts additionally normalise the simulated samples (MC backgrounds and signal) to the
integrated luminosity; their output is what the signal systematics and the plots read.

Each step is driven entirely from the command line or from a small settings block at the top of a
script. **Nothing in the C++ ever has to be edited to change a configuration.**

---

## 0. Disclaimer

I wrote all the code in this project myself. However, this README and the comments in the scripts were generated using Claude.

---

## 1. Method in a nutshell

The prediction relies on an **ABCD method** in the two-dimensional plane (`Fpixel`, `pT`), which are
uncorrelated for SM background but strongly correlated for a heavy signal.

```
   pT
    |            |             |
    |     C      |      D'     |     D        D  = search region (blinded above 300 GeV)
    |    (CR)    |     (VR)    |  (blind)     D' = validation region
 70 |------------|-------------|--------      C  = control region: the templates
    |            |             |                   are extracted here
    |     A      |      B'     |     B
 55 |____________|_____________|________
   0.3          0.8           0.9      1.0    Fpixel
```

**Normalisation.** The yield in the target region is obtained from the three other regions:

```
N_D = N_B * N_C / N_A
```

**Shape.** The mass spectrum is built by convolving two templates extracted in the control region:
the momentum template `(1/p, η)` and the ionisation template `(η, Ih)`. Because `p`, `pT` and `η`
are correlated, the convolution is performed **in bins of `η`**. For every bin triplet
`(i, j, k) = (η, 1/p, Ih)` the mass and its weight are

```
M_ijk = p_ji * sqrt( (Ih_ik - C) / K )
w_ijk = ( n_ji / Σ_m n_mi ) * h_ik
```

Only the `1/p` template is normalised to unit area; the absolute scale of the prediction is imposed
afterwards by the ABCD factor above.

**η reweighting.** The `η` distribution of the control region does not match the one of the target
region, so the template is reweighted by `η_B / η_A` before the convolution.

**Template tails.** To suppress statistical fluctuations where the statistics run out, the tail of
the `Ih` template is replaced by a fitted function, and the low tail of the `1/p` template (i.e. the
*high*-momentum region) likewise. A fit is used for a given bin only when it is sparsely populated
(< 100 entries) and sits beyond the switching threshold; elsewhere the raw bin contents are used.
When a fit replaces the data, the bin is sub-divided into 5 slices so that the strongly varying
mass mapping is sampled finely enough.

| Template | Shape | Fit range | Used below/above |
| --- | --- | --- | --- |
| `Ih`, `useOldIhFit = 0` | Gaussian | `[1.1 × peak, 6]` MeV/cm | above `1.1 × peak` |
| `Ih`, `useOldIhFit = 1` | legacy shape | `[3, 6]` MeV/cm | above `3.5` MeV/cm |
| `1/p`, `useOld1oPFit = 1` | `[0]·([1] + erf((log x − [2])/[3]))` | `[0, f × peak]`, `f` scanned over `0.9 … 0.5` until the fit converges | below `0.2 × f × peak` |
| `1/p`, `useOld1oPFit = 0` | `0.5·(e^(ax²+bx) + e^(−ax²−bx)) − 1` | `[0, 0.6 × peak]` | below `0.2 × f × peak` |

The `0.2` factor is a deliberate safety margin: the fit is only trusted well inside its fitted
range, never up to its upper bound. The 1/p fit is seeded from a pre-fit of the full control-region
spectrum, which is what makes the per-toy fits converge reliably. A failing fit is **not** fatal: a
diagnostic is printed (`Bad fit Ih …` / `Bad fit 1/p …`) and that `η` slice falls back to the raw
template.

**Statistical uncertainty.** The whole prediction is repeated for `nPE` pseudo-experiments (default
200), run in parallel with `ROOT::TProcessExecutor`. In each toy every template bin — and the yields
entering the ABCD normalisation — is resampled from a Poisson distribution. The bin-by-bin mean over
the toys is the central value, the RMS is the statistical uncertainty. Seeds are derived from the
toy index (`seed = 1 + workerID`), so a run is **bit-for-bit reproducible**.

**Mass calibration.** `M = p · sqrt((Ih − C)/K)`, with `(K, C)` selected from the `sampleType` key
in `GetDeDxCalib()` (`Regions.h`):

| `sampleType` | K | C |
| --- | --- | --- |
| `data2017` | 2.54 | 3.14 |
| `data2018` | 2.55 | 3.14 |
| `data2024` | 2.8202 | 2.9784 |
| `mc2017` | 2.48 | 3.19 |
| `mc2018` | 2.49 | 3.19 |
| `mc2024` | 2.83894 | 3.01756 |
| `ttbar2024` | 2.83894 | 3.01756 |

An unknown key falls back to the 2024 data constants **with a warning only**, so always check that
`sampleType` matches the dataset: reading MC with the data calibration shifts the whole mass
spectrum silently.

---

## 2. Repository content

| File | Role |
| --- | --- |
| `RescaleBKG.py` | **Step 1b.** Scales the histograms of the simulated background samples (QCD, ttbar, W+jets) to the integrated luminosity and writes them to `<sample>_weighted.root`. |
| `RescaleSignal.py` | **Step 1b.** Same for the signal samples (gluino, stop, stau). |
| `LaunchBkgPred.py` | **Step 2 driver.** Writes one config file per systematic variation, runs the macro on it, and files the outputs into a structured directory tree. |
| `BkgPrediction.C` | Main ROOT macro: parses the config file, resolves the calibration and the histogram names, loads the ABCD regions, runs the estimate for the VR and/or the SR, writes one output file. |
| `Regions.h` | `Region` class (histogram container), mass formula, dE/dx calibration table, `Fpixel`-correlation corrections, blinding, and the core `fillPredMass()` convolution. |
| `CommonFunctions.h` | Histogram I/O helpers, η rebinning tables, `|η|` folding, `loadHistograms()`, the pseudo-experiment machinery, the η reweighting, and `bckgEstimate()` — the top-level driver of the toys. |
| `systBckg.py` | **Step 3a.** Combines the systematic variations of the background prediction into a single uncertainty file. |
| `systSignal.py` | **Step 3b.** Computes the systematic uncertainties of the gluino signal samples. |
| `ShowPlots.py` | **Step 3c driver.** Resolves paths and calls the plotting macro, once per η range. |
| `PlottingMacro.py` | Mass-spectrum plotter: observation vs. prediction, uncertainty band, signal overlays, ratio and pull panels. |
| `tdrstyle.py` | Standard CMS plotting style. |

---

## 3. Setup

The C++ is **plain ROOT** — no CMSSW dependency remains (`Region::initHisto` is templated on the
directory type, so the header compiles both inside step 1 and standalone here). Any ROOT 6
installation with a C++17 compiler works; sourcing a CMSSW area is only convenient because the
input paths live there:

```bash
export SCRAM_ARCH=<arch matching the release>
cmsrel CMSSW_15_0_13_patch1
cd CMSSW_15_0_13_patch1/src/
cmsenv
```

Cloning requires an SSH key associated with your GitHub account
(see [connecting-to-github-with-ssh-key](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent)):

```bash
git clone -b master git@github.com:gcoulon/massSpectrum_bckgPrediction.git massSpectrum_bckgPrediction
cd massSpectrum_bckgPrediction
```

`DebugFit/` is created automatically when `saveFits = 1`, so nothing has to be prepared by hand.

The two rescaling scripts (§5) import `tqdm` and `numpy` on top of PyROOT: both must be importable
for the scripts to start, even though neither is actually used.

---

## 4. Overall workflow

```
   step1 (TupleAnalysis): one .root per sample, with the per-region histograms
                        |
                        +---- simulated samples only ----+
                        |                                v
                        |             [1b] RescaleBKG.py / RescaleSignal.py
                        |                  -> <sample>_weighted.root, next to the input
                        |                     (read by [3b] and [3c])
                        v
  [2]  LaunchBkgPred.py
            -> configFile_<label>.txt   (one per variation)
            -> root -l -q -b 'BkgPrediction.C+("configFile_<label>.txt")'
            -> one .root per variation, filed under
               <outputDir>/<sampleType>_V<version>__<region>_<labelDir>/<etaRange>/
                        |
            +-----------+-------------------------+
            |                                     |
            v                                     v
  [3a] systBckg.py                        [3b] systSignal.py
       -> <indir>/SystCombined/                -> systSignal/Gluino_<m>_.../
          sysTotBinned_<eta>_<region>.root        sysTotBinned_signal.root
            |
            v
  [3c] ShowPlots.py -> PlottingMacro.py
       -> <indir>/<eta>/Plots_<region>/*.pdf|.root|.C
```

---

## 5. Step 1b — Normalising the simulated samples to the luminosity

```bash
python3 RescaleBKG.py        # simulated backgrounds: QCD, ttbar, W+jets
python3 RescaleSignal.py     # signal: gluino, stop, stau
```

The step 1 histograms of a simulated sample hold raw event counts. These two scripts write, for
every input file, a copy in which each histogram is scaled to the yield expected in the 2024
dataset:

```
w = L_int × σ / N_preTrig
```

| Term | Taken from |
| --- | --- |
| `L_int` | `intLumi` = 108 950 pb⁻¹ (2024) |
| `σ` | `crossSectionArray`, in pb, looked up by the base name of the file |
| `N_preTrig` | first bin of the `EventCutflow` histogram of that file, i.e. the number of events before the trigger |

Data are never rescaled, so this step is **not** needed to run the prediction on data. It is needed
by everything that reads a simulated sample:

| Reader | `_weighted.root` files it opens |
| --- | --- |
| `PlottingMacro.py`, always | The three gluino overlays (2000, 2400, 2600 GeV). They are opened unconditionally, so even a data plot fails without them. |
| `PlottingMacro.py`, MC mode | W+jets, ttbar dileptonic, ttbar semileptonic and QCD, which make up the stack and the "observation". |
| `systSignal.py` | Every gluino mass point, unless `--raw` is given. |

### 5.1 What a pass does

For each file of the active list (`BackgroundSamples_HISTO` — the variable keeps that name in
`RescaleSignal.py` too, where it holds signal files):

1. the file is skipped, with a one-line message, if it does not exist, has no entry in
   `crossSectionArray`, has no `EventCutflow`, or has an empty first cut-flow bin;
2. the weight is computed and printed (`<file> -> nEvetsPreTrig = … -> weight = …`);
3. `<name>_weighted.root` is written **next to the input**, opened `RECREATE`. It contains
   - every top-level histogram scaled by `w` (all `TH1`-derived classes, `TH2` included), except
     those whose name contains one of the `SelectionNames` tags;
   - every other top-level object, copied unchanged;
   - two `TParameter<double>`, `weight` and `lumi`, recording the normalisation that was applied.

The input file is opened read-only and is never modified. A skipped file is not an error: the
script exits with status 0 either way, so the printout is the only record of what was actually
processed.

### 5.2 Settings

Neither script takes command-line options; everything is set in the block at the top.

| Variable | Meaning |
| --- | --- |
| `version` | Version string inserted into **every** file name, both in the cross-section table and in the sample list (`…_V<version>….root`). |
| `intLumi` | Integrated luminosity, in pb⁻¹. |
| `crossSectionArray` | `{file base name: σ in pb}`. |
| `pathBKG_*`, `pathSig*` | Step 1 output directory of each sample family. |
| `BackgroundSamples_HISTO` | The files processed by this pass. The complete list is kept in the triple-quoted block just above it; the active list is the subset to run. |
| `SelectionNames` | Sub-strings of histogram names to leave unscaled (e.g. `"Nm1"`, `"EventCutflow"`). Empty by default: everything is scaled. |

Sample families, and the file names the scripts expect:

| Script | Family | File name | Points |
| --- | --- | --- | --- |
| `RescaleBKG.py` | QCD (`mu`) | `QCD2024_mu_pt<range>_V<version><i>.root` | 12 pT bins, from 15 to above 1000 GeV, `<i>` = 1…12 |
| | ttbar | `TTbar2024_V<version>.root`, `TTbarSemiLep2024_V<version>.root` | dileptonic, semileptonic |
| | W+jets, binned | `Wjets2024_<1J\|2J>_pt<range>_V<version><i>.root` | 2 × 5 pT bins, from 40 to above 600 GeV, `<i>` = 1…10 |
| | W+jets, inclusive | `WjetMuNu2024_V<version>.root` | — |
| `RescaleSignal.py` | gluino | `Gluino_Run3_MET_<m>_V<version>.root` | 9 points, 1000 → 2600 GeV |
| | gluino (MadGraph) | `Gluino_Run3_MET_madgraph_<m>_V<version>.root` | 10 points, 1100 → 2600 GeV |
| | stop (MadGraph) | `Stop_Run3_MET_madgraph_<m>_V<version>.root` | 12 points, 700 → 2600 GeV |
| | stau | `Stau_Run3_MET_<m>_V<version>.root` | 11 points, 247 → 1599 GeV |

In the binned samples the bin index `<i>` is appended directly to the version, with no separator:
`version = "16p3"` gives `…_pt15to20_V16p31.root`, …, `…_pt1000_V16p312.root`.

**One family per pass.** `version` is a single variable, but the families were not produced at the
same version: the files currently read by `PlottingMacro.py` are QCD `16p3`, ttbar dileptonic
`15p10`, ttbar semileptonic `22p1`, W+jets `14p14` and gluino `19p12` (`VSIGNAL`), with stau `20p0`
and stop `21p0` in its commented-out alternatives. So, for each family: set `version`, put the
matching block in `BackgroundSamples_HISTO`, run, and repeat. A file whose name does not match is
simply reported as `File not found` and skipped. As committed, `RescaleBKG.py` is set up for QCD
(`16p3`) and `RescaleSignal.py` for the staus (`20p0`).

### 5.3 Things to know

- **QCD has to be merged by hand.** `RescaleBKG.py` writes one weighted file per pT bin
  (`QCD2024_mu_pt15to20_V16p31_weighted.root`, …), while `PlottingMacro.py` reads a single
  `QCD2024_mu_V16p3_weighted.root`. No script of this repository builds it: add the weighted bins
  together once they exist, e.g.

  ```bash
  hadd QCD2024_mu_V16p3_weighted.root QCD2024_mu_pt*_V16p3*_weighted.root
  ```

  Always merge the **weighted** files: each bin has its own weight, so a sum of the raw files
  cannot be rescaled afterwards.
- **The weighted copies go stale.** They are independent files: re-running step 1 on a simulated
  sample does not refresh them, and an outdated copy is picked up downstream without any warning.
  Re-run the rescaling after every step 1 production.
- **`EventCutflow` is scaled as well** with the default empty `SelectionNames`: in a weighted file
  its first bin reads `L_int × σ`, no longer the number of generated events. Add `"EventCutflow"`
  to `SelectionNames` to keep the raw cut-flow.
- **Statistical errors.** `RescaleSignal.py` explicitly enables `Sumw2` and the normal error option
  before scaling, so a bin filled with `N` unweighted entries carries `w·√N`. `RescaleBKG.py` calls
  `Scale` directly and relies on ROOT 6 enabling `Sumw2` by itself inside `TH1::Scale`.
- **Luminosity.** The samples are weighted to 108.95 fb⁻¹, whereas `PlottingMacro.py` rescales the
  signal overlays by `lumi / 109.0`, i.e. assumes 109.0 fb⁻¹: a 0.05 % mismatch, negligible but to
  be kept in mind if either number changes. Only the signal overlays are rescaled when an era
  (`F`, `G`) is selected; the MC stack always corresponds to the full 2024 luminosity.
- **Only the MadGraph gluinos are read downstream.** `systSignal.py` and `PlottingMacro.py` both
  point at `Gluino_Run3_MET_madgraph_<m>_V<vsignal>_weighted.root`. The other gluino set, the stops
  and the staus are weighted but not consumed by any script here (`PlottingMacro.py` only keeps two
  commented-out lines to overlay a stau or a stop point instead of a gluino one).
- The complete sample list kept as a comment in `RescaleBKG.py` lacks a comma before its last entry
  (`WjetMuNu2024`), so it cannot be un-commented as it stands: copy the block of the family you
  need instead.

---

## 6. Step 2 — Running the background prediction

```bash
python LaunchBkgPred.py                    # nominal only (fast path)
python LaunchBkgPred.py --all              # full systematic sweep (14 runs)
python LaunchBkgPred.py --only binEtaUp,fitIhDown
```

`--only` wins over `--all`. An unknown label aborts immediately and prints the valid ones.
The full sweep takes several hours — launch it inside a `screen` / `tmux` session.

### 6.1 What the launcher does

For every (dataset, variation) pair it:

1. merges three layers of settings into one flat dict —
   `common` (session settings) → `nominal` (binning and fit defaults) → the variation's own
   overrides. The merge order means a key present in both `common` and `nominal` would be won by
   `nominal`; **keep the two sets disjoint**;
2. writes that dict as a `key = value` file, `configFile_<label>.txt`;
3. runs `root -l -q -b 'BkgPrediction.C+("configFile_<label>.txt")'`;
4. collects the freshly written `.root` files and moves them into the destination directory.

On success the config file is deleted; on failure it is kept so the run can be replayed verbatim.
The macro echoes the whole parsed configuration into its log, which is then the only record of what
a given output was produced with.

### 6.2 Settings

Everything below is edited in the `common` dict at the top of `LaunchBkgPred.py`, except
`datasetList`, `outputDir`, `etaRangeName`, `sampleTypeName` and `labelDir`, which sit just above
it.

| Key | Values | Meaning |
| --- | --- | --- |
| `sample` | full path **without** `.root` | Input file. Set from `datasetList`. |
| `sampleType` | `data2017` … `mc2024`, `ttbar2024` | Selects `(K, C)`. Carried per dataset, so data and MC cannot be mixed up. `ttbar2024` has the `mc2024` constants: what it changes is the name of the working directory (`ttbar2024_V<version>__…`), used for the ttbar-only closure test of §9. |
| `label` | free string, **must be unique** | Appended to the output file name. An empty label is a hard error. |
| `nPE` | integer > 0 | Pseudo-experiments. Runtime scales linearly; use 20 for a quick test, 200 for a result. |
| `rebin` | 0 \| 1 | Master switch. At 0 the step 1 histograms are used as they are and `rebinEta/Ih/Mom` are ignored. |
| `rebinEta` | 2 \| 4 \| 8 | Fine / nominal / coarse **variable** η binning table. Any other value triggers a plain uniform `Rebin2D` by that factor. |
| `rebinIh` | integer | Merging factor on the `Ih` axis. Nominal 4. |
| `rebinMom` | 1 \| 2 \| 4 | Code, remapped internally to an actual merging factor of 4 / 6 / 8 on the `1/p` axis. Nominal 2. |
| `fitIh`, `fitMom` | 0 \| 1 \| 2 | 1 = nominal, 2 = +1σ, 0 = −1σ on the fitted template, propagated with `TF1::IntegralError` from the MINUIT covariance matrix. |
| `useFit` | 0 \| 1 | 0 = pure bin-by-bin convolution, no tail fit at all. |
| `corrTemplateIh`, `corrTemplate1oP` | 0 \| 1 | Apply the linear `Fpixel`-correlation correction to the corresponding template. |
| `etaRange` | `Eta1` \| `Eta1_2p4` \| `Eta2p4` \| `Eta1p2_2p2` \| `Eta1p2_2p4` | Must exist on the step 1 side. |
| `eopCut` | `""` \| `0p1` | step 1 E/p cut. |
| `sigmaPtCut` | `""` \| `0p5` | step 1 σ(pT)/pT cut. |
| `useOldIhFit` | 0 \| 1 | Shape of the `Ih` tail fit (see the table in §1). |
| `useOld1oPFit` | 0 \| 1 | Shape of the `1/p` tail fit. **Careful:** this does not appear in the output name, so two runs differing only by this setting overwrite each other. |
| `saveFits` | 0 \| 1 | Dump every fit into `DebugFit/`. At 1 this is `nPE` files **per systematic** (200 × 14 = 2800 with `--all`). |
| `takeAbsEta` | 0 \| 1 | Fold the templates onto `|η|` before the estimate. |
| `runVR` | 0 \| 1 | Validation region, `0.8 < Fpixel ≤ 0.9`, unblinded. |
| `runSR` | 0 \| 1 | Search region, `0.9 < Fpixel ≤ 1.0`, blinded above 300 GeV. |
| `nWorkers` | integer | Forked processes for the toys. Must stay ≤ the number of available cores. |


`runVR` and `runSR` write into the **same** output file, with the `_8fp9` / `_9fp10` suffixes on the
histogram names. Set both to 1 to get them together; running two separate passes does not work,
because the output is opened `RECREATE` and the second pass overwrites the first.

### 6.3 The systematic variations

The `config` list holds one entry per variation, listing only the keys that change. The first entry
is the nominal and is the one that runs without `--all`.

| Label | Override | Source |
| --- | --- | --- |
| `nominal` | — | reference |
| `binEtaUp` / `binEtaDown` | `rebinEta = 2 / 8` | η binning granularity |
| `binIhUp` / `binIhDown` | `rebinIh = 2 / 8` | `Ih` binning granularity |
| `binMomUp` / `binMomDown` | `rebinMom = 1 / 4` | `1/p` binning granularity |
| `fitIhUp` / `fitIhDown` | `fitIh = 2 / 0` | ±1σ on the `Ih` tail fit |
| `fitMomUp` / `fitMomDown` | `fitMom = 2 / 0` | ±1σ on the `1/p` tail fit |
| `noFit` | `useFit = 0` | cross-check: contribution of the fits (not combined by `systBckg.py`, see §7) |
| `corrTemplateIh` | `corrTemplateIh = 1` | `Fpixel` correlation of the `Ih` template |
| `corrTemplate1oP` | `corrTemplate1oP = 1` | `Fpixel` correlation of the `1/p` template |

The nominal binning and fit values live in the `nominal` dict: an entry in `config` only overrides
what it mentions, everything else falls back there.

### 6.4 Output

One ROOT file per configuration, written **next to the input sample** (the `sample` key is a full
path, so the output name is one too):

```
<dataset path>_<etaRange>[_SigmaPtoverPt_<x>][_EoP_<y>]_<label>.root
```

The launcher then moves it to

```
<outputDir>/<sampleType>_V<version>__<region>_<labelDir>/<etaRange>/
```

where `<version>` is what follows `_V` in the dataset base name (`…_V12p35` → `12p35`) and
`<labelDir>` gains a `SigmaPtoverPt_<x>_EoP_<y>_` prefix when both cuts are active. Everything that
distinguishes one production from another is in the path, so outputs never silently overwrite each
other.

Main content of the file, with `<r>` = `8fp9` or `9fp10`:

| Histogram | Content |
| --- | --- |
| `mass_predBC_<r>` | **Predicted mass spectrum** — bin-by-bin mean of the toys, RMS as the error. |
| `mass_obs_<r>` | Observed mass spectrum (zeroed above 300 GeV in the SR). |
| `mass_predBCR_<r>` | Cumulative ratio observation / prediction (`∫` from bin *i* to the end). |
| `pred_mass_eta_mean_<r>` | Predicted mass vs. η, and its η projection. |
| `h_norm_<r>` | Distribution of the ABCD closure over the toys (predicted / observed yield); should peak at 1. |
| `toys_<r>` | Canvas with every individual toy superimposed. |
| `ih_eta_mean_<r>`, `ih_VR_<r>` | Toy-mean `Ih` template vs. what is observed in D. |
| `oP_eta_mean_<r>`, `eta_VR_<r>` | Same for the reweighted `1/p` template. |
| `mass_C_<r>` | Control-region mass spectrum rescaled to the ABCD normalisation: shows what the prediction would look like without the convolution. |
| `mass_eta_C_<r>`, `mass_eta_D_<r>` | Observed mass vs. η in C and D, plus their η projections. |
| `eta_1oP_*`, `ih_eta_*` | The raw, unsmeared input templates, so the inputs of a run can always be re-inspected from the output alone. |

On success the very last line printed is `Done: <output>.root`. That sentinel is the reliable
success test: `TFile::Open(…, "RECREATE")` creates the file on disk *before* the macro has any
chance to fail, so the presence of an output file proves nothing. The launcher currently checks the
ROOT exit code plus the appearance of fresh `.root` files; grepping the log for `Done:` is the
stronger check when a run looks suspicious.

---

## 7. Step 3a — Background systematics

```bash
python3 systBckg.py
python3 systBckg.py --etas Eta1,Eta1_2p4,Eta2p4
python3 systBckg.py --region 9fp10
python3 systBckg.py --only Eta,Ih          # subset of sources
python3 systBckg.py --nominal-only         # statistical term only
python3 systBckg.py --no-plots
```

Set `BASE`, `DATASET`, `SAMPLETYPE`, `SUFFIX`, `CUTS`, `REGION` and `ERA` in the settings block at
the top; every one of them can also be overridden from the command line.

The script opens the nominal file **and every variation file** produced at step 2, rebins them onto
the analysis mass binning, folds under/overflow in, normalises each to unit area (so the comparison
is a pure **shape** comparison — the ABCD normalisation is common to all variations and cancels),
and computes for each source

```
syst = 100 × max( |1 − varUp/nominal| , |1 − varDown/nominal| )    [%]
```

both bin by bin (`*_binned`) and on the cumulative integral. One-sided sources are symmetrised. A
two-sided source missing one of its two files is **skipped entirely**, since keeping the surviving
side would silently halve the envelope.

| Key | Variations | In the total |
| --- | --- | --- |
| `Stat` | RMS of the toys of the nominal | yes |
| `Eta` | `binEtaUp` / `binEtaDown` | yes |
| `Ih` | `binIhUp` / `binIhDown` | yes |
| `P` | `binMomUp` / `binMomDown` | yes |
| `FitIh` | `fitIhUp` / `fitIhDown` | yes |
| `FitP` | `fitMomUp` / `fitMomDown` | yes |
| `CorrIh` | `corrTemplateIh` (one-sided) | yes |
| `Corr1oP` | `corrTemplate1oP` (one-sided) | yes |

The keys of this table, `Stat` excepted, are the values accepted by `--only`. The total is the
quadratic sum of the sources flagged `inTotal` — all of them at present. Note that the statistical
term is included, so the key `systTotalBinned` is really a *total* uncertainty, not a purely
systematic one; keep that in mind when combining it downstream.

The `noFit` variation is still produced by the launcher (§6.3) but is no longer in the
`SYSTEMATICS` table of `systBckg.py`: it is neither plotted nor summed here. To look at it, plot
that variation directly with `python ShowPlots.py --label noFit`.

Outputs, under `<indir>/<systdir>/` (`systdir` defaults to `SystCombined`):

- `sysTotBinned_<eta>_<region>.root` — `mass_predBC_nominal`, `Stat`, `Stat_binned`, one `<key>` and
  one `<key>_binned` per source, `systTotal` and **`systTotalBinned`** (the key read by the plotting
  macro — do not rename it);
- `individualSyst/<eta>/plot_<key>.pdf` — nominal vs. up vs. down, with both ratio forms;
- `pdf/`, `root/`, `Cfile/summary_binned_syst_<stem>_<eta>.*` — all sources on one plot.

---

## 8. Step 3b — Signal systematics

```bash
python3 systSignal.py
python3 systSignal.py --etas Eta1,Eta2p4 --region 9fp10
python3 systSignal.py --masses 2000,2400,2600
python3 systSignal.py --raw                # original step 1 files, not rescaled to the luminosity
python3 systSignal.py --no-plots
```

Structurally parallel to `systBckg.py`, with three differences:

1. **Input layout.** All the variations of a mass point sit in the *same* step 1 file, distinguished
   by a suffix on the histogram name — they are not re-derived here:
   `METanalysis_TestPUppiMETCut<cuts>_<eta>_<region>_SignalMass_<variation>`.
2. **No normalisation.** The spectra keep their absolute yields, since the flat systematics are
   multiplicative on yields.
3. **No statistical term.** Only the systematic sources are summed.

The file of a mass point is `<IDIR>/Gluino_Run3_MET_madgraph_<mass>_V<vsignal>_weighted.root`, i.e.
the luminosity-normalised copy written by `RescaleSignal.py` (§5); `IDIR`, `VSIGNAL` and `WEIGHTED`
are set at the top of the script. `--raw` reads the original step 1 file instead (same name
without `_weighted`). The relative uncertainties are the same either way, since the weight is
common to the nominal and to its variations; only the stored `mass_nominal` changes scale.

| Source | Variation suffixes | Value |
| --- | --- | --- |
| `K` | `KUp` / `KDown` | from step 1 |
| `C` | `CUp` / `CDown` | from step 1 |
| `PU` | `PUUp` / `PUDown` | from step 1 |
| `Trigger` | `TriggerSFUp` / `TriggerSFDown` | from step 1 |
| `Jet` | `JetUp` / `JetDown` | from step 1 |
| `lumi` | — | flat **1.4 %** |
| `Fpix` | — | flat **1.6 %** |

Mass points: 1100 → 2600 GeV. Outputs, under `--odir` (default `systSignal/`):

- `Gluino_<mass>_<region><cuts>_<eta>/sysTotBinned_signal.root` — `mass_nominal`, one `syst_<key>`
  per source, and `systTotalBinned` (same key name as on the background side, so both can be
  consumed identically);
- `Gluino_<mass>_…/summary_Gluino_<mass>_<eta>.pdf` and one up/down plot per source;
- `sysTot_allMasses_<eta>.pdf` — total uncertainty for the highlighted points (2000, 2400, 2600 GeV).

---

## 9. Step 3c — Mass-spectrum plots

```bash
python ShowPlots.py
python ShowPlots.py --etas Eta1,Eta1_2p4,Eta2p4
python ShowPlots.py --label binEtaUp
python ShowPlots.py --syst
python ShowPlots.py --dry-run              # print the commands without running them
```

`ShowPlots.py` does no physics: it rebuilds the step 2 output path for the requested systematic
label, optionally locates the combined-systematics file, and shells out to `PlottingMacro.py` once
per η range. Set `DATASET`, `SAMPLETYPE`, `SUFFIX`, `CUTS`, `VSIGNAL`, `REGION`, `YEAR`, `ERA` and
`ISTTBAR` in the block at the top; all of them can be overridden from the command line.

Without `--syst` the band is the statistical uncertainty of the prediction alone. With `--syst` the
combined file written by `systBckg.py` is required; if it is missing the region is **skipped**
rather than silently falling back to the nominal-only plot, which would look identical but mean
something different.

**MC mode.** `ShowPlots.py` passes `--isMC True` as soon as `SAMPLETYPE` starts with `mc` or
`ttbar`. The "observation" is then the sum of the simulated processes, read from the
luminosity-normalised files of §5, which lets the closure of the method be tested where the truth
is known.

**ttbar-only closure test.** Run step 2 on the dileptonic ttbar sample with `sampleType = ttbar2024`,
then plot with `SAMPLETYPE = "ttbar2024"`, `DATASET` pointing at the same ttbar file, **and**
`--ttbar` (or `ISTTBAR = True`). The two switches are independent: the sample type alone does not
restrict the observation to ttbar. In this mode the prediction is **normalised to the integral of
the observation** before being drawn. As committed, step 2 reads the raw step 1 file
(`TTbar2024_V15p10`) whereas the observation comes from its `_weighted` copy, so the two are not
on the same scale to begin with; the comparison is therefore shape-only and says nothing about the
ABCD normalisation.

> **Careful.** That normalisation is applied whenever `--isTTbar True` reaches the macro, whatever
> the sample: `--ttbar` on a data or full-MC working directory silently scales the prediction to
> the observation.

| `PlottingMacro.py` option | Meaning |
| --- | --- |
| `--ifile` | Step 2 output file. |
| `--cuts` | step 1 selection fragment, used to rebuild the histogram names. |
| `--ofile` | Prefix of the output files. |
| `--region` | `8fp9` (VR) or `9fp10` (SR). The SR is blinded above 300 GeV in the *drawn* histograms. |
| `--odir` | Output directory (created if missing). |
| `--nom` | `True` = statistical band only; `False` = read `--systfile`. |
| `--systfile` | `sysTotBinned_<eta>_<region>.root`, key `systTotalBinned`. |
| `--eta` | η range; drives the label drawn on the plot. |
| `--isMC` | `True` = the "observation" is the stack of the simulated processes. |
| `--isTTbar` | MC only: restrict the "observation" to the dileptonic ttbar, and normalise the prediction to it. |
| `--vsignal`, `--year`, `--era` | Signal version, year, and era (`""` = whole year, or `F` / `G`). |

The figure has four pads: the spectra (log y, with the prediction band, the observed points, the
optional MC stack and the gluino overlays), an optional cumulative-ratio pad, the bin-by-bin ratio
`obs/pred`, and the pull `(N_obs − N_pred)/σ`. Three formats are written —
`.pdf`, `.root` and `.C` — under `<indir>/<eta>/Plots_<region>/`, with every switch encoded in the
file name so two configurations cannot overwrite each other.

> **Note.** `PlottingMacro.py` still contains hard-coded paths to the gluino samples (`Gluino_V19`)
> and to the MC samples used in the stack (W+jets, ttbar dileptonic and semileptonic, QCD). All of
> them are `_weighted.root` files produced at step 1b (§5). Update the paths when the sample
> versions change.

---

## 10. Expected input histograms (step 1)

`loadHistograms()` reads, for each region name `region{A,B,C,D}_<FpixelRange><Ext>`:

| Name | Axes |
| --- | --- |
| `eta_1oP_<region>` | X = `10⁴/p` [GeV⁻¹], Y = `η` |
| `ih_eta_<region>` | X = `η`, Y = `Ih` [MeV/cm] |
| `mass_<region>` | mass |
| `mass_eta_<region>` | X = mass, Y = `η` |

with

```
Ext = "_METanalysis_TestPUppiMETCut" [+ "_SigmaPtoverPt_<x>"] [+ "_EoP_<y>"] + "_" + <etaRange>
```

so a full name reads e.g. `eta_1oP_regionA_3fp8_METanalysis_TestPUppiMETCut_Eta2p4`. The leading
`_METanalysis_TestPUppiMETCut` is the step 1 selection tag and must be updated here (and in
`PlottingMacro.py` / `systSignal.py`) whenever step 1 changes its naming.

`<FpixelRange>` is `3fp8` (0.3 < Fpixel ≤ 0.8), `3fp9`, `8fp9` or `9fp10`. The validation region
needs `regionA_3fp8`, `regionB_8fp9`, `regionC_3fp8`, `regionD_8fp9`; the search region needs
`regionA_3fp9`, `regionB_9fp10`, `regionC_3fp9`, `regionD_9fp10`.

`TH1F`/`TH2F` inputs are converted explicitly to `TH1D`/`TH2D` on read-out (`GetAsTH1D` /
`GetAsTH2D`): a direct cast is undefined behaviour and produces garbage silently.

The η rebinning uses variable-width edge tables selected from the region name. **The order of the
substring tests matters** in both `PickEtaBinning()` and the correction functions `corrIh` /
`corr1oP`: `Eta1_2p4` and `Eta1p2_2p2` also contain `Eta1`, so the most specific names are tested
first and reordering the blocks silently makes branches unreachable.

---

## 11. Path conventions

The layout is built independently in three places and nothing enforces that they agree — a mismatch
shows up as a "file not found" list rather than as a wrong plot. The three must produce the same
string:

| Script | Directory | Notes |
| --- | --- | --- |
| `LaunchBkgPred.py` | `<outputDir>/<sampleType>_V<version>__<region>_<labelDir>/<etaRange>/` | `labelDir` receives the `SigmaPtoverPt_…_EoP_…_` prefix automatically when both cuts are set. |
| `ShowPlots.py` | `<BASE>/<SAMPLETYPE>_V<version>__<REGION><CUTS>_<SUFFIX>/` | `CUTS` is inserted for you; `SUFFIX` is the bare suffix (`v2`). |
| `systBckg.py` | `<BASE>/<SAMPLETYPE>_V<version>__<REGION>_<SUFFIX>/` | `CUTS` is **not** inserted here: `SUFFIX` must already include it (`SigmaPtoverPt_0p5_EoP_0p1_v2`). |

File names inside that directory:

```
<eta>/<dataset stem>_<eta><CUTS>_<label>.root                 prediction (step 2)
<systdir>/sysTotBinned_<eta>_<region>.root                    combined systematics (step 3a)
<eta>/Plots_<region>/<ofile>_region<region>_<year>…_<eta>.pdf final plot (step 3c)
```

The analysis mass binning (36 edges, from 0 to 4000 GeV, widening towards the tail) is duplicated in
`rebinHisto()` (`CommonFunctions.h`), `REBINNING` (`systBckg.py`, `systSignal.py`) and `rebinning`
(`PlottingMacro.py`). **The four must stay identical**, otherwise the uncertainty band does not line
up with the points.

---

## 12. Practical notes

- **Do not combine `TProcessExecutor` with `EnableImplicitMT`** — the CPU is over-subscribed
  (`nWorkers × nThreads`). Likewise, shell-level stdout redirection (`> log 2>&1`) breaks the forked
  workers; capture the output through Python instead if you need a log file.
- `nWorkers` and `nPE` drive both the runtime and the memory usage. Reduce them on a shared machine.
- A crashed worker silently drops its result, which would bias the mean and under-estimate the
  spread; `bckgEstimate()` checks the count explicitly and fails the whole configuration rather than
  writing a partial result.
- `saveFits = 1` writes one file per toy in `DebugFit/`; with `--all` that is `nPE × 14` files.
  Useful for one variation, unmanageable for the full sweep.
- `systErr_` in `Regions.h` is set to `0` on purpose, so that the ratio stored by `ratioIntegral()`
  carries the statistical uncertainty only while the systematics are being measured. Setting it
  non-zero would double-count the spread between variations.
- The region directory name is deduced from `runVR` / `runSR` by an `elif` chain: with both flags at
  1 the directory is labelled `8fp9` even though the file also contains the SR histograms.
- `Region` copies are **shallow** (raw pointers, implicit copy constructor). This is used
  deliberately in `bckgEstimate()`, where a local copy has its template pointers reassigned to
  per-toy histograms and cleared afterwards. Never let two `Region`s delete the same pointer.
- The convolution loops include the **overflow** bin on purpose (dropping it would lose the
  highest-mass entries) and exclude the underflow (unphysical values).