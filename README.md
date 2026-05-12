# Two-Stage HMREM Model for Relational Event Data

This repository contains simulation, case-study, and reporting scripts for a
two-stage Hidden Markov Model - Relational Event Model (HMMREM/HMREM) workflow.
The code is designed to evaluate how well latent behavioral states can be
recovered from relational event sequences and how those state assignments affect
second-stage REM estimates.

## Overview

The project combines two modeling steps:

1. **Hidden Markov Model (HMM)**: infers latent states from inter-event times.
2. **Relational Event Model (REM)**: estimates event-generation mechanisms using
   observed event covariates and, where relevant, inferred state information.

The current analysis includes a scenario simulation study, Apollo 13 case-study
scripts, uncertainty diagnostics, overlap diagnostics, and Excel/plot export
utilities.

## Current Analysis Design

The main simulation script is `HMMREM_scenarios.R`. It now runs a systematic
scenario grid over:

- Scenarios: `Easy`, `Medium`, `Hard`, and `ExtremeHard`
- Event sequence lengths: `T = 3000` and `T = 6000`
- Actor counts: `N = 10`, `15`, and `20`
- Replications: `R = 100`
- HMM state counts: `1:4` when `RUN_HMM_REM <- TRUE`

The 2-state HMM is always fitted because it is used for the core recovery,
confidence, and overlap diagnostics. Extended HMM model comparison can be
controlled with:

```r
RUN_HMM_REM <- TRUE
```

## Main Files

| File | Purpose |
| --- | --- |
| `HMMREM_scenarios.R` | Main simulation pipeline for scenario generation, HMM fitting, state recovery, REM fitting, diagnostics, plots, and saved result bundles. |
| `export_results_to_excel.R` | Loads `hmmrem_analysis_bundle.rds` or runs the analysis if needed, then exports summary tables to `HMMREM_results.xlsx`. |
| `create_plots_and_excel_from_saved_results.R` | Rebuilds plots and Excel output from saved `.rds` results without rerunning the full simulation. |
| `Apollo2StateCase_analysis.R` | Apollo 13 analysis with HMM model selection, confidence diagnostics, REM/HMMREM comparison, plots, and Excel export. |
| `run_hmmrem_pipeline_server.R` | Server-oriented runner that checks packages and runs the full scenario pipeline. |
| `2stateHMREM.Rmd` | Earlier 2-state HMREM simulation and analysis notebook. |
| `3stateHMREM.Rmd` | Earlier 3-state HMREM simulation and analysis notebook. |
| `Apollo2StateCase.Rmd` | Apollo 13 case-study notebook. |
| `predictive_performance.Rmd` | In-sample and out-of-sample REM predictive performance helpers. |
| `HMMREM_analysis.md` | Technical documentation for the scenario-analysis methodology. |

## Diagnostics Added in Recent Updates

The current scripts include several additions beyond the original README:

- State-assignment confidence from HMM posterior probabilities.
- Uncertainty rates using `tau1 = 0.9` and `tau2 = 0.8`.
- State overlap diagnostics for the 2-state HMM:
  - histogram overlap (`OVL_hist`)
  - kernel-density overlap (`OVL_kde`)
  - Wasserstein distance (`W1`)
  - parametric overlap / divergence summaries where stable
- REM uncertainty calibration comparing empirical coefficient SDs with reported
  REM standard errors.
- Saved density data and regenerated density plots for each configuration.
- Saved analysis bundle support via `hmmrem_analysis_bundle.rds`.

## Outputs

Typical generated outputs include:

| Output | Description |
| --- | --- |
| `hmmrem_analysis_bundle.rds` | Saved R objects needed to recreate summaries, plots, and Excel output. |
| `HMMREM_results.xlsx` | Excel workbook containing summary, scenario parameters, diagnostics, HMM results, confidence metrics, state overlap, and REM uncertainty sheets. |
| `HMMREM_scenario_comparison.png` | BIC/model-selection comparison plot. |
| `HMMREM_confidence_comparison.png` | State-assignment confidence and uncertainty comparison plot. |
| `HMMREM_overlap_comparison.png` | State-overlap diagnostic plot. |
| `density_data/*.rds` | Per-configuration density data used for state-separation plots. |
| `plots/density/*.png` | Per-configuration inter-event-time density plots by inferred state. |
| `Apollo_HMMREM_Results.xlsx` | Apollo 13 analysis workbook. |

## Running the Simulation Pipeline

Run the full scenario analysis from the repository root:

```r
source("HMMREM_scenarios.R")
```

Or from the command line:

```bash
Rscript HMMREM_scenarios.R
```

To export the latest available results to Excel:

```r
source("export_results_to_excel.R")
```

To rebuild plots and Excel files from saved results without rerunning the full
simulation:

```bash
Rscript create_plots_and_excel_from_saved_results.R
```

You can also pass a specific saved bundle:

```bash
Rscript create_plots_and_excel_from_saved_results.R path/to/hmmrem_analysis_bundle.rds
```

## Dependencies

The scripts use the following R packages:

- `dplyr`
- `ggplot2`
- `plotly`
- `momentuHMM`
- `remulate`
- `tidyr`
- `knitr`
- `remify`
- `remstats`
- `remstimate`
- `openxlsx`
- `readxl`
- `readr`
- `zoo`
- `fastDummies`
- `abind`

Install the common dependencies with:

```r
install.packages(c(
  "dplyr", "ggplot2", "plotly", "momentuHMM", "remulate", "tidyr",
  "knitr", "openxlsx", "readxl", "readr", "zoo", "fastDummies", "abind"
))
```

The REM packages `remify`, `remstats`, `remstimate`, and `remulate` may require
installation from their current R package sources if they are not available from
your default CRAN setup.

## Reproducibility Notes

- The scenario pipeline sets `set.seed(1234)`.
- `hmmrem_analysis_bundle.rds` is the preferred handoff file for rebuilding
  reports and plots.
- If the bundle is missing, `export_results_to_excel.R` will attempt to run
  `HMMREM_scenarios.R` before exporting.
- Full simulation runs can take time because the design covers multiple
  scenarios, actor counts, sequence lengths, replications, and HMM state counts.

## Project Status

The repository currently tracks both `origin/master` and `origin/main`.
Recent updates focused on the scenario-analysis scripts, Excel export, REM
uncertainty calculations, state-overlap diagnostics, and regenerated plot formats.
