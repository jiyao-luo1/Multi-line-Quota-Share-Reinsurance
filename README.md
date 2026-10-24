# Multi-line Quota-Share Reinsurance

R code accompanying **Diversification Benefits of Multi-line Quota-Share Reinsurance: Data-Driven Decision**, by Yichun Chi, Tsz Chai Fung, Jiyao Luo, and Liang Peng.

The repository contains the supplied implementations for the simulation studies, normal Q-Q diagnostics, empirical quota-share analysis, and descriptive cross-coverage dependence summaries. The code is organized for local execution. The European motor-insurance dataset is proprietary and subject to a non-disclosure agreement; it is not distributed here.

## Requirements

R and the packages `nleqslv` and `QRM` (with their automatically installed dependencies) are required for the simulation workflows. The empirical optimization scripts require `nleqslv`. The descriptive dependence script uses base R only. Parallel workflows use R's built-in `parallel` package.

From the repository root, install missing packages with:

```sh
Rscript code/install_dependencies.R
```

The source was checked locally with R 4.5.2. Record `sessionInfo()` when running a study, since package versions and numerical solvers can affect results.

## Code map

| Script | Purpose |
| --- | --- |
| `code/Simu_True.R` | Compute large-sample reference estimates and simulation diagnostic plots. |
| `code/Simu.R` | Run Monte Carlo estimates and random-weight bootstrap draws for one simulation configuration. |
| `code/Summary_Simu.R` | Summarize coverage from matching Monte Carlo output and reference values. |
| `code/Simu_QQ_Rev5.R` | Generate or redraw normal Q-Q plots and upper-tail diagnostics. |
| `code/real_eur1_R1.R` | Empirical Case 1: estimate ceded proportions and diversification measures over sensitivity settings. |
| `code/real_eur2_R1.R` | Empirical Case 2, using the Case 1 numerical functions. |
| `code/Data_Dependence_Rev8.R` | Export aggregate correlations and finite-threshold conditional exceedance summaries. |

Generated files belong under `outputs/`, which is excluded from version control. Keep any authorized confidential inputs in `private/` or elsewhere outside the repository. [Data requirements](DATA_REQUIREMENTS.md) specify the input schemas.

## Simulation studies

Run commands from the repository root. The simulation workflow is:

1. Generate reference values with `Simu_True.R`.
2. Generate Monte Carlo and bootstrap draws with `Simu.R` for the same configuration.
3. Run `Summary_Simu.R` with matching reference and simulation settings.

Full settings are computationally intensive. Use the smaller examples below to check the workflow before running the full study. These examples check execution; they are not estimates of the paper's numerical results.

```sh
REINSURANCE_CASE=1 REINSURANCE_SIGMA=1 REINSURANCE_SCALE=1 REINSURANCE_SHAPE=3.5 \
  REINSURANCE_TRUE_N=100 REINSURANCE_TRUE_REPS=3 Rscript code/Simu_True.R

REINSURANCE_CASE=1 REINSURANCE_SIGMA=1 REINSURANCE_SCALE=1 REINSURANCE_SHAPE=3.5 \
  REINSURANCE_N=100 REINSURANCE_REPS=2 REINSURANCE_BOOT=40 Rscript code/Simu.R

REINSURANCE_CASE=1 REINSURANCE_SIGMA=1 REINSURANCE_SCALE=1 REINSURANCE_SHAPE=3.5 \
  REINSURANCE_N=100 REINSURANCE_BOOT=40 Rscript code/Summary_Simu.R
```

For a full run, use the settings in [Portability notes](PORTABILITY_NOTES.md), with matching case, scale, shape, and sample size across the three steps. `Simu_True.R` defaults to 100,000 observations and 10,000 replications for each of 16 settings; select one setting with `REINSURANCE_JOB_ID` or explicit scenario variables. `Simu.R` defaults to 1,000 observations, 1,000 Monte Carlo replications, and 2,000 bootstrap draws. Set `REINSURANCE_N=3000` for the second sample size. The summary must use the same `REINSURANCE_BOOT` value as its simulation inputs. Use `REINSURANCE_OUT` to choose separate output directories when comparing runs, and set `REINSURANCE_SIM_DIR` and `REINSURANCE_TRUE_DIR` when summarizing nondefault locations.

For normal Q-Q diagnostics, a small example is:

```sh
Rscript code/Simu_QQ_Rev5.R --n=100 --reps=20 --cases=1 --shapes=3.5 --scale=1 --workers=1 --out=outputs/qq_smoke
```

The Q-Q script defaults to `n=3000`, `reps=2000`, cases 1-4, both shape parameters 3.5 and 6, scale configuration 1, and four workers. Specify `--scale=2` in a separate output directory for the unequal-scale configuration. Add `--plot-only` to redraw complete saved replicates with the same configuration. Choose a fresh output directory when changing any configuration.

## Empirical analysis

Authorized users can run the two empirical analyses after providing the appropriate private input:

```sh
REINSURANCE_DATA_FILE=/absolute/path/to/data_eur.Rda Rscript code/real_eur1_R1.R
REINSURANCE_DATA_FILE=/absolute/path/to/data_eur.Rda Rscript code/real_eur2_R1.R
```

The defaults are 1,000 bootstrap draws, up to 14 workers, and 51 sensitivity-grid points. The environment variables `REINSURANCE_B`, `REINSURANCE_NCORE`, `REINSURANCE_NSIM`, and `REINSURANCE_OUTPUT_DIR` control these settings. Each case writes its analysis object to `outputs/empirical/` by default. The supplied scripts compute estimates and intervals; separate scripts for recreating the published empirical figure layouts were not included in the source materials.

For descriptive dependence summaries, the source data use the second schema documented in [Data requirements](DATA_REQUIREMENTS.md):

```sh
Rscript code/Data_Dependence_Rev8.R /absolute/path/to/data_eur_new.Rda outputs/dependence
```

Only aggregate summaries are written. The conditional exceedance probabilities are finite-threshold descriptive measures, rather than asymptotic tail-dependence coefficients.

## Implementation and verification scope

The public copies retain the supplied estimation formulas, random-generation formulas, empirical row-shuffling seed, and numerical defaults. Changes to paths and run configuration are described in [Portability notes](PORTABILITY_NOTES.md). The simulation generator and exposure-target expressions are retained explicitly so the supplied numerical implementation can be inspected; they should not be treated as a separately corrected implementation of the theoretical model.

The local checks cover R parsing and reduced synthetic execution. The full Monte Carlo/bootstrap study and the proprietary-data results have not been rerun for this release. Historical empirical split/hybrid estimates and the full set of manuscript tables and figure layouts have not been independently verified by these checks.

## Reference

Chi, Y., Fung, T. C., Luo, J., and Peng, L. *Diversification Benefits of Multi-line Quota-Share Reinsurance: Data-Driven Decision*. Manuscript.
