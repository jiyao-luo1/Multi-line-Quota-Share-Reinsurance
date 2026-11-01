# Multi-line Quota-Share Reinsurance

R code for **Diversification Benefits of Multi-line Quota-Share Reinsurance: Data-Driven Decision**, by Yichun Chi, Tsz Chai Fung, Jiyao Luo, and Liang Peng.

## Setup

Requires R and the packages `snowfall`, `rlecuyer`, `copula`, `nleqslv`, `mvtnorm`, and `QRM`. Install them from the repository root:

```sh
Rscript code/install_dependencies.R
```

## Files

| Script | Purpose |
| --- | --- |
| `Simu_True.R` | Reference values, density plots, and normal Q-Q plots. |
| `Simu.R` | Monte Carlo simulation and random-weight bootstrap. |
| `Summary_Simu.R` | Coverage deviations from nominal levels. |
| `real_eur1_R1.R`, `real_eur2_R1.R` | Empirical estimates and sensitivity analyses for Cases 1 and 2. |
| `Data_Dependence_Rev8.R` | Aggregate correlations and conditional exceedance summaries. |

All scripts are in `code/`.

## Simulations

The 16 configurations combine Cases 1–4, scale vectors (1000, 1000, 1000) and (800, 1000, 1200), and Pareto shapes 3.5 and 6. All use `Sigma=1`, a Student-t copula with 3 degrees of freedom and off-diagonal correlation 0.25, and population joint ceded proportion 0.45.

Configure the source files before running:

1. `Simu_True.R` runs all 16 configurations with `n=100000`, `N.sim=10000`, and `Boot=0`.
2. In `Simu.R`, uncomment exactly one scenario line, such as `Case=1;Sigma=1;Scale=1;p.shape=3.5;`. Defaults are `n=1000`, `N.sim=1000`, and `Boot=2000`.
3. In `Summary_Simu.R`, update both local paths (`setwd(...)` and `Myname.TV`) and match the scenario, `myn`, and `myboot` to the generated files. Its supplied selection is n=3000, Case 4, Sigma 1, Scale 2, shape 6, with 2000 bootstrap draws.

Run in this order from the repository root:

```sh
Rscript code/Simu_True.R
Rscript code/Simu.R
Rscript code/Summary_Simu.R
```

Reference CSVs, density/Q-Q PDFs, and simulation files are written to the current working directory. The summary prints coverage minus the nominal level for three estimators at 90% and 95%, using at most 1000 replications. Full simulation settings are computationally intensive.

## Empirical analysis

The proprietary European motor-insurance data are subject to a non-disclosure agreement and are not included. Supply authorized local inputs with numeric, finite, nonnegative losses, including zero-loss policies.

| Analysis | Input file and saved object | Required columns |
| --- | --- | --- |
| Case analyses | `code/data_eur.Rda`, containing `final_df` | `"1"`, `"2"`, `"3"` |
| Dependence summaries | A local Rda file containing data frame `data_eur` | `CLAIM_TPL`, `CLAIM_CAR`, `CLAIM_UMC` |

Columns represent third-party liability, car damage, and uninsured motorist losses, respectively. Keep both case scripts and their input in the same directory; Case 2 sources the Case 1 functions.

```sh
Rscript code/real_eur1_R1.R
Rscript code/real_eur2_R1.R
Rscript code/Data_Dependence_Rev8.R /absolute/path/data_eur_new.Rda /absolute/path/dependence_results
```

The case analyses save their results beside the scripts. Optional environment controls are `REINSURANCE_B=1000`, `REINSURANCE_NCORE=14`, and `REINSURANCE_NSIM=51` by default; the sensitivity-grid size must be odd.

For dependence summaries, pass the input and output paths explicitly as shown. Outputs contain aggregate correlations, finite-threshold conditional exceedance probabilities, a LaTeX table, and computation notes.
