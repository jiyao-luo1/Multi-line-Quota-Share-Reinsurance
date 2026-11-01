# Multi-line Quota-Share Reinsurance

R code accompanying **Diversification Benefits of Multi-line Quota-Share Reinsurance: Data-Driven Decision**, by Yichun Chi, Tsz Chai Fung, Jiyao Luo, and Liang Peng.

This repository provides six scientific R scripts supplied by the authors for simulation studies, empirical quota-share analysis, and descriptive cross-coverage dependence summaries. The proprietary European motor-insurance dataset is subject to a non-disclosure agreement and is not distributed here.

## Requirements

Install R and the packages imported by the simulation scripts: `snowfall`, `rlecuyer`, `copula`, `nleqslv`, `mvtnorm`, and `QRM`. From the repository root:

```sh
Rscript code/install_dependencies.R
```

The empirical case analyses require `nleqslv` and use R's built-in `parallel` package. The dependence script uses base R. Save `sessionInfo()` with local computations to record the software environment.

## Code map

| Script | Purpose |
| --- | --- |
| `code/Simu_True.R` | Large-sample reference estimates, density plots, and built-in normal Q-Q plots for 16 configurations. |
| `code/Simu.R` | Monte Carlo estimates and random-weight bootstrap draws for one selected configuration. |
| `code/Summary_Simu.R` | Coverage deviations from nominal levels using matching simulation and reference files. |
| `code/real_eur1_R1.R` | Empirical Case 1 estimates, bootstrap draws, and sensitivity analysis. |
| `code/real_eur2_R1.R` | Empirical Case 2, using the Case 1 functions from the same directory. |
| `code/Data_Dependence_Rev8.R` | Aggregate correlations and finite-threshold conditional exceedance summaries. |

Normal Q-Q plots are generated within `Simu_True.R`; a separate Q-Q script is no longer required.

## Simulation workflow

The simulation scripts use settings written in the source files. They do not support the previous simulation `REINSURANCE_*` environment variables. Make any run-specific changes in a local working copy.

1. **Generate reference values.** `Simu_True.R` runs all 16 configurations by default: Cases 1-4, two scale vectors, and shapes 3.5 and 6, with `Sigma=1`. Each configuration uses `n=100000`, `N.sim=10000`, and `Boot=0`. It writes reference CSVs, density PDFs, and three normal Q-Q PDFs per configuration to the current working directory. These defaults are computationally intensive.

   ```sh
   Rscript code/Simu_True.R
   ```

2. **Choose a simulation configuration before running.** All scenario assignments near the end of `Simu.R` are commented out. Uncomment exactly one line, for example:

   ```r
   Case=1;Sigma=1;Scale=1;p.shape=3.5;
   ```

   The defaults are `n=1000`, `N.sim=1000`, and `Boot=2000`. The supplied execution loop runs one batch, with `job_id=1`. After selecting the configuration:

   ```sh
   Rscript code/Simu.R
   ```

   The result is saved as `Simu_n..._Case..._Sigma..._Scale..._Shape..._1.Rdata` in the current working directory.

3. **Configure the coverage summary.** In `Summary_Simu.R`, replace the two author-specific locations: the directory in `setwd(...)` and the reference-CSV prefix in `Myname.TV`. Set `myn`, `mycase`, `mysigma`, `myscale`, and `myshape` to match the available results. Its supplied selection is `n=3000`, Case 4, Sigma 1, Scale 2, shape 6; this differs from the simulation sample-size default of `n=1000`. Its `myboot=2000` must also match the generated draws.

   ```sh
   Rscript code/Summary_Simu.R
   ```

   The console reports six coverage-minus-nominal deviations: the three estimators at 90%, followed by the three estimators at 95%. Negative values indicate coverage below the nominal level.

Both simulation generators use the uniform output of `rcopula.t` directly in the Pareto transformation and set the population joint ceded proportion to `0.45`. See [Execution notes](PORTABILITY_NOTES.md) for the exact target formula and file conventions.

## Empirical analysis

Authorized users should place the private file `data_eur.Rda` beside both empirical case scripts, normally at `code/data_eur.Rda`. Keep the two scripts together, then run:

```sh
Rscript code/real_eur1_R1.R
Rscript code/real_eur2_R1.R
```

Each case saves its `.Rda` result in the same directory as the scripts. The supported run controls are `REINSURANCE_B` (default `1000`), `REINSURANCE_NCORE` (`14`), and `REINSURANCE_NSIM` (`51`, odd). The earlier data-file and output-directory environment overrides are not supported by these supplied files.

For dependence summaries, provide the authorized input and output paths explicitly:

```sh
Rscript code/Data_Dependence_Rev8.R /absolute/local/path/data_eur_new.Rda /absolute/local/path/dependence_results
```

Without arguments, this script uses `Code/data_eur_new.Rda` and `Empirical_Rev8`, relative to the current working directory. The capitalized default `Code/` differs from this repository's `code/` directory. [Data requirements](DATA_REQUIREMENTS.md) describe the two input schemas. No private input data are included.

## Verification scope

The six scientific R files in this update are supplied unchanged. Release checks cover R syntax parsing and file-integrity verification. The new full reference calculations, Monte Carlo/bootstrap experiments, and proprietary-data analyses have not been rerun for this update; these checks do not establish exact reproduction of every manuscript table or figure.

## Reference

Chi, Y., Fung, T. C., Luo, J., and Peng, L. *Diversification Benefits of Multi-line Quota-Share Reinsurance: Data-Driven Decision*. Manuscript.
