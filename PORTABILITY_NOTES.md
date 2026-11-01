# Execution notes and verification scope

The six scientific R files are published exactly as supplied by the authors. Configure settings and local paths in a working copy before running a study. The dependency installer and repository documentation are separate from those scientific files.

## Simulation design

Both `Simu.R` and `Simu_True.R` use `d=3`, `alpha=0.95`, and a Student-t copula with `df.nu=3`. The supplied `Sigma=1` matrix has off-diagonal entries `0.25`. The two scale vectors are `(1000, 1000, 1000)` and `(800, 1000, 1200)`; the shape values are `3.5` and `6`.

The generator uses the copula output directly:

```r
U = rcopula.t(n, df=df.nu, p.sigma)
U_unif = U
X[,s] = p.scale[s] * ((1-U_unif[,s])^(-1/p.shape)-1)
```

The exposure target is defined by:

```r
mu_true = p.scale/(p.shape-1)
mu0 = 0.45*tilde_beta_scalar*sum(mu_true)-sum(beta_vec*mu_true)
```

This sets the population joint ceded proportion to `0.45`. The estimated joint proportion uses sample means and need not equal `0.45` in a finite sample. The supplied smoothing exponent remains `-0.33`.

## Running and locating results

| Script | Supplied configuration and execution | Output or required local adjustment |
| --- | --- | --- |
| `Simu.R` | `n=1000`, `N.sim=1000`, `Boot=2000`; every scenario assignment is commented out. Select exactly one scenario before running. The active loop runs `job_id=1`. | Writes one `Simu_n..._Case..._Sigma..._Scale..._Shape..._1.Rdata` file to the current working directory. |
| `Simu_True.R` | `n=100000`, `N.sim=10000`, `Boot=0`; the active loop runs `job_id=1:16` to select all 16 configurations. | Writes `True_Case..._Sigma..._Scale..._Shape....csv`, a density PDF, and `_QQ1.pdf`, `_QQ2.pdf`, `_QQ3.pdf` files per configuration to the current working directory. |
| `Summary_Simu.R` | Selects `myn=3000`, `mycase=4`, `mysigma=1`, `myscale=2`, `myshape=6`, and `myboot=2000`. | Edit both author-specific paths and match the selected settings to existing reference and simulation files. Prints coverage-minus-nominal deviations to the console. |
| `real_eur1_R1.R` | Case 1; defaults: 1,000 bootstrap draws, up to 14 workers, 51 sensitivity-grid points. | Reads `data_eur.Rda` and saves `real_eur1_R1.Rda` in the script directory. |
| `real_eur2_R1.R` | Case 2; the same run-control defaults. Sources `real_eur1_R1.R` from the same directory. | Reads `data_eur.Rda` and saves `real_eur2_R1.Rda` in the script directory. |
| `Data_Dependence_Rev8.R` | Accepts an input file and output directory as positional arguments. | Defaults to `Code/data_eur_new.Rda` and `Empirical_Rev8`, relative to the working directory. |

The simulation scripts import `snowfall`, `rlecuyer`, `copula`, `nleqslv`, `mvtnorm`, and `QRM`. Their active outer loops are sequential; the supplied snowfall launch commands are commented out. Full runs can be computationally substantial. Normal Q-Q generation is built into `Simu_True.R`.

### Matching summary inputs

The two local locations in `Summary_Simu.R` are:

- `setwd("/Users/lpeng/2026/Quota_ASTIN_Rev2/Simu_Data")`, which selects the simulation-file directory;
- the `~/2026/Quota_ASTIN_Rev2/Simu_True/` prefix in `Myname.TV`, which selects the reference CSV.

Replace both locations for the local machine. The selected case, Sigma, scale, and shape must agree across the reference and simulation files. `myn` must match the simulation sample size: the supplied summary uses `3000`, while `Simu.R` defaults to `1000`. `myboot` must match `Boot`; the supplied summary uses the fixed order statistics for 2,000 bootstrap draws.

The summary combines files matching the selected filename pattern, removes rows with missing values, and uses at most 1,000 retained rows. It uses the first reference value in the selected CSV for the three estimators. The six printed deviations are empirical coverage minus `(0.90, 0.90, 0.90, 0.95, 0.95, 0.95)`. The script does not export a coverage CSV.

### Supported controls

The simulation and summary scripts do not read the previous simulation `REINSURANCE_*` environment variables; change their source settings and local paths when needed. Their output location follows the working directory, including the summary's explicit `setwd(...)`.

The empirical case scripts support `REINSURANCE_B`, `REINSURANCE_NCORE`, and `REINSURANCE_NSIM`; `REINSURANCE_NSIM` must be odd. `REINSURANCE_SKIP_MAIN=1` is available for sourcing functions without running the main analysis. These files do not implement `REINSURANCE_DATA_FILE` or `REINSURANCE_OUTPUT_DIR`. The dependence script uses positional paths rather than environment overrides.

## Verification scope

This update was checked by parsing the R sources and verifying file integrity. No new full reference, Monte Carlo/bootstrap, or confidential-data computations were performed for the update. The repository provides the supplied code and usage information without asserting exact reproduction of all numerical results or figure layouts.
