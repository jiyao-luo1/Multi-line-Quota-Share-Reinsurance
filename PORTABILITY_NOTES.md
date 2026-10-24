# Portability changes and reproducibility scope

The four public simulation scripts are copies of the supplied scripts in `Code/`.
Their estimating equations, nonlinear solvers, smoothing formulas, bootstrap
weights, and historical random generator have been preserved.

## Changes to execution

- `Simu.R` now has an active default scenario (Case 1, Sigma 1, Scale 1,
  shape 3.5); the original file had every scenario assignment commented out.
  Its paper-sized defaults remain `n=1000`, 1,000 replications, and 2,000
  bootstrap replicates. Use `REINSURANCE_N=3000` for the second sample size.
- Unused direct imports of `snowfall`, `rlecuyer`, `copula`, and `mvtnorm`
  were removed from `Simu.R` and `Simu_True.R`. Their required direct packages
  are `QRM` and `nleqslv`; R installs package dependencies automatically.
- Input settings are configurable through environment variables. Relative
  paths are resolved from the repository root; absolute paths are accepted.
  The default output directories are `outputs/simulation`, `outputs/true`,
  `outputs/summary`, and `outputs/qq`, respectively. `REINSURANCE_OUT`
  overrides the destination for each invocation.
- `Summary_Simu.R` uses `REINSURANCE_SIM_DIR` and `REINSURANCE_TRUE_DIR`
  instead of author-specific absolute paths, checks missing/incompatible
  outputs, and writes a coverage CSV and session information in addition
  to its original console summary. It retains the original interval order
  statistics, the original file aggregation order, and removal of rows with
  missing values. `REINSURANCE_BOOT` must match the simulation output and
  be a multiple of 40 for those interval order statistics.
- `Simu_QQ_Rev5.R` retains its existing command-line options and resumable
  checkpoints; explicit command-line options take precedence over environment
  defaults. Windows uses one worker because fork-based `mclapply` is not
  available there. Other platforms retain the original worker behavior.
- `Simu.R` remains sequential within a job. Its `REINSURANCE_JOB_ID`
  selects the original seed/batch. `Simu_True.R` runs all 16 original scenario
  branches by default; `REINSURANCE_JOB_ID=1` selects just one. Any explicit
  scenario override selects branch 1 unless a job ID is also supplied.
  The original fixed truth seed `102025` is retained across branches.

## Environment settings

| Script | Settings and defaults |
| --- | --- |
| `Simu.R` | `REINSURANCE_N=1000`, `REINSURANCE_REPS=1000`, `REINSURANCE_BOOT=2000`, `REINSURANCE_CASE=1`, `REINSURANCE_SIGMA=1`, `REINSURANCE_SCALE=1`, `REINSURANCE_SHAPE=3.5`, `REINSURANCE_JOB_ID=1` |
| `Simu_True.R` | `REINSURANCE_TRUE_N=100000`, `REINSURANCE_TRUE_REPS=10000`; optional `REINSURANCE_JOB_ID` (1--16) and the same case/sigma/scale/shape overrides |
| `Summary_Simu.R` | `REINSURANCE_N=3000`, `REINSURANCE_CASE=4`, `REINSURANCE_SIGMA=1`, `REINSURANCE_SCALE=2`, `REINSURANCE_SHAPE=6`, `REINSURANCE_BOOT=2000`, `REINSURANCE_SUMMARY_LIMIT=1000`, optional input-directory overrides above |
| `Simu_QQ_Rev5.R` | `REINSURANCE_N=3000`, `REINSURANCE_REPS=2000`, `REINSURANCE_SEED=20261021`, `REINSURANCE_CASES=1,2,3,4`, `REINSURANCE_SHAPES=3.5,6`, `REINSURANCE_SCALE=1`, `REINSURANCE_WORKERS=4`; singular case/shape and `REINSURANCE_NCORE` are also accepted as fallbacks |

Only small smoke runs were executed while preparing these public copies.
The full truth calculation and coverage experiment are computationally
substantial and were not rerun. The truth CSVs and original simulation
replicate files are not present in the supplied working folder; run the
truth and simulation scripts before requesting a coverage summary.

## Historical implementation details retained

All three simulation generators first call `QRM::rcopula.t` and then apply
`pt(..., df=3)` to its output. `QRM::rcopula.t` already returns uniform
copula variates, so the second transform changes the margins; the resulting
sample is not Pareto with the distribution stated in the manuscript. This
behavior was verified from the installed `QRM` function and is retained to
avoid silently changing the supplied implementation or historical results.

The supplied simulation scripts also set
`mu0 = 0.2 * 0.1 * sum(scales) / (shape - 1)` while keeping
`beta_vec = c(0.05, 0.05, 0.05)`. This differs from the manuscript's stated
profitability-target formula, which subtracts `beta_vec` times the population
means. The original simulation smoothing exponent `-0.33` is also retained,
rather than replacing it with exactly `-1/3`.

These implementation-to-model differences require scientific reconciliation
before claiming that a new full simulation run validates the stated model.
The portability edits do not resolve or conceal those differences and do
not establish exact reproduction of every historical table or figure.
