# Private empirical data requirements

The empirical dataset is proprietary and subject to a non-disclosure agreement. It is not distributed in this repository. Keep all private data locally and do not commit or upload individual records, identifiers, or real `.Rda` files. The simulation scripts do not require these private empirical data.

The supplied empirical scripts use two distinct local input schemas. The object names and column names below are requirements of the corresponding scripts; these files are not interchangeable.

| Scripts | Default local file | Required saved object | Required loss columns, in order |
| --- | --- | --- | --- |
| `code/real_eur1_R1.R`, `code/real_eur2_R1.R` | `private/data_eur.Rda` | `final_df`, a data frame | `"1"`, `"2"`, `"3"` |
| `code/Data_Dependence_Rev8.R` | `private/data_eur_new.Rda` | `data_eur`, a data frame | `CLAIM_TPL`, `CLAIM_CAR`, `CLAIM_UMC` |

Each row represents one complete policy-level loss vector. The three columns represent third-party liability, car damage, and uninsured motorist losses, respectively. The columns must be numeric, finite, and nonnegative, retaining policies with zero losses. Additional columns are ignored. The case-analysis loader requires at least four rows for the split calculation; meaningful estimation also requires nondegenerate losses and sufficient tail observations. The dependence loader requires at least two rows. Constant columns or empty exceedance sets can produce undefined descriptive statistics.

The supplied case-analysis scripts preserve the original row shuffle with seed `9` before splitting the observations. They also preserve the original estimation, smoothing, bootstrap procedures, pricing inputs, and default settings. Retain the supplied row ordering and shuffle procedure when comparing an authorized run with the original workflow. Public code availability does not imply that a different extract, row ordering, or software environment will reproduce the paper's numerical tables exactly.

## Running the empirical analyses

Run these commands from the repository root after placing the authorized files at their default local locations:

```sh
Rscript code/real_eur1_R1.R
Rscript code/real_eur2_R1.R
Rscript code/Data_Dependence_Rev8.R
```

The case-analysis scripts require the R package `nleqslv`; the dependence script uses base R. Case 2 sources the shared implementation in `code/real_eur1_R1.R`, so keep both scripts together.

The case-analysis defaults are `REINSURANCE_B=1000`, `REINSURANCE_NCORE=14`, and `REINSURANCE_NSIM=51`. Available processor counts cap the requested number of workers. The two case analyses use their original bootstrap seeds, `91000` and `92000`. `REINSURANCE_NSIM` must be odd. These production runs can be computationally expensive.

To use a different authorized local file or output location, set `REINSURANCE_DATA_FILE` and `REINSURANCE_OUTPUT_DIR` for the corresponding command. For example:

```sh
REINSURANCE_DATA_FILE=/absolute/local/path/data_eur.Rda \
REINSURANCE_OUTPUT_DIR=/absolute/local/path/empirical-results \
Rscript code/real_eur1_R1.R
```

Default paths are resolved relative to the repository, independently of the working directory used to launch the script. Explicit environment-variable or positional paths follow normal R path resolution; absolute paths avoid ambiguity. The default case-analysis results are saved as `outputs/empirical/real_eur1_R1.Rda` and `outputs/empirical/real_eur2_R1.Rda`. These are analysis outputs, not copies of the input records.

The dependence script also retains its optional positional arguments, which take precedence over the environment variables:

```sh
Rscript code/Data_Dependence_Rev8.R \
  /absolute/local/path/data_eur_new.Rda \
  /absolute/local/path/dependence-results
```

Its default output directory is `outputs/dependence/`. It writes aggregate Pearson and Spearman correlations, finite-threshold conditional exceedance summaries, a LaTeX table, and computation notes. Its exceedance quantities at `q = 0.95` and `0.99` are descriptive finite-threshold conditional probabilities, not asymptotic tail-dependence coefficients. No individual records are exported.

The repository does not provide the private input files or assert exact reproduction of all paper tables. An authorized user must supply the corresponding local data extract in the required schema.
