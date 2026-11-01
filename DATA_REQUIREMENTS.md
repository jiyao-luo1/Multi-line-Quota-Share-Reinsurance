# Private empirical data requirements

The European motor-insurance dataset is proprietary and subject to a non-disclosure agreement. No individual records or original private data files are distributed in this repository. The simulation scripts generate their own samples and do not require these empirical inputs.

## Input schemas

The empirical scripts use two distinct saved-object schemas. The files cannot be substituted for one another without preparing the required object locally under authorized access.

| Scripts | Required saved object | Required loss columns, in order | Input location |
| --- | --- | --- | --- |
| `code/real_eur1_R1.R`, `code/real_eur2_R1.R` | `final_df`, with named loss columns | `"1"`, `"2"`, `"3"` | `data_eur.Rda` beside the two case scripts, normally `code/data_eur.Rda`. |
| `code/Data_Dependence_Rev8.R` | `data_eur`, a data frame | `CLAIM_TPL`, `CLAIM_CAR`, `CLAIM_UMC` | First positional argument; otherwise `Code/data_eur_new.Rda`, relative to the working directory. |

Each row should contain one complete policy-level loss vector. The columns represent third-party liability, car damage, and uninsured motorist losses in that order. Supply numeric, finite, nonnegative losses, retaining zero-loss policies. Additional columns are ignored. The case-analysis workflow requires enough nondegenerate observations for splitting and tail estimation. The dependence script explicitly requires at least two rows; constant columns or empty exceedance sets can produce undefined descriptive statistics.

## Case analyses

Place the authorized `data_eur.Rda` file in the same directory as both case scripts. Keep those scripts together because Case 2 sources the Case 1 functions. From the repository root:

```sh
Rscript code/real_eur1_R1.R
Rscript code/real_eur2_R1.R
```

Case 1 saves `real_eur1_R1.Rda`, containing `sim.list.1`; Case 2 saves `real_eur2_R1.Rda`, containing `sim.list.2`. Both outputs are written beside the scripts, regardless of the launch working directory. The case analyses require `nleqslv` and use the built-in `parallel` package for bootstrap work.

The supported controls are `REINSURANCE_B=1000`, `REINSURANCE_NCORE=14`, and `REINSURANCE_NSIM=51` by default. The grid size must be odd, and the worker count is capped by available processors. These supplied files do not support data-file or output-directory environment overrides.

The loader shuffles complete rows with seed `9` before the split analysis. The bootstrap seeds are `91000` and `92000` for Cases 1 and 2. Preserve the authorized input extract, its row ordering, and the supplied shuffle procedure when comparing runs. Code availability alone does not establish exact reproduction of the manuscript tables from a different data extract or software environment.

## Descriptive dependence summaries

Supply explicit local paths to avoid relying on the original default directory names:

```sh
Rscript code/Data_Dependence_Rev8.R /absolute/local/path/data_eur_new.Rda /absolute/local/path/dependence_results
```

The first argument selects the private input; the second selects the output directory. If omitted, the defaults are `Code/data_eur_new.Rda` and `Empirical_Rev8`, relative to the working directory. The capitalized `Code/` default is distinct from this repository's `code/` directory. The script does not read data-path or output-path environment variables.

The dependence script uses base R and writes `dependence_summary.csv`, `tail_exceedance_counts.csv`, `dependence_table.tex`, and `computation_notes.txt`. These contain aggregate summaries, with no individual records exported. Correlations use all supplied rows, including zero losses. Exceedance quantities at `q=0.95` and `0.99` are descriptive finite-threshold conditional probabilities, not asymptotic tail-dependence coefficients.

Keep all private inputs local. The repository's exclusions cover private data and generated result files; inspect any separate submission or archive before sharing it. The confidential-data analyses were not rerun as part of this update.
