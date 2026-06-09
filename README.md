# Q-methodology MCEV workflow

This repository contains `qmcev.r`, an R script that runs a project-specific Q-methodology analysis workflow using the [`qmethod`](https://cran.r-project.org/package=qmethod) R package.

The script is configured for a fixed MCEV Q-sort design:

- **49 statements**.
- **8 participant Q-sorts**, named `P1` through `P8`.
- A **forced distribution** from `-6` to `6` with counts `1, 2, 3, 4, 5, 6, 7, 6, 5, 4, 3, 2, 1`.
- Q-method factor solutions from **1 factor through 4 factors**.
- **PCA** extraction, **Pearson** correlations, and **varimax** rotation for multi-factor solutions.

## What `qmcev.r` does

`qmcev.r` is an end-to-end batch workflow for validating Q-sort data, running multiple Q-method factor solutions, and exporting analysis products for review.

At a high level, it:

1. Reads an input CSV named `Qsort_results.csv` from the repository root.
2. Validates that the CSV has the required metadata columns and expected participant columns.
3. Validates that the dataset contains exactly 49 statements.
4. Validates each participant's Q-sort against the expected forced distribution.
5. Converts participant columns into a numeric Q-sort matrix.
6. Writes validation reports and the prepared analysis matrix.
7. Computes a participant correlation matrix.
8. Runs `qmethod` solutions for 1-, 2-, 3-, and 4-factor models.
9. Exports factor loadings and defining-sort classifications for each solution.
10. Exports factor arrays with statement metadata, factor z-scores, and forced-distribution scores.
11. Exports distinguishing and consensus statement tables for multi-factor solutions when available from `qmethod`.
12. Exports solution diagnostics and a comparison table.
13. Generates a strict LaTeX solution-comparison report.
14. Saves each raw `qmethod` result object as an `.rds` file for reproducibility and follow-up analysis.

## Repository contents

| Path | Description |
| --- | --- |
| `qmcev.r` | Main executable R workflow. |
| `README.md` | Project documentation and run instructions. |
| `LICENSE` | Repository license. |

The input data file, `Qsort_results.csv`, is expected to be supplied by the user and is not included in this repository.

## Input requirements

Place a UTF-8 CSV file named `Qsort_results.csv` in the repository root before running the script.

### Required columns

The CSV must include these metadata columns:

| Column | Description |
| --- | --- |
| `statement_id` | Unique statement identifier. Must not be blank or duplicated. |
| `statement_text` | Full text of the Q statement. Must not be blank. |
| `category` | Statement category or grouping. Must not be blank. |

The CSV must also include exactly these participant columns:

```text
P1, P2, P3, P4, P5, P6, P7, P8
```

### Required shape

- The file must contain **49 rows**, one row per Q statement.
- Each participant column must contain numeric integer Q-sort scores.
- Scores must be one of the allowed values from `-6` through `6`.
- Each participant must match the forced distribution below.

### Forced distribution

| Score | Expected count per participant |
| ---: | ---: |
| -6 | 1 |
| -5 | 2 |
| -4 | 3 |
| -3 | 4 |
| -2 | 5 |
| -1 | 6 |
| 0 | 7 |
| 1 | 6 |
| 2 | 5 |
| 3 | 4 |
| 4 | 3 |
| 5 | 2 |
| 6 | 1 |

### Minimal CSV header example

```csv
statement_id,statement_text,category,P1,P2,P3,P4,P5,P6,P7,P8
S01,"Example statement text",Category A,-6,0,1,2,-1,3,4,5
```

The example above shows the required columns, but it is not a complete valid input file because a real input must contain all 49 statements and each participant column must satisfy the forced distribution.

## Outputs

All outputs are written under the `qresults/` directory. The script creates the directory tree automatically.

```text
qresults/
  validation/
  correlations/
  factor_solutions/
  factor_arrays/
  distinguishing_consensus/
  reports/
  rds/
```

### Validation outputs

| Output | Description |
| --- | --- |
| `qresults/validation/input_validation_report.csv` | Machine-readable validation report. |
| `qresults/validation/input_validation_report.txt` | Human-readable validation report. |
| `qresults/validation/qsort_distribution_check.csv` | Expected vs. observed forced-distribution counts by score and participant. |
| `qresults/validation/participant_distribution_status.csv` | Pass/fail forced-distribution status for each participant. |
| `qresults/validation/statement_metadata.csv` | Extracted statement metadata. |
| `qresults/validation/q_matrix_numeric.csv` | Numeric Q-sort matrix used for analysis. |

If validation fails, the script writes the validation outputs that can be produced safely and then stops before running the factor analysis.

### Correlation output

| Output | Description |
| --- | --- |
| `qresults/correlations/participant_correlation_matrix.csv` | Pearson correlation matrix among participant Q-sorts. |

### Factor-solution outputs

For each solution label `1F`, `2F`, `3F`, and `4F`, the script writes:

| Output pattern | Description |
| --- | --- |
| `qresults/factor_solutions/qmethod_summary_<solution>.txt` | Captured printed `qmethod` output and summary. |
| `qresults/factor_solutions/loadings_<solution>.csv` | Factor loadings and factor-flag columns. |
| `qresults/factor_solutions/defining_sorts_<solution>.csv` | Classification of each participant as defining, confounded, or non-defining. |
| `qresults/factor_solutions/factor_characteristics_<solution>.txt` | Factor characteristics when available from `qmethod`. |

The script uses `2.58 / sqrt(49)` as the loading-significance threshold when it needs to fall back to loading-based flagging.

### Factor-array outputs

For each solution, the script writes:

| Output pattern | Description |
| --- | --- |
| `qresults/factor_arrays/factor_arrays_<solution>.csv` | Statement metadata plus factor z-scores and factor scores. |
| `qresults/factor_arrays/factor_array_<solution>_F<factor>_ranked.csv` | One ranked statement table per factor. |

If `qmethod` does not provide normalized scores directly, the script assigns forced-distribution scores from factor z-scores.

### Distinguishing and consensus outputs

For multi-factor solutions only (`2F`, `3F`, and `4F`), the script writes distinguishing/consensus statement tables when those can be extracted from `qmethod`:

| Output pattern | Description |
| --- | --- |
| `qresults/distinguishing_consensus/distinguishing_consensus_<solution>.csv` | Distinguishing and consensus statement table, usually with statement metadata attached. |
| `qresults/distinguishing_consensus/distinguishing_consensus_<solution>_not_available.txt` | Written if the table cannot be extracted. |

For the 1-factor solution, distinguishing/consensus analysis is not applicable, so the script writes:

```text
qresults/distinguishing_consensus/distinguishing_consensus_1F_not_applicable.txt
```

### Report and reproducibility outputs

| Output | Description |
| --- | --- |
| `qresults/reports/solution_comparison_table.csv` | Comparison table for 1F through 4F solutions, including variance, defining-sort counts, distinguishing/consensus counts, and diagnostic warning flags. |
| `qresults/reports/solution_comparison_report.tex` | LaTeX report summarizing validation results, the forced distribution, participant status, solution comparison, interpretive cautions, and generated files. |
| `qresults/rds/qmethod_result_<solution>.rds` | Raw serialized `qmethod` result object for each factor solution. |

## Dependencies

### Runtime dependencies

- R.
- The R package `qmethod`.

The script uses several base/recommended R functions from packages that ship with R, including `utils` and `stats`, but it only explicitly checks for and loads `qmethod`.

### Optional dependencies

- A LaTeX distribution, such as TeX Live or TinyTeX, if you want to compile `qresults/reports/solution_comparison_report.tex` to PDF.
- `micromamba`, if you want to create an isolated environment using the instructions below.

## Install with micromamba

The following instructions create an isolated environment containing R and the required R package.

### 1. Install micromamba

If `micromamba` is not already installed, follow the official installation instructions for your operating system:

<https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html>

After installation, open a new shell or initialize your shell as instructed by micromamba.

### 2. Create the environment

From anywhere on your system, create an environment named `qmcev`:

```bash
micromamba create -n qmcev -c conda-forge r-base r-qmethod
```

If your conda-forge mirror does not provide `r-qmethod`, create the environment with R first and install `qmethod` from CRAN inside R:

```bash
micromamba create -n qmcev -c conda-forge r-base
micromamba activate qmcev
Rscript -e 'install.packages("qmethod", repos = "https://cloud.r-project.org")'
```

### 3. Activate the environment

```bash
micromamba activate qmcev
```

### 4. Verify the installation

```bash
Rscript -e 'packageVersion("qmethod")'
```

You should see the installed `qmethod` package version printed without an error.

## Alternative install without micromamba

If you already have R installed, install the only required R package from CRAN:

```r
install.packages("qmethod", repos = "https://cloud.r-project.org")
```

Then verify it from a shell:

```bash
Rscript -e 'packageVersion("qmethod")'
```

## How to run

1. Clone or download this repository.
2. Put your input file at the repository root and name it exactly:

   ```text
   Qsort_results.csv
   ```

3. Activate your environment:

   ```bash
   micromamba activate qmcev
   ```

4. Run the script from the repository root:

   ```bash
   Rscript qmcev.r
   ```

5. Review the outputs under:

   ```text
   qresults/
   ```

The script prints a completion message when it finishes successfully, including the paths to the solution comparison CSV and the LaTeX report.

## Configuration

The script is currently configured directly in the `USER SETTINGS` section near the top of `qmcev.r`.

Important settings include:

| Setting | Current value | Meaning |
| --- | --- | --- |
| `input_file` | `Qsort_results.csv` | CSV file read from the working directory. |
| `output_dir` | `qresults` | Output directory. |
| `expected_n_statements` | `49` | Required number of Q statements. |
| `expected_participants` | `P1` through `P8` | Required participant column names. |
| `factor_solutions` | `1:4` | Factor solutions to run. |
| `extraction_method` | `PCA` | Extraction method passed to `qmethod` when supported by the installed package version. |
| `rotation_method_multi_factor` | `varimax` | Rotation used for solutions with more than one factor. |
| `correlation_method` | `pearson` | Participant correlation method. |

If you change the study design, update these settings and the forced-distribution vector before running the script.

## Troubleshooting

### `Input file not found: Qsort_results.csv`

Run the script from the repository root or update `input_file` in `qmcev.r` to point to the correct CSV path.

### `Missing required package(s): qmethod`

Install `qmethod` in the R environment you are using:

```bash
Rscript -e 'install.packages("qmethod", repos = "https://cloud.r-project.org")'
```

If using micromamba, first activate the environment:

```bash
micromamba activate qmcev
```

### Validation fails

Open these files first:

```text
qresults/validation/input_validation_report.txt
qresults/validation/qsort_distribution_check.csv
qresults/validation/participant_distribution_status.csv
```

Common causes are missing columns, fewer or more than 49 rows, blank statement metadata, non-numeric participant cells, invalid score values, or participant columns that do not match the expected forced distribution.

### LaTeX report does not compile to PDF

The script only writes the `.tex` file. To compile it to PDF, install a LaTeX distribution and run a LaTeX engine manually, for example:

```bash
pdflatex -output-directory qresults/reports qresults/reports/solution_comparison_report.tex
```

## Notes on interpretation

The generated solution comparison table is a diagnostic aid. Do not choose the retained factor solution based only on explained variance. Also consider:

- The number of defining Q-sorts per factor.
- Whether any factors have zero or only one defining sort.
- Whether any sorts are confounded or non-defining.
- Whether factor arrays are coherent and substantively interpretable.
- Whether distinguishing and consensus statements support a meaningful interpretation.

