## 👩‍💻 Author
Dr. Abhishikha Sharma
# Data wrangling with dplyr (biological example)

Single-tutorial R script that walks through common [**dplyr**](https://dplyr.tidyverse.org/) patterns using one assay-level biological tibble (`bio`) plus small lookup tables for joins and row binding.

## Requirements

- [R](https://www.r-project.org/) **4.1+** (native pipe `|>`)
- Package **dplyr**

Install dplyr in R:

```r
install.packages("dplyr")
```

## Quick start

Clone or download this repository, then from the project directory:

```bash
Rscript dplyr_analysis_example.R
```

Or open `dplyr_analysis_example.R` in RStudio and run it line by line.

## Files

| File | Description |
|------|-------------|
| [`dplyr_analysis_example.R`](dplyr_analysis_example.R) | Full script: defines `bio`, comments every command, prints results. |
| `README.md` | This page: setup, data dictionary, topic index. |

## Example data

### `bio` — assay-level measurements

Each row is one measurement (assay) for a subject, tissue, and gene.

| Column | Type | Description |
|--------|------|-------------|
| `assay_id` | integer | Unique assay identifier |
| `subject_id` | character | Subject ID (`S01`, `S02`, `S03`) |
| `tissue` | character | Sample tissue (`Liver`, `Muscle`, …) |
| `gene` | character | Target gene (e.g. `GAPDH`, `IL6`) |
| `raw_counts` | numeric | Read counts (toy values) |
| `tpm` | numeric | TPM-style normalized expression (toy values) |
| `collection_date` | date | Collection date |

The script also defines:

- **`subject_meta`** — one row per subject: `cohort`, `species` (used for joins).
- **`extra_assay`** — one extra assay row (used with `bind_rows()`).

> **Note:** Combined metrics such as `expression_load = raw_counts * tpm / 1e6` are for **teaching dplyr only**, not a real biological normalization.

## dplyr topics covered

The script is organized in the same order as the sections below. Each line in the R file has a short comment explaining the verb or pattern.

| Topic | Verbs / patterns |
|--------|-------------------|
| Inspect | `glimpse()`, `print()`, `nrow()`, `ncol()`, `names()` |
| Subset rows | `filter()`, `slice()`, `slice_head()`, `slice_tail()`, `slice_max()`, `slice_min()`, `slice_sample()` |
| Unique rows | `distinct()` |
| Sort | `arrange()`, `desc()` |
| Columns | `select()`, `rename()`, `rename_with()`, `relocate()`, tidyselect (`starts_with()`, …) |
| New columns | `mutate()`, `transmute()`, `if_else()`, `case_when()` |
| By row | `rowwise()`, `c_across()`, `ungroup()` |
| Groups | `group_by()`, `summarise()`, `reframe()`, `count()`, `add_count()` |
| Many columns | `across()`, `where()` |
| Combine tables | `left_join()`, `right_join()`, `inner_join()`, `full_join()`, `semi_join()`, `anti_join()` |
| Stack rows | `bind_rows()` |
| Extract vector | `pull()` |


