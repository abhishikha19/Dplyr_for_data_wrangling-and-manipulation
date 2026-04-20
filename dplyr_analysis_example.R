# =============================================================================
# Single R script: dplyr analysis using ONE biological example dataframe `bio`.
# Run: Rscript dplyr_analysis_example.R   |   Needs: install.packages("dplyr")
# =============================================================================

library(dplyr)

# Create the example dataframe once: per-assay rows (subject, tissue, gene, counts, TPM, date).
# Every command below starts from `bio`.
bio <- tibble(
  assay_id = 1:12,
  subject_id = rep(c("S01", "S02", "S01", "S03"), c(3, 4, 2, 3)),
  tissue = rep(c("Liver", "Muscle", "Liver"), c(5, 4, 3)),
  gene = rep(c("GAPDH", "IL6", "GAPDH"), each = 4),
  raw_counts = c(200, 120, 350, 500, 180, 420, 90, 210, 600, 100, 220, 310),
  tpm = c(1200, 800, 1100, 2500, 900, 2100, 750, 1900, 500, 850, 700, 1800),
  collection_date = as.Date("2026-01-01") + c(0, 1, 5, 0, 2, 3, 7, 8, 1, 4, 6, 9)
)

# glimpse — quick column types and sample values (tibble print helper).
glimpse(bio)

# print — show the table (tibbles truncate wide output).
print(bio)

# nrow / ncol / names — base R size and column names (often used with dplyr pipelines).
nrow(bio)
ncol(bio)
names(bio)

# --- Rows: filter, slice family, distinct ------------------------------------

# filter — keep rows where the condition is TRUE (e.g. only liver samples).
bio |> filter(tissue == "Liver")

# filter — comma between conditions means AND.
bio |> filter(subject_id == "S01", tissue == "Liver")

# filter — | means OR (either gene).
bio |> filter(gene == "GAPDH" | gene == "IL6")

# filter — %in% matches any element in a vector (several subjects).
bio |> filter(subject_id %in% c("S01", "S03"))

# filter — numeric comparisons (high read depth).
bio |> filter(raw_counts >= 300)

# filter — dates (samples collected in the first week).
bio |> filter(collection_date <= as.Date("2026-01-07"))

# filter — drop missing (no NAs here; pattern is common with real lab data).
bio |> filter(!is.na(tissue))

# slice — rows by position (current row order).
bio |> slice(1:5)

# slice_head — first n rows.
bio |> slice_head(n = 3)

# slice_tail — last n rows.
bio |> slice_tail(n = 3)

# slice_max — rows with largest raw_counts (ties kept by default).
bio |> slice_max(order_by = raw_counts, n = 3)

# slice_min — rows with smallest raw_counts.
bio |> slice_min(order_by = raw_counts, n = 3)

# slice_sample — random n assays (set seed for reproducibility).
set.seed(1)
bio |> slice_sample(n = 4)

# distinct — unique rows across all columns.
bio |> distinct()

# distinct — unique combinations of subject and tissue (other columns dropped).
bio |> distinct(subject_id, tissue)

# --- Sort -------------------------------------------------------------------

# arrange — sort by collection date.
bio |> arrange(collection_date)

# arrange + desc — tissue A–Z, then raw_counts high to low within tissue.
bio |> arrange(tissue, desc(raw_counts))

# --- Columns: select, rename, relocate --------------------------------------

# select — keep only key assay columns.
bio |> select(subject_id, gene, raw_counts, tpm)

# select — tidyselect helper: columns whose names start with "collection".
bio |> select(assay_id, starts_with("collection"))

# select — rename while selecting (new_name = old_column).
bio |> select(assay_id, patient = subject_id, target_gene = gene)

# rename — change names without dropping columns.
bio |> rename(tpm_estimate = tpm)

# rename_with — transform many names with a function (e.g. uppercase).
bio |> rename_with(toupper)

# relocate — move columns (here: put collection_date first).
bio |> relocate(collection_date, .before = everything())

# --- New columns: mutate, transmute, row-wise helpers ------------------------

# mutate — add a toy combined metric (for dplyr practice; not a real biological model).
bio |> mutate(expression_load = raw_counts * tpm / 1e6)

# mutate — several columns; later expressions can use earlier new columns.
bio |> mutate(
  expression_load = raw_counts * tpm / 1e6,
  high_depth = raw_counts >= 300,
  year = as.integer(format(collection_date, "%Y"))
)

# transmute — like mutate but only keep listed / computed columns.
bio |> transmute(assay_id, expression_load = raw_counts * tpm / 1e6)

# mutate + if_else — two-way vectorised branch (e.g. high vs low TPM).
bio |> mutate(tpm_level = if_else(tpm >= 1500, "high", "low"))

# mutate + case_when — many-way branch on read depth (first match wins).
bio |> mutate(
  depth_band = case_when(
    raw_counts < 150 ~ "low",
    raw_counts <= 350 ~ "medium",
    TRUE ~ "high"
  )
)

# mutate + row_number after group_by — rank assays within each tissue by counts.
bio |> group_by(tissue) |> mutate(rank_counts_in_tissue = row_number(desc(raw_counts))) |> ungroup()

# rowwise + c_across — per-row product of raw_counts and tpm (same as vectorised multiply here).
bio |>
  rowwise() |>
  mutate(count_tpm_product = prod(c_across(all_of(c("raw_counts", "tpm"))))) |>
  ungroup()

# --- Groups: summarise, reframe, count --------------------------------------

# group_by + summarise — one summary row per subject (collapse assays).
bio |>
  group_by(subject_id) |>
  summarise(
    n_assays = n(),
    total_counts = sum(raw_counts),
    mean_tpm = mean(tpm),
    .groups = "drop"
  )

# group_by + summarise — multiple grouping columns (subject × gene).
bio |>
  group_by(subject_id, gene) |>
  summarise(mean_counts = mean(raw_counts), .groups = "drop")

# summarise — no group_by: whole table is one group (one row of summaries).
bio |>
  summarise(
    n_assays = n(),
    mean_raw_counts = mean(raw_counts),
    median_tpm = median(tpm)
  )

# reframe — can return more than one row per group (dplyr 1.1+); list each assay per subject.
bio |>
  group_by(subject_id) |>
  reframe(assay_id, expression_load = raw_counts * tpm / 1e6)

# count — frequency table: how many assays per tissue.
bio |> count(tissue, name = "n_assays", sort = TRUE)

# add_count — append how many assays exist per tissue to each row.
bio |> add_count(tissue, name = "assays_in_tissue")

# group_by + mutate — fraction of that tissue’s total counts on each row.
bio |>
  group_by(tissue) |>
  mutate(pct_of_tissue_counts = raw_counts / sum(raw_counts)) |>
  ungroup()

# group_by + slice_max — assay with highest expression_load per subject.
bio |>
  group_by(subject_id) |>
  slice_max(order_by = raw_counts * tpm, n = 1, with_ties = FALSE) |>
  ungroup()

# --- across (summarise / mutate over many columns) --------------------------

# summarise + across — mean of every numeric column (assay_id included; often you’d exclude IDs).
bio |> summarise(across(where(is.numeric), mean, .names = "mean_{.col}"))

# --- Joins: analyse `bio` with a small subject metadata lookup --------------

subject_meta <- tibble(
  subject_id = c("S01", "S02", "S03"),
  cohort = c("Treatment", "Control", "Treatment"),
  species = rep("Mouse", 3)
)

# left_join — keep every assay; add cohort/species where subject matches.
bio |> left_join(subject_meta, by = "subject_id")

# right_join — keep every subject in metadata; assays without matching subject drop.
bio |> right_join(subject_meta, by = "subject_id")

# inner_join — only assays whose subject appears in metadata.
bio |> inner_join(subject_meta, by = "subject_id")

# full_join — all keys from both sides.
bio |> full_join(subject_meta, by = "subject_id")

# semi_join — keep assays for subjects present in metadata (no extra columns from y).
bio |> semi_join(subject_meta, by = "subject_id")

# anti_join — assays whose subject is NOT in metadata (QA / orphan detection).
bio |> anti_join(subject_meta, by = "subject_id")

# --- Bind rows (same schema as bio: one extra assay row) --------------------

extra_assay <- tibble(
  assay_id = 13L,
  subject_id = "S04",
  tissue = "Brain",
  gene = "IL6",
  raw_counts = 150L,
  tpm = 950,
  collection_date = as.Date("2026-01-20")
)

# bind_rows — stack compatible assay tables (same columns ideal).
bind_rows(bio, extra_assay)

# --- pull — extract one column as a vector ----------------------------------

bio |> pull(raw_counts)

# End: all pipelines above read from the same `bio` object defined at the top.
