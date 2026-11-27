#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

# =====================================================================
# Mixed Oral-proxy: study-level Train / External split
#  - Set == 'Oral_proxy' only
#  - 60% of studies in Train and 40% in External (rounded)
#  - seek ~60:40 samples, both outcome classes, and all 4 sites in both sets
#  - split selection uses metadata ONLY; never uses taxonomy or model performance
#  - outputs named Train / External, no duplicate Test files
# =====================================================================

# 0. Parameters --------------------------------------------------------
ROOT_DIR  <- Sys.getenv("OMWI_ROOT", getwd())
TAX_FILE  <- file.path(ROOT_DIR, "data/all_oral", "taxonomy.tsv")
META_FILE <- file.path(ROOT_DIR, "data/all_oral", "all_oral_metadata.csv")

# Fresh folder: avoid confusing old meta_test and otu_test files with new files.
OUT_DIR   <- file.path(ROOT_DIR, "results/preparation/mixed_split")

SET_KEEP              <- "Oral_proxy"
SEED                  <- 704L
RATIO_TRAIN           <- 0.60
SAMPLE_RATIO_TOL      <- 0.05      # require Train sample share between 55% and 65%
N_CANDIDATES          <- 30000L    # repeat metadata-only constrained random splits
MIN_SITE_N_PER_SPLIT  <- 10L       # each observed oral site has >=10 samples in both sets
MIN_CLASS_N_PER_SPLIT <- 10L       # healthy and disease each >=10 in both sets

# 1. Metadata ----------------------------------------------------------
if (!file.exists(META_FILE)) stop("Missing metadata: ", META_FILE)
if (!file.exists(TAX_FILE)) stop("Missing taxonomy: ", TAX_FILE)

cat("[1/5] Read metadata: ", META_FILE, "\n", sep = "")
meta <- fread(META_FILE, check.names = FALSE)
required <- c("sample_id", "study_name", "Set", "body_site", "disease")
missing_cols <- setdiff(required, names(meta))
if (length(missing_cols)) stop("Missing metadata columns: ", paste(missing_cols, collapse = ", "))

for (nm in required) {
  set(meta, j = nm, value = trimws(as.character(meta[[nm]])))
}

n_before <- nrow(meta)
meta <- meta[!is.na(Set) & tolower(Set) == tolower(SET_KEEP)]
cat(sprintf("  Set == %s: %d / %d samples\n", SET_KEEP, nrow(meta), n_before))
if (!nrow(meta)) stop("No Oral_proxy samples after filtering Set.")

for (nm in c("sample_id", "study_name", "body_site", "disease")) {
  if (anyNA(meta[[nm]]) || any(meta[[nm]] == "")) {
    stop("Missing/blank values in metadata column: ", nm)
  }
}
if (anyDuplicated(meta$sample_id)) {
  stop("Duplicate sample_id values among Oral_proxy samples.")
}

meta[, disease := tolower(disease)]
if (!setequal(unique(meta$disease), c("healthy", "disease"))) {
  stop("'disease' must contain only healthy and disease after filtering. Observed: ",
       paste(sort(unique(meta$disease)), collapse = ", "))
}

system_col <- intersect(c("system", "System"), names(meta))
if (length(system_col)) {
  meta[, system_use := trimws(as.character(get(system_col[1L])))]
  meta[system_use %in% c("", "NA", "Na", "N/A", "nan", "NULL"),
       system_use := NA_character_]
}
if (!("system_use" %in% names(meta))) {
  meta[, system_use := NA_character_]
}
if ("group" %in% names(meta)) {
  meta[, group_norm := tolower(trimws(as.character(group)))]
  meta[group_norm %in% c("control", "healthy", "no disease", "normal"),
       group_norm := "control"]
}
cat("  Site distribution:\n")
print(meta[, .N, by = body_site][order(-N)])

# 2. Match taxonomy ----------------------------------------------------
cat("\n[2/5] Read taxonomy: ", TAX_FILE, "\n", sep = "")
tax_dt <- fread(TAX_FILE, check.names = FALSE)
if (!("ID" %in% names(tax_dt))) stop("taxonomy.tsv requires an ID column.")
if (anyDuplicated(names(tax_dt))) stop("Duplicate taxonomy column names.")
if (anyNA(tax_dt$ID) || anyDuplicated(tax_dt$ID)) {
  stop("Missing or duplicate taxonomy feature IDs.")
}

tax_samples <- setdiff(names(tax_dt), "ID")
common <- intersect(meta$sample_id, tax_samples)
cat(sprintf("  Oral_proxy metadata: %d; taxonomy samples: %d; overlap: %d\n",
            nrow(meta), length(tax_samples), length(common)))
if (!length(common)) stop("No overlapping sample_id values between metadata and taxonomy.")

missing_tax <- setdiff(meta$sample_id, common)
meta <- meta[sample_id %in% common]
if (length(missing_tax)) {
  warning(length(missing_tax), " Oral_proxy samples absent from taxonomy and omitted.")
}

# 3. Split by entire studies -------------------------------------------
cat("\n[3/5] Constrained study-level Train / External split ...\n")
studies <- sort(unique(meta$study_name))
n_studies <- length(studies)
if (n_studies < 2L) stop("Fewer than 2 studies remain.")
n_train_studies <- min(n_studies - 1L,
                       max(1L, as.integer(round(n_studies * RATIO_TRAIN))))

site_mat <- as.matrix(xtabs(~ study_name + body_site, data = as.data.frame(meta)))
class_mat <- as.matrix(xtabs(~ study_name + disease, data = as.data.frame(meta)))
site_mat <- site_mat[studies, , drop = FALSE]
class_mat <- class_mat[studies, , drop = FALSE]

n_by_study <- rowSums(site_mat)
n_all <- nrow(meta)
n_by_site <- colSums(site_mat)
n_by_class <- colSums(class_mat)
full_site_frac <- n_by_site / n_all
full_healthy_frac <- n_by_class[["healthy"]] / n_all

set.seed(SEED)
best_score <- Inf
best_train_idx <- NULL
valid_n <- 0L

for (i in seq_len(N_CANDIDATES)) {
  tr_idx <- sample.int(n_studies, n_train_studies, replace = FALSE)
  tr_n <- sum(n_by_study[tr_idx])
  ext_n <- n_all - tr_n
  train_frac <- tr_n / n_all

  if (abs(train_frac - RATIO_TRAIN) > SAMPLE_RATIO_TOL) next

  tr_site <- colSums(site_mat[tr_idx, , drop = FALSE])
  ext_site <- n_by_site - tr_site
  if (any(tr_site < MIN_SITE_N_PER_SPLIT) ||
      any(ext_site < MIN_SITE_N_PER_SPLIT)) next

  tr_class <- colSums(class_mat[tr_idx, , drop = FALSE])
  ext_class <- n_by_class - tr_class
  if (any(tr_class < MIN_CLASS_N_PER_SPLIT) ||
      any(ext_class < MIN_CLASS_N_PER_SPLIT)) next

  valid_n <- valid_n + 1L
  tr_healthy_frac <- tr_class[["healthy"]] / tr_n
  ext_healthy_frac <- ext_class[["healthy"]] / ext_n
  tr_site_frac <- tr_site / tr_n
  ext_site_frac <- ext_site / ext_n

  # Pre-specified metadata-only balance criteria; not an AUC-based selection.
  score <- 8 * abs(train_frac - RATIO_TRAIN) +
           3 * (abs(tr_healthy_frac - full_healthy_frac) +
                abs(ext_healthy_frac - full_healthy_frac)) +
           sum(abs(tr_site_frac - full_site_frac)) +
           sum(abs(ext_site_frac - full_site_frac))

  if (score < best_score) {
    best_score <- score
    best_train_idx <- tr_idx
  }
}

if (is.null(best_train_idx)) {
  stop("No feasible split in ", N_CANDIDATES,
       " trials. Increase N_CANDIDATES or relax metadata constraints; " ,
       "do NOT choose splits based on model AUC.")
}

TRAIN_studies <- studies[best_train_idx]
EXTERNAL_studies <- setdiff(studies, TRAIN_studies)
meta_train <- meta[study_name %in% TRAIN_studies]
meta_external <- meta[study_name %in% EXTERNAL_studies]

stopifnot(
  !length(intersect(TRAIN_studies, EXTERNAL_studies)),
  !length(intersect(meta_train$sample_id, meta_external$sample_id)),
  nrow(meta_train) + nrow(meta_external) == nrow(meta),
  uniqueN(meta_train$study_name) == length(TRAIN_studies),
  uniqueN(meta_external$study_name) == length(EXTERNAL_studies)
)

cat(sprintf("  Valid candidate splits: %d / %d\n", valid_n, N_CANDIDATES))
cat(sprintf("  Train:    %d studies, %d samples (%.1f%%)\n",
            length(TRAIN_studies), nrow(meta_train),
            100 * nrow(meta_train) / n_all))
cat(sprintf("  External: %d studies, %d samples (%.1f%%)\n",
            length(EXTERNAL_studies), nrow(meta_external),
            100 * nrow(meta_external) / n_all))

# 4. Write files -------------------------------------------------------
cat("\n[4/5] Write output files ...\n")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Prevent old legacy files from being mistaken for new groups.
legacy <- file.path(OUT_DIR,
                    c("meta_test.csv", "otu_test_species.tsv", "test_studies.txt"))
if (any(file.exists(legacy))) {
  stop("OUT_DIR contains old test aliases; choose a clean directory: ", OUT_DIR)
}

if (length(missing_tax)) {
  fwrite(data.table(sample_id = missing_tax),
         file.path(OUT_DIR, "oral_proxy_missing_in_taxonomy.csv"))
}

fwrite(data.table(study_name = TRAIN_studies),
       file.path(OUT_DIR, "train_studies.txt"), col.names = FALSE)
fwrite(data.table(study_name = EXTERNAL_studies),
       file.path(OUT_DIR, "external_studies.txt"), col.names = FALSE)

fwrite(meta, file.path(OUT_DIR, "mix_sited_metadata.csv"))
fwrite(meta_train, file.path(OUT_DIR, "meta_train.csv"))
fwrite(meta_external, file.path(OUT_DIR, "meta_external.csv"))

train_cols <- c("ID", meta_train$sample_id)
external_cols <- c("ID", meta_external$sample_id)
tax_train <- copy(tax_dt[, ..train_cols])
tax_external <- copy(tax_dt[, ..external_cols])
setnames(tax_train, "ID", "Taxon")
setnames(tax_external, "ID", "Taxon")
fwrite(tax_train, file.path(OUT_DIR, "otu_train_species.tsv"), sep = "\t")
fwrite(tax_external, file.path(OUT_DIR, "otu_external_species.tsv"), sep = "\t")

study_summary <- meta[, .(
  n_samples = .N,
  n_healthy = sum(disease == "healthy"),
  n_disease = sum(disease == "disease"),
  body_sites = paste(sort(unique(body_site)), collapse = "; ")
), by = study_name]
study_summary[, split := fifelse(study_name %in% TRAIN_studies,
                                "train", "external")]
setcolorder(study_summary,
            c("study_name", "split", "n_samples", "n_healthy", "n_disease", "body_sites"))
setorder(study_summary, split, -n_samples, study_name)
fwrite(study_summary, file.path(OUT_DIR, "study_split_summary.csv"))

meta_counts <- rbindlist(list(
  meta_train[, .(split = "train", body_site, disease, system_use)],
  meta_external[, .(split = "external", body_site, disease, system_use)]
))

fwrite(meta_counts[, .(n = .N), by = .(split, body_site, disease)][order(split, body_site, disease)],
       file.path(OUT_DIR, "split_site_disease_summary.csv"))
fwrite(meta_counts[, .(n = .N), by = .(split, system_use)][order(split, system_use)],
       file.path(OUT_DIR, "split_system_summary.csv"))

fwrite(data.table(
  parameter = c("seed", "ratio_train_studies", "ratio_train_samples_target",
                "sample_ratio_tolerance", "n_candidates", "n_valid_candidates",
                "min_site_n_per_split", "min_class_n_per_split", "selection_score"),
  value = as.character(c(SEED, RATIO_TRAIN, RATIO_TRAIN,
                         SAMPLE_RATIO_TOL, N_CANDIDATES, valid_n,
                         MIN_SITE_N_PER_SPLIT, MIN_CLASS_N_PER_SPLIT,
                         signif(best_score, 6)))
), file.path(OUT_DIR, "split_parameters.csv"))

# 5. Verify ------------------------------------------------------------
cat("\n[5/5] Verify split\n")
print(data.table(
  split = c("train", "external", "total"),
  n_studies = c(length(TRAIN_studies), length(EXTERNAL_studies), n_studies),
  n_samples = c(nrow(meta_train), nrow(meta_external), n_all)
))
for (nm in c("disease", "body_site", "system_use", "group_norm")) {
  if (nm %in% names(meta)) {
    cat("\n", nm, " x split:\n", sep = "")
    print(table(
      split = c(rep("train", nrow(meta_train)),
                rep("external", nrow(meta_external))),
      value = c(meta_train[[nm]], meta_external[[nm]]),
      useNA = "ifany"
    ))
  }
}
cat("\nDone. Outputs in: ", OUT_DIR, "\n", sep = "")
