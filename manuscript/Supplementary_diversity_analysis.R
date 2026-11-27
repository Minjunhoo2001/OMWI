#!/usr/bin/env Rscript
# ============================================================================
# Supplementary Fig. 2: original 12 panels + 2D PCoA + study-adjusted PERMANOVA
# Reference: Supplementaryfig2(2).pdf, 4 rows x 3 columns, A4 portrait.
#
# IMPORTANT:
# - This is NOT the Raw-A/LOSO stability analysis (no glmnet or RDS folds needed).
# - Calculations use original NON-NEGATIVE species relative-abundance profiles.
# - PCoA1/2/3 are RECOMPUTED from species abundance via Bray-Curtis + cmdscale.
# - No previously saved PCoA coordinates are used.
# - The PDF cannot recover sample-level abundances; numeric reproduction requires
#   the same original unfiltered species profiles and sample exclusion criteria.
# - P values: two-sided Wilcoxon rank-sum; d: Cliff's delta (Healthy - Non-healthy).
# - New: separate 2D PCoA visualization and Bray-Curtis PERMANOVA,
#   adjusted for study, with phenotype permutations restricted within study.
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(vegan)
  library(ggplot2)
  library(patchwork)
})
options(stringsAsFactors = FALSE)

# 0. PATHS ------------------------------------------------------------------
ROOT_DIR <- Sys.getenv("OMWI_ROOT", getwd())
RAW_DIR  <- file.path(ROOT_DIR, "data/abundance")
META_DIR <- file.path(ROOT_DIR, "data/metadata")
OUT_DIR  <- file.path(ROOT_DIR, "Supplementaryfigure", "ab")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Optional exact files. Change these only if automatic lookup cannot locate
# an UNFILTERED species abundance matrix or correctly split metadata.
# Use forward slashes in Windows paths. NA means choose from the naming rules.
TRAIN_SPECIES_FILE <- file.path(Sys.getenv("OMWI_ROOT", getwd()), "data/abundance", "train_species.tsv")
TEST_SPECIES_FILE  <- file.path(Sys.getenv("OMWI_ROOT", getwd()), "data/abundance", "external_species.tsv")
TRAIN_META_FILE    <- file.path(META_DIR, "meta_train.csv")
TEST_META_FILE     <- file.path(META_DIR, "meta_external.csv")

# PCoA recomputation choices. No old PCoA coordinate files are accepted.
# combined: one Bray-Curtis ordination on training + external samples, then
#           split for panels c/d. This keeps axes comparable across cohorts.
# separate: independently recompute ordination for each cohort (axes NOT comparable).
PCOA_SCOPE <- "combined"        # change to "separate" only if original used that method
BRAY_PCOA_ADD_CONSTANT <- FALSE # standard cmdscale without Euclidean correction
# Eigenvector signs are mathematically arbitrary. Set a data-only sign convention,
# making reruns more reproducible without using phenotype labels.
ORIENT_AXES_DETERMINISTIC <- TRUE

# Expected original panel counts, used for validation rather than for labels.
# Reference sample counts from SupplementaryFig2_recomputed_PCoA.pdf (latest).
# These are QC-only: do not force the underlying data to match them.
EXPECTED_ALPHA <- list(Training = c(Healthy = 574L, `Non-healthy` = 503L),
                       `External validation` = c(Healthy = 385L, `Non-healthy` = 331L))
EXPECTED_PCOA  <- list(Training = c(Healthy = 574L, `Non-healthy` = 503L),
                       `External validation` = c(Healthy = 385L, `Non-healthy` = 330L))

REQUIRE_REFERENCE_ALPHA_COUNTS <- FALSE   # warn on differences, still produce QC and figure
# If no original PCoA table is supplied, we warn about mismatched PCoA n.
FONT_FAMILY <- "Arial"                   # available on Windows
COL_HEALTH <- "#4DBBD5"
COL_NONHEALTH <- "#D95F59"
FIG_WIDTH <- 8.27                         # A4 portrait, inches
FIG_HEIGHT <- 11.69
FIG_DPI <- 600
SEED_JITTER <- 704
# New analysis settings: at least one study must contain BOTH groups for
# a meaningfully identifiable study-adjusted group effect.
PERMANOVA_NPERM <- 999L
PERMANOVA_SEED <- 20261121L
PCOA_2D_WIDTH <- 7.3
PCOA_2D_HEIGHT <- 8.5

# 1. FILE DISCOVERY ---------------------------------------------------------
assert_dir <- function(path) {
  if (!dir.exists(path)) stop("Directory not found: ", path, call. = FALSE)
}
assert_dir(RAW_DIR)
assert_dir(META_DIR)

choose_input <- function(explicit, directory, candidates, pattern, label) {
  if (!is.na(explicit) && nzchar(explicit)) {
    if (!file.exists(explicit)) stop(label, " not found: ", explicit, call. = FALSE)
    return(normalizePath(explicit, winslash = "/"))
  }
  # Prioritize exact filenames before regex; never choose from two fuzzy matches.
  for (nm in candidates) {
    path <- file.path(directory, nm)
    if (file.exists(path)) return(normalizePath(path, winslash = "/"))
  }
  all <- list.files(directory, full.names = TRUE, recursive = TRUE,
                    pattern = "\\.(csv|tsv|txt)$", ignore.case = TRUE)
  hits <- all[grepl(pattern, basename(all), ignore.case = TRUE, perl = TRUE)]
  if (length(hits) != 1L) {
    stop(label, ": unable to identify exactly one input file in ", directory,
         "; found ", length(hits), " candidate(s).\n",
         if (length(hits)) paste(hits, collapse = "\n") else "(none)",
         "\nSet the corresponding *_FILE variable at the top of this script.",
         call. = FALSE)
  }
  normalizePath(hits[1], winslash = "/")
}

train_ab_file <- choose_input(
  TRAIN_SPECIES_FILE, RAW_DIR,
  c("train_species.raw.tsv", "train_species.raw.csv",
    "train_species.tsv", "train_species.csv",
    "train_species.biom.tsv", "train_species.biom.csv"),
  "^(train.*species|species.*train).*\\.(tsv|csv|txt)$", "Training species table")
test_ab_file <- choose_input(
  TEST_SPECIES_FILE, RAW_DIR,
  c("test_species.raw.tsv", "test_species.raw.csv",
    "test_species.tsv", "test_species.csv",
    "test_species.biom.tsv", "test_species.biom.csv"),
  "^(test.*species|species.*test).*\\.(tsv|csv|txt)$", "Test species table")
train_meta_file <- choose_input(
  TRAIN_META_FILE, META_DIR,
  c("train_metadata.csv", "train_metadata.tsv", "train_meta.csv",
    "metadata_train.csv", "metadata_train.tsv", "train_metadata_mixed64.csv"),
  "^(train.*meta|meta.*train).*\\.(csv|tsv|txt)$", "Training metadata")
test_meta_file <- choose_input(
  TEST_META_FILE, META_DIR,
  c("test_metadata.csv", "test_metadata.tsv", "test_meta.csv",
    "metadata_test.csv", "metadata_test.tsv", "test_metadata_mixed64.csv"),
  "^(test.*meta|meta.*test).*\\.(csv|tsv|txt)$", "Test metadata")

paths <- data.frame(
  role = c("training_species", "test_species", "training_metadata", "test_metadata"),
  path = c(train_ab_file, test_ab_file, train_meta_file, test_meta_file))
print(paths, row.names = FALSE)
fwrite(paths, file.path(OUT_DIR, "SuppFig2_input_manifest.csv"))

# 2. METADATA AND RAW MATRIX -----------------------------------------------
find_col <- function(df, candidates, label) {
  hit <- intersect(candidates, names(df))
  if (!length(hit)) stop("Cannot find ", label, " column. Available: ",
                         paste(names(df), collapse = ", "), call. = FALSE)
  hit[1]
}

read_metadata <- function(file) {
  x <- as.data.frame(fread(file, check.names = FALSE), check.names = FALSE)
  sid <- find_col(x, c("sample_id", "SampleID", "sample_name", "Sample_Name",
                       "sample", "Sample", "sampleid", "run_accession", "Run"),
                  "sample ID")
  grp <- find_col(x, c("disease", "group", "Group", "phenotype",
                       "health_status", "HealthStatus", "status"), "health status")
  ids <- trimws(as.character(x[[sid]]))
  raw <- trimws(as.character(x[[grp]]))
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids))
    stop("Missing/duplicated sample IDs in metadata: ", file, call. = FALSE)
  raw_lower <- tolower(raw)
  if (!any(raw_lower %in% c("healthy", "health", "control", "ctrl", "normal")))
    stop("No explicit Healthy records in ", file,
         ". Check label mapping; numeric 0/1 is intentionally NOT guessed.",
         call. = FALSE)
  if (anyNA(raw) || any(!nzchar(raw)))
    stop("Missing phenotype labels in ", file, call. = FALSE)
  group <- ifelse(raw_lower %in% c("healthy", "health", "control", "ctrl", "normal"),
                  "Healthy", "Non-healthy")
  # Optional for the plots, essential for study-adjusted PERMANOVA.
  # Do not guess a study identifier from a sample ID or use Training/External as a study.
  study_candidates <- c("study_name", "study_use", "study", "Study",
                        "study_id", "Study_ID", "StudyID", "studyName",
                        "dataset_name", "dataset_id", "Dataset")
  study_hit <- intersect(study_candidates, names(x))
  if (!length(study_hit)) {
    warning("No study origin column found in ", file,
            ". 12-panel / 2D figures will run, but study-adjusted PERMANOVA is skipped.",
            call. = FALSE)
    study <- rep(NA_character_, length(ids))
  } else {
    study <- trimws(as.character(x[[study_hit[1L]]]))
    study[is.na(study) | !nzchar(study) | study %in% c("NA", "N/A", "Unknown")] <- NA_character_
  }
  data.frame(sample_id = ids, group = group, study = study,
             stringsAsFactors = FALSE)
}

read_abundance <- function(file, sample_ids) {
  dat <- as.data.frame(fread(file, check.names = FALSE), check.names = FALSE)
  if (ncol(dat) < 2L) stop("Not a valid species abundance table: ", file)
  first <- trimws(as.character(dat[[1L]]))
  col_hits <- intersect(sample_ids, names(dat)[-1L])
  row_hits <- intersect(sample_ids, first)
  if (length(col_hits) > length(row_hits) && length(col_hits) >= 5L) {
    # MetaPhlAn-like: first column = taxon, other columns = samples
    taxa <- first
    if (anyNA(taxa) || any(!nzchar(taxa)) || anyDuplicated(taxa))
      stop("Invalid / duplicate feature names in ", file)
    num <- dat[, col_hits, drop = FALSE]
    x <- t(as.matrix(data.frame(lapply(num, as.numeric), check.names = FALSE)))
    rownames(x) <- col_hits
    colnames(x) <- taxa
  } else if (length(row_hits) >= 5L && length(row_hits) >= length(col_hits)) {
    # Samples x taxa
    ids <- first
    if (anyDuplicated(ids)) stop("Duplicate sample IDs in ", file)
    taxa <- names(dat)[-1L]
    num <- dat[, -1L, drop = FALSE]
    x <- as.matrix(data.frame(lapply(num, as.numeric), check.names = FALSE))
    rownames(x) <- ids
    colnames(x) <- taxa
    x <- x[row_hits, , drop = FALSE]
  } else {
    stop("Cannot establish sample/feature table orientation in ", file,
         ". Metadata matches columns: ", length(col_hits),
         ", first-column rows: ", length(row_hits), call. = FALSE)
  }
  storage.mode(x) <- "double"
  if (anyNA(x) || any(!is.finite(x))) stop("NA/Inf abundance in ", file)
  if (any(x < 0)) stop("Negative abundances found in ", file,
                       ". This is likely CLR/batch-corrected data. Use RAW, non-negative abundances.")
  if (ncol(x) < 100) warning("Only ", ncol(x), " taxa in ", file,
                             "; richness will not match an unfiltered original figure.")
  x
}

prepare_split <- function(ab_file, meta_file, split_name) {
  md <- read_metadata(meta_file)
  x <- read_abundance(ab_file, md$sample_id)
  matched <- md$sample_id[md$sample_id %in% rownames(x)]
  omitted <- setdiff(md$sample_id, matched)
  if (length(omitted)) warning(split_name, ": ", length(omitted),
                               " metadata samples missing from species table")
  x <- x[matched, , drop = FALSE]
  md <- md[match(matched, md$sample_id), , drop = FALSE]
  if (!identical(rownames(x), md$sample_id)) stop("Sample order mismatch")
  totals <- rowSums(x)
  if (any(totals == 0)) warning(split_name, ": ", sum(totals == 0),
                                " all-zero samples: included in alpha, excluded in PCoA")
  nonzero <- totals > 0
  x[nonzero, ] <- sweep(x[nonzero, , drop = FALSE], 1, totals[nonzero], "/")
  # vegan's diversity returns 0 for all-zero samples; specnumber returns 0.
  shannon <- vegan::diversity(x, index = "shannon")
  simpson <- vegan::diversity(x, index = "simpson")
  richness <- vegan::specnumber(x)
  alpha <- data.frame(
    set = split_name, sample_id = md$sample_id, group = md$group,
    Shannon = shannon, Simpson = simpson, Richness = richness,
    stringsAsFactors = FALSE)
  list(meta = md, x = x, nonzero = nonzero, alpha = alpha)
}

training <- prepare_split(train_ab_file, train_meta_file, "Training")
external <- prepare_split(test_ab_file, test_meta_file, "External validation")
alpha_all <- rbind(training$alpha, external$alpha)
fwrite(alpha_all, file.path(OUT_DIR, "SuppFig2_alpha_per_sample.csv"))

count_groups <- function(df) {
  z <- table(factor(df$group, levels = c("Healthy", "Non-healthy")))
  setNames(as.integer(z), names(z))
}
for (nm in names(EXPECTED_ALPHA)) {
  observed <- count_groups(alpha_all[alpha_all$set == nm, , drop = FALSE])
  target <- EXPECTED_ALPHA[[nm]]
  if (!identical(observed, target)) {
    txt <- paste0(nm, ": alpha counts observed ", paste(observed, collapse = "/"),
                  ", original ", paste(target, collapse = "/"), ".")
    if (REQUIRE_REFERENCE_ALPHA_COUNTS) stop(txt, " Check input files.")
    warning(txt)
  }
}

# 3. PCoA: RECOMPUTE from species abundances, never load old coordinates -------
# The two split abundance tables may contain different taxon columns. Align
# them to their union (missing taxa = 0); never silently drop split-specific taxa.
common_feature_space <- function(a, b) {
  if (is.null(colnames(a)) || is.null(colnames(b)))
    stop("Both matrices must have taxon names")
  taxa <- sort(union(colnames(a), colnames(b)))
  build <- function(x) {
    z <- matrix(0, nrow = nrow(x), ncol = length(taxa),
                dimnames = list(rownames(x), taxa))
    z[, colnames(x)] <- x
    z
  }
  list(train = build(a), external = build(b), taxa = taxa)
}

# Read preparatory objects from alpha stage; discard only all-zero samples for
# Bray-Curtis, which is undefined for two empty profiles. No phenotype-based
# sample filtering is performed.
training_nonzero <- training$meta[training$nonzero, , drop = FALSE]
external_nonzero <- external$meta[external$nonzero, , drop = FALSE]
training_nonzero$set <- "Training"
external_nonzero$set <- "External validation"
metadata_pcoa <- rbind(training_nonzero, external_nonzero)
if (anyDuplicated(metadata_pcoa$sample_id))
  stop("Sample IDs overlap between training and external splits; cannot recompute PCoA safely.")

aligned <- common_feature_space(
  training$x[training$nonzero, , drop = FALSE],
  external$x[external$nonzero, , drop = FALSE]
)
if (length(aligned$taxa) < 3L)
  stop("PCoA requires at least three species features")

removed_zero <- rbind(
  data.frame(set = rep("Training", sum(!training$nonzero)),
             sample_id = training$meta$sample_id[!training$nonzero]),
  data.frame(set = rep("External validation", sum(!external$nonzero)),
             sample_id = external$meta$sample_id[!external$nonzero])
)
if (nrow(removed_zero)) {
  removed_zero$reason <- "All-zero species abundance; Bray-Curtis undefined"
} else {
  removed_zero$reason <- character(0)
}
fwrite(removed_zero, file.path(OUT_DIR, "SuppFig2_PCoA_excluded_samples.csv"))

calculate_bray_pcoa <- function(ab_mat, scope_label) {
  # Deterministic sample ordering; input values are already row-normalized.
  ab_mat <- ab_mat[order(rownames(ab_mat)), , drop = FALSE]
  ab_mat <- ab_mat[, colSums(ab_mat > 0) > 0, drop = FALSE]
  if (nrow(ab_mat) < 5L || ncol(ab_mat) < 3L)
    stop("Not enough samples/taxa for PCoA in ", scope_label)
  if (any(rowSums(ab_mat) == 0))
    stop("Unexpected all-zero sample reached Bray-Curtis: ", scope_label)
  
  # Full recomputation, not a replot of old PCoA coordinates.
  bc <- vegan::vegdist(ab_mat, method = "bray")
  ord <- stats::cmdscale(bc, k = 3L, eig = TRUE,
                         add = BRAY_PCOA_ADD_CONSTANT)
  if (is.null(dim(ord$points)) || ncol(ord$points) < 3L)
    stop("Unable to estimate three PCoA axes in ", scope_label)
  coords <- as.data.frame(ord$points[, 1:3, drop = FALSE])
  names(coords) <- c("PCoA1", "PCoA2", "PCoA3")
  if (ORIENT_AXES_DETERMINISTIC) {
    for (nm in names(coords)) {
      farthest <- which.max(abs(coords[[nm]]))
      if (coords[[nm]][farthest] < 0) coords[[nm]] <- -coords[[nm]]
    }
  }
  coords$sample_id <- rownames(ab_mat)
  
  eigen <- as.numeric(ord$eig)
  denom <- sum(eigen[eigen > 0])
  eigen_report <- data.frame(
    scope = scope_label,
    axis = paste0("PCoA", 1:3),
    eigenvalue = eigen[1:3],
    percent_positive_eigenvalues = 100 * eigen[1:3] / denom,
    negative_eigenvalues_n = sum(eigen < -1e-10),
    n_samples = nrow(ab_mat),
    n_taxa_nonzero = ncol(ab_mat),
    stringsAsFactors = FALSE
  )
  list(coords = coords, eigen = eigen_report)
}

if (!PCOA_SCOPE %in% c("combined", "separate"))
  stop("PCOA_SCOPE must be 'combined' or 'separate'")

if (PCOA_SCOPE == "combined") {
  # A single ordination is appropriate if axes from panels c and d are meant
  # to have the same meaning. It is done BEFORE splitting for plotting.
  all_x <- rbind(aligned$train, aligned$external)
  pcoa_fit <- calculate_bray_pcoa(all_x, "Combined training + external")
  pcoa_points <- pcoa_fit$coords
  eigen_report <- pcoa_fit$eigen
} else {
  fit_tr <- calculate_bray_pcoa(aligned$train, "Training")
  fit_ex <- calculate_bray_pcoa(aligned$external, "External validation")
  pcoa_points <- rbind(fit_tr$coords, fit_ex$coords)
  eigen_report <- rbind(fit_tr$eigen, fit_ex$eigen)
}

row_idx <- match(pcoa_points$sample_id, metadata_pcoa$sample_id)
if (anyNA(row_idx)) stop("PCoA samples could not be matched to metadata")
pc_all <- data.frame(
  set = metadata_pcoa$set[row_idx],
  sample_id = pcoa_points$sample_id,
  group = metadata_pcoa$group[row_idx],
  study = metadata_pcoa$study[row_idx],
  PCoA1 = pcoa_points$PCoA1,
  PCoA2 = pcoa_points$PCoA2,
  PCoA3 = pcoa_points$PCoA3,
  stringsAsFactors = FALSE
)
pc_all <- pc_all[order(match(pc_all$set, c("Training", "External validation")),
                       pc_all$sample_id), , drop = FALSE]
fwrite(pc_all, file.path(OUT_DIR, "SuppFig2_PCoA_coordinates_RECALCULATED.csv"))
fwrite(eigen_report, file.path(OUT_DIR, "SuppFig2_PCoA_eigenvalues.csv"))

for (nm in names(EXPECTED_PCOA)) {
  observed <- count_groups(pc_all[pc_all$set == nm, , drop = FALSE])
  target <- EXPECTED_PCOA[[nm]]
  if (!identical(observed, target)) {
    warning(nm, ": PCoA counts recalculated ", paste(observed, collapse = "/"),
            ", original plot ", paste(target, collapse = "/"),
            ". Do not remove an extra sample merely to match the plot;",
            " inspect SuppFig2_PCoA_excluded_samples.csv and raw profiles.")
  }
}
message("PCoA RECOMPUTED from Bray-Curtis; scope = ", PCOA_SCOPE,
        "; retained samples = ", nrow(pc_all),
        "; species features = ", length(aligned$taxa))

# 3b. STUDY-ADJUSTED BRAY-CURTIS PERMANOVA ---------------------------------
# Testing is conducted SEPARATELY in Training and External validation.
# Unlike the ordination visualization, PERMANOVA is computed on the full
# Bray-Curtis distance matrix within each set, not on PCoA axes.
# Study is a fixed nuisance effect: group is tested conditional on study.
# Permutations are restricted within the original study; this matters because
# cross-study exchangeability cannot be assumed in multi-cohort data.
# If study origin is missing or no within-study group contrast is available,
# DO NOT silently replace this with a naive pooled PERMANOVA.

study_counts <- as.data.frame(table(
  set = factor(metadata_pcoa$set, levels = c("Training", "External validation")),
  study = factor(metadata_pcoa$study, exclude = NULL),
  group = factor(metadata_pcoa$group, levels = c("Healthy", "Non-healthy")),
  useNA = "ifany"), stringsAsFactors = FALSE)
study_counts <- study_counts[study_counts$Freq > 0L, , drop = FALSE]
fwrite(study_counts, file.path(OUT_DIR, "SuppFig2_study_group_counts.csv"))

permanova_one <- function(cohort, ab_mat, md) {
  # Match PCoA-eligible (nonzero) samples; do not reuse ordination coordinates.
  md <- md[match(rownames(ab_mat), md$sample_id), , drop = FALSE]
  if (anyNA(md$sample_id) || !identical(md$sample_id, rownames(ab_mat)))
    stop("PERMANOVA sample / metadata mismatch in ", cohort)
  total_n <- nrow(md)
  valid <- !is.na(md$study) & nzchar(md$study)
  missing_n <- sum(!valid)
  md <- md[valid, , drop = FALSE]
  x <- ab_mat[valid, , drop = FALSE]
  report <- data.frame(
    set = cohort, status = NA_character_, n_PCoA = total_n,
    n_study_missing = missing_n, n_tested = nrow(md),
    n_studies = length(unique(md$study)), n_mixed_studies = NA_integer_,
    n_samples_in_mixed_studies = NA_integer_,
    n_permutations = PERMANOVA_NPERM,
    group_df = NA_real_, group_R2 = NA_real_, group_F = NA_real_,
    group_p = NA_real_, stringsAsFactors = FALSE)
  if (missing_n) warning(cohort, ": ", missing_n,
                         " PCoA samples have missing study origin; omitted from adjusted PERMANOVA")
  if (nrow(md) < 5L || length(unique(md$study)) < 2L ||
      length(unique(md$group)) < 2L) {
    report$status <- "SKIPPED: inadequate study/group data"
    return(report)
  }
  tab <- table(md$study, md$group)
  if (!all(c("Healthy", "Non-healthy") %in% colnames(tab))) {
    report$status <- "SKIPPED: one phenotype only"
    return(report)
  }
  mixed <- rownames(tab)[rowSums(tab > 0L) == 2L]
  report$n_mixed_studies <- length(mixed)
  report$n_samples_in_mixed_studies <- sum(md$study %in% mixed)
  if (!length(mixed)) {
    report$status <- "SKIPPED: phenotype fully confounded by study"
    warning(cohort, ": no study includes both groups; adjusted group effect cannot be tested")
    return(report)
  }
  md$study <- factor(md$study)
  md$group <- factor(md$group, levels = c("Healthy", "Non-healthy"))
  # Study-adjusted model; permutation blocks preserve original study identities.
  result <- tryCatch({
    bc <- vegan::vegdist(x, method = "bray")
    perm <- permute::how(nperm = PERMANOVA_NPERM, blocks = md$study)
    fit <- vegan::adonis2(bc ~ study + group, data = md,
                          permutations = perm, by = "margin")
    f <- as.data.frame(fit)
    row <- match("group", rownames(f))
    if (is.na(row)) stop("Cannot find group term in adonis2 marginal table")
    report$group_df <- as.numeric(f[row, "Df"])
    report$group_R2 <- as.numeric(f[row, "R2"])
    report$group_F <- as.numeric(f[row, "F"])
    report$group_p <- as.numeric(f[row, "Pr(>F)"])
    report$status <- "OK: group adjusted for study; permutations blocked by study"
    report
  }, error = function(e) {
    report$status <- paste0("FAILED: ", conditionMessage(e))
    warning(cohort, ": study-adjusted PERMANOVA error: ", conditionMessage(e))
    report
  })
  result
}

set.seed(PERMANOVA_SEED)
permanova_out <- rbind(
  permanova_one("Training", aligned$train,
                training_nonzero),
  permanova_one("External validation", aligned$external,
                external_nonzero)
)
fwrite(permanova_out, file.path(OUT_DIR, "SuppFig2_PERMANOVA_study_adjusted.csv"))
print(permanova_out, row.names = FALSE)

# 4. STATS: Wilcoxon P & Cliff's delta --------------------------------------
stat_one <- function(df, metric) {
  sub <- df[is.finite(df[[metric]]) & df$group %in% c("Healthy", "Non-healthy"), , drop = FALSE]
  h <- sub[[metric]][sub$group == "Healthy"]
  n <- sub[[metric]][sub$group == "Non-healthy"]
  if (!length(h) || !length(n)) stop("Missing either phenotype for ", metric)
  w <- suppressWarnings(stats::wilcox.test(h, n, alternative = "two.sided",
                                           exact = FALSE, correct = TRUE))
  # Wilcoxon's W statistic is the Mann-Whitney U for the first sample.
  d <- 2 * as.numeric(w$statistic) / (length(h) * length(n)) - 1
  data.frame(metric = metric, p_value = unname(w$p.value),
             cliffs_delta = d, n_healthy = length(h), n_nonhealthy = length(n))
}

stat_rows <- list()
for (cohort in c("Training", "External validation")) {
  for (metric in c("Shannon", "Simpson", "Richness")) {
    s <- stat_one(alpha_all[alpha_all$set == cohort, , drop = FALSE], metric)
    s$set <- cohort; s$source <- "alpha"; stat_rows[[length(stat_rows) + 1L]] <- s
  }
  for (metric in c("PCoA1", "PCoA2", "PCoA3")) {
    s <- stat_one(pc_all[pc_all$set == cohort, , drop = FALSE], metric)
    s$set <- cohort; s$source <- "PCoA"; stat_rows[[length(stat_rows) + 1L]] <- s
  }
}
stats_all <- do.call(rbind, stat_rows)
stats_all <- stats_all[, c("set", "source", "metric", "n_healthy", "n_nonhealthy",
                           "p_value", "cliffs_delta")]
fwrite(stats_all, file.path(OUT_DIR, "SuppFig2_12panel_stats.csv"))
print(stats_all, row.names = FALSE)

# 5. FIGURE THEME -----------------------------------------------------------
fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("< 0.001")
  if (p < 0.01) return("< 0.01")
  if (p < 0.05) return("< 0.05")
  paste0("= ", sprintf("%.3f", p))
}

metric_yname <- c(Shannon = "Shannon index", Simpson = "Simpson index",
                  Richness = "Species richness", PCoA1 = "PCoA1 coordinate",
                  PCoA2 = "PCoA2 coordinate", PCoA3 = "PCoA3 coordinate")

plot_panel <- function(df, cohort, metric) {
  sub <- df[df$set == cohort & is.finite(df[[metric]]), , drop = FALSE]
  sub$group <- factor(sub$group, levels = c("Healthy", "Non-healthy"))
  sub$value <- sub[[metric]]
  st <- stat_one(sub, "value")
  lab <- paste0("P ", fmt_p(st$p_value), "\nd = ", sprintf("%.2f", st$cliffs_delta))
  xlabs <- c(paste0("Healthy\n(n=", st$n_healthy, ")"),
             paste0("Non-healthy\n(n=", st$n_nonhealthy, ")"))
  vmin <- min(sub$value)
  vmax <- max(sub$value)
  vspan <- max(vmax - vmin, 1e-4)
  # Put statistics in a fixed, unclipped upper band without covering violins.
  vtop <- vmax + 0.13 * vspan
  vbottom <- if (metric %in% c("Shannon", "Simpson", "Richness")) {
    0
  } else vmin - 0.06 * vspan
  set.seed(SEED_JITTER)
  ggplot(sub, aes(x = group, y = value, fill = group, colour = group)) +
    geom_violin(trim = FALSE, width = 0.80, alpha = 0.32, linewidth = 0.70) +
    geom_point(position = position_jitter(width = 0.11, height = 0),
               alpha = 0.20, size = 0.25, stroke = 0) +
    geom_boxplot(aes(group = group), width = 0.20, outlier.shape = NA,
                 fill = "white", colour = "#777777", alpha = 0.70,
                 linewidth = 0.55, fatten = 1.0) +
    annotate("text", x = 1.5, y = vtop, label = lab,
             family = FONT_FAMILY, fontface = "italic", size = 2.60,
             hjust = 0.5, vjust = 1.0, lineheight = 0.95) +
    scale_fill_manual(values = c("Healthy" = COL_HEALTH,
                                 "Non-healthy" = COL_NONHEALTH)) +
    scale_colour_manual(values = c("Healthy" = COL_HEALTH,
                                   "Non-healthy" = COL_NONHEALTH)) +
    scale_x_discrete(labels = xlabs, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02)),
                       breaks = scales::breaks_pretty(n = 4)) +
    coord_cartesian(ylim = c(vbottom, vtop + 0.11 * vspan), clip = "off") +
    labs(title = metric, x = NULL, y = metric_yname[[metric]]) +
    theme_bw(base_family = FONT_FAMILY, base_size = 7.0) +
    theme(
      legend.position = "none",
      panel.border = element_rect(fill = NA, colour = "#565656", linewidth = 0.60),
      panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(colour = "#D8D8D8", linetype = "dashed", linewidth = 0.42),
      axis.ticks.x = element_blank(),
      axis.text = element_text(colour = "#111111", size = 7.0),
      axis.text.x = element_text(size = 6.9, lineheight = 0.90,
                                 margin = margin(t = 3)),
      axis.title.y = element_text(size = 7.5, margin = margin(r = 2.5)),
      axis.title.x = element_blank(),
      plot.title = element_text(face = "bold", colour = "black", size = 8.5,
                                hjust = 0.5, margin = margin(b = 4)),
      plot.margin = margin(t = 2, r = 5.5, b = 4, l = 2)
    )
}

# Narrow label column (a-d, left-rotated cohort name), same layout every row.
row_caption <- function(letter, cohort) {
  caption_label <- if (cohort == "Training") "Training" else "External validation"
  ggplot() +
    annotate("text", x = 0.17, y = 0.96, label = letter,
             family = FONT_FAMILY, fontface = "bold", size = 3.8,
             hjust = 0, vjust = 1) +
    annotate("text", x = 0.47, y = 0.48, label = caption_label,
             angle = 90, family = FONT_FAMILY, fontface = "bold", size = 2.9) +
    xlim(0, 1) + ylim(0, 1) + theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
}

make_row <- function(letter, cohort, source, metrics) {
  dat <- if (source == "alpha") alpha_all else pc_all
  row_caption(letter, cohort) +
    plot_panel(dat, cohort, metrics[1]) +
    plot_panel(dat, cohort, metrics[2]) +
    plot_panel(dat, cohort, metrics[3]) +
    plot_layout(ncol = 4, widths = c(0.16, 1, 1, 1))
}

row_a <- make_row("a", "Training", "alpha", c("Shannon", "Simpson", "Richness"))
row_b <- make_row("b", "External validation", "alpha", c("Shannon", "Simpson", "Richness"))
row_c <- make_row("c", "Training", "PCoA", c("PCoA1", "PCoA2", "PCoA3"))
row_d <- make_row("d", "External validation", "PCoA", c("PCoA1", "PCoA2", "PCoA3"))
fig <- (row_a / row_b / row_c / row_d) +
  plot_layout(heights = c(1, 1, 1, 1)) +
  plot_annotation(
    theme = theme(
      plot.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(1, 1, 1, 1)
    )
  )

print(fig)

pdf_path <- file.path(OUT_DIR, "SupplementaryFig2_recomputed_PCoA.pdf")
png_path <- file.path(OUT_DIR, "SupplementaryFig2_recomputed_PCoA.png")
svg_path <- file.path(OUT_DIR, "SupplementaryFig2_recomputed_PCoA.svg")

# cairo_pdf works with Arial on Windows and embeds/scales text better than base PDF.
ggsave(pdf_path, fig, width = FIG_WIDTH, height = FIG_HEIGHT,
       units = "in", device = grDevices::cairo_pdf, bg = "white", limitsize = FALSE)
ggsave(png_path, fig, width = FIG_WIDTH, height = FIG_HEIGHT,
       units = "in", dpi = FIG_DPI, bg = "white", limitsize = FALSE)
# SVG is useful for subsequent Adobe Illustrator edits (requires svglite).
if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(svg_path, fig, width = FIG_WIDTH, height = FIG_HEIGHT,
         units = "in", device = svglite::svglite, bg = "white")
}
message("\nFinished. Output directory: ", OUT_DIR)
message("PDF: ", pdf_path)
message("PNG: ", png_path)
message("Statistics: ", file.path(OUT_DIR, "SuppFig2_12panel_stats.csv"))
message("Fresh PCoA: ", file.path(OUT_DIR, "SuppFig2_PCoA_coordinates_RECALCULATED.csv"))

# 6. ADDITIONAL 2D PCoA (PCoA1 vs PCoA2) -------------------------------------
# This figure does NOT replace or modify the original 12-panel Fig. 2.
# It visualizes the first two axes of the SAME recomputed Bray-Curtis PCoA,
# not a new ordination (therefore no possibility of mismatched coordinates).

pc_plot <- pc_all
pc_plot$set <- factor(pc_plot$set, levels = c("Training", "External validation"))
# Continued panel numbering after the original four rows (a-d).
pc_plot$facet_label <- factor(pc_plot$set,
                              levels = c("Training", "External validation"),
                              labels = c("e   Training", "f   External validation"))
pc_plot$group <- factor(pc_plot$group, levels = c("Healthy", "Non-healthy"))
if (PCOA_SCOPE == "combined") {
  axis_rows <- eigen_report[eigen_report$scope == "Combined training + external", , drop = FALSE]
  axis_x <- axis_rows$percent_positive_eigenvalues[match("PCoA1", axis_rows$axis)]
  axis_y <- axis_rows$percent_positive_eigenvalues[match("PCoA2", axis_rows$axis)]
  x_label <- sprintf("PCoA1 (%.1f%%)", axis_x)
  y_label <- sprintf("PCoA2 (%.1f%%)", axis_y)
  facet_scales <- "fixed"
} else {
  # Separate ordinations have different eigenvalue percentages and axes.
  x_label <- "PCoA1 coordinate (separate ordinations)"
  y_label <- "PCoA2 coordinate (separate ordinations)"
  facet_scales <- "free"
}

pcoa_scatter <- ggplot(pc_plot, aes(x = PCoA1, y = PCoA2, colour = group)) +
  geom_hline(yintercept = 0, colour = "grey80", linetype = "dotted", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey80", linetype = "dotted", linewidth = 0.3) +
  geom_point(size = 1.05, alpha = 0.46, stroke = 0) +
  facet_wrap(~facet_label, nrow = 2, scales = facet_scales) +
  scale_colour_manual(values = c("Healthy" = COL_HEALTH,
                                 "Non-healthy" = COL_NONHEALTH), drop = FALSE) +
  labs(x = x_label, y = y_label, colour = NULL,
       title = "Bray–Curtis PCoA of oral microbial communities",
       subtitle = paste0("Each point represents one sample; ordination scope: ", PCOA_SCOPE)) +
  theme_bw(base_family = FONT_FAMILY, base_size = 10) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "grey55", fill = NA),
    strip.background = element_rect(fill = "grey95", colour = "grey55"),
    strip.text = element_text(face = "bold", size = 11),
    legend.position = "bottom",
    legend.text = element_text(size = 9),
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.subtitle = element_text(size = 8, hjust = 0.5),
    axis.title = element_text(size = 10),
    plot.margin = margin(10, 12, 8, 12)
  )
# Coordinate aspect ratio is meaningful for a shared (combined) ordination;
# avoid coord_fixed with free facets in the separate-ordination option.
if (PCOA_SCOPE == "combined") {
  pcoa_scatter <- pcoa_scatter + coord_fixed(ratio = 1)
}

pcoa_2d_pdf <- file.path(OUT_DIR, "SupplementaryFig2_PCoA1_vs_PCoA2_2D.pdf")
pcoa_2d_png <- file.path(OUT_DIR, "SupplementaryFig2_PCoA1_vs_PCoA2_2D.png")
pcoa_2d_svg <- file.path(OUT_DIR, "SupplementaryFig2_PCoA1_vs_PCoA2_2D.svg")

ggsave(pcoa_2d_pdf, pcoa_scatter, width = PCOA_2D_WIDTH,
       height = PCOA_2D_HEIGHT, units = "in", device = grDevices::cairo_pdf,
       bg = "white", limitsize = FALSE)
ggsave(pcoa_2d_png, pcoa_scatter, width = PCOA_2D_WIDTH,
       height = PCOA_2D_HEIGHT, units = "in", dpi = FIG_DPI,
       bg = "white", limitsize = FALSE)
if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(pcoa_2d_svg, pcoa_scatter, width = PCOA_2D_WIDTH,
         height = PCOA_2D_HEIGHT, units = "in", device = svglite::svglite,
         bg = "white", limitsize = FALSE)
}

# Convenient optional 2-page PDF: original violin panels on page 1,
# 2D Bray-Curtis PCoA on page 2. Original 1-page PDF is still preserved.
extended_pdf <- file.path(OUT_DIR, "SupplementaryFig2_extended_2pages.pdf")
grDevices::cairo_pdf(extended_pdf, width = FIG_WIDTH, height = FIG_HEIGHT,
                     onefile = TRUE, bg = "white")
tryCatch({
  print(fig)
  print(pcoa_scatter)
}, finally = grDevices::dev.off())

message("Additional 2D PCoA PDF: ", pcoa_2d_pdf)
message("Study-adjusted PERMANOVA CSV: ",
        file.path(OUT_DIR, "SuppFig2_PERMANOVA_study_adjusted.csv"))
message("Optional 2-page combined PDF: ", extended_pdf)
