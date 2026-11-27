#!/usr/bin/env Rscript
# ============================================================================
# OMWI | Supplementary Fig. 1 | TRUE leave-one-study-out coefficient heatmap
# Reuse the 17 locally regenerated LOSO fold selections and tuned alpha/gamma; DO NOT redo inner CV or stability rank.
# Refit one glmnet path per held-out study to recover signed coefficients.
# Strictly verify reconstructed held-out predictions against saved RDS results.
# Authoring specification: raw A, panel-only TSS + CLR, no MMUPHin.
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(glmnet)
  library(ggplot2)
})
options(stringsAsFactors = FALSE)

# 0. Windows paths: edit ROOT_DIR only if you moved the project.
ROOT_DIR <- Sys.getenv("OMWI_ROOT", getwd())
OUT_DIR <- file.path(ROOT_DIR, "Supplementaryfigure", "bac")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# These files are NOT contained in the saved stability ZIP. They must be present
# You can set each directory explicitly if your raw A inputs live elsewhere.
RAW_DIR_OVERRIDE <- file.path(ROOT_DIR, "data/abundance")
META_DIR_OVERRIDE <- file.path(ROOT_DIR, "data/metadata")
MODEL_DIR_OVERRIDE <- file.path(ROOT_DIR, "results/model_refit")

# Validate that re-fitted predictions agree with the original LOSO RDS results.
MAX_PREDICTION_DIFFERENCE <- 1e-4
COEF_ZERO_TOL <- 1e-12
TOP_PER_DIRECTION <- 25L   # 50-taxa layout like the old supplementary figure
COL_HEALTH <- "#4DBBD5"
COL_NONHEALTH <- "#D95F59"
COL_GREY <- "#D0D0D0"

# IMPORTANT: match source code model_development.R exactly.
GAMMAS <- exp(seq(log(1), log(0.025), length.out = 16L))
CLR_PC_PERCENT <- 1e-6
MIN_BATCH_N <- 2L

# ---- Inputs (no silent substitution of filtered vs original raw data) ------
choose_dir <- function(explicit, subfolder) {
  if (!is.na(explicit) && nzchar(explicit)) {
    if (!dir.exists(explicit)) stop("Configured directory missing: ", explicit)
    return(explicit)
  }
  options <- c(file.path(ROOT_DIR, subfolder),
               file.path(ROOT_DIR, "oral_proxy", subfolder))
  found <- options[dir.exists(options)]
  if (!length(found)) {
    stop("Input directory not found: ", subfolder, "\nSearched:\n  ",
         paste(options, collapse = "\n  "),
         "\nSet RAW_DIR_OVERRIDE / META_DIR_OVERRIDE / MODEL_DIR_OVERRIDE as appropriate.")
  }
  normalizePath(found[1L], winslash = "/", mustWork = TRUE)
}

MODEL_DIR <- choose_dir(MODEL_DIR_OVERRIDE, "data/model_outputs")
RAW_DIR <- choose_dir(RAW_DIR_OVERRIDE, "data/abundance")
META_DIR <- choose_dir(META_DIR_OVERRIDE, "data/metadata")

META_FILE <- file.path(META_DIR, "meta_train.csv")
TAXA_FILE <- file.path(RAW_DIR, "train_species.tsv")
FINAL_COEF_FILE <- file.path(MODEL_DIR, "final", "species_all_coefficients.csv")
for (f in c(META_FILE, TAXA_FILE, FINAL_COEF_FILE))
  if (!file.exists(f)) stop("Required file missing: ", f,
                             "\nDo not substitute a different model version's data.")
fold_files <- sort(list.files(file.path(MODEL_DIR, "loso_folds"),
                        pattern = "^study_[0-9]+[.]rds$", full.names = TRUE))
if (!length(fold_files)) stop("No saved LOSO fold RDS objects found.")

# ---- Exact raw-A abundance reader; samples x taxa --------------------------
meta <- as.data.frame(data.table::fread(META_FILE, check.names = FALSE))
for (v in c("sample_id", "study_name", "disease"))
  if (!v %in% names(meta)) stop("Missing meta_train column: ", v)
meta$sample_id <- trimws(as.character(meta$sample_id))
meta$study_name <- trimws(as.character(meta$study_name))
meta$disease <- tolower(trimws(as.character(meta$disease)))
if (anyNA(meta[, c("sample_id", "study_name", "disease")]) ||
    anyDuplicated(meta$sample_id) ||
    any(!nzchar(meta$sample_id)) || any(!nzchar(meta$study_name)))
  stop("Invalid/duplicated training sample or study identifiers")
meta$y <- as.integer(meta$disease == "healthy")
if (length(unique(meta$y)) != 2L) stop("Both outcome classes are required")

raw <- data.table::fread(TAXA_FILE, data.table = FALSE, check.names = FALSE)
if (ncol(raw) < 2L) stop("Invalid train_species.tsv")
names(raw)[1L] <- "Taxon"
raw$Taxon <- trimws(as.character(raw$Taxon))
names(raw)[-1L] <- trimws(names(raw)[-1L])
if (anyNA(raw$Taxon) || any(!nzchar(raw$Taxon)) || anyDuplicated(names(raw)[-1L]))
  stop("Duplicate/missing raw taxa/sample identifiers")
missing_samples <- setdiff(meta$sample_id, names(raw)[-1L])
if (length(missing_samples))
  stop("Samples missing from train_species.tsv: ", paste(head(missing_samples), collapse = ", "))
xa <- as.matrix(raw[, meta$sample_id, drop = FALSE])
storage.mode(xa) <- "double"
rownames(xa) <- raw$Taxon
if (anyNA(xa) || any(!is.finite(xa)) || any(xa < 0))
  stop("Raw matrix contains NA, infinite or negative abundances")
if (anyDuplicated(rownames(xa))) xa <- rowsum(xa, group = rownames(xa), reorder = FALSE)
x <- t(xa)
rm(raw, xa)
if (!identical(rownames(x), meta$sample_id)) stop("X / meta sample order mismatch")

# ---- Exact preprocessing and original glmnet path --------------------------
tss_raw_A <- function(z) {
  den <- rowSums(z)
  if (any(!is.finite(den))) stop("Nonfinite panel sum")
  ans <- z
  good <- den > 0
  if (any(good)) ans[good, ] <- sweep(z[good, , drop = FALSE], 1L, den[good], "/")
  if (any(!good)) ans[!good, ] <- 0
  ans
}
clr_raw_A <- function(z) {
  if (any(!is.finite(z)) || any(z < 0)) stop("CLR requires nonnegative finite abundance")
  logz <- log(z * 100 + CLR_PC_PERCENT)
  sweep(logz, 1L, rowMeans(logz), FUN = "-")
}
drop_singleton <- function(xtr, mtr) {
  if (!identical(rownames(xtr), mtr$sample_id)) stop("Fold sample ordering mismatch")
  tb <- table(mtr$study_name)
  keep <- mtr$study_name %in% names(tb[tb >= MIN_BATCH_N])
  if (sum(keep) < 10L || length(unique(mtr$y[keep])) != 2L)
    stop("Invalid fold after singleton-study exclusion")
  list(x = xtr[keep, , drop = FALSE], m = mtr[keep, , drop = FALSE],
       dropped = sum(!keep))
}

refit_saved_fold <- function(study, saved) {
  tr <- which(meta$study_name != study)
  te <- which(meta$study_name == study)
  if (!length(te)) stop("Study absent from original train metadata: ", study)
  if (is.null(saved$panel) || is.null(saved$chosen) ||
      !all(c("K", "alpha", "gamma") %in% names(saved$chosen)))
    stop("RDS missing trained panel/alpha/gamma: ", study)
  panel <- as.character(saved$panel)
  if (anyDuplicated(panel) || length(panel) != as.integer(saved$chosen$K[1L]))
    stop("RDS panel length/K mismatch: ", study)
  if (!all(panel %in% colnames(x)))
    stop("Original raw abundance file lacks some chosen panel taxa: ", study)
  if (!setequal(as.character(saved$predictions$sample_id), meta$sample_id[te]))
    stop("Held-out sample IDs disagree with archived LOSO result: ", study)
  d <- drop_singleton(x[tr, , drop = FALSE], meta[tr, , drop = FALSE])
  ztr <- clr_raw_A(tss_raw_A(d$x[, panel, drop = FALSE]))
  zte <- clr_raw_A(tss_raw_A(x[te, panel, drop = FALSE]))
  alpha <- as.numeric(saved$chosen$alpha[1L])
  gamma <- as.numeric(saved$chosen$gamma[1L])
  if (!is.finite(alpha) || !is.finite(gamma)) stop("Invalid saved parameters: ", study)
  gi <- which.min(abs(GAMMAS - gamma))
  if (abs(GAMMAS[gi] - gamma) > 1e-8) stop("Gamma grid differs from original script: ", study)

  # Identical to fit_path() in your model_development.R.
  prelim <- glmnet::glmnet(ztr, d$m$y, family = "binomial", alpha = alpha,
                           standardize = TRUE, nlambda = 5L)
  lams <- as.numeric(prelim$lambda[1L] * GAMMAS)
  fit <- glmnet::glmnet(ztr, d$m$y, family = "binomial", alpha = alpha,
                        standardize = TRUE, lambda = lams)
  available <- lams >= min(fit$lambda) * (1 - 1e-8)
  if (!available[gi]) stop("Selected lambda unavailable in refit path: ", study)
  chosen_lambda <- lams[gi]
  pred_prob <- as.numeric(predict(fit, newx = zte, s = chosen_lambda, type = "response"))
  pred_logit <- as.numeric(predict(fit, newx = zte, s = chosen_lambda, type = "link"))
  r <- saved$predictions
  saved_order <- match(meta$sample_id[te], as.character(r$sample_id))
  p_diff <- max(abs(pred_prob - as.numeric(r$prob_healthy[saved_order])))
  logit_diff <- max(abs(pred_logit - as.numeric(r$OMHI[saved_order])))
  if (!is.finite(p_diff) || !is.finite(logit_diff)) stop("Missing/nonfinite parity values: ", study)
  audit <- data.frame(study_name = study, n_holdout = length(te), K = length(panel),
                      alpha = alpha, gamma = gamma, lambda = chosen_lambda,
                      max_prob_diff = p_diff, max_logit_diff = logit_diff,
                      n_train_used = nrow(ztr), dropped_train = d$dropped,
                      stringsAsFactors = FALSE)
  if (p_diff > MAX_PREDICTION_DIFFERENCE || logit_diff > MAX_PREDICTION_DIFFERENCE)
    stop(sprintf("Prediction parity FAILED for %s (prob %.6g, logit %.6g).\nCheck raw input files and glmnet version.",
                 study, p_diff, logit_diff))
  if ("n_train" %in% names(r) && any(r$n_train != nrow(ztr)))
    stop("Original RDS training count differs: ", study)
  if ("dropped_train" %in% names(r) && any(r$dropped_train != d$dropped))
    stop("Original RDS dropped-train count differs: ", study)

  cc <- as.matrix(stats::coef(fit, s = chosen_lambda))
  values <- setNames(as.numeric(cc[-1L, 1L]), rownames(cc)[-1L])
  # A zero coefficient includes panel taxa that entered K but were shrunk to 0.
  list(coef = values, panel = panel, audit = audit)
}

# ---- Read saved LOSO folds and reconstruct exactly their coefficients -------
expected_n <- length(unique(meta$study_name))
if (length(fold_files) != expected_n)
  stop(sprintf("LOSO archive has %d fold RDS, metadata has %d studies.\nDo not silently draw partial results.",
               length(fold_files), expected_n))

results <- vector("list", length(fold_files))
names(results) <- basename(fold_files)
fold_studies <- character(length(fold_files))
audit_rows <- vector("list", length(fold_files))
for (j in seq_along(fold_files)) {
  obj <- readRDS(fold_files[j])
  st <- unique(as.character(obj$predictions$study_name))
  if (length(st) != 1L || is.na(st) || !nzchar(st))
    stop("RDS does not contain a unique study: ", fold_files[j])
  fold_studies[j] <- st
  if (anyDuplicated(fold_studies[seq_len(j)])) stop("Duplicate LOSO held-out study: ", st)
  message(sprintf("Refitting saved LOSO %02d/%02d: %s", j, length(fold_files), st))
  one <- refit_saved_fold(st, obj)
  results[[j]] <- one
  audit_rows[[j]] <- one$audit
}
if (!setequal(fold_studies, unique(meta$study_name))) stop("LOSO study coverage is incomplete")
# Order by the original study name, as in the old supplementary plot.
ix <- order(fold_studies)
fold_studies <- fold_studies[ix]
results <- results[ix]
audit <- data.table::rbindlist(audit_rows, fill = TRUE)
audit <- audit[match(fold_studies, audit$study_name), ]
data.table::fwrite(audit, file.path(OUT_DIR, "SuppFig3_LOSO_refit_prediction_parity.csv"))

# Read FROZEN final coefficients (do not refit the final model or mix versions).
final <- as.data.frame(data.table::fread(FINAL_COEF_FILE, check.names = FALSE))
if (!all(c("feature", "coefficient") %in% names(final)))
  stop("Frozen coefficient table must include feature, coefficient")
final <- final[final$feature != "(Intercept)", c("feature", "coefficient")]
final$coefficient <- as.numeric(final$coefficient)
if (anyDuplicated(final$feature) || anyNA(final$coefficient)) stop("Invalid final coefficients")
final <- final[abs(final$coefficient) > COEF_ZERO_TOL, , drop = FALSE]
if (nrow(final) < 2L) stop("Too few final nonzero species")

full_coef <- data.frame(feature = final$feature, final_coefficient = final$coefficient,
                        stringsAsFactors = FALSE)
for (j in seq_along(results)) {
  st <- fold_studies[j]
  co <- results[[j]]$coef
  full_coef[[st]] <- as.numeric(co[match(final$feature, names(co))])
  full_coef[[st]][is.na(full_coef[[st]])] <- 0
}
full_coef$`Final model` <- full_coef$final_coefficient
# Save the EXACT numeric per-study coefficients, not only sign categories.
data.table::fwrite(full_coef, file.path(OUT_DIR, "SuppFig3_LOSO_species_coefficients.csv"))

# Per-taxon direction stability using 17 real refitted models.
coef_mat <- as.matrix(full_coef[, fold_studies, drop = FALSE])
ref_sign <- sign(full_coef$final_coefficient)
selected <- abs(coef_mat) > COEF_ZERO_TOL
same_sign <- selected & sweep(sign(coef_mat), 1L, ref_sign, `==`)
flip_sign <- selected & !same_sign
frequency <- data.frame(feature = full_coef$feature,
  final_coefficient = full_coef$final_coefficient,
  n_loso_models = length(fold_studies),
  n_selected = rowSums(selected),
  selection_frequency = rowMeans(selected),
  n_same_sign = rowSums(same_sign),
  n_opposite_sign = rowSums(flip_sign),
  sign_agreement_given_selected = ifelse(rowSums(selected) > 0,
                                         rowSums(same_sign) / rowSums(selected), NA_real_),
  stringsAsFactors = FALSE)
data.table::fwrite(frequency, file.path(OUT_DIR, "SuppFig3_LOSO_taxon_stability_summary.csv"))

# Signed status from each ACTUAL LOSO refit (NOT final-coefficient imputation).
plot_long <- do.call(rbind, lapply(c(fold_studies, "Final model"), function(st) {
  val <- full_coef[[st]]
  data.frame(feature = full_coef$feature, fold = st,
             coefficient = val,
             status = ifelse(val > COEF_ZERO_TOL, "Health-associated",
                      ifelse(val < -COEF_ZERO_TOL, "Non-healthy-associated", "Not selected")),
             stringsAsFactors = FALSE)
}))
data.table::fwrite(plot_long, file.path(OUT_DIR, "SuppFig3_LOSO_heatmap_long.csv"))

# Arrange rows to match the old figure: negative at TOP, positive at BOTTOM.
# Put strongest positive at bottom and strongest negative at top.
choose_order <- function(df) {
  p <- df$feature[df$final_coefficient > 0]
  n <- df$feature[df$final_coefficient < 0]
  p <- p[order(-df$final_coefficient[match(p, df$feature)])]
  n <- n[order(-df$final_coefficient[match(n, df$feature)])]
  c(p, n)   # factor levels start from bottom of y-axis
}
fold_levels <- c(fold_studies, "Final model")
palette <- c("Health-associated" = COL_HEALTH,
             "Non-healthy-associated" = COL_NONHEALTH,
             "Not selected" = COL_GREY)

save_plot <- function(taxa, stem, width, height, text_size) {
  z <- plot_long[plot_long$feature %in% taxa, , drop = FALSE]
  z$feature <- factor(z$feature, levels = choose_order(full_coef[full_coef$feature %in% taxa, ]))
  z$fold <- factor(z$fold, levels = fold_levels)
  z$status <- factor(z$status, levels = names(palette))
  p <- ggplot(z, aes(x = fold, y = feature, fill = status)) +
    geom_tile(width = 0.98, height = 0.98, color = "white", linewidth = 0.22) +
    scale_fill_manual(values = palette, drop = FALSE, name = NULL) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(x = "Excluded training study", y = NULL) +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme_classic(base_size = 10) +
    theme(panel.grid = element_blank(), axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_text(angle = 65, hjust = 1, vjust = 1,
                                     size = 8, color = "black"),
          axis.text.y = element_text(size = text_size, color = "black"),
          axis.title.x = element_text(size = 11, margin = margin(t = 9)),
          legend.position = "top", legend.direction = "horizontal",
          legend.text = element_text(size = 9),
          legend.key.width = grid::unit(0.35, "cm"),
          legend.margin = margin(b = 5),
          plot.margin = margin(5, 7, 6, 5))
  ggsave(file.path(OUT_DIR, paste0(stem, ".pdf")), p,
         width = width, height = height, units = "in", device = "pdf",
         useDingbats = FALSE, bg = "white")
  ggsave(file.path(OUT_DIR, paste0(stem, ".png")), p,
         width = width, height = height, units = "in", dpi = 600,
         limitsize = FALSE, bg = "white")
  invisible(p)
}

# Primary faithful layout: 25 top negative + 25 top positive final coefficients.
pos <- full_coef[full_coef$final_coefficient > 0, , drop = FALSE]
neg <- full_coef[full_coef$final_coefficient < 0, , drop = FALSE]
pos <- pos[order(-pos$final_coefficient), , drop = FALSE]
neg <- neg[order(neg$final_coefficient), , drop = FALSE]
top_taxa <- c(head(pos$feature, TOP_PER_DIRECTION),
              head(neg$feature, TOP_PER_DIRECTION))
save_plot(top_taxa, "SuppFig3_LOSO_50taxa_like_original", 10.5, 10.2, 8.4)

# Exhaustive supplemental audit: ALL final nonzero species.
all_taxa <- full_coef$feature
save_plot(all_taxa, "SuppFig3_LOSO_all_final_nonzero_taxa", 11, 15.2, 7.2)

message("DONE. Strictly verified ", nrow(audit), " LOSO refits against saved predictions.")
message("Final model nonzero taxa: ", nrow(final),
        "; figure compact taxa: ", length(top_taxa))
message("Output: ", OUT_DIR)
