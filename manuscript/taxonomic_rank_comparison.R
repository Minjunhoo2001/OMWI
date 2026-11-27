#!/usr/bin/env Rscript
# OMWI | Six-rank RAW taxonomic benchmark, NO batch correction of any kind.
# Six independent inputs: phylum/class/order/family/genus/species; no SGB.
# NO external cohort or external outcome is loaded by this script.
#
# Modes: preflight | cv REPEAT FOLD | loso STUDY_INDEX | aggregate | final-fit
# Outer: 10 repeated stratified 10-fold CV + study-wise LOSO.
# Inner: stratified 5-fold CV to tune alpha/prevalence/lambda_max ratio.
# Fit feature filtering, TSS and CLR entirely within training partitions.
# only for direct apples-to-apples compatibility; this is NOT batch correction.

suppressPackageStartupMessages({
  library(data.table)
  library(glmnet)
  library(pROC)
})
options(stringsAsFactors = FALSE)

ROOT <- Sys.getenv("OMWI_ROOT", getwd())
RAW_DIR <- Sys.getenv("OMWI_RAW_DIR", file.path(ROOT, "data", "abundance"))
META_DIR <- Sys.getenv("OMWI_META_DIR", file.path(ROOT, "data", "metadata"))
RANKS <- c("phylum", "class", "order", "family", "genus", "species")
RANK <- tolower(trimws(Sys.getenv("OMWI_RANK", "")))
if (!RANK %in% RANKS) stop("Set OMWI_RANK to one of: ", paste(RANKS, collapse=", "))
BENCH_ROOT <- Sys.getenv("OMWI_RANK_RAW_OUT", file.path(ROOT, "results/rank_benchmark"))
OUT_DIR <- file.path(BENCH_ROOT, RANK)
protected <- c(file.path(ROOT, "05.raw_A_species"), file.path(ROOT, "05.nested_A_species"))
if (normalizePath(OUT_DIR, mustWork=FALSE) %in%
    vapply(protected, normalizePath, character(1L), mustWork=FALSE))
  stop("Refusing to write benchmark outputs inside an existing model directory")
SEED <- 704L
N_REPEATS <- 10L
N_OUTER <- 10L
N_INNER <- 5L
ALPHAS <- seq(0.1, 1.0, by = 0.1)
PREVS <- c(0.01, 0.05, 0.10, 0.20, 0.30)
# lambda is selected as a RATIO to each fit's lambda_max; avoids estimating
# a common absolute lambda from samples held out in each INNER fold.
GAMMAS <- exp(seq(log(1), log(0.02), length.out = 24L))
CLR_PC_PERCENT <- 1e-6
MIN_TAXA <- 2L
MIN_BATCH_N <- 2L

for (p in c("cv_folds", "loso_folds", "summary", "final")) {
  dir.create(file.path(OUT_DIR, p), recursive = TRUE, showWarnings = FALSE)
}

save_csv <- function(x, path) {
  if (file.exists(path)) stop("Refusing overwrite: ", path)
  tmp <- paste0(path, ".", Sys.getpid(), ".tmp")
  data.table::fwrite(x, tmp)
  if (!file.rename(tmp, path)) stop("Cannot rename output to ", path)
}
save_rds <- function(x, path) {
  if (file.exists(path)) stop("Refusing overwrite: ", path)
  tmp <- paste0(path, ".", Sys.getpid(), ".tmp")
  saveRDS(x, tmp, compress = "gzip")
  if (!file.rename(tmp, path)) stop("Cannot rename output to ", path)
}
read_raw_taxa <- function(path, sample_ids, missing_taxa = NULL) {
  if (!file.exists(path)) stop("Missing raw taxa table: ", path)
  d <- data.table::fread(path, data.table = FALSE, check.names = FALSE)
  if (ncol(d) < 2L) stop("Invalid taxa table: ", path)
  names(d)[1L] <- "Taxon"
  d$Taxon <- trimws(as.character(d$Taxon))
  if (anyNA(d$Taxon) || any(!nzchar(d$Taxon))) stop("Invalid taxa IDs: ", path)
  names(d) <- c("Taxon", trimws(names(d)[-1L]))
  if (anyDuplicated(names(d)[-1L])) stop("Duplicate raw sample column names")
  absent_samples <- setdiff(sample_ids, names(d))
  if (length(absent_samples)) {
    sample_ids_first_column <- sum(sample_ids %in% d$Taxon)
    stop("TRAIN INPUT MISMATCH [rank=", RANK, "]\n",
         "File read: ", normalizePath(path, mustWork = FALSE), "\n",
         "Metadata: ", file.path(META_DIR, "meta_train.csv"), "\n",
         "Metadata samples: ", length(sample_ids),
         "; sample columns in abundance file: ", ncol(d) - 1L, "\n",
         "Missing sample IDs: ", length(absent_samples), "\n",
         "First missing IDs: ", paste(head(absent_samples, 8L), collapse = ", "), "\n",
         "Found in first column (possible transposed input): ", sample_ids_first_column, "\n",
         "Check OMWI_RAW_DIR and file header; no samples have been removed. ",
         "Run: bash submit_rank_raw.sh check-inputs")
  }
  x <- as.matrix(d[, sample_ids, drop = FALSE])
  storage.mode(x) <- "double"
  rownames(x) <- d$Taxon
  if (anyNA(x) || any(!is.finite(x)) || any(x < 0)) stop("Invalid abundance values")
  if (anyDuplicated(rownames(x))) {
    x <- rowsum(x, group = rownames(x), reorder = FALSE)
  }
  if (!is.null(missing_taxa)) {
    # During external prediction, absent taxa are treated as undetected (0).
    added <- setdiff(missing_taxa, rownames(x))
    if (length(added)) {
      z <- matrix(0, nrow = length(added), ncol = ncol(x),
                  dimnames = list(added, colnames(x)))
      x <- rbind(x, z)
      warning(length(added), " expected species not listed in external input; set to zero")
    }
    x <- x[missing_taxa, , drop = FALSE]
  }
  x <- t(x) # samples x taxa
  if (!identical(rownames(x), sample_ids)) stop("Sample ordering error")
  x
}
read_training <- function() {
  path <- file.path(META_DIR, "meta_train.csv")
  if (!file.exists(path)) stop("Missing training metadata: ", path)
  m <- data.table::fread(path, data.table = FALSE)
  for (v in c("sample_id", "study_name", "disease")) {
    if (!v %in% names(m)) stop("Missing training metadata variable: ", v)
  }
  m$sample_id <- trimws(as.character(m$sample_id))
  m$study_name <- trimws(as.character(m$study_name))
  m$disease <- trimws(tolower(as.character(m$disease)))
  if (anyNA(m[, c("sample_id", "study_name", "disease")]) ||
      any(!nzchar(m$sample_id)) || any(!nzchar(m$study_name)) ||
      any(!nzchar(m$disease)) || anyDuplicated(m$sample_id)) {
    stop("Missing / duplicate IDs or labels in training metadata")
  }
  m$y <- as.integer(m$disease == "healthy")
  m$health_bin <- factor(ifelse(m$y == 1L, "Healthy", "Nonhealthy"),
                         levels = c("Healthy", "Nonhealthy"))
  if (length(unique(m$y)) < 2L) stop("Both classes required")
  cat("[INPUT] rank=", RANK, " raw_dir=", normalizePath(RAW_DIR, mustWork = FALSE),
      " metadata=", normalizePath(path, mustWork = FALSE), "\n", sep = "")
  x <- read_raw_taxa(file.path(RAW_DIR, paste0("train_", RANK, ".tsv")), m$sample_id)
  list(X = x, meta = m)
}

make_stratified <- function(y, k, seed) {
  k <- min(as.integer(k), sum(y == 0L), sum(y == 1L))
  if (k < 2L) stop("Not enough samples per class for inner CV")
  set.seed(seed)
  f <- integer(length(y))
  for (cl in c(0L, 1L)) {
    ids <- which(y == cl)
    f[ids] <- sample(rep(seq_len(k), length.out = length(ids)))
  }
  f
}

safe_auc <- function(y, p) {
  if (length(unique(y)) < 2L || anyNA(p)) return(NA_real_)
  as.numeric(pROC::roc(response = y, predictor = p,
                       levels = c(0, 1), direction = "<", quiet = TRUE)$auc)
}
get_metrics <- function(y, p, link = NULL) {
  pred <- as.integer(p >= 0.5)
  tpr <- if (any(y == 1L)) mean(pred[y == 1L] == 1L) else NA_real_
  tnr <- if (any(y == 0L)) mean(pred[y == 0L] == 0L) else NA_real_
  data.frame(n = length(y), n_healthy = sum(y == 1L), n_nonhealthy = sum(y == 0L),
             auc = safe_auc(y, p), balanced_accuracy = mean(c(tpr, tnr)),
             sensitivity = tpr, specificity = tnr)
}
clr <- function(x) {
  if (any(x < 0) || any(!is.finite(x))) stop("CLR received invalid values")
  z <- log(x * 100 + CLR_PC_PERCENT) # percentage units, same on both sides
  sweep(z, 1, rowMeans(z), "-")
}

# TSS is computed for EACH sample independently. A zero held-out sample
# (among training-selected available taxa) is preserved as all-zero.
tss <- function(x, allow_zero = FALSE) {
  rs <- rowSums(x)
  if (any(!is.finite(rs))) stop("Invalid row sums")
  if (!allow_zero && any(rs <= 0)) stop("Zero-sum training sample")
  z <- x
  valid <- rs > 0
  z[valid, ] <- sweep(x[valid, , drop = FALSE], 1, rs[valid], "/")
  if (any(!valid)) z[!valid, ] <- 0
  z
}

# Raw-only preprocessing. No batch estimation, correction, or combat.
# Preserve original raw arm's singleton-study exclusion and TSS/CLR policy.
# Coarse ranks may have fewer than 10 variable taxa, so minimum is 2.
# Held-out samples never influence any fitted transformation.
preprocess_train <- function(Xtrain, meta_train, Xholdout = NULL) {
  if (!identical(rownames(Xtrain), meta_train$sample_id)) stop("Training alignment mismatch")
  if (!is.null(Xholdout) && !identical(colnames(Xtrain), colnames(Xholdout))) {
    stop("Train/holdout taxa columns must match")
  }
  counts <- table(meta_train$study_name)
  use <- meta_train$study_name %in% names(counts[counts >= MIN_BATCH_N])
  omitted <- sum(!use)
  if (omitted) message("Matching the corrected arm: dropped ", omitted,
                       " singleton-study training sample(s)")
  m <- meta_train[use, , drop = FALSE]
  a <- Xtrain[use, , drop = FALSE]
  if (nrow(a) < 10L || length(unique(m$y)) < 2L)
    stop("Too few train samples / classes")
  variances <- apply(a, 2, stats::var)
  pre_taxa <- names(variances)[is.finite(variances) & variances > 0]
  if (length(pre_taxa) < MIN_TAXA) stop("Too few varying taxa (need >=2)")
  raw_a <- a[, pre_taxa, drop = FALSE]
  raw_h <- if (!is.null(Xholdout)) Xholdout[, pre_taxa, drop = FALSE] else NULL
  prev <- colMeans(raw_a > 0)
  train_tss <- tss(raw_a)
  hold_tss <- if (!is.null(raw_h)) tss(raw_h, allow_zero = TRUE) else NULL
  # Raw TSS values only (not adjusted).
  list(adj = train_tss, held_raw = hold_tss, prev = prev,
       meta_used = m, pre_taxa = pre_taxa, n_dropped = omitted)
}

get_xy <- function(pre, prev_cut, skip_if_small = FALSE) {
  variation <- apply(pre$adj, 2, stats::var)
  keep <- names(pre$prev)[pre$prev >= prev_cut &
                            is.finite(variation) & variation > 0]
  if (length(keep) < MIN_TAXA) {
    if (skip_if_small) return(NULL)
    stop("Too few features at prevalence ", prev_cut)
  }
  list(x = clr(pre$adj[, keep, drop = FALSE]),
       x_hold = if (!is.null(pre$held_raw))
         clr(pre$held_raw[, keep, drop = FALSE]) else NULL,
       y = pre$meta_used$y, taxa = keep)
}

# Fit regularization path. Gamma is a dimensionless ratio to lambda_max,
# estimated from this model's TRAINING fold only.
fit_path <- function(x, y, alpha) {
  # A preliminary TRAIN-only fit identifies lambda_max. The supplied full
  # lambda sequence then guarantees the identical gamma grid in every split.
  prelim <- glmnet::glmnet(x, y, family = "binomial", alpha = alpha,
                           standardize = TRUE, nlambda = 5L)
  s <- as.numeric(prelim$lambda[1L] * GAMMAS)
  f <- glmnet::glmnet(x, y, family = "binomial", alpha = alpha,
                      standardize = TRUE, lambda = s)
  if (length(f$lambda) != length(s) ||
      max(abs(f$lambda - s)) > 1e-8 * max(s)) {
    stop("glmnet did not fit the requested lambda grid")
  }
  list(fit = f, s = s, valid = rep(TRUE, length(s)))
}

# Inner 5-fold: preprocess each inner TRAIN fold; held-out part is never
# included in feature filtering, TSS/CLR feature selection or lambda_max.
# The same 5 preprocessing outputs are reused across all 50 combinations.
tune_inner <- function(X, meta, seed) {
  foldid <- make_stratified(meta$y, N_INNER, seed)
  k <- max(foldid)
  pre <- lapply(seq_len(k), function(j) {
    te <- which(foldid == j); tr <- which(foldid != j)
    preprocess_train(X[tr, , drop = FALSE], meta[tr, , drop = FALSE],
                     X[te, , drop = FALSE])
  })
  candidates <- vector("list", length(PREVS) * length(ALPHAS))
  ix <- 0L
  for (pc in PREVS) {
    xs <- lapply(pre, function(z) get_xy(z, pc, skip_if_small = TRUE))
    # Very coarse ranks may have too few taxa at high prevalence thresholds.
    # Skip such hyperparameter settings in ALL folds, never drop a fold.
    if (any(vapply(xs, is.null, logical(1)))) next
    for (alpha in ALPHAS) {
      ix <- ix + 1L
      scores <- matrix(NA_real_, nrow = k, ncol = length(GAMMAS))
      for (j in seq_len(k)) {
        p <- fit_path(xs[[j]]$x, xs[[j]]$y, alpha)
        valid <- which(p$valid)
        if (!length(valid)) next
        prob <- as.matrix(predict(p$fit, newx = xs[[j]]$x_hold,
                                  s = p$s[valid], type = "response"))
        y_held <- meta$y[which(foldid == j)]
        for (h in seq_along(valid)) {
          scores[j, valid[h]] <- safe_auc(y_held, prob[, h])
        }
      }
      complete <- colSums(is.finite(scores)) == k
      if (!any(complete)) next
      means <- colMeans(scores)
      means[!complete] <- -Inf
      best_g <- which.max(means)
      # Analogue of lambda.1se: largest penalty with mean AUC within 1 SE
      se <- stats::sd(scores[, best_g]) / sqrt(k)
      if (!is.finite(se)) se <- 0
      acceptable <- which(complete & means >= means[best_g] - se)
      chosen_g <- acceptable[which.max(GAMMAS[acceptable])]
      candidates[[ix]] <- data.frame(alpha = alpha, prevalence = pc,
                                      gamma = GAMMAS[chosen_g],
                                      inner_auc = means[chosen_g],
                                      best_inner_auc = means[best_g],
                                      inner_se_at_best = se,
                                      inner_folds = k)
    }
  }
  candidates <- Filter(Negate(is.null), candidates)
  if (!length(candidates)) stop("No valid nested tuning candidates for ", RANK)
  c <- do.call(rbind, candidates)
  # alpha / prevalence selected exclusively on inner OOF scores; ties parsimonious.
  c <- c[order(-c$inner_auc, -c$prevalence, -c$gamma, -c$alpha), , drop = FALSE]
  rownames(c) <- NULL
  c[1L, , drop = FALSE]
}

fit_selected <- function(Xtr, mtr, Xte, selected) {
  pre <- preprocess_train(Xtr, mtr, Xte)
  xy <- get_xy(pre, selected$prevalence)
  p <- fit_path(xy$x, xy$y, selected$alpha)
  wanted <- p$fit$lambda[1L] * selected$gamma
  if (wanted < min(p$fit$lambda) * (1 - 1e-8)) {
    stop("Selected gamma outside full-refit model path")
  }
  prob <- if (!is.null(Xte)) as.numeric(predict(p$fit, newx = xy$x_hold,
                                                s = wanted, type = "response")) else NULL
  link <- if (!is.null(Xte)) as.numeric(predict(p$fit, newx = xy$x_hold,
                                                s = wanted, type = "link")) else NULL
  train_prob <- as.numeric(predict(p$fit, newx = xy$x, s = wanted,
                                    type = "response"))
  train_link <- as.numeric(predict(p$fit, newx = xy$x, s = wanted,
                                    type = "link"))
  list(prob = prob, link = link, model = p$fit, lambda = wanted,
       selected_taxa = xy$taxa, pre_taxa = pre$pre_taxa,
       train_sample_ids = pre$meta_used$sample_id,
       train_prob = train_prob, train_link = train_link,
       n_train = nrow(xy$x), n_train_dropped = pre$n_dropped,
       n_features = length(xy$taxa), corrected = pre$adj)
}

run_outer_cv <- function(repeat_id, fold_id) {
  if (!repeat_id %in% seq_len(N_REPEATS) || !fold_id %in% seq_len(N_OUTER))
    stop("Invalid repeat/fold ID")
  dest <- file.path(OUT_DIR, "cv_folds", sprintf("r%02d_f%02d.csv", repeat_id, fold_id))
  if (file.exists(dest)) stop("Already exists: ", dest)
  d <- read_training(); X <- d$X; m <- d$meta
  folds <- make_stratified(m$y, N_OUTER, SEED + repeat_id * 1000L)
  te <- which(folds == fold_id); tr <- which(folds != fold_id)
  if (!length(te)) stop("Empty outer test fold")
  choice <- tune_inner(X[tr, , drop = FALSE], m[tr, , drop = FALSE],
                       SEED + repeat_id * 1000L + fold_id * 23L)
  fit <- fit_selected(X[tr, , drop = FALSE], m[tr, , drop = FALSE],
                      X[te, , drop = FALSE], choice)
  out <- data.frame(sample_id = m$sample_id[te], disease = m$disease[te],
                    y_true = m$y[te], study_name = m$study_name[te],
                    prob_healthy = fit$prob, OMHI = fit$link,
                    repeat_id = repeat_id, fold_id = fold_id,
                    alpha = choice$alpha, prev_cutoff = choice$prevalence,
                    gamma = choice$gamma, inner_auc = choice$inner_auc,
                    lambda_final_fold = fit$lambda,
                    n_train = fit$n_train, dropped_train = fit$n_train_dropped,
                    n_model_input_taxa = fit$n_features,
                    rank = RANK)
  save_csv(out, dest)
  cat("Completed repeated outer CV ", repeat_id, "/", fold_id,
      " n_holdout=", nrow(out), "\n", sep = "")
}

run_outer_loso <- function(study_idx) {
  d <- read_training(); X <- d$X; m <- d$meta
  studies <- sort(unique(m$study_name))
  if (!study_idx %in% seq_along(studies)) stop("Study index should be 1..", length(studies))
  target <- file.path(OUT_DIR, "loso_folds", sprintf("study_%02d.csv", study_idx))
  if (file.exists(target)) stop("Already exists: ", target)
  st <- studies[study_idx]
  te <- which(m$study_name == st); tr <- which(m$study_name != st)
  choice <- tune_inner(X[tr, , drop = FALSE], m[tr, , drop = FALSE],
                       SEED + 40000L + study_idx * 23L)
  fit <- fit_selected(X[tr, , drop = FALSE], m[tr, , drop = FALSE],
                      X[te, , drop = FALSE], choice)
  out <- data.frame(sample_id = m$sample_id[te], disease = m$disease[te],
                    y_true = m$y[te], study_name = m$study_name[te],
                    prob_healthy = fit$prob, OMHI = fit$link,
                    study_index = study_idx, alpha = choice$alpha,
                    prev_cutoff = choice$prevalence, gamma = choice$gamma,
                    inner_auc = choice$inner_auc,
                    n_train = fit$n_train, dropped_train = fit$n_train_dropped,
                    n_model_input_taxa = fit$n_features,
                    rank = RANK)
  save_csv(out, target)
  cat("Completed LOSO ", study_idx, "/", length(studies),
      ": ", st, " n=", nrow(out), "\n", sep = "")
}

aggregate <- function() {
  d <- read_training(); m <- d$meta
  cv_files <- unlist(lapply(seq_len(N_REPEATS), function(r) {
    file.path(OUT_DIR, "cv_folds", sprintf("r%02d_f%02d.csv", r, seq_len(N_OUTER)))
  }))
  studies <- sort(unique(m$study_name))
  loso_files <- file.path(OUT_DIR, "loso_folds", sprintf("study_%02d.csv", seq_along(studies)))
  needed <- c(cv_files, loso_files)
  if (any(!file.exists(needed)))
    stop("Missing fold output(s): ", paste(head(basename(needed[!file.exists(needed)]), 12L), collapse = ", "))
  cv <- data.table::rbindlist(lapply(cv_files, data.table::fread), fill = TRUE)
  loso <- data.table::rbindlist(lapply(loso_files, data.table::fread), fill = TRUE)
  if (nrow(cv) != nrow(m) * N_REPEATS ||
      anyDuplicated(paste(cv$repeat_id, cv$sample_id)) ||
      anyNA(cv$prob_healthy)) stop("Incomplete/duplicated repeated CV predictions")
  if (nrow(loso) != nrow(m) || anyDuplicated(loso$sample_id) || anyNA(loso$prob_healthy))
    stop("Incomplete/duplicated LOSO predictions")
  if (!setequal(loso$sample_id, m$sample_id)) stop("LOSO sample mismatch")
  save_csv(cv, file.path(OUT_DIR, "summary", "nested_repeated10fold_OOF.csv"))
  save_csv(loso, file.path(OUT_DIR, "summary", "nested_LOSO_OOF.csv"))
  each_repeat <- do.call(rbind, lapply(seq_len(N_REPEATS), function(r) {
    z <- cv[cv$repeat_id == r, ]; a <- get_metrics(z$y_true, z$prob_healthy)
    a$repeat_id <- r; a
  }))
  save_csv(each_repeat, file.path(OUT_DIR, "summary", "nested_repeated10fold_each.csv"))
  mean_rep <- data.frame(auc_mean = mean(each_repeat$auc),
                         auc_sd = sd(each_repeat$auc),
                         BA_mean = mean(each_repeat$balanced_accuracy),
                         BA_sd = sd(each_repeat$balanced_accuracy),
                         n_repeats = N_REPEATS)
  save_csv(mean_rep, file.path(OUT_DIR, "summary", "nested_repeated10fold_mean.csv"))
  each_study <- do.call(rbind, lapply(split(as.data.frame(loso), loso$study_name), function(z) {
    a <- get_metrics(z$y_true, z$prob_healthy)
    a$study_name <- z$study_name[1L]
    a
  }))
  save_csv(each_study, file.path(OUT_DIR, "summary", "nested_LOSO_each_study.csv"))
  pooled <- get_metrics(loso$y_true, loso$prob_healthy)
  eligible <- each_study[is.finite(each_study$auc), , drop = FALSE]
  pooled$mean_study_auc <- if (nrow(eligible)) mean(eligible$auc) else NA_real_
  pooled$weighted_study_auc <- if (nrow(eligible))
    weighted.mean(eligible$auc, w = eligible$n) else NA_real_
  pooled$n_studies_with_auc <- nrow(eligible)
  save_csv(pooled, file.path(OUT_DIR, "summary", "nested_LOSO_summary.csv"))
  print(mean_rep); print(pooled)
}

final_fit <- function() {
  dest <- file.path(OUT_DIR, "final", paste0(RANK, "_raw_nested_model.rds"))
  if (file.exists(dest)) stop("Final model exists; select a new output dir")
  d <- read_training()
  chosen <- tune_inner(d$X, d$meta, SEED + 90000L)
  # Uses exactly the preprocessing policy used in each outer training partition.
  obj <- fit_selected(d$X, d$meta, Xte = NULL, chosen)
  model <- list(fit = obj$model, lambda = obj$lambda, level = RANK,
                alpha = chosen$alpha, prevalence = chosen$prevalence,
                gamma = chosen$gamma, chosen = chosen,
                model_taxa = obj$selected_taxa,
                pre_filter_taxa = obj$pre_taxa,
                clr_pseudocount_percentage = CLR_PC_PERCENT,
                batch_correction = "none",
                benchmark_only = TRUE,
                control_covariate = NA_character_,
                batch_variable = NA_character_)
  save_rds(model, dest)
  save_csv(chosen, file.path(OUT_DIR, "final", paste0(RANK, "_selected_parameters.csv")))
  coefs <- as.matrix(coef(obj$model, s = obj$lambda))
  co <- data.frame(feature = rownames(coefs), coefficient = as.numeric(coefs[, 1L]))
  save_csv(co, file.path(OUT_DIR, "final", paste0(RANK, "_all_coefficients.csv")))
  scores <- d$meta[match(obj$train_sample_ids, d$meta$sample_id), , drop = FALSE]
  scores$OMHI <- obj$train_link
  scores$prob_healthy <- obj$train_prob
  scores$dataset <- "Training (in-sample; NOT CV)"
  save_csv(scores, file.path(OUT_DIR, "final", paste0(RANK, "_train_in_sample_scores.csv")))
  # No adjusted matrix is generated for the raw-control arm.
  cat("Final raw benchmark: ", RANK, ": ", length(obj$selected_taxa),
      " input taxa; nonzero=", sum(co$coefficient[co$feature != "(Intercept)"] != 0),
      "; alpha=", chosen$alpha, "; prev=", chosen$prevalence,
      "; gamma=", chosen$gamma, "\n", sep = "")
}


preflight <- function() {
  d <- read_training()
  X <- d$X; m <- d$meta
  cat("Training samples: ", nrow(X), "\nStudies: ",
      length(unique(m$study_name)), "\nRaw ", RANK, " taxa: ", ncol(X),
      "\nHealthy: ", sum(m$y == 1L), " Nonhealthy: ", sum(m$y == 0L), "\n", sep = "")
  print(table(m$study_name, m$health_bin))
  cat("Testing one fold of uncorrected TSS/CLR preprocessing...\n")
  fold <- make_stratified(m$y, N_OUTER, SEED + 1000L)
  j <- which(fold == 1L); i <- which(fold != 1L)
  a <- preprocess_train(X[i, , drop = FALSE], m[i, , drop = FALSE],
                        X[j, , drop = FALSE])
  b <- get_xy(a, PREVS[1L], skip_if_small = TRUE)
  if (is.null(b)) stop("No valid features even at 1% prevalence: ", RANK)
  cat("Successful fold: train ", nrow(b$x), ", heldout ", nrow(b$x_hold),
      ", features ", ncol(b$x), "\n", sep = "")
  invisible(TRUE)
}

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) {
  stop("Usage: OMWI_RANK=phylum|class|order|family|genus|species Rscript 05_OMWI_rank_raw_nested.R preflight|cv REPEAT FOLD|loso STUDY_INDEX|aggregate|final-fit")
}
mode <- args[1L]
arg_int <- function(i) {
  if (length(args) < i) stop("Missing numeric command argument ", i)
  z <- suppressWarnings(as.integer(args[i]))
  if (!is.finite(z)) stop("Invalid numeric command argument ", i)
  z
}
if (mode == "preflight") {
  preflight()
} else if (mode == "cv") {
  run_outer_cv(arg_int(2L), arg_int(3L))
} else if (mode == "loso") {
  run_outer_loso(arg_int(2L))
} else if (mode == "aggregate") {
  aggregate()
} else if (mode == "final-fit") {
  final_fit()
} else {
  stop("Unknown mode: ", mode)
}
