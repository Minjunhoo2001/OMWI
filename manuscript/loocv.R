#!/usr/bin/env Rscript
# OMWI raw A | TRAIN-ONLY resampling stability ranking + compact-panel refitting
# LOOCV-only rerun: fixed elastic-net alpha=0.1; reselect K and gamma inside each held-out split.
# Inner: 5-fold. Stability ranking is relearnt on EACH inner training partition.
# The compact panel is self-contained: TSS and CLR only use the K panel taxa.
# No MMUPHin, no external outcomes during selection.
suppressPackageStartupMessages({
  library(data.table)
  library(glmnet)
  library(pROC)
})
options(stringsAsFactors = FALSE)

ROOT <- Sys.getenv("OMWI_ROOT", getwd())
RAW_DIR <- Sys.getenv("OMWI_RAW_DIR", file.path(ROOT, "data/abundance"))
META_DIR <- Sys.getenv("OMWI_META_DIR", file.path(ROOT, "data/metadata"))
OUT_DIR <- Sys.getenv("OMWI_LOOCV_OUT", file.path(ROOT, "results/loocv"))
SEED <- 704L
N_REPEATS <- 10L
N_OUTER <- 10L
N_INNER <- 5L
MIN_BATCH_N <- 2L  # preserve raw A exclusion of singleton-study training batches
CLR_PC_PERCENT <- 1e-6
# Exploratory, fixed before running this package; NOT a formal noninferiority claim.
AUC_LOSS_BUDGET <- 0.02
K_GRID <- c(20L,30L,50L,80L,100L,120L,150L)
ALPHAS <- 0.1
GAMMAS <- exp(seq(log(1), log(0.025), length.out = 16L))
# Fixed ranking learner: lasso, lambda / lambda_max = 0.15.
RANK_ALPHA <- 1.0
RANK_GAMMA <- 0.15
RANK_PREV <- 0.1
RANK_B <- 40L
SUBSAMPLE_FRACTION <- 0.8

if (!is.finite(AUC_LOSS_BUDGET) || AUC_LOSS_BUDGET < 0 || AUC_LOSS_BUDGET > 0.10)
  stop("AUC loss budget must be in [0, 0.10]")
if (!length(K_GRID) || anyNA(K_GRID) || any(K_GRID < 2L) || anyDuplicated(K_GRID))
  stop("Invalid or duplicated K_GRID")
K_GRID <- sort(K_GRID)
if (is.na(RANK_B) || RANK_B < 10L) stop("Use >=10 stability resamples")
if (normalizePath(OUT_DIR, mustWork = FALSE) %in% normalizePath(c(ROOT, RAW_DIR, META_DIR,
                                  file.path(ROOT, "05.raw_A_species")), mustWork = FALSE))
  stop("Refusing to overwrite original or input directory")
if (!startsWith(normalizePath(OUT_DIR, winslash = "/", mustWork = FALSE),
                paste0(normalizePath(ROOT, winslash = "/", mustWork = FALSE), "/")))
  stop("Output directory must be a separate child of OMWI_ROOT")
if (basename(normalizePath(OUT_DIR, mustWork=FALSE)) %in% c("data/model_outputs", "05.raw_A_species")) stop("Do not use an existing model output directory for this rerun")
if (dir.exists(file.path(OUT_DIR,"final")) || dir.exists(file.path(OUT_DIR,"cv_folds"))) stop("Use a dedicated LOOCV-only output directory")
if (!dir.exists(RAW_DIR) || !dir.exists(META_DIR)) stop("Missing raw data or metadata directory")
for (sub in c("loocv_folds", "summary"))
  dir.create(file.path(OUT_DIR, sub), showWarnings = FALSE, recursive = TRUE)
CONFIG_STAMP <- paste("loocv_alpha01_prev10_v2", ALPHAS, SEED, N_INNER, normalizePath(RAW_DIR, mustWork = FALSE),
                      normalizePath(META_DIR, mustWork = FALSE), RANK_B,
                      paste(K_GRID, collapse = ";"), AUC_LOSS_BUDGET,
                      RANK_ALPHA, RANK_GAMMA, RANK_PREV,
                      SUBSAMPLE_FRACTION, CLR_PC_PERCENT,
                      paste(tools::md5sum(c(file.path(RAW_DIR,"train_species.tsv"),file.path(META_DIR,"meta_train.csv"))),collapse=";"), sep = "|")

save_once <- function(x, path, rds = FALSE) {
  if (file.exists(path)) stop("Refusing to overwrite: ", path)
  tmp <- paste0(path, ".", Sys.getpid(), ".tmp")
  if (rds) saveRDS(x, tmp, compress = "gzip") else data.table::fwrite(x, tmp)
  if (!file.rename(tmp, path)) stop("Failed to save atomically: ", path)
}
read_taxa <- function(path, ids, panel_taxa = NULL) {
  if (!file.exists(path)) stop("Missing taxa file: ", path)
  dat <- data.table::fread(path, data.table = FALSE, check.names = FALSE)
  if (ncol(dat) < 2L) stop("Empty taxa file")
  names(dat)[1L] <- "Taxon"
  dat$Taxon <- trimws(as.character(dat$Taxon))
  names(dat)[-1L] <- trimws(names(dat)[-1L])
  if (anyNA(dat$Taxon) || any(!nzchar(dat$Taxon)) || anyDuplicated(names(dat)[-1L]))
    stop("Missing/invalid taxa IDs or duplicate sample IDs")
  absent_ids <- setdiff(ids, names(dat)[-1L])
  if (length(absent_ids)) stop("Sample IDs absent from taxa table: ", paste(head(absent_ids), collapse=","))
  a <- as.matrix(dat[, ids, drop = FALSE]); storage.mode(a) <- "double"
  rownames(a) <- dat$Taxon
  if (anyNA(a) || any(!is.finite(a)) || any(a < 0)) stop("Non-finite/negative taxa abundance")
  if (anyDuplicated(rownames(a))) a <- rowsum(a, group = rownames(a), reorder = FALSE)
  if (!is.null(panel_taxa)) {
    missing <- setdiff(panel_taxa, rownames(a))
    if (length(missing)) {
      extra <- matrix(0, nrow = length(missing), ncol = ncol(a),
                      dimnames = list(missing, colnames(a)))
      a <- rbind(a, extra)
      warning(length(missing), " panel species missing from external taxa table; set to zero")
    }
    a <- a[panel_taxa, , drop = FALSE]
  }
  x <- t(a)
  if (!identical(rownames(x), ids)) stop("Misaligned taxa and sample IDs")
  x
}
read_training <- function() {
  mp <- file.path(META_DIR, "meta_train.csv")
  if (!file.exists(mp)) stop("Training metadata missing")
  m <- as.data.frame(data.table::fread(mp))
  for (v in c("sample_id", "study_name", "disease"))
    if (!v %in% names(m)) stop("Missing metadata column: ", v)
  m$sample_id <- trimws(as.character(m$sample_id))
  m$study_name <- trimws(as.character(m$study_name))
  m$disease <- tolower(trimws(as.character(m$disease)))
  if (anyNA(m[, c("sample_id", "study_name", "disease")]) ||
      any(!nzchar(m$sample_id)) || any(!nzchar(m$study_name)) ||
      any(!nzchar(m$disease)) || anyDuplicated(m$sample_id))
    stop("Training metadata has invalid identifiers or labels")
  m$y <- as.integer(m$disease == "healthy")
  if (length(unique(m$y)) != 2L) stop("Both outcome classes required")
  x <- read_taxa(file.path(RAW_DIR, "train_species.tsv"), m$sample_id)
  list(x=x, m=m)
}
make_stratified <- function(y, k, seed) {
  k <- min(as.integer(k), sum(y == 0L), sum(y == 1L))
  if (k < 2L) stop("Fewer than two folds possible")
  set.seed(seed); f <- integer(length(y))
  for (cl in c(0L, 1L)) {
    ix <- which(y == cl)
    f[ix] <- sample(rep(seq_len(k), length.out = length(ix)))
  }
  f
}
auc_safe <- function(y, pred) {
  if (length(unique(y)) < 2L || anyNA(pred)) return(NA_real_)
  as.numeric(pROC::roc(y, pred, levels=c(0,1), direction="<", quiet=TRUE)$auc)
}
metrics <- function(y, pred) {
  classes <- as.integer(pred >= 0.5)
  sensitivity <- if (any(y == 1L)) mean(classes[y == 1L] == 1L) else NA_real_
  specificity <- if (any(y == 0L)) mean(classes[y == 0L] == 0L) else NA_real_
  data.frame(n=length(y), n_healthy=sum(y==1L), n_nonhealthy=sum(y==0L),
             auc=auc_safe(y,pred), balanced_accuracy=if (is.finite(sensitivity) && is.finite(specificity))
                  mean(c(sensitivity,specificity)) else NA_real_,
             sensitivity=sensitivity, specificity=specificity)
}
clr <- function(x) {
  if (any(!is.finite(x)) || any(x < 0)) stop("CLR requires finite nonnegative inputs")
  zz <- log(x * 100 + CLR_PC_PERCENT)
  sweep(zz, 1L, rowMeans(zz), FUN="-")
}
tss <- function(x) {
  sums <- rowSums(x)
  if (any(!is.finite(sums))) stop("Non-finite TSS sum")
  out <- x
  good <- sums > 0
  if (any(good)) out[good, ] <- sweep(x[good,,drop=FALSE], 1L, sums[good], FUN="/")
  if (any(!good)) out[!good, ] <- 0
  out
}
drop_singleton_train <- function(x, m) {
  if (!identical(rownames(x), m$sample_id)) stop("Training X/meta order mismatch")
  ct <- table(m$study_name)
  ok <- m$study_name %in% names(ct[ct >= MIN_BATCH_N])
  if (sum(ok) < 10L || length(unique(m$y[ok])) < 2L)
    stop("Insufficient samples/classes after singleton-study exclusion")
  list(x=x[ok,,drop=FALSE], m=m[ok,,drop=FALSE], dropped=sum(!ok))
}
# Reproducible study x outcome stratified SUBSAMPLING. Sampling happens only within a training fold.
subsample_ids <- function(m, fraction, seed) {
  set.seed(seed)
  strata <- interaction(m$study_name, m$y, drop=TRUE, lex.order=TRUE)
  ix <- split(seq_len(nrow(m)), strata)
  sort(unlist(lapply(ix, function(v) {
    n <- length(v)
    take <- min(n, max(if (n <= 2L) n else 2L, floor(n * fraction)))
    sample(v, size=take, replace=FALSE)
  }), use.names=FALSE))
}
# Fixed lasso base learner for selection frequency; NOT tuned on held-out labels.
rank_once <- function(x, m) {
  d <- drop_singleton_train(x, m)
  vari <- apply(d$x, 2, stats::var)
  prev <- colMeans(d$x > 0)
  available <- names(vari)[is.finite(vari) & vari>0 & prev >= RANK_PREV]
  if (length(available) < 2L) stop("Too few species for ranking")
  # The ranking learner mirrors raw A's TSS/CLR feature calculation.
  z <- clr(tss(d$x[,available,drop=FALSE]))
  prelim <- glmnet::glmnet(z, d$m$y, family="binomial", alpha=RANK_ALPHA,
                           standardize=TRUE, nlambda=5L)
  lambda <- as.numeric(prelim$lambda[1L] * RANK_GAMMA)
  fit <- glmnet::glmnet(z, d$m$y, family="binomial", alpha=RANK_ALPHA,
                        standardize=TRUE, lambda=lambda)
  cc <- as.matrix(stats::coef(fit))
  stats::setNames(as.numeric(cc[-1L, 1L]), rownames(cc)[-1L])
}
# Rank ALL varying species. Never look at validation labels for the ranking.
stability_rank <- function(x, m, seed) {
  d <- drop_singleton_train(x, m)
  vari <- apply(d$x, 2L, stats::var)
  prev <- colMeans(d$x > 0)
  universe <- names(vari)[is.finite(vari) & vari > 0 & prev >= RANK_PREV]
  if (length(universe) < min(K_GRID)) stop("Not enough species for smallest candidate K")
  n_sel <- stats::setNames(integer(length(universe)), universe)
  n_pos <- n_sel; n_neg <- n_sel
  sum_abs <- stats::setNames(numeric(length(universe)), universe)
  failures <- 0L
  for (b in seq_len(RANK_B)) {
    idx <- subsample_ids(d$m, SUBSAMPLE_FRACTION, seed + b * 101L)
    cc <- tryCatch(rank_once(d$x[idx,,drop=FALSE], d$m[idx,,drop=FALSE]),
                   error=function(e) { warning("Ranking resample ", b, ": ", conditionMessage(e)); NULL })
    if (is.null(cc)) { failures <- failures+1L; next }
    cc <- cc[names(cc) %in% universe]
    active <- names(cc)[is.finite(cc) & abs(cc) > 1e-12]
    n_sel[active] <- n_sel[active] + 1L
    pos <- active[cc[active] > 0]
    neg <- active[cc[active] < 0]
    n_pos[pos] <- n_pos[pos] + 1L
    n_neg[neg] <- n_neg[neg] + 1L
    sum_abs[active] <- sum_abs[active] + abs(cc[active])
  }
  succeeded <- RANK_B - failures
  if (succeeded < ceiling(0.8 * RANK_B))
    stop("Stability sampling failed too often: ", succeeded, "/", RANK_B)
  sign_consistency <- ifelse(n_sel > 0L, pmax(n_pos, n_neg) / pmax(1L,n_sel), 0)
  result <- data.frame(feature=universe, selected_n=as.integer(n_sel[universe]),
                       selection_frequency=as.numeric(n_sel[universe]/succeeded),
                       sign_consistency=as.numeric(sign_consistency[universe]),
                       mean_abs_coefficient=as.numeric(sum_abs[universe]/succeeded),
                       raw_prevalence=as.numeric(prev[universe]),
                       n_resamples_successful=succeeded, n_resamples_failed=failures)
  result <- result[order(-result$selection_frequency, -result$sign_consistency,
                         -result$mean_abs_coefficient, -result$raw_prevalence,
                         result$feature), , drop=FALSE]
  rownames(result) <- NULL
  result$stability_rank <- seq_len(nrow(result))
  result
}
# Refitted model ONLY depends on its K panel taxa (including TSS denominator and CLR center).
prepare_panel <- function(xtr, mtr, xhold, panel) {
  d <- drop_singleton_train(xtr, mtr)
  if (length(panel) < 2L || anyDuplicated(panel) || !all(panel %in% colnames(xtr)))
    stop("Invalid panel feature list")
  train_zero <- sum(rowSums(d$x[,panel,drop=FALSE]) == 0)
  hold_zero <- if (!is.null(xhold)) sum(rowSums(xhold[,panel,drop=FALSE]) == 0) else NA_integer_
  list(x=clr(tss(d$x[,panel,drop=FALSE])), y=d$m$y,
       xhold=if (!is.null(xhold)) clr(tss(xhold[,panel,drop=FALSE])) else NULL,
       train_zero=train_zero, hold_zero=hold_zero, dropped=d$dropped,
       train_ids=d$m$sample_id)
}
fit_path <- function(x, y, alpha) {
  prelim <- glmnet::glmnet(x, y, family="binomial", alpha=alpha,
                           standardize=TRUE, nlambda=5L)
  lams <- as.numeric(prelim$lambda[1L] * GAMMAS)
  fit <- glmnet::glmnet(x, y, family="binomial", alpha=alpha,
                        standardize=TRUE, lambda=lams)
  available <- lams >= min(fit$lambda) * (1 - 1e-8)
  list(fit=fit, lambdas=lams, available=available)
}
fit_final_panel <- function(xtr, mtr, xhold, panel, alpha, gamma) {
  xx <- prepare_panel(xtr, mtr, xhold, panel)
  path <- fit_path(xx$x, xx$y, alpha)
  gi <- which.min(abs(GAMMAS - gamma))
  if (!path$available[gi]) stop("Chosen gamma absent from full-refit lambda path")
  lambda <- path$lambdas[gi]
  pred <- if (!is.null(xx$xhold)) as.numeric(predict(path$fit, newx=xx$xhold,
                                                    s=lambda, type="response")) else NULL
  logit <- if (!is.null(xx$xhold)) as.numeric(predict(path$fit, newx=xx$xhold,
                                                     s=lambda, type="link")) else NULL
  cc <- as.matrix(stats::coef(path$fit, s=lambda))
  list(fit=path$fit, lambda=lambda, prob=pred, logit=logit,
       n_nonzero=sum(abs(cc[-1L, 1L]) > 1e-12),
       n_zero_train=xx$train_zero, n_zero_hold=xx$hold_zero,
       n_train=length(xx$y), dropped=xx$dropped, train_ids=xx$train_ids,
       train_prob=as.numeric(predict(path$fit,newx=xx$x,s=lambda,type="response")),
       train_link=as.numeric(predict(path$fit,newx=xx$x,s=lambda,type="link")))
}
# Critical nesting: each INNER fold ranks species only on its inner TRAINING portion.
tune_inner <- function(x, m, seed, return_grid=FALSE) {
  foldid <- make_stratified(m$y, N_INNER, seed)
  k <- max(foldid)
  scores <- list()
  for (j in seq_len(k)) {
    tr <- which(foldid != j); te <- which(foldid == j)
    xt <- x[tr,,drop=FALSE]; mt <- m[tr,,drop=FALSE]
    xh <- x[te,,drop=FALSE]
    ranking <- stability_rank(xt, mt, seed + j * 100003L)
    for (K in K_GRID[K_GRID <= nrow(ranking)]) {
      panel <- ranking$feature[seq_len(K)]
      xx <- prepare_panel(xt, mt, xh, panel)
      for (alpha in ALPHAS) {
        path <- fit_path(xx$x, xx$y, alpha)
        use <- which(path$available)
        if (!length(use)) next
        pred <- as.matrix(predict(path$fit, newx=xx$xhold,
                                  s=path$lambdas[use], type="response"))
        if (ncol(pred) != length(use)) stop("Unexpected prediction-column count")
        for (h in seq_along(use)) {
          scores[[length(scores)+1L]] <- data.frame(
            inner_fold=j, K=K, alpha=alpha, gamma=GAMMAS[use[h]],
            auc=auc_safe(m$y[te],pred[,h]))
        }
      }
    }
  }
  if (!length(scores)) stop("No valid inner-CV candidate")
  tab <- data.table::rbindlist(scores)
  grid <- tab[, .(n_inner_folds=.N, mean_auc=if (all(is.finite(auc))) mean(auc) else NA_real_,
                   auc_sd=if (all(is.finite(auc))) stats::sd(auc) else NA_real_),
              by=.(K,alpha,gamma)]
  grid$auc_se <- grid$auc_sd / sqrt(k)
  grid <- as.data.frame(grid)
  grid <- grid[grid$n_inner_folds == k & is.finite(grid$mean_auc), , drop=FALSE]
  if (!nrow(grid)) stop("No candidate evaluated on all inner folds")
  best <- max(grid$mean_auc)
  cutoff <- best - AUC_LOSS_BUDGET
  grid$within_budget <- grid$mean_auc >= cutoff - 1e-12
  eligible <- grid[grid$within_budget,,drop=FALSE]
  eligible <- eligible[order(eligible$K, -eligible$mean_auc,
                             -eligible$gamma, -eligible$alpha),,drop=FALSE]
  chosen <- eligible[1L,,drop=FALSE]
  chosen$best_inner_auc <- best
  chosen$acceptable_inner_auc <- cutoff
  chosen$auc_loss_budget <- AUC_LOSS_BUDGET
  chosen$n_eligible <- nrow(eligible)
  chosen$selection_rule <- "inner_auc_loss_budget_then_min_panel_K"
  rownames(chosen) <- NULL
  if (return_grid) return(list(chosen=chosen, grid=grid))
  chosen
}
make_fold_result <- function(x, m, tr, te, seed, extra) {
  choice <- tune_inner(x[tr,,drop=FALSE],m[tr,,drop=FALSE],seed)
  ranking <- stability_rank(x[tr,,drop=FALSE],m[tr,,drop=FALSE],seed+7000003L)
  if (nrow(ranking) < choice$K) stop("Chosen K larger than eligible taxa")
  panel <- ranking$feature[seq_len(choice$K)]
  fit <- fit_final_panel(x[tr,,drop=FALSE],m[tr,,drop=FALSE],x[te,,drop=FALSE],
                         panel,choice$alpha,choice$gamma)
  pred <- data.frame(sample_id=m$sample_id[te], disease=m$disease[te],
                     y_true=m$y[te], study_name=m$study_name[te],
                     prob_healthy=fit$prob, OMHI=fit$logit,
                     K_panel=length(panel), n_nonzero=fit$n_nonzero,
                     alpha=choice$alpha, gamma=choice$gamma,
                     inner_auc=choice$mean_auc, best_inner_auc=choice$best_inner_auc,
                     n_zero_panel_hold_total=fit$n_zero_hold, n_train=fit$n_train,
                     dropped_train=fit$dropped)
  for (nm in names(extra)) pred[[nm]] <- extra[[nm]]
  list(config=CONFIG_STAMP, predictions=pred, panel=panel,
       panel_rank=ranking[seq_len(choice$K),,drop=FALSE], chosen=choice)
}
run_loocv <- function(sample_index) {
  d <- read_training()
  n <- nrow(d$m)
  if (length(sample_index) != 1L || is.na(sample_index) ||
      !(sample_index %in% seq_len(n))) stop("Invalid LOOCV sample index")
  dest <- file.path(OUT_DIR, "loocv_folds", sprintf("sample_%05d.rds", sample_index))
  if (file.exists(dest)) {
    prior <- readRDS(dest)
    if (!identical(prior$config, CONFIG_STAMP) || nrow(prior$predictions) != 1L ||
        !identical(as.character(prior$predictions$sample_id), d$m$sample_id[sample_index]) ||
        !isTRUE(all.equal(as.numeric(prior$predictions$alpha), ALPHAS))) stop("Existing fold incompatible: ",dest)
    cat("Verified existing fold; skipping: ",dest,"\n",sep=""); return(invisible(NULL))
  }
  te <- sample_index
  tr <- setdiff(seq_len(n), te)
  res <- make_fold_result(d$x, d$m, tr, te,
                          SEED + 100000L + sample_index * 23L,
                          list(loo_id=sample_index, cv_scheme="LOOCV"))
  if (nrow(res$predictions) != 1L ||
      !identical(as.character(res$predictions$sample_id),
                 as.character(d$m$sample_id[sample_index])))
    stop("LOOCV prediction must correspond to the single held-out sample")
  save_once(res, dest, rds=TRUE)
  cat("Completed LOOCV ", sample_index, "/", n,
      " sample=", res$predictions$sample_id,
      " panel K=", length(res$panel), "\n", sep="")
}

# Aggregate LOOCV independently from existing CV/LOSO/final-model results,
# so adding LOOCV never requires rerunning or overwriting earlier analyses.
aggregate_loocv <- function() {
  d <- read_training()
  m <- d$m
  n <- nrow(m)
  paths <- file.path(OUT_DIR, "loocv_folds",
                     sprintf("sample_%05d.rds", seq_len(n)))
  missing <- paths[!file.exists(paths)]
  if (length(missing)) stop("Incomplete LOOCV: ", length(missing),
                            " missing, e.g. ", paste(head(basename(missing)),collapse=","))
  results <- lapply(paths, readRDS)
  if (!all(vapply(results, function(z) identical(z$config, CONFIG_STAMP), logical(1L))))
    stop("LOOCV folds were generated under inconsistent configurations")
  if (any(vapply(results, function(z) nrow(z$predictions) != 1L, logical(1L))))
    stop("Each LOOCV fold must contain exactly one held-out prediction")
  preds <- as.data.frame(data.table::rbindlist(lapply(results, `[[`, "predictions"), fill=TRUE))
  if (nrow(preds) != n || anyNA(preds$alpha) || any(preds$alpha != ALPHAS) || anyNA(preds$loo_id) ||
      !identical(as.integer(preds$loo_id), seq_len(n)) ||
      !identical(as.character(preds$sample_id), as.character(m$sample_id)) ||
      !identical(as.character(preds$study_name), as.character(m$study_name)) ||
      !identical(as.integer(preds$y_true), as.integer(m$y)) ||
      anyDuplicated(preds$sample_id) || anyNA(preds$prob_healthy) ||
      any(!is.finite(preds$prob_healthy)) ||
      any(preds$prob_healthy < 0 | preds$prob_healthy > 1))
    stop("LOOCV indices, held-out IDs, outcomes or probabilities do not align")

  summary <- metrics(preds$y_true, preds$prob_healthy)
  summary$cv_scheme <- "nested_LOOCV"
  summary$n_outer_models <- n
  summary$median_K_panel <- stats::median(preds$K_panel)
  summary$mean_K_panel <- mean(preds$K_panel)
  summary$median_n_nonzero <- stats::median(preds$n_nonzero)
  summary$mean_n_nonzero <- mean(preds$n_nonzero)
  summary$classification_threshold <- 0.5
  # AUC and BA are computed ONCE from the pooled n held-out predictions,
  # not averaged across single-sample folds (where AUC is undefined).
  save_once(preds, file.path(OUT_DIR, "summary", "nested_LOOCV_OOF.csv"))
  save_once(summary, file.path(OUT_DIR, "summary", "nested_LOOCV_summary.csv"))
  complexity <- preds[, c("loo_id", "sample_id", "study_name", "y_true",
                          "K_panel", "n_nonzero", "alpha", "gamma",
                          "inner_auc", "best_inner_auc", "n_train",
                          "dropped_train", "n_zero_panel_hold_total"), drop=FALSE]
  save_once(complexity, file.path(OUT_DIR, "summary", "nested_LOOCV_model_complexity.csv"))
  panel_taxa <- unlist(lapply(results, `[[`, "panel"), use.names=FALSE)
  freq <- as.data.frame(table(panel_taxa), stringsAsFactors=FALSE)
  names(freq) <- c("feature", "n_LOOCV_panel_selections")
  freq$frequency_across_LOOCV_models <- freq$n_LOOCV_panel_selections / n
  freq <- freq[order(-freq$frequency_across_LOOCV_models, freq$feature), , drop=FALSE]
  save_once(freq, file.path(OUT_DIR, "summary", "nested_LOOCV_panel_frequency.csv"))
  print(summary)
  cat("Saved complete nested LOOCV (", n, " held-out samples).\n",sep="")
}

preflight <- function() {
  d <- read_training()
  if (nrow(d$m)!=1077L || length(unique(d$m$study_name))!=17L) stop("Expected original 1,077-sample/17-study training set")
  cat("Refit alpha=",ALPHAS,"; ranking prevalence=",RANK_PREV,"\n",sep="")
  cat("INPUT RAW: ",RAW_DIR,"\nINPUT META: ",META_DIR,"\nOUTPUT: ",OUT_DIR,
      "\nTRAIN n=",nrow(d$x)," studies=",length(unique(d$m$study_name)),
      " taxa=",ncol(d$x),"\n",sep="")
  cat("K grid: ",paste(K_GRID,collapse=","),"; stability B=",RANK_B,
      "; inner AUC loss budget=",AUC_LOSS_BUDGET,"\n",sep="")
  f <- make_stratified(d$m$y,N_OUTER,SEED+1000L)
  idx <- which(f!=1L)
  s <- subsample_ids(d$m[idx,,drop=FALSE],SUBSAMPLE_FRACTION,SEED+4L)
  co <- rank_once(d$x[idx[s],,drop=FALSE],d$m[idx[s],,drop=FALSE])
  if (!length(co)) stop("Ranking smoke failed")
  cat("Ranking fit works, coefficients=",length(co),". Full smoke next.\n",sep="")
}
args <- commandArgs(trailingOnly=TRUE)
if (!length(args)) stop("Usage: preflight | loocv INDEX | aggregate-loocv")
mode <- args[1L]
if (mode == "preflight") preflight() else if (mode == "loocv") {
  if (length(args)!=2L) stop("Provide a sample index")
  run_loocv(suppressWarnings(as.integer(args[2L])))
} else if (mode == "aggregate-loocv") aggregate_loocv() else stop("Unknown mode: ",mode)
