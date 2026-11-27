#!/usr/bin/env Rscript
# ============================================================================
# OMWI | Age/sex-adjusted logistic ROC | Nature-style, two panels
# a Training cohort; b Independent external validation cohort
#   Unadjusted:      nonhealthy ~ OMWI
#   Age/sex-adjusted: nonhealthy ~ OMWI + age_num + sex
# Both logistic models are fitted ONLY in the training cohort, then applied
# unchanged to the external validation cohort (no external refitting).
# Comparison is performed on the SAME age/sex-complete samples in each cohort.
# No matching, no retraining of the original OMWI microbial model.
# ============================================================================

suppressPackageStartupMessages({
  library(pROC)
  library(ggplot2)
})

# ------------------------------ 1. PATHS ------------------------------------
ROOT_DIR     <- Sys.getenv("OMWI_ROOT", getwd())
MODEL_FOLDER <- "data/model_outputs"

# Explicit frozen-model inputs; never search alternative model versions.
TRAIN_CSV <- file.path(ROOT_DIR, MODEL_FOLDER, "final", "species_train_in_sample_scores.csv")
EXTERNAL_CSV <- file.path(ROOT_DIR, MODEL_FOLDER, "final", "external_independent_predictions.csv")
EXPECTED_TRAIN_N <- 1077L
EXPECTED_EXTERNAL_N <- 716L

OUT_DIR <- file.path(ROOT_DIR, "Supplementaryfigure", "bias",
                     "age_sex_separate_ROC_UTF8")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

FONT_FAMILY <- "sans"  # Change to "Arial" if available.
COL_RAW     <- "#747474"  # Muted gray (unadjusted)
COL_ADJ     <- "#008F82"  # Teal (age/sex-adjusted)

# -------------------------- 2. READ FIXED INPUTS ----------------------------
read_input <- function(path) {
  if (!file.exists(path)) stop("Missing input: ", path)
  withCallingHandlers(
    utils::read.csv(path, fileEncoding = "UTF-8-BOM",
                    stringsAsFactors = FALSE, check.names = FALSE),
    warning = function(w) stop("CSV read warning: ", path, ": ", conditionMessage(w), call. = FALSE)
  )
}

# -------------------------- 2b. SAMPLE COUNT AUDIT --------------------------
# Diagnostics ONLY: does not modify model, ROC, or sample-selection rules.
# In particular, numeric gender coding is NOT guessed or silently recoded.
QC_DIR <- file.path(OUT_DIR, "sample_count_audit")
dir.create(QC_DIR, recursive = TRUE, showWarnings = FALSE)

write_qc <- function(df, name) {
  utils::write.csv(df, file.path(QC_DIR, name), row.names = FALSE,
                   fileEncoding = "UTF-8", na = "")
}

run_sample_audit <- function(df, label, source_path) {
  require_cols <- c("sample_id", "age_num", "gender", "disease")
  missing_cols <- setdiff(require_cols, names(df))
  if (length(missing_cols)) stop("Audit ", label, " missing required columns: ",
                                 paste(missing_cols, collapse = ", "))
  present_scores <- intersect(c("OMHI", "OMWI", "OMWI_score"), names(df))
  if (!length(present_scores)) stop("Audit ", label, ": no OMWI score column")
  score_col <- present_scores[1L] # EXACTLY the original script's priority
  if (length(present_scores) > 1L) {
    warning(label, " has MULTIPLE score columns: ",
            paste(present_scores, collapse = ", "),
            "; current ROC uses '" , score_col, "' (original script priority).")
  }
  
  norm <- function(x) tolower(trimws(as.character(x)))
  phen <- norm(df$disease)
  sex  <- norm(df$gender)
  valid_y <- phen %in% c("disease", "diseased", "nonhealthy", "non-healthy",
                         "case", "1", "healthy", "health", "control", "ctrl", "0")
  valid_sex <- sex %in% c("male", "m", "1", "female", "f", "2")
  sc <- suppressWarnings(as.numeric(df[[score_col]]))
  ag <- suppressWarnings(as.numeric(df$age_num))
  valid_score <- is.finite(sc)
  valid_age <- is.finite(ag)
  adult <- valid_age & ag >= 18
  keep <- valid_y & valid_score & adult & valid_sex
  n <- nrow(df)
  n_at_step <- c(
    n,
    sum(valid_y),
    sum(valid_y & valid_score),
    sum(valid_y & valid_score & valid_age),
    sum(valid_y & valid_score & adult),
    sum(keep)
  )
  audit_steps <- data.frame(
    Cohort = label,
    Step = c("00 CSV original rows", "01 recognizable health labels",
             "02 + finite selected OMWI score", "03 + numeric nonmissing age",
             "04 + age >= 18", "05 + recognized sex (FINAL N)"),
    N_remaining = n_at_step,
    N_lost_at_step = c(0L, -diff(n_at_step)),
    Source_file = normalizePath(source_path, winslash = "/", mustWork = TRUE),
    Score_column_used = score_col,
    stringsAsFactors = FALSE
  )
  write_qc(audit_steps, paste0("01_", label, "_stepwise_counts.csv"))
  
  lost <- data.frame(
    Cohort = label, Metric = c("Original CSV rows", "Accepted all conditions",
                               "Invalid / unrecognized disease (overlaps allowed)",
                               "Missing / nonnumeric / nonfinite score (overlaps allowed)",
                               "Missing / nonnumeric / nonfinite age (overlaps allowed)",
                               "Age under 18 among numeric age (overlaps allowed)",
                               "Invalid / unrecognized sex (overlaps allowed)"),
    N = c(n, sum(keep), sum(!valid_y), sum(!valid_score), sum(!valid_age),
          sum(valid_age & ag < 18), sum(!valid_sex)),
    stringsAsFactors = FALSE
  )
  write_qc(lost, paste0("02_", label, "_invalid_field_counts.csv"))
  
  reasons <- vapply(seq_len(n), function(i) {
    r <- character()
    if (!valid_y[i])     r <- c(r, "unrecognized_disease")
    if (!valid_score[i]) r <- c(r, "missing_or_nonnumeric_score")
    if (!valid_age[i])   r <- c(r, "missing_or_nonnumeric_age")
    if (valid_age[i] && ag[i] < 18) r <- c(r, "age_under_18")
    if (!valid_sex[i])   r <- c(r, "unrecognized_gender")
    if (!length(r)) "INCLUDED" else paste(r, collapse = "; ")
  }, character(1))
  study_col <- intersect(c("study", "study_id", "study_name", "Study", "STUDY",
                           "dataset", "dataset_id", "project", "study_accession"), names(df))
  study_values <- if (length(study_col)) as.character(df[[study_col[1L]]]) else rep(NA_character_, n)
  samples <- data.frame(
    sample_id = as.character(df$sample_id), study = study_values,
    disease_raw = as.character(df$disease),
    score_raw = as.character(df[[score_col]]),
    age_raw = as.character(df$age_num),
    gender_raw = as.character(df$gender),
    recognized_disease = valid_y, finite_score = valid_score,
    numeric_age = valid_age, age_at_least_18 = adult,
    recognized_gender = valid_sex, included_in_ROC = keep,
    exclusion_reasons = reasons, stringsAsFactors = FALSE
  )
  write_qc(samples, paste0("03_", label, "_every_sample_exclusion_reasons.csv"))
  
  if (length(study_col)) {
    studies <- as.character(df[[study_col[1L]]])
    studies[is.na(studies) | !nzchar(trimws(studies))] <- "<MISSING>"
    study_levels <- unique(studies)
    perstudy <- do.call(rbind, lapply(study_levels, function(z) {
      ii <- studies == z
      data.frame(Study = z, Original_N = sum(ii), Included_N = sum(ii & keep),
                 Excluded_N = sum(ii & !keep),
                 Missing_age_N = sum(ii & !valid_age),
                 Invalid_gender_N = sum(ii & !valid_sex),
                 Invalid_disease_N = sum(ii & !valid_y),
                 Invalid_score_N = sum(ii & !valid_score),
                 stringsAsFactors = FALSE)
    }))
    perstudy <- perstudy[order(-perstudy$Excluded_N, -perstudy$Original_N), ]
    write_qc(perstudy, paste0("04_", label, "_by_study.csv"))
  }
  
  raw_levels <- do.call(rbind, lapply(c("gender", "disease"), function(column) {
    x <- as.character(df[[column]])
    x[is.na(x)] <- "<NA>"
    z <- sort(table(x), decreasing = TRUE)
    data.frame(Cohort = label, Field = column, Raw_value = names(z),
               N = as.integer(z), stringsAsFactors = FALSE)
  }))
  write_qc(raw_levels, paste0("05_", label, "_raw_label_frequencies.csv"))
  
  message("\n========== ", label, " INPUT AUDIT ==========")
  message("SOURCE = ", normalizePath(source_path, winslash = "/", mustWork = TRUE))
  message("SCORE COLUMN = ", score_col)
  print(audit_steps[, c("Step", "N_remaining", "N_lost_at_step")], row.names = FALSE)
  message("Independent invalid-field counts (overlap; DO NOT sum these counts):")
  print(lost[, c("Metric", "N")], row.names = FALSE)
  message("Sex labels as imported:")
  print(subset(raw_levels, Field == "gender"), row.names = FALSE)
  message("Disease labels as imported:")
  print(subset(raw_levels, Field == "disease"), row.names = FALSE)
  if (label == "Training") {
    if (n <= 200L) {
      warning("Training source CSV has only ", n,
              " rows BEFORE filtering. Check whether the correct FULL training prediction file was selected.")
    } else if (sum(keep) <= 200L) {
      warning("Training CSV has ", n, " original rows but only ", sum(keep),
              " pass the complete-case rules. Inspect stepwise counts and raw labels.")
    }
  }
  invisible(list(n_original = n, n_final = sum(keep), score_col = score_col))
}


# ------------------------ 3. DEFINE VARIABLES -------------------------------
prepare_data <- function(df, dataset_label) {
  needed <- c("sample_id", "age_num", "gender", "disease")
  missing <- setdiff(needed, names(df))
  if (length(missing)) stop(dataset_label, " CSV is missing: ",
                            paste(missing, collapse = ", "))
  score_col <- intersect(c("OMHI", "OMWI", "OMWI_score"), names(df))
  if (length(score_col) == 0L) {
    stop(dataset_label, " CSV lacks OMWI score (expected OMHI or OMWI).")
  }
  score_col <- score_col[1L]
  
  ids <- trimws(as.character(df$sample_id))
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop(dataset_label, ": missing or duplicate sample_id values.")
  }
  
  disease <- tolower(trimws(as.character(df$disease)))
  y <- ifelse(disease %in% c("disease", "diseased", "nonhealthy",
                             "non-healthy", "case", "1"), 1L,
              ifelse(disease %in% c("healthy", "health", "control",
                                    "ctrl", "0"), 0L, NA_integer_))
  sex_chr <- tolower(trimws(as.character(df$gender)))
  sex_chr <- ifelse(sex_chr %in% c("male", "m", "1"), "Male",
                    ifelse(sex_chr %in% c("female", "f", "2"), "Female", NA_character_))
  
  d <- data.frame(sample_id = ids, nonhealthy = y,
                  OMWI = suppressWarnings(as.numeric(df[[score_col]])),
                  age_num = suppressWarnings(as.numeric(df$age_num)),
                  sex = factor(sex_chr, levels = c("Male", "Female")))
  ok <- !is.na(d$nonhealthy) & is.finite(d$OMWI) &
    is.finite(d$age_num) & d$age_num >= 18 & !is.na(d$sex)
  n_original <- nrow(d)
  d <- d[ok, , drop = FALSE]
  rownames(d) <- NULL
  
  if (length(unique(d$nonhealthy)) != 2L) {
    stop(dataset_label, " must contain both healthy and non-healthy samples.")
  }
  message(dataset_label, ": complete cases = ", nrow(d), " / ", n_original,
          "; healthy = ", sum(d$nonhealthy == 0L),
          "; non-healthy = ", sum(d$nonhealthy == 1L))
  d
}

train_path <- TRAIN_CSV
external_path <- EXTERNAL_CSV
if (normalizePath(train_path) == normalizePath(external_path)) {
  stop("Training and external CSV paths must differ.")
}
# Read each selected file ONCE; audit the exact data used for the figure.
train_raw <- read_input(train_path)
external_raw <- read_input(external_path)
check_input <- function(df, expected, label) {
  if (nrow(df) != expected) stop(label, ": expected ", expected, " rows, got ", nrow(df))
  if (!"sample_id" %in% names(df)) stop(label, ": missing sample_id")
  ids <- trimws(as.character(df$sample_id))
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) stop(label, ": invalid or duplicate IDs")
  message(label, ": verified ", nrow(df), " unique samples")
}
check_input(train_raw, EXPECTED_TRAIN_N, "Training")
check_input(external_raw, EXPECTED_EXTERNAL_N, "External")
if (length(intersect(trimws(train_raw$sample_id), trimws(external_raw$sample_id)))) {
  stop("Training/external overlap before complete-case filtering")
}

train_qc <- run_sample_audit(train_raw, "Training", train_path)
external_qc <- run_sample_audit(external_raw, "External", external_path)
train <- prepare_data(train_raw, "Training")
external <- prepare_data(external_raw, "External validation")
stopifnot(nrow(train) == train_qc$n_final,
          nrow(external) == external_qc$n_final)
message("
Audit CSVs saved to: ", QC_DIR, "
")

shared_ids <- intersect(train$sample_id, external$sample_id)
if (length(shared_ids)) {
  stop("Training/external sample_id overlap (n = ", length(shared_ids),
       "). Check dataset independence.")
}
if (length(unique(train$sex)) < 2L) {
  stop("Both sexes are required in training to estimate the sex coefficient.")
}

# ---------------------- 4. TRAIN TWO LOGISTIC MODELS ------------------------
# The original OMWI microbial index is fixed; the two GLMs below only evaluate
# the contribution of age/sex to disease discrimination.
fit_raw <- glm(nonhealthy ~ OMWI, data = train, family = binomial())
fit_adj <- glm(nonhealthy ~ OMWI + age_num + sex,
               data = train, family = binomial())
if (!fit_raw$converged || !fit_adj$converged || any(!is.finite(coef(fit_raw))) || any(!is.finite(coef(fit_adj)))) {
  stop("A logistic coefficient could not be estimated. Check data separation/collinearity.")
}

# IMPORTANT: use the SAME fitted training coefficients in both cohorts.
train$pred_raw    <- predict(fit_raw, newdata = train, type = "response")
train$pred_adj    <- predict(fit_adj, newdata = train, type = "response")
external$pred_raw <- predict(fit_raw, newdata = external, type = "response")
external$pred_adj <- predict(fit_adj, newdata = external, type = "response")

# --------------------------- 5. AUC + 95% CI --------------------------------
make_roc <- function(d, model_name, cohort) {
  predictor <- if (model_name == "Unadjusted") d$pred_raw else d$pred_adj
  if (anyNA(predictor) || any(!is.finite(predictor))) {
    stop(cohort, ": non-finite predicted probability.")
  }
  roc_obj <- pROC::roc(response = d$nonhealthy, predictor = predictor,
                       levels = c(0, 1), direction = "<", quiet = TRUE)
  ci <- as.numeric(pROC::ci.auc(roc_obj, method = "delong"))
  summary <- data.frame(Cohort = cohort, Model = model_name,
                        N = nrow(d), N_healthy = sum(d$nonhealthy == 0L),
                        N_nonhealthy = sum(d$nonhealthy == 1L),
                        AUC = as.numeric(pROC::auc(roc_obj)),
                        CI_low = ci[1L], CI_high = ci[3L],
                        stringsAsFactors = FALSE)
  points <- as.data.frame(pROC::coords(roc_obj, x = "all",
                                       input = "threshold",
                                       ret = c("specificity", "sensitivity"),
                                       transpose = FALSE))
  curve <- data.frame(Cohort = cohort, Model = model_name,
                      FPR = 1 - points$specificity,
                      TPR = points$sensitivity)
  curve <- curve[is.finite(curve$FPR) & is.finite(curve$TPR), ]
  curve <- curve[order(curve$FPR, curve$TPR), , drop = FALSE]
  list(summary = summary, curve = curve)
}

results <- list(make_roc(train, "Unadjusted", "Training"),
                make_roc(train, "Age/sex-adjusted", "Training"),
                make_roc(external, "Unadjusted", "External"),
                make_roc(external, "Age/sex-adjusted", "External"))
auc_table <- do.call(rbind, lapply(results, `[[`, "summary"))
roc_df <- do.call(rbind, lapply(results, `[[`, "curve"))
write.csv(auc_table, file.path(OUT_DIR, "ROC_auc_95CI_summary.csv"), row.names = FALSE)
write.csv(roc_df, file.path(OUT_DIR, "ROC_coordinates.csv"), row.names = FALSE)

# --------------------------- 6. NATURE-STYLE FIGURE --------------------------
lab_training <- paste0("a   Training cohort (n = ", nrow(train), ")")
lab_external <- paste0("b   External validation (n = ", nrow(external), ")")
roc_df$Panel <- factor(ifelse(roc_df$Cohort == "Training",
                              lab_training, lab_external),
                       levels = c(lab_training, lab_external))
roc_df$Model <- factor(roc_df$Model,
                       levels = c("Unadjusted", "Age/sex-adjusted"))

labels <- auc_table
labels$Panel <- factor(ifelse(labels$Cohort == "Training",
                              lab_training, lab_external),
                       levels = c(lab_training, lab_external))
labels$Model <- factor(labels$Model,
                       levels = c("Unadjusted", "Age/sex-adjusted"))
labels$y <- ifelse(labels$Model == "Unadjusted", 0.235, 0.145)
labels$x <- 0.965
labels$text <- sprintf("%s  AUC = %.4f", as.character(labels$Model), labels$AUC)

fig <- ggplot(roc_df, aes(x = FPR, y = TPR, colour = Model, group = Model)) +
  geom_abline(slope = 1, intercept = 0, colour = "#BBBBBB",
              linetype = "dashed", linewidth = 0.55) +
  geom_step(direction = "hv", linewidth = 0.95) +
  geom_text(data = labels, aes(x = x, y = y, label = text, colour = Model),
            inherit.aes = FALSE, hjust = 1, vjust = 0.5,
            size = 3.25, family = FONT_FAMILY, show.legend = FALSE) +
  facet_wrap(~ Panel, nrow = 1) +
  scale_colour_manual(values = c("Unadjusted" = COL_RAW,
                                 "Age/sex-adjusted" = COL_ADJ)) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                     expand = expansion(mult = 0)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                     expand = expansion(mult = 0)) +
  coord_fixed(ratio = 1) +
  labs(x = "1 - Specificity", y = "Sensitivity") +
  theme_classic(base_size = 11, base_family = FONT_FAMILY) +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.75),
    panel.grid = element_blank(),
    axis.line = element_blank(),
    axis.text = element_text(size = 9, colour = "#111111"),
    axis.title = element_text(size = 10.5, colour = "#111111"),
    axis.ticks = element_line(colour = "#111111", linewidth = 0.45),
    axis.ticks.length = grid::unit(2.5, "pt"),
    strip.background = element_blank(),
    strip.text = element_text(size = 11, face = "bold", hjust = 0,
                              colour = "#111111", margin = margin(b = 8)),
    legend.position = "none",  # direct AUC labels already identify each curve
    panel.spacing = grid::unit(0.38, "in"),
    plot.margin = margin(8, 12, 9, 8)
  )

base_out <- file.path(OUT_DIR, "Fig_AgeSex_Adjusted_ROC_AB_Nature")
ggsave(paste0(base_out, ".png"), fig, width = 9.0, height = 4.55,
       units = "in", dpi = 600, bg = "white")
if (capabilities("cairo")) {
  ggsave(paste0(base_out, ".pdf"), fig, width = 9.0, height = 4.55,
         units = "in", device = cairo_pdf, bg = "white")
} else {
  ggsave(paste0(base_out, ".pdf"), fig, width = 9.0, height = 4.55,
         units = "in", device = "pdf", bg = "white")
}
ggsave(paste0(base_out, ".tiff"), fig, width = 9.0, height = 4.55,
       units = "in", dpi = 600, compression = "lzw", bg = "white")

cat("\n===== Training and external ROC summary =====\n")
print(auc_table, row.names = FALSE, digits = 4)
cat("\nOutputs saved in: ", OUT_DIR, "\n", sep = "")
cat("NOTE: adjusted external AUC comes from training-fitted coefficients, ",
    "not external refitting.\n", sep = "")
cat("NOTE: training AUC is apparent performance (logistic models evaluated on their fitting data).\n")

# --------------------- 7. SINGLE-COVARIATE COMPARISONS ---------------------
# Primary comparisons use the identical age/sex-complete subset as section 4.
# Secondary comparisons use all eligible samples with the relevant covariate.
# No change to the frozen microbial OMWI model.
prepare_single <- function(df, label, variable) {
  score_col <- intersect(c("OMHI", "OMWI", "OMWI_score"), names(df))[1L]
  disease <- tolower(trimws(as.character(df$disease)))
  gender <- tolower(trimws(as.character(df$gender)))
  d <- data.frame(sample_id = trimws(as.character(df$sample_id)),
                  nonhealthy = ifelse(disease %in% c("disease", "diseased", "nonhealthy", "non-healthy", "case", "1"), 1L,
                                      ifelse(disease %in% c("healthy", "health", "control", "ctrl", "0"), 0L, NA_integer_)),
                  OMWI = suppressWarnings(as.numeric(df[[score_col]])),
                  age_num = suppressWarnings(as.numeric(df$age_num)),
                  sex = factor(ifelse(gender %in% c("male", "m", "1"), "Male",
                                      ifelse(gender %in% c("female", "f", "2"), "Female", NA_character_)), levels = c("Male", "Female")))
  ok <- !is.na(d$nonhealthy) & is.finite(d$OMWI)
  if (variable == "age_num") ok <- ok & is.finite(d$age_num) & d$age_num >= 18
  # Sex-only analysis does not require known age, but excludes known minors.
  if (variable == "sex") ok <- ok & !is.na(d$sex) & (!is.finite(d$age_num) | d$age_num >= 18)
  d[ok, , drop = FALSE]
}

run_comparison <- function(tr, ex, variable, analysis_name, adjusted_label, adjusted_color) {
  if (length(unique(tr$nonhealthy)) != 2L || length(unique(ex$nonhealthy)) != 2L) stop("Both classes required: ", analysis_name)
  raw_fit <- glm(nonhealthy ~ OMWI, data = tr, family = binomial(), na.action = na.fail)
  adj_fit <- glm(reformulate(c("OMWI", variable), response = "nonhealthy"),
                 data = tr, family = binomial(), na.action = na.fail)
  for (fit in list(raw_fit, adj_fit)) {
    if (!fit$converged || any(!is.finite(coef(fit)))) stop("Invalid coefficients: ", analysis_name)
  }
  pred <- function(d) {
    d$pred_raw <- predict(raw_fit, newdata = d, type = "response")
    d$pred_adj <- predict(adj_fit, newdata = d, type = "response")
    d
  }
  tr <- pred(tr); ex <- pred(ex)
  rr <- list(make_roc(tr, "Unadjusted", "Training"), make_roc(tr, adjusted_label, "Training"),
             make_roc(ex, "Unadjusted", "External"), make_roc(ex, adjusted_label, "External"))
  tab <- do.call(rbind, lapply(rr, `[[`, "summary"))
  tab$Analysis <- analysis_name
  curves <- do.call(rbind, lapply(rr, `[[`, "curve"))
  labs_panel <- c(paste0("a   Training cohort (n = ", nrow(tr), ")"),
                  paste0("b   External validation (n = ", nrow(ex), ")"))
  panel <- function(x) factor(ifelse(x == "Training", labs_panel[1], labs_panel[2]), levels = labs_panel)
  curves$Panel <- panel(curves$Cohort)
  ann <- tab; ann$Panel <- panel(ann$Cohort)
  ann$x <- 0.965; ann$y <- ifelse(ann$Model == "Unadjusted", 0.235, 0.145)
  ann$text <- sprintf("%s  AUC = %.4f", ann$Model, ann$AUC)
  cols <- setNames(c(COL_RAW, adjusted_color), c("Unadjusted", adjusted_label))
  # Reuse the original figure theme; replace only data, annotation and color scale.
  p <- ggplot(curves, aes(FPR, TPR, colour = Model, group = Model)) +
    geom_abline(slope = 1, intercept = 0, colour = "#BBBBBB", linetype = "dashed", linewidth = 0.55) +
    geom_step(direction = "hv", linewidth = 0.95) +
    geom_text(data = ann, aes(x = x, y = y, label = text, colour = Model), inherit.aes = FALSE,
              hjust = 1, size = 3.25, family = FONT_FAMILY, show.legend = FALSE) +
    facet_wrap(~Panel, nrow = 1) + scale_colour_manual(values = cols) +
    scale_x_continuous(limits = c(0,1), breaks = seq(0,1,.25), expand = expansion(mult = 0)) +
    scale_y_continuous(limits = c(0,1), breaks = seq(0,1,.25), expand = expansion(mult = 0)) +
    coord_fixed() + labs(x = "1 - Specificity", y = "Sensitivity") + fig$theme
  stem <- file.path(OUT_DIR, paste0("Fig_", analysis_name))
  ggsave(paste0(stem,".png"), p, width = 9, height = 4.55, dpi = 600, bg = "white")
  ggsave(paste0(stem,".pdf"), p, width = 9, height = 4.55, device = if (capabilities("cairo")) cairo_pdf else "pdf", bg = "white")
  ggsave(paste0(stem,".tiff"), p, width = 9, height = 4.55, dpi = 600, compression = "lzw", bg = "white")
  write.csv(curves, paste0(stem,"_coordinates.csv"), row.names = FALSE)
  write.csv(rbind(transform(tr, Cohort = "Training"), transform(ex, Cohort = "External")),
            paste0(stem,"_predictions.csv"), row.names = FALSE)
  co <- function(f, label) data.frame(Model = label, Term = names(coef(f)), Coefficient = unname(coef(f)))
  write.csv(rbind(co(raw_fit,"Unadjusted"), co(adj_fit,adjusted_label)), paste0(stem,"_coefficients.csv"), row.names = FALSE)
  tab
}

all_results <- list(
  transform(auc_table, Analysis = "AgeSex_common"),
  run_comparison(train, external, "age_num", "Age_common", "Age-adjusted", "#3C5488"),
  run_comparison(train, external, "sex", "Sex_common", "Sex-adjusted", "#E64B35"),
  run_comparison(prepare_single(train_raw,"Training","age_num"), prepare_single(external_raw,"External","age_num"),
                 "age_num", "Age_available", "Age-adjusted", "#3C5488"),
  run_comparison(prepare_single(train_raw,"Training","sex"), prepare_single(external_raw,"External","sex"),
                 "sex", "Sex_available", "Sex-adjusted", "#E64B35")
)
all_auc <- do.call(rbind, all_results)
all_auc$Delta_AUC <- ave(all_auc$AUC, all_auc$Analysis, all_auc$Cohort,
                         FUN = function(x) x - x[1L])
write.csv(all_auc, file.path(OUT_DIR,"All_adjustment_AUC_summary.csv"), row.names = FALSE)
print(all_auc, row.names = FALSE, digits = 6)
