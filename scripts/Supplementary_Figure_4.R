#!/usr/bin/env Rscript

# =========================================================
# Supplementary figure: robustness of OMWI across age and sex
#
# Panels (matched to the submitted manuscript legend):
#   a Age-group OMWI distributions
#   b Healthy versus disease OMWI within age groups
#   c Sex-specific OMWI distributions
#   d Healthy versus disease OMWI within sex groups
#
# Notes:
#   1) Smoking, BMI and OMWI < 0 OR panels are not included.
#   2) Age groups are 18-44, 45-59 and >=60 years.
#   3) ROC analyses are retained in the exported source-data tables but are not
#      displayed because they are not described in the submitted figure legend.
#   4) The patchwork ampersand operator is not used anywhere in figure assembly.
# =========================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)
  library(patchwork)
  library(pROC)
  library(broom)
})

options(stringsAsFactors = FALSE)

# =========================================================
# 0. Paths and analysis parameters
# =========================================================
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx)) return(default)
  if (idx == length(args)) stop("Missing value for ", flag)
  args[[idx + 1L]]
}

META_TRAIN_FILE  <- get_arg("--meta-train")
META_TEST_FILE   <- get_arg("--meta-test")
SCORE_TRAIN_FILE <- get_arg("--score-train")
SCORE_TEST_FILE  <- get_arg("--score-test")
OUT_DIR          <- get_arg("--out-dir", "results/Supplementary_Figure_4")

required_args <- c(META_TRAIN_FILE, META_TEST_FILE, SCORE_TRAIN_FILE, SCORE_TEST_FILE)
if (any(vapply(required_args, is.null, logical(1)))) {
  stop(
    "Usage: Rscript Supplementary_Figure_4.R ",
    "--meta-train meta_train.csv --meta-test meta_test.csv ",
    "--score-train omwi_train.csv --score-test omwi_test.csv ",
    "--out-dir results/Supplementary_Figure_4"
  )
}

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT_DIR, "individual_panels"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT_DIR, "tables"), showWarnings = FALSE, recursive = TRUE)

# "all" uses train and test together, matching the previous host-factor figure.
# "test" uses only the independent external validation set and is more conservative
# when the figure is interpreted as model-performance evidence.
ANALYSIS_SCOPE <- "all"

MIN_PER_CLASS_FOR_VIOLIN <- 5
MIN_PER_CLASS_FOR_AUC    <- 10
ROC_BOOT_N               <- 2000
SET_SEED                  <- 704

# =========================================================
# 1. Colors and plotting theme, matched to Fig. 4
# =========================================================
COL_1 <- "#4DBBD5"
COL_2 <- "#FEE5D9"
COL_3 <- "#FCAE91"
COL_4 <- "#FB6A4A"
COL_5 <- "#CB181D"

COL_HEALTH  <- COL_1
COL_DISEASE <- COL_5
COL_MALE    <- COL_1
COL_FEMALE  <- "#3C5488"

AGE_COLORS <- c(
  "18-44" = COL_1,
  "45-59" = COL_3,
  ">=60"   = COL_5
)

SEX_COLORS <- c(
  "Male"   = COL_MALE,
  "Female" = COL_FEMALE
)

HD_COLORS <- c(
  "Healthy" = COL_HEALTH,
  "Disease" = COL_DISEASE
)

theme_cell <- function(base_size = 10.5) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text = element_text(color = "black", size = 9),
      axis.title = element_text(color = "black", size = 10),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.6),
      strip.text = element_text(face = "bold", color = "black", size = 9.5),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
      legend.title = element_text(face = "bold", color = "black", size = 9.5),
      legend.text = element_text(color = "black", size = 8.5),
      plot.title = element_text(face = "bold", color = "black", hjust = 0.5),
      plot.tag = element_text(face = "bold", color = "black", size = 12),
      plot.tag.position = c(0.01, 0.99)
    )
}

format_p_label <- function(p, prefix = "P") {
  if (is.na(p)) return(paste0(prefix, " = NA"))
  if (p < 0.001) return(paste0(prefix, " < 0.001"))
  paste0(prefix, " = ", sprintf("%.3f", p))
}

p_to_star <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE ~ "ns"
  )
}

# Convert mixed UTF-8 / GBK / GB18030 metadata text to valid UTF-8.
# Conversion is done element by element because one metadata column may
# contain strings written with different encodings.
safe_text <- function(x) {
  x_chr <- as.character(x)
  out <- rep(NA_character_, length(x_chr))

  for (i in seq_along(x_chr)) {
    xi <- x_chr[i]
    if (is.na(xi)) next

    is_valid_utf8 <- tryCatch(
      isTRUE(validUTF8(xi)),
      error = function(e) FALSE
    )

    if (is_valid_utf8) {
      yi <- enc2utf8(xi)
    } else {
      yi <- suppressWarnings(
        iconv(xi, from = "GB18030", to = "UTF-8", sub = "byte")
      )

      if (is.na(yi)) {
        yi <- suppressWarnings(
          iconv(xi, from = "CP936", to = "UTF-8", sub = "byte")
        )
      }

      if (is.na(yi)) {
        yi <- suppressWarnings(
          iconv(xi, from = "latin1", to = "UTF-8", sub = "byte")
        )
      }
    }

    out[i] <- yi
  }

  out <- gsub("\u00A0", " ", out, fixed = TRUE)
  trimws(out)
}

# Strict numeric conversion after repairing text encoding.
# Values such as ">=18", ">=18", "18-44" or "unknown" are intentionally
# returned as NA because they are not exact individual-level ages.
safe_numeric <- function(x) {
  if (is.numeric(x)) return(as.numeric(x))

  x_fix <- safe_text(x)
  x_fix <- gsub(",", "", x_fix, fixed = TRUE)
  x_fix <- gsub("[[:space:]]+", "", x_fix, perl = TRUE)
  x_fix <- gsub("(years?|yrs?|yearold|岁)$", "", x_fix,
                ignore.case = TRUE, perl = TRUE)

  missing_text <- is.na(x_fix) | x_fix == "" |
    tolower(x_fix) %in% c("na", "n/a", "nan", "null", "none", "unknown", "missing")

  numeric_pattern <- "^[+-]?(([0-9]+\\.?[0-9]*)|(\\.[0-9]+))([eE][+-]?[0-9]+)?$"
  is_exact_numeric <- !missing_text & grepl(numeric_pattern, x_fix, perl = TRUE)

  out <- rep(NA_real_, length(x_fix))
  out[is_exact_numeric] <- suppressWarnings(as.numeric(x_fix[is_exact_numeric]))
  out
}

pick_first_existing_col <- function(df, candidates) {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

# =========================================================
# 2. Read metadata and OMWI scores
# =========================================================
required_files <- c(META_TRAIN_FILE, META_TEST_FILE, SCORE_TRAIN_FILE, SCORE_TEST_FILE)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing input files:\n", paste(missing_files, collapse = "\n"))
}

meta_train <- fread(META_TRAIN_FILE,  data.table = FALSE, check.names = FALSE)
meta_test  <- fread(META_TEST_FILE,   data.table = FALSE, check.names = FALSE)
omwi_train <- fread(SCORE_TRAIN_FILE, data.table = FALSE, check.names = FALSE)
omwi_test  <- fread(SCORE_TEST_FILE,  data.table = FALSE, check.names = FALSE)

meta_train$dataset_split <- "train"
meta_test$dataset_split  <- "test"
omwi_train$dataset_split <- "train"
omwi_test$dataset_split  <- "test"

meta_all <- bind_rows(meta_train, meta_test)
omwi_all <- bind_rows(omwi_train, omwi_test)

if (!"sample_id" %in% colnames(meta_all)) {
  stop("The metadata files do not contain sample_id.")
}
if (!"sample_id_use" %in% colnames(omwi_all)) {
  stop("The OMWI files do not contain sample_id_use.")
}
if (!"OMWI" %in% colnames(omwi_all)) {
  stop("The OMWI files do not contain OMWI.")
}

omwi_use <- omwi_all %>%
  transmute(
    sample_id = as.character(sample_id_use),
    dataset_split = as.character(dataset_split),
    OMWI = safe_numeric(OMWI)
  )

# Avoid duplicated split columns after joining.
meta_use <- meta_all %>%
  mutate(
    sample_id = as.character(sample_id),
    dataset_split = as.character(dataset_split)
  )

dat <- meta_use %>%
  left_join(
    omwi_use %>% select(sample_id, dataset_split, OMWI),
    by = c("sample_id", "dataset_split")
  )

if (all(is.na(dat$OMWI))) {
  stop("OMWI is entirely missing after joining metadata and score files.")
}

# =========================================================
# 3. Standardize outcome, age, sex and study variables
# =========================================================
if (!"disease" %in% colnames(dat)) {
  stop("The metadata do not contain the disease column.")
}

raw_disease <- tolower(safe_text(dat$disease))
dat$health_status <- dplyr::case_when(
  raw_disease %in% c("healthy", "health", "control", "ctrl", "0") ~ "Healthy",
  !is.na(raw_disease) & raw_disease != "" ~ "Disease",
  TRUE ~ NA_character_
)
dat$health_status <- factor(dat$health_status, levels = c("Healthy", "Disease"))
dat$disease01 <- ifelse(dat$health_status == "Disease", 1, ifelse(dat$health_status == "Healthy", 0, NA))

age_col <- pick_first_existing_col(dat, c("age", "Age", "AGE", "age_years", "Age_years"))
if (is.na(age_col)) stop("No age column was found.")

age_text_fixed <- safe_text(dat[[age_col]])
dat$age_num <- safe_numeric(age_text_fixed)

age_invalid <- unique(age_text_fixed[
  !is.na(age_text_fixed) & age_text_fixed != "" & is.na(dat$age_num)
])
if (length(age_invalid) > 0) {
  message(
    "Age entries excluded because they are not exact numeric ages (n = ",
    sum(!is.na(age_text_fixed) & age_text_fixed != "" & is.na(dat$age_num)),
    "). Examples: ",
    paste(utils::head(age_invalid, 10), collapse = "; ")
  )
}
dat$age_group <- dplyr::case_when(
  !is.na(dat$age_num) & dat$age_num >= 18 & dat$age_num <= 44 ~ "18-44",
  !is.na(dat$age_num) & dat$age_num >= 45 & dat$age_num <= 59 ~ "45-59",
  !is.na(dat$age_num) & dat$age_num >= 60 ~ ">=60",
  TRUE ~ NA_character_
)
dat$age_group <- factor(dat$age_group, levels = c("18-44", "45-59", ">=60"), ordered = TRUE)
dat$age_order <- as.numeric(dat$age_group)

sex_col <- pick_first_existing_col(dat, c("gender", "Gender", "sex", "Sex"))
if (is.na(sex_col)) stop("No sex or gender column was found.")

raw_sex <- tolower(safe_text(dat[[sex_col]]))
dat$sex_group <- dplyr::case_when(
  raw_sex %in% c("male", "m", "man", "boy", "1", "男") ~ "Male",
  raw_sex %in% c("female", "f", "woman", "girl", "2", "女") ~ "Female",
  TRUE ~ NA_character_
)
dat$sex_group <- factor(dat$sex_group, levels = c("Male", "Female"))

study_col <- pick_first_existing_col(
  dat,
  c("study_name", "Study", "study", "study_id", "StudyID", "project", "Project", "cohort", "Cohort", "PMID")
)

if (!is.na(study_col)) {
  dat$study_use <- factor(as.character(dat[[study_col]]))
} else {
  dat$study_use <- factor("single_study")
}

if (ANALYSIS_SCOPE == "test") {
  dat_use <- dat %>% filter(dataset_split == "test")
} else if (ANALYSIS_SCOPE == "all") {
  dat_use <- dat
} else {
  stop("ANALYSIS_SCOPE must be either 'all' or 'test'.")
}

dat_use <- dat_use %>%
  filter(!is.na(OMWI), !is.na(health_status))

cat("\nAnalysis scope:", ANALYSIS_SCOPE, "\n")
cat("Samples with OMWI and outcome:", nrow(dat_use), "\n")

# =========================================================
# 4. Statistical helper functions
# =========================================================
get_pairwise_wilcox <- function(df, group_var) {
  groups <- levels(droplevels(factor(df[[group_var]])))
  if (length(groups) < 2) return(data.frame())

  pairs <- combn(groups, 2, simplify = FALSE)
  out <- lapply(pairs, function(pair_now) {
    tmp <- df %>% filter(.data[[group_var]] %in% pair_now)
    tmp$group_tmp <- factor(tmp[[group_var]])
    p_now <- tryCatch(
      wilcox.test(OMWI ~ group_tmp, data = tmp, exact = FALSE)$p.value,
      error = function(e) NA_real_
    )
    data.frame(
      group1 = pair_now[1],
      group2 = pair_now[2],
      p = p_now,
      stringsAsFactors = FALSE
    )
  })

  bind_rows(out) %>%
    mutate(
      p_adj = p.adjust(p, method = "BH"),
      p_label = p_to_star(p_adj)
    )
}

get_stratified_wilcox <- function(df, strata_var) {
  stat_df <- df %>%
    filter(!is.na(.data[[strata_var]]), !is.na(health_status), !is.na(OMWI)) %>%
    group_by(stratum = .data[[strata_var]]) %>%
    group_modify(function(.x, .y) {
      class_n <- table(.x$health_status)
      has_both <- all(c("Healthy", "Disease") %in% names(class_n))
      enough_n <- has_both && all(class_n[c("Healthy", "Disease")] >= MIN_PER_CLASS_FOR_VIOLIN)

      if (!enough_n) {
        return(tibble(
          n_healthy = ifelse("Healthy" %in% names(class_n), as.numeric(class_n["Healthy"]), 0),
          n_disease = ifelse("Disease" %in% names(class_n), as.numeric(class_n["Disease"]), 0),
          p = NA_real_
        ))
      }

      p_now <- tryCatch(
        wilcox.test(OMWI ~ health_status, data = .x, exact = FALSE)$p.value,
        error = function(e) NA_real_
      )

      tibble(
        n_healthy = as.numeric(class_n["Healthy"]),
        n_disease = as.numeric(class_n["Disease"]),
        p = p_now
      )
    }) %>%
    ungroup() %>%
    mutate(
      p_adj = p.adjust(p, method = "BH"),
      p_label = p_to_star(p_adj)
    )

  stat_df
}

fit_interaction_lrt <- function(df, type = c("age", "sex"), adjusted = TRUE) {
  type <- match.arg(type)

  if (type == "age") {
    model_df <- df %>%
      filter(!is.na(disease01), !is.na(OMWI), !is.na(age_order)) %>%
      mutate(age_order = as.numeric(age_order))

    adjustment_terms <- character(0)
    if (adjusted && n_distinct(model_df$sex_group, na.rm = TRUE) > 1) {
      adjustment_terms <- c(adjustment_terms, "sex_group")
    }
    if (adjusted && n_distinct(model_df$study_use, na.rm = TRUE) > 1) {
      adjustment_terms <- c(adjustment_terms, "study_use")
    }

    reduced_terms <- c("OMWI", "age_order", adjustment_terms)
    full_terms <- c("OMWI * age_order", adjustment_terms)
  } else {
    model_df <- df %>%
      filter(!is.na(disease01), !is.na(OMWI), !is.na(sex_group))

    adjustment_terms <- character(0)
    if (adjusted && sum(!is.na(model_df$age_num)) >= 20) {
      adjustment_terms <- c(adjustment_terms, "age_num")
    }
    if (adjusted && n_distinct(model_df$study_use, na.rm = TRUE) > 1) {
      adjustment_terms <- c(adjustment_terms, "study_use")
    }

    reduced_terms <- c("OMWI", "sex_group", adjustment_terms)
    full_terms <- c("OMWI * sex_group", adjustment_terms)
  }

  reduced_formula <- as.formula(paste("disease01 ~", paste(reduced_terms, collapse = " + ")))
  full_formula <- as.formula(paste("disease01 ~", paste(full_terms, collapse = " + ")))

  fit_reduced <- tryCatch(
    suppressWarnings(glm(reduced_formula, data = model_df, family = binomial())),
    error = function(e) NULL
  )
  fit_full <- tryCatch(
    suppressWarnings(glm(full_formula, data = model_df, family = binomial())),
    error = function(e) NULL
  )

  if (is.null(fit_reduced) || is.null(fit_full)) {
    return(tibble(
      analysis = type,
      adjusted = adjusted,
      n = nrow(model_df),
      p_interaction = NA_real_,
      model_reduced = paste(deparse(reduced_formula), collapse = ""),
      model_full = paste(deparse(full_formula), collapse = "")
    ))
  }

  lrt <- tryCatch(anova(fit_reduced, fit_full, test = "LRT"), error = function(e) NULL)
  p_int <- if (is.null(lrt)) NA_real_ else as.numeric(lrt$`Pr(>Chi)`[2])

  tibble(
    analysis = type,
    adjusted = adjusted,
    n = nrow(model_df),
    p_interaction = p_int,
    model_reduced = paste(deparse(reduced_formula), collapse = ""),
    model_full = paste(deparse(full_formula), collapse = "")
  )
}

calc_stratified_roc <- function(df, strata_var, strata_levels, colors_map) {
  roc_objects <- list()
  roc_curve_list <- list()
  roc_stat_list <- list()

  for (stratum_now in strata_levels) {
    tmp <- df %>%
      filter(
        as.character(.data[[strata_var]]) == stratum_now,
        !is.na(disease01),
        !is.na(OMWI)
      )

    class_n <- table(tmp$disease01)
    n_healthy <- ifelse("0" %in% names(class_n), as.numeric(class_n["0"]), 0)
    n_disease <- ifelse("1" %in% names(class_n), as.numeric(class_n["1"]), 0)

    if (n_healthy < MIN_PER_CLASS_FOR_AUC || n_disease < MIN_PER_CLASS_FOR_AUC) {
      roc_stat_list[[stratum_now]] <- tibble(
        stratum = stratum_now,
        n_healthy = n_healthy,
        n_disease = n_disease,
        auc = NA_real_,
        ci_low = NA_real_,
        ci_high = NA_real_,
        sensitivity_at_OMWI0 = NA_real_,
        specificity_at_OMWI0 = NA_real_,
        balanced_accuracy_at_OMWI0 = NA_real_
      )
      next
    }

    roc_obj <- tryCatch(
      pROC::roc(
        response = tmp$disease01,
        predictor = -tmp$OMWI,
        levels = c(0, 1),
        direction = "<",
        quiet = TRUE
      ),
      error = function(e) NULL
    )

    if (is.null(roc_obj)) next

    set.seed(SET_SEED)
    auc_ci <- tryCatch(
      pROC::ci.auc(roc_obj, method = "bootstrap", boot.n = ROC_BOOT_N),
      error = function(e) c(NA_real_, NA_real_, NA_real_)
    )

    roc_df <- tibble(
      specificity = as.numeric(roc_obj$specificities),
      sensitivity = as.numeric(roc_obj$sensitivities),
      stratum = stratum_now,
      fpr = 1 - specificity
    ) %>%
      arrange(fpr, sensitivity)

    roc_objects[[stratum_now]] <- roc_obj
    roc_curve_list[[stratum_now]] <- roc_df
    pred_disease_at_zero <- ifelse(tmp$OMWI < 0, 1, 0)
    sensitivity_at_zero <- mean(pred_disease_at_zero[tmp$disease01 == 1] == 1, na.rm = TRUE)
    specificity_at_zero <- mean(pred_disease_at_zero[tmp$disease01 == 0] == 0, na.rm = TRUE)
    balanced_accuracy_at_zero <- mean(c(sensitivity_at_zero, specificity_at_zero), na.rm = TRUE)

    roc_stat_list[[stratum_now]] <- tibble(
      stratum = stratum_now,
      n_healthy = n_healthy,
      n_disease = n_disease,
      auc = as.numeric(pROC::auc(roc_obj)),
      ci_low = as.numeric(auc_ci[1]),
      ci_high = as.numeric(auc_ci[3]),
      sensitivity_at_OMWI0 = sensitivity_at_zero,
      specificity_at_OMWI0 = specificity_at_zero,
      balanced_accuracy_at_OMWI0 = balanced_accuracy_at_zero
    )
  }

  stats_df <- bind_rows(roc_stat_list) %>%
    mutate(
      stratum = factor(stratum, levels = strata_levels),
      color = unname(colors_map[as.character(stratum)]),
      legend_label = ifelse(
        is.na(auc),
        paste0(as.character(stratum), ": insufficient n"),
        paste0(
          as.character(stratum), ": AUC ", sprintf("%.2f", auc),
          " (", sprintf("%.2f", ci_low), "-", sprintf("%.2f", ci_high), ")"
        )
      )
    )

  if (length(roc_curve_list) == 0) {
    curve_df <- tibble(
      specificity = numeric(0),
      sensitivity = numeric(0),
      stratum = factor(character(0), levels = strata_levels),
      fpr = numeric(0),
      legend_label = character(0)
    )
  } else {
    curve_df <- bind_rows(roc_curve_list) %>%
      mutate(
        stratum = factor(stratum, levels = strata_levels),
        legend_label = stats_df$legend_label[match(as.character(stratum), as.character(stats_df$stratum))]
      )
  }

  list(
    roc_objects = roc_objects,
    curve_df = curve_df,
    stats_df = stats_df
  )
}

compare_auc_delong <- function(roc_objects) {
  valid_names <- names(roc_objects)
  if (length(valid_names) < 2) return(tibble())

  pairs <- combn(valid_names, 2, simplify = FALSE)
  out <- lapply(pairs, function(pair_now) {
    test_now <- tryCatch(
      pROC::roc.test(
        roc_objects[[pair_now[1]]],
        roc_objects[[pair_now[2]]],
        method = "delong",
        paired = FALSE
      ),
      error = function(e) NULL
    )

    p_now <- if (is.null(test_now)) NA_real_ else as.numeric(test_now$p.value)

    tibble(
      group1 = pair_now[1],
      group2 = pair_now[2],
      p = p_now
    )
  })

  bind_rows(out) %>%
    mutate(p_adj = p.adjust(p, method = "BH"))
}

# =========================================================
# 5. Plot functions
# =========================================================
plot_overall_distribution <- function(df, group_var, group_levels, colors_map,
                                      x_title, p_label, tag_label,
                                      rotate_x = FALSE, draw_median_line = FALSE,
                                      pairwise_df = NULL) {
  sub <- df %>%
    filter(!is.na(.data[[group_var]]), !is.na(OMWI)) %>%
    mutate(group_plot = factor(as.character(.data[[group_var]]), levels = group_levels)) %>%
    filter(!is.na(group_plot))

  if (n_distinct(sub$group_plot) < 2) return(NULL)

  median_df <- sub %>%
    group_by(group_plot) %>%
    summarise(median_OMWI = median(OMWI, na.rm = TRUE), .groups = "drop")

  p <- ggplot(sub, aes(x = group_plot, y = OMWI, fill = group_plot, color = group_plot)) +
    geom_violin(width = 0.88, trim = FALSE, linewidth = 0.65, alpha = 0.65) +
    geom_jitter(width = 0.06, size = 0.9, alpha = 0.35, show.legend = FALSE) +
    geom_boxplot(
      width = 0.20,
      outlier.shape = NA,
      fill = "transparent",
      color = "black",
      linewidth = 0.45
    ) +
    geom_point(
      data = median_df,
      aes(x = group_plot, y = median_OMWI),
      inherit.aes = FALSE,
      shape = 21,
      fill = "white",
      color = "black",
      stroke = 0.8,
      size = 2.5
    ) +
    scale_fill_manual(values = colors_map, drop = FALSE) +
    scale_color_manual(values = colors_map, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.32))) +
    annotate(
      "text",
      x = Inf,
      y = Inf,
      label = p_label,
      hjust = 1.04,
      vjust = 1.25,
      size = 3.2
    ) +
    labs(x = x_title, y = "OMWI", tag = tag_label) +
    theme_cell(10.5) +
    theme(
      legend.position = "none",
      axis.text.x = if (rotate_x) {
        element_text(angle = 25, hjust = 1, vjust = 1)
      } else {
        element_text(angle = 0, hjust = 0.5, vjust = 0.5)
      },
      plot.margin = margin(6, 7, 5.5, 5.5)
    )

  if (draw_median_line) {
    p <- p +
      geom_line(
        data = median_df,
        aes(x = group_plot, y = median_OMWI, group = 1),
        inherit.aes = FALSE,
        color = "black",
        linewidth = 0.75
      )
  }

  if (!is.null(pairwise_df) && nrow(pairwise_df) > 0) {
    pw_plot <- pairwise_df %>%
      filter(!is.na(p_adj), p_label != "ns", p_label != "") %>%
      mutate(
        x1 = match(group1, group_levels),
        x2 = match(group2, group_levels)
      ) %>%
      filter(!is.na(x1), !is.na(x2)) %>%
      arrange(x1, x2)

    if (nrow(pw_plot) > 0) {
      y_max <- max(sub$OMWI, na.rm = TRUE)
      y_min <- min(sub$OMWI, na.rm = TRUE)
      y_rng <- y_max - y_min
      if (!is.finite(y_rng) || y_rng <= 0) y_rng <- 1

      pw_plot <- pw_plot %>%
        mutate(
          y_bracket = y_max + (0.07 + 0.09 * (row_number() - 1)) * y_rng,
          y_label = y_bracket + 0.025 * y_rng,
          x_mid = (x1 + x2) / 2
        )

      p <- p +
        geom_segment(
          data = pw_plot,
          aes(x = x1, xend = x2, y = y_bracket, yend = y_bracket),
          inherit.aes = FALSE,
          linewidth = 0.4,
          color = "black"
        ) +
        geom_segment(
          data = pw_plot,
          aes(x = x1, xend = x1, y = y_bracket, yend = y_bracket - 0.02 * y_rng),
          inherit.aes = FALSE,
          linewidth = 0.4,
          color = "black"
        ) +
        geom_segment(
          data = pw_plot,
          aes(x = x2, xend = x2, y = y_bracket, yend = y_bracket - 0.02 * y_rng),
          inherit.aes = FALSE,
          linewidth = 0.4,
          color = "black"
        ) +
        geom_text(
          data = pw_plot,
          aes(x = x_mid, y = y_label, label = p_label),
          inherit.aes = FALSE,
          size = 3.4,
          color = "black"
        )
    }
  }

  p
}

plot_stratified_violin <- function(df, strata_var, strata_levels, tag_label) {
  sub <- df %>%
    filter(!is.na(.data[[strata_var]]), !is.na(health_status), !is.na(OMWI)) %>%
    mutate(
      stratum = factor(as.character(.data[[strata_var]]), levels = strata_levels),
      health_status = factor(health_status, levels = c("Healthy", "Disease"))
    ) %>%
    filter(!is.na(stratum))

  stat_df <- get_stratified_wilcox(sub, "stratum")
  keep_strata <- stat_df %>%
    filter(n_healthy >= MIN_PER_CLASS_FOR_VIOLIN, n_disease >= MIN_PER_CLASS_FOR_VIOLIN) %>%
    pull(stratum) %>%
    as.character()

  sub <- sub %>%
    filter(as.character(stratum) %in% keep_strata) %>%
    mutate(stratum = droplevels(stratum))

  stat_df <- stat_df %>%
    filter(as.character(stratum) %in% keep_strata)

  if (nrow(sub) == 0) return(list(plot = NULL, stats = stat_df))

  y_stats <- sub %>%
    group_by(stratum) %>%
    summarise(
      y_max = max(OMWI, na.rm = TRUE),
      y_min = min(OMWI, na.rm = TRUE),
      y_rng = y_max - y_min,
      .groups = "drop"
    ) %>%
    mutate(y_rng = ifelse(y_rng <= 0, 1, y_rng))

  stat_plot <- stat_df %>%
    left_join(y_stats, by = "stratum") %>%
    mutate(
      # Use numeric x positions so that the annotation layers do not
      # inherit x = health_status from the main ggplot data.
      x1 = 1,
      x2 = 2,
      x_mid = 1.5,
      y_bracket = y_max + 0.13 * y_rng,
      y_tip = y_bracket - 0.025 * y_rng,
      y_label = y_bracket + 0.025 * y_rng,
      p_label_plot = ifelse(p_label == "ns", "ns", p_label)
    )

  p <- ggplot(sub, aes(x = health_status, y = OMWI, fill = health_status, color = health_status)) +
    geom_violin(width = 0.88, trim = FALSE, linewidth = 0.65, alpha = 0.65) +
    geom_jitter(width = 0.06, size = 0.8, alpha = 0.30, show.legend = FALSE) +
    geom_boxplot(
      width = 0.20,
      outlier.shape = NA,
      fill = "transparent",
      color = "black",
      linewidth = 0.45
    ) +
    # Draw the significance bracket manually. This avoids compatibility
    # problems in older ggpubr versions where stat_pvalue_manual may
    # inherit the parent health_status aesthetic.
    geom_segment(
      data = stat_plot,
      aes(x = x1, xend = x2, y = y_bracket, yend = y_bracket),
      inherit.aes = FALSE,
      linewidth = 0.4,
      color = "black"
    ) +
    geom_segment(
      data = stat_plot,
      aes(x = x1, xend = x1, y = y_bracket, yend = y_tip),
      inherit.aes = FALSE,
      linewidth = 0.4,
      color = "black"
    ) +
    geom_segment(
      data = stat_plot,
      aes(x = x2, xend = x2, y = y_bracket, yend = y_tip),
      inherit.aes = FALSE,
      linewidth = 0.4,
      color = "black"
    ) +
    geom_text(
      data = stat_plot,
      aes(x = x_mid, y = y_label, label = p_label_plot),
      inherit.aes = FALSE,
      size = 3.3,
      color = "black"
    ) +
    facet_grid(. ~ stratum, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = HD_COLORS, drop = FALSE) +
    scale_color_manual(values = HD_COLORS, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
    labs(x = NULL, y = "OMWI", tag = tag_label) +
    theme_cell(10.5) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1),
      panel.spacing.x = grid::unit(0.25, "lines"),
      plot.margin = margin(6, 6, 5.5, 5.5)
    )

  list(plot = p, stats = stat_df)
}

plot_roc_panel <- function(roc_result, colors_map, interaction_label, tag_label) {
  curve_df <- roc_result$curve_df
  stats_df <- roc_result$stats_df %>% filter(!is.na(auc))

  if (nrow(curve_df) == 0) return(NULL)

  legend_values <- unname(colors_map[as.character(stats_df$stratum)])
  names(legend_values) <- stats_df$legend_label

  curve_df <- curve_df %>%
    mutate(legend_label = factor(legend_label, levels = stats_df$legend_label))

  ggplot(curve_df, aes(x = fpr, y = sensitivity, color = legend_label)) +
    geom_step(linewidth = 1.0, direction = "vh") +
    geom_abline(intercept = 0, slope = 1, linetype = 2, color = "grey70", linewidth = 0.6) +
    annotate(
      "text",
      x = 0.98,
      y = 0.04,
      label = interaction_label,
      hjust = 1,
      vjust = 0,
      size = 3.1
    ) +
    scale_color_manual(values = legend_values, drop = FALSE) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    labs(
      x = "1 - Specificity",
      y = "Sensitivity",
      color = NULL,
      tag = tag_label
    ) +
    theme_cell(10.5) +
    theme(
      legend.position = c(0.67, 0.27),
      legend.background = element_rect(fill = scales::alpha("white", 0.82), color = NA),
      legend.key.height = grid::unit(0.45, "lines"),
      legend.spacing.y = grid::unit(0.05, "lines"),
      panel.grid = element_blank(),
      plot.margin = margin(6, 8, 5.5, 5.5)
    )
}

# =========================================================
# 6. Overall age and sex tests
# =========================================================
age_overall_df <- dat_use %>%
  filter(!is.na(age_group), !is.na(OMWI))

age_trend_test <- tryCatch(
  suppressWarnings(cor.test(as.numeric(age_overall_df$age_group), age_overall_df$OMWI, method = "spearman", exact = FALSE)),
  error = function(e) NULL
)
age_trend_p <- if (is.null(age_trend_test)) NA_real_ else age_trend_test$p.value

sex_overall_df <- dat_use %>%
  filter(!is.na(sex_group), !is.na(OMWI))

sex_wilcox <- tryCatch(
  wilcox.test(OMWI ~ sex_group, data = sex_overall_df, exact = FALSE),
  error = function(e) NULL
)
sex_p <- if (is.null(sex_wilcox)) NA_real_ else sex_wilcox$p.value
# =========================================================
# Overall tests shown in panels a and d
# =========================================================
overall_host_test_summary <- bind_rows(
  tibble(
    variable = "Age category",
    comparison = "18-44 vs 45-59 vs >=60",
    test = "Spearman rank correlation using ordered age groups",
    n = nrow(age_overall_df),
    statistic_name = "Spearman rho",
    statistic = if (is.null(age_trend_test)) {
      NA_real_
    } else {
      as.numeric(age_trend_test$estimate)
    },
    p_value = age_trend_p,
    figure_label = format_p_label(age_trend_p, "Trend P")
  ),
  tibble(
    variable = "Gender",
    comparison = "Male vs Female",
    test = "Two-sided Wilcoxon rank-sum test",
    n = nrow(sex_overall_df),
    statistic_name = "Wilcoxon W",
    statistic = if (is.null(sex_wilcox)) {
      NA_real_
    } else {
      as.numeric(sex_wilcox$statistic)
    },
    p_value = sex_p,
    figure_label = format_p_label(sex_p, "P")
  )
)

print(overall_host_test_summary)
age_pairwise <- get_pairwise_wilcox(age_overall_df, "age_group")
sex_pairwise <- get_pairwise_wilcox(sex_overall_df, "sex_group")

# =========================================================
# 7. Stratified ROC analyses and heterogeneity tests
# =========================================================
age_levels <- c("18-44", "45-59", ">=60")
sex_levels <- c("Male", "Female")

age_roc <- calc_stratified_roc(
  df = dat_use,
  strata_var = "age_group",
  strata_levels = age_levels,
  colors_map = AGE_COLORS
)

sex_roc <- calc_stratified_roc(
  df = dat_use,
  strata_var = "sex_group",
  strata_levels = sex_levels,
  colors_map = SEX_COLORS
)

age_delong <- compare_auc_delong(age_roc$roc_objects) %>% mutate(analysis = "age")
sex_delong <- compare_auc_delong(sex_roc$roc_objects) %>% mutate(analysis = "sex")

interaction_tests <- bind_rows(
  fit_interaction_lrt(dat_use, type = "age", adjusted = FALSE),
  fit_interaction_lrt(dat_use, type = "age", adjusted = TRUE),
  fit_interaction_lrt(dat_use, type = "sex", adjusted = FALSE),
  fit_interaction_lrt(dat_use, type = "sex", adjusted = TRUE)
)

age_p_interaction <- interaction_tests %>%
  filter(analysis == "age", adjusted) %>%
  pull(p_interaction)
if (length(age_p_interaction) == 0 || is.na(age_p_interaction[1])) {
  age_p_interaction <- interaction_tests %>%
    filter(analysis == "age", !adjusted) %>%
    pull(p_interaction)
}
age_p_interaction <- age_p_interaction[1]

sex_p_interaction <- interaction_tests %>%
  filter(analysis == "sex", adjusted) %>%
  pull(p_interaction)
if (length(sex_p_interaction) == 0 || is.na(sex_p_interaction[1])) {
  sex_p_interaction <- interaction_tests %>%
    filter(analysis == "sex", !adjusted) %>%
    pull(p_interaction)
}
sex_p_interaction <- sex_p_interaction[1]

age_interaction_label <- format_p_label(age_p_interaction, "Interaction trend P")
sex_interaction_label <- format_p_label(sex_p_interaction, "Interaction P")

# =========================================================
# 7.1 Subgroup sample-count tables
# Define these before figure rendering so they remain available even if
# a graphics device or plot layer fails later.
# =========================================================
age_counts <- dat_use %>%
  filter(!is.na(age_group), !is.na(health_status), !is.na(OMWI)) %>%
  count(age_group, health_status, name = "n") %>%
  complete(
    age_group = factor(age_levels, levels = age_levels, ordered = TRUE),
    health_status = factor(c("Healthy", "Disease"), levels = c("Healthy", "Disease")),
    fill = list(n = 0)
  )

sex_counts <- dat_use %>%
  filter(!is.na(sex_group), !is.na(health_status), !is.na(OMWI)) %>%
  count(sex_group, health_status, name = "n") %>%
  complete(
    sex_group = factor(sex_levels, levels = sex_levels),
    health_status = factor(c("Healthy", "Disease"), levels = c("Healthy", "Disease")),
    fill = list(n = 0)
  )

# =========================================================
# 8. Build panels
# =========================================================
p_age_overall <- plot_overall_distribution(
  df = dat_use,
  group_var = "age_group",
  group_levels = age_levels,
  colors_map = AGE_COLORS,
  x_title = "Age category",
  p_label = format_p_label(age_trend_p, "Trend P"),
  tag_label = "a",
  rotate_x = TRUE,
  draw_median_line = TRUE,
  pairwise_df = age_pairwise
)

age_stratified <- plot_stratified_violin(
  df = dat_use,
  strata_var = "age_group",
  strata_levels = age_levels,
  tag_label = "b"
)
p_age_stratified <- age_stratified$plot

p_age_roc <- plot_roc_panel(
  roc_result = age_roc,
  colors_map = AGE_COLORS,
  interaction_label = age_interaction_label,
  tag_label = "c"
)

p_sex_overall <- plot_overall_distribution(
  df = dat_use,
  group_var = "sex_group",
  group_levels = sex_levels,
  colors_map = SEX_COLORS,
  x_title = "Gender",
  p_label = format_p_label(sex_p, "P"),
  tag_label = "c",
  rotate_x = FALSE,
  draw_median_line = FALSE,
  pairwise_df = sex_pairwise
)

sex_stratified <- plot_stratified_violin(
  df = dat_use,
  strata_var = "sex_group",
  strata_levels = sex_levels,
  tag_label = "d"
)
p_sex_stratified <- sex_stratified$plot

p_sex_roc <- plot_roc_panel(
  roc_result = sex_roc,
  colors_map = SEX_COLORS,
  interaction_label = sex_interaction_label,
  tag_label = "f"
)

plot_list <- list(
  p_age_overall,
  p_age_stratified,
  p_sex_overall,
  p_sex_stratified
)

if (any(vapply(plot_list, is.null, logical(1)))) {
  stop("At least one panel could not be generated. Check subgroup sample counts and metadata columns.")
}

# Figure assembly deliberately uses wrap_plots and does not use the patchwork ampersand operator.
figure_all <- patchwork::wrap_plots(
  plots = plot_list,
  ncol = 2,
  widths = c(0.92, 1.32),
  heights = c(1, 1)
)

# =========================================================
# 9. Export figure and individual panels
# =========================================================
ggsave(
  filename = file.path(OUT_DIR, "Supplementary_Figure_4.pdf"),
  plot = figure_all,
  width = 7.2,
  height = 6.4,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = file.path(OUT_DIR, "Supplementary_Figure_4.png"),
  plot = figure_all,
  width = 7.2,
  height = 6.4,
  units = "in",
  dpi = 600,
  bg = "white"
)

panel_names <- c(
  "a_age_overall",
  "b_age_stratified",
  "c_gender_overall",
  "d_gender_stratified"
)

panel_widths <- c(4.0, 5.1, 4.0, 4.4)

for (i in seq_along(plot_list)) {
  ggsave(
    filename = file.path(OUT_DIR, "individual_panels", paste0(panel_names[i], ".pdf")),
    plot = plot_list[[i]],
    width = panel_widths[i],
    height = 3.4,
    units = "in",
    device = cairo_pdf
  )
  ggsave(
    filename = file.path(OUT_DIR, "individual_panels", paste0(panel_names[i], ".png")),
    plot = plot_list[[i]],
    width = panel_widths[i],
    height = 3.4,
    units = "in",
    dpi = 600,
    bg = "white"
  )
}

# =========================================================
# 10. Export statistical tables
# =========================================================
fwrite(age_counts, file.path(OUT_DIR, "tables", "age_group_sample_counts.csv"))
fwrite(sex_counts, file.path(OUT_DIR, "tables", "gender_sample_counts.csv"))
fwrite(age_pairwise, file.path(OUT_DIR, "tables", "age_overall_pairwise_wilcoxon.csv"))
fwrite(sex_pairwise, file.path(OUT_DIR, "tables", "gender_overall_wilcoxon.csv"))
fwrite(age_stratified$stats, file.path(OUT_DIR, "tables", "age_stratified_healthy_disease_wilcoxon.csv"))
fwrite(sex_stratified$stats, file.path(OUT_DIR, "tables", "gender_stratified_healthy_disease_wilcoxon.csv"))
fwrite(age_roc$stats_df, file.path(OUT_DIR, "tables", "age_stratified_auc.csv"))
fwrite(sex_roc$stats_df, file.path(OUT_DIR, "tables", "gender_stratified_auc.csv"))
fwrite(age_delong, file.path(OUT_DIR, "tables", "age_auc_pairwise_delong.csv"))
fwrite(sex_delong, file.path(OUT_DIR, "tables", "gender_auc_pairwise_delong.csv"))
fwrite(interaction_tests, file.path(OUT_DIR, "tables", "omwi_age_gender_interaction_tests.csv"))
fwrite(
  overall_host_test_summary,
  file.path(
    OUT_DIR,
    "tables",
    "overall_age_gender_tests_shown_in_figure.csv"
  )
)
capture.output(sessionInfo(), file = file.path(OUT_DIR, "sessionInfo.txt"))
cat("\n[DONE] Output directory: ", OUT_DIR, "\n", sep = "")
cat("Age discrimination test: ", age_interaction_label, "\n", sep = "")
cat("Sex discrimination test: ", sex_interaction_label, "\n", sep = "")
