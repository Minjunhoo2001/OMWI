#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(pROC)
  library(purrr)
})

options(stringsAsFactors = FALSE)

# =========================================================
# 0. Path settings
# =========================================================
ROOT_DIR <- "/share/home/HeMinjun/metagenomic/Oral_proxy_v2/04.lasso64"

COEF_DIR      <- file.path(ROOT_DIR, "coefficients")
OMHI_DIR      <- file.path(ROOT_DIR, "omhi_scores")
CV_OOF_DIR    <- file.path(ROOT_DIR, "cv_oof")
LOOCV_OOF_DIR <- file.path(ROOT_DIR, "loocv_oof")

OUT_DIR <- "/share/home/HeMinjun/metagenomic/GitHub/results/Figure2"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

LEVEL <- "species"
BOOT_N <- 2000

# =========================================================
# 1. color and theme
# =========================================================
COL_HEALTH  <- "#4DBBD5"
COL_DISEASE <- "#E64B35"

COL_LOOCV <- "#7F7F7F"
COL_CV    <- "#F4A261"
COL_TEST  <- "#66C2A5"

theme_fig <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", linewidth = 0.2),
      axis.line = element_line(color = "black", linewidth = 0.2),
      axis.ticks = element_line(color = "black", linewidth = 0.2),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      plot.title = element_blank(),
      legend.title = element_blank(),
      legend.background = element_rect(fill = "white", colour = "black", linewidth = 0.4),
      legend.key = element_rect(fill = "white", colour = NA),
      strip.background = element_blank(),
      strip.text = element_text(color = "black", face = "bold")
    )
}

# =========================================================
# 2. General functions
# =========================================================
alpha_tag <- function(a) {
  paste0("alpha", gsub("\\.", "", sprintf("%.2f", a)))
}

pick_file <- function(dir, patterns, must = TRUE) {
  if (!dir.exists(dir)) {
    if (must) stop("The folder does not exist: ", dir)
    return(NULL)
  }
  
  for (pat in patterns) {
    fs <- list.files(
      dir,
      pattern = pat,
      full.names = TRUE,
      ignore.case = TRUE
    )
    if (length(fs) > 0) return(fs[1])
  }
  
  if (must) {
    stop(
      "未找到文件。目录: ", dir,
      "\nCandidate pattern:\n",
      paste(patterns, collapse = "\n")
    )
  }
  
  NULL
}

pick_col <- function(df, candidates, required = TRUE) {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) == 0) {
    if (required) {
      stop("Column name not found. Candidate columns are: ", paste(candidates, collapse = ", "))
    } else {
      return(NULL)
    }
  }
  hit[1]
}

check_required_cols <- function(df, cols, label = "data") {
  miss <- setdiff(cols, colnames(df))
  if (length(miss) > 0) {
    stop(label, " Required columns are missing: ", paste(miss, collapse = ", "))
  }
}

std_binary_label <- function(x) {
  x0 <- trimws(tolower(as.character(x)))
  num <- suppressWarnings(as.numeric(x0))
  
  out <- rep(NA_real_, length(x0))
  
  out[!is.na(num) & num == 1] <- 1
  out[!is.na(num) & num == 0] <- 0
  
  out[x0 %in% c("healthy", "health", "control", "ctrl", "normal")] <- 1
  out[x0 %in% c("disease", "case", "patient", "unhealthy", "oral disease", "nonhealthy", "non-healthy")] <- 0
  
  out
}

get_health01 <- function(df, set_name) {
  if (str_detect(set_name, "OOF") && "y_true" %in% colnames(df)) {
    return(std_binary_label(df$y_true))
  }
  
  if ("y_healthy" %in% colnames(df)) {
    return(std_binary_label(df$y_healthy))
  }
  
  if ("health01" %in% colnames(df)) {
    return(std_binary_label(df$health01))
  }
  
  if ("disease" %in% colnames(df)) {
    return(std_binary_label(df$disease))
  }
  
  if ("label" %in% colnames(df)) {
    return(std_binary_label(df$label))
  }
  
  stop(set_name, " The label could not be found for y_true / y_healthy / disease / label.")
}

read_score_one <- function(file, set_name, need_omhi = TRUE) {
  x <- fread(file, data.table = FALSE, check.names = FALSE)
  x$.source_file <- basename(file)
  x$set <- set_name
  
  check_required_cols(x, "prob_healthy", label = set_name)
  
  if (need_omhi) {
    check_required_cols(x, "OMHI", label = set_name)
  }
  
  x$prob_healthy <- as.numeric(x$prob_healthy)
  
  if ("OMHI" %in% colnames(x)) {
    x$OMHI <- as.numeric(x$OMHI)
  } else {
    x$OMHI <- NA_real_
  }
  
  x$health01 <- get_health01(x, set_name)
  
  x <- x %>%
    mutate(
      Disease = factor(
        ifelse(health01 == 1, "Healthy", "Disease"),
        levels = c("Healthy", "Disease")
      )
    ) %>%
    filter(
      !is.na(prob_healthy),
      !is.na(health01),
      !is.na(Disease)
    )
  
  if (need_omhi) {
    x <- x %>% filter(!is.na(OMHI))
  }
  
  x
}

safe_auc <- function(y, score) {
  keep <- is.finite(y) & is.finite(score)
  y <- y[keep]
  score <- score[keep]
  
  if (length(unique(y)) < 2) return(NA_real_)
  
  as.numeric(
    pROC::roc(
      response = y,
      predictor = score,
      levels = c(0, 1),
      direction = "<",
      quiet = TRUE
    )$auc
  )
}

safe_auc_ci <- function(y, score, conf.level = 0.95) {
  keep <- is.finite(y) & is.finite(score)
  y <- y[keep]
  score <- score[keep]
  
  if (length(unique(y)) < 2) {
    return(data.frame(
      auc = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_
    ))
  }
  
  roc_obj <- pROC::roc(
    response = y,
    predictor = score,
    levels = c(0, 1),
    direction = "<",
    quiet = TRUE
  )
  
  auc_val <- as.numeric(pROC::auc(roc_obj))
  ci_val <- as.numeric(pROC::ci.auc(roc_obj, conf.level = conf.level))
  
  data.frame(
    auc = auc_val,
    ci_low = ci_val[1],
    ci_high = ci_val[3]
  )
}

safe_roc_with_ci_df <- function(y, score, label, conf.level = 0.95, boot.n = 2000) {
  keep <- is.finite(y) & is.finite(score)
  y <- y[keep]
  score <- score[keep]
  
  if (length(unique(y)) < 2) {
    stop("ROC calculation failed:", label, "There are fewer than two categories")
  }
  
  roc_obj <- pROC::roc(
    response = y,
    predictor = score,
    levels = c(0, 1),
    direction = "<",
    quiet = TRUE
  )
  
  roc_df <- data.frame(
    fpr = 1 - roc_obj$specificities,
    tpr = roc_obj$sensitivities,
    set = label
  ) %>%
    arrange(fpr, tpr)
  
  spec_grid <- seq(0, 1, by = 0.01)
  
  ci_mat <- tryCatch(
    pROC::ci.se(
      roc_obj,
      specificities = spec_grid,
      conf.level = conf.level,
      boot.n = boot.n,
      progress = "none"
    ),
    error = function(e) NULL
  )
  
  if (is.null(ci_mat)) {
    ci_df <- data.frame(
      specificity = numeric(0),
      tpr_low = numeric(0),
      tpr_mid = numeric(0),
      tpr_high = numeric(0),
      fpr = numeric(0),
      set = character(0)
    )
  } else {
    ci_df <- data.frame(
      specificity = as.numeric(rownames(ci_mat)),
      tpr_low  = ci_mat[, 1],
      tpr_mid  = ci_mat[, 2],
      tpr_high = ci_mat[, 3]
    ) %>%
      mutate(
        fpr = 1 - specificity,
        set = label
      ) %>%
      arrange(fpr)
  }
  
  list(
    roc = roc_df,
    ci = ci_df,
    roc_obj = roc_obj
  )
}

format_p <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "P = NA",
    p < 0.001 ~ "P < 0.001",
    TRUE ~ paste0("P = ", signif(p, 2))
  )
}

calc_cliff_delta <- function(x_healthy, x_disease) {
  x <- x_healthy[is.finite(x_healthy)]
  y <- x_disease[is.finite(x_disease)]
  
  if (length(x) == 0 || length(y) == 0) return(NA_real_)
  
  wt <- suppressWarnings(wilcox.test(x, y, exact = FALSE))
  W <- as.numeric(wt$statistic)
  U <- W - length(x) * (length(x) + 1) / 2
  
  2 * U / (length(x) * length(y)) - 1
}

calc_balanced_accuracy <- function(y_true, score) {
  pred_healthy <- ifelse(score > 0, 1, 0)
  
  healthy_ccr <- if (sum(y_true == 1, na.rm = TRUE) > 0) {
    mean(pred_healthy[y_true == 1] == 1, na.rm = TRUE)
  } else {
    NA_real_
  }
  
  disease_ccr <- if (sum(y_true == 0, na.rm = TRUE) > 0) {
    mean(pred_healthy[y_true == 0] == 0, na.rm = TRUE)
  } else {
    NA_real_
  }
  
  mean(c(healthy_ccr, disease_ccr), na.rm = TRUE)
}

save_panel <- function(p, prefix, width, height) {
  ggsave(
    filename = file.path(OUT_DIR, paste0(prefix, ".png")),
    plot = p,
    width = width,
    height = height,
    dpi = 600,
    bg = "white"
  )
  
  ggsave(
    filename = file.path(OUT_DIR, paste0(prefix, ".pdf")),
    plot = p,
    width = width,
    height = height,
    bg = "white"
  )
}

# =========================================================
# 3. Find the best species model parameters
# =========================================================
best_file <- file.path(ROOT_DIR, "LASSO_best_model_by_level_full.csv")

if (file.exists(best_file)) {
  best_tab <- fread(best_file, data.table = FALSE)
  best_species <- best_tab %>% filter(level == LEVEL)
  
  if (nrow(best_species) != 1) {
    stop("best_table 里 species 层级不是唯一一行，请检查: ", best_file)
  }
  
  alpha_val <- best_species$alpha[1]
  prev_val  <- best_species$prev_cutoff[1]
} else {
  warning("未找到 best_table，默认使用 alpha = 0.5")
  alpha_val <- 0.5
  prev_val <- NA_real_
}

atag <- alpha_tag(alpha_val)

cat("\nBest species model:\n")
cat("alpha =", alpha_val, "\n")
cat("prev_cutoff =", prev_val, "\n")
cat("alpha tag =", atag, "\n")

# =========================================================
# 4. Read data
# =========================================================
cat("\n================ Read OMHI / OOF score files ================\n")

train_file <- pick_file(
  OMHI_DIR,
  c(
    paste0("^OMHI_train_", LEVEL, "_", atag, "\\.csv$"),
    paste0("OMHI_train_", LEVEL, ".*", atag),
    paste0("OMHI_train_", LEVEL)
  ),
  must = TRUE
)

test_file <- pick_file(
  OMHI_DIR,
  c(
    paste0("^OMHI_test_", LEVEL, "_", atag, "\\.csv$"),
    paste0("OMHI_test_", LEVEL, ".*", atag),
    paste0("OMHI_test_", LEVEL)
  ),
  must = TRUE
)

loocv_file <- pick_file(
  LOOCV_OOF_DIR,
  c(
    paste0("^OOF_LOOCV_", LEVEL, "_", atag, "_best_model\\.csv$"),
    paste0("OOF_LOOCV_", LEVEL, ".*", atag),
    paste0("OOF_LOOCV_", LEVEL)
  ),
  must = TRUE
)

cv10_file <- pick_file(
  CV_OOF_DIR,
  c(
    paste0("^OOF_10fold_repeated_", LEVEL, "_", atag, "_best_model\\.csv$"),
    paste0("^OOF_10fold_", LEVEL, "_", atag, "_best_model\\.csv$"),
    paste0("OOF_10fold_repeated_", LEVEL, ".*", atag),
    paste0("OOF_10fold_", LEVEL, ".*", atag),
    paste0("OOF_10fold.*", LEVEL)
  ),
  must = TRUE
)

cat("Training file: ", basename(train_file), "\n")
cat("External file: ", basename(test_file), "\n")
cat("LOOCV file: ", basename(loocv_file), "\n")
cat("10-fold file: ", basename(cv10_file), "\n")

train_df <- read_score_one(train_file, "Training", need_omhi = TRUE)
test_df  <- read_score_one(test_file, "External validation", need_omhi = TRUE)
loocv_df <- read_score_one(loocv_file, "LOOCV OOF", need_omhi = TRUE)
cv10_df  <- read_score_one(cv10_file, "10-fold CV OOF", need_omhi = TRUE)

score_df_all <- bind_rows(
  train_df,
  loocv_df,
  cv10_df,
  test_df
) %>%
  mutate(
    set = factor(
      set,
      levels = c("Training", "LOOCV OOF", "10-fold CV OOF", "External validation")
    )
  )

cat("\n检查 set × Disease 样本数：\n")
print(table(score_df_all$set, score_df_all$Disease, useNA = "ifany"))

# =========================================================
# 5. Fig2A: Lollipop
# =========================================================
cat("\n================ Fig2A: coefficient lollipop ================\n")

coef_file <- pick_file(
  COEF_DIR,
  c(
    paste0("^coef_", LEVEL, "_", atag, "_best_model\\.csv$"),
    paste0("coef_", LEVEL, ".*", atag),
    paste0("coef_", LEVEL)
  ),
  must = TRUE
)

coef_df <- fread(coef_file, data.table = FALSE, check.names = FALSE)

feature_col <- pick_col(coef_df, c("feature", "taxon", "taxa", "species", "name"))
coef_col    <- pick_col(coef_df, c("coefficient", "coef", "estimate", "beta", "weight"))

coef_plot_df <- coef_df %>%
  mutate(
    feature = as.character(.data[[feature_col]]),
    feature = gsub("^s__", "", feature),
    feature = gsub("_", " ", feature),
    coef = as.numeric(.data[[coef_col]])
  ) %>%
  filter(!is.na(coef), coef != 0) %>%
  mutate(
    direction = ifelse(coef > 0, "Health-associated", "Disease-associated"),
    direction = factor(
      direction,
      levels = c("Health-associated", "Disease-associated")
    )
  ) %>%
  arrange(coef)

coef_plot_df$feature <- make.unique(coef_plot_df$feature)

coef_plot_df <- coef_plot_df %>%
  mutate(feature = factor(feature, levels = feature))

p_A <- ggplot(coef_plot_df, aes(x = coef, y = feature)) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey60",
    linewidth = 0.55
  ) +
  geom_segment(
    aes(x = 0, xend = coef, y = feature, yend = feature, color = direction),
    linewidth = 0.75
  ) +
  geom_point(
    aes(color = direction),
    size = 2.6
  ) +
  scale_color_manual(
    values = c(
      "Health-associated" = COL_HEALTH,
      "Disease-associated" = COL_DISEASE
    )
  ) +
  labs(
    x = "Coefficient",
    y = NULL,
    color = NULL
  ) +
  theme_fig(10) +
  theme(panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = c(0.72, 0.18),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(size = 12, face = "italic"),
    axis.text.x = element_text(size = 13),
    axis.title.x = element_text(size = 14)
  )

print(p_A)

save_panel(
  p_A,
  "Fig2A_species_coefficients_lollipop",
  width = 6.0,
  height = 10
)

write.csv(
  coef_plot_df,
  file.path(OUT_DIR, "Fig2A_species_coefficients_lollipop_data.csv"),
  row.names = FALSE
)

# =========================================================
# 6. Fig2B：ROC curve ：LOOCV OOF + 10-fold CV OOF + External validation
# =========================================================
cat("\n================ Fig2B: ROC curves ================\n")

auc_loocv_info <- safe_auc_ci(loocv_df$health01, loocv_df$prob_healthy)
auc_cv_info    <- safe_auc_ci(cv10_df$health01,  cv10_df$prob_healthy)
auc_test_info  <- safe_auc_ci(test_df$health01,  test_df$prob_healthy)

auc_loocv <- auc_loocv_info$auc
auc_cv    <- auc_cv_info$auc
auc_test  <- auc_test_info$auc

if ("repeat_id" %in% colnames(cv10_df)) {
  fold_auc_df <- cv10_df %>%
    group_by(repeat_id, fold) %>%
    summarise(
      fold_auc = safe_auc(health01, prob_healthy),
      n = n(),
      .groups = "drop"
    )
} else if ("fold" %in% colnames(cv10_df)) {
  fold_auc_df <- cv10_df %>%
    group_by(fold) %>%
    summarise(
      fold_auc = safe_auc(health01, prob_healthy),
      n = n(),
      .groups = "drop"
    )
} else {
  fold_auc_df <- data.frame(
    fold_auc = safe_auc(cv10_df$health01, cv10_df$prob_healthy),
    n = nrow(cv10_df)
  )
}

mean_cv_auc <- mean(fold_auc_df$fold_auc, na.rm = TRUE)

roc_loocv_obj <- safe_roc_with_ci_df(
  y = loocv_df$health01,
  score = loocv_df$prob_healthy,
  label = "LOOCV OOF",
  boot.n = BOOT_N
)

roc_cv_obj <- safe_roc_with_ci_df(
  y = cv10_df$health01,
  score = cv10_df$prob_healthy,
  label = "10-fold CV OOF",
  boot.n = BOOT_N
)

roc_test_obj <- safe_roc_with_ci_df(
  y = test_df$health01,
  score = test_df$prob_healthy,
  label = "External validation",
  boot.n = BOOT_N
)

roc_all <- bind_rows(
  roc_loocv_obj$roc,
  roc_cv_obj$roc,
  roc_test_obj$roc
) %>%
  mutate(
    set = factor(
      set,
      levels = c("LOOCV OOF", "10-fold CV OOF", "External validation")
    )
  )

ci_all <- bind_rows(
  roc_loocv_obj$ci,
  roc_cv_obj$ci,
  roc_test_obj$ci
) %>%
  mutate(
    set = factor(
      set,
      levels = c("LOOCV OOF", "10-fold CV OOF", "External validation")
    )
  )

legend_df <- data.frame(
  set = factor(
    c("LOOCV OOF", "10-fold CV OOF", "External validation"),
    levels = c("LOOCV OOF", "10-fold CV OOF", "External validation")
  ),
  set_label = c(
    paste0("LOOCV OOF (AUC = ", sprintf("%.2f", auc_loocv), ")"),
    paste0("10-fold CV OOF (AUC = ", sprintf("%.2f", auc_cv), ")"),
    paste0("External validation (AUC = ", sprintf("%.2f", auc_test), ")")
  )
)

roc_cols <- c(
  "LOOCV OOF" = COL_LOOCV,
  "10-fold CV OOF" = COL_CV,
  "External validation" = COL_TEST
)

p_B <- ggplot() +
  geom_abline(
    intercept = 0,
    slope = 1,
    color = "black",
    linewidth = 0.8
  ) +
  geom_ribbon(
    data = ci_all,
    aes(x = fpr, ymin = tpr_low, ymax = tpr_high, fill = set),
    alpha = 0.10,
    inherit.aes = FALSE
  ) +
  geom_line(
    data = roc_all,
    aes(x = fpr, y = tpr, color = set),
    linewidth = 1.2
  ) +
  scale_color_manual(
    values = roc_cols,
    labels = legend_df$set_label
  ) +
  scale_fill_manual(
    values = roc_cols,
    guide = "none"
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    expand = c(0.01, 0.01)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    expand = c(0.01, 0.01)
  ) +
  coord_equal() +
  labs(
    x = "False positive rate (1−Specificity)",
    y = "True positive rate (Sensitivity)"
  ) +
  theme_fig(11) +
  theme(
    legend.position = c(0.70, 0.23),
    legend.text = element_text(size = 8.6),
    legend.key.height = unit(0.45, "lines"),
    legend.key.width = unit(1.0, "lines")
  )

print(p_B)

save_panel(
  p_B,
  "Fig2B_species_ROC_3lines",
  width = 4.4,
  height = 4.2
)

auc_ci_tab <- bind_rows(
  data.frame(set = "LOOCV OOF", auc_loocv_info),
  data.frame(set = "10-fold CV OOF", auc_cv_info),
  data.frame(set = "External validation", auc_test_info)
)

write.csv(
  auc_ci_tab,
  file.path(OUT_DIR, "Fig2B_species_ROC_AUC_95CI.csv"),
  row.names = FALSE
)

write.csv(
  fold_auc_df,
  file.path(OUT_DIR, "Fig2B_species_10fold_foldwise_auc.csv"),
  row.names = FALSE
)

# =========================================================
# 7. Fig2C：OMHI Violin
# =========================================================
cat("\n================ Fig2C: OMHI violin plot ================\n")

violin_dat <- bind_rows(train_df, test_df) %>%
  mutate(
    set = factor(
      set,
      levels = c("Training", "External validation")
    )
  ) %>%
  group_by(set, Disease) %>%
  mutate(
    Disease_n = paste0(as.character(Disease), "\n(n=", n(), ")")
  ) %>%
  ungroup()

x_levels <- violin_dat %>%
  distinct(set, Disease, Disease_n) %>%
  arrange(set, Disease) %>%
  pull(Disease_n) %>%
  unique()

violin_dat <- violin_dat %>%
  mutate(Disease_n = factor(Disease_n, levels = x_levels))

stat_C <- violin_dat %>%
  group_by(set) %>%
  summarise(
    n_healthy = sum(Disease == "Healthy", na.rm = TRUE),
    n_disease = sum(Disease == "Disease", na.rm = TRUE),
    
    p_value = tryCatch(
      wilcox.test(
        OMHI[Disease == "Healthy"],
        OMHI[Disease == "Disease"],
        exact = FALSE
      )$p.value,
      error = function(e) NA_real_
    ),
    
    cliff_delta = calc_cliff_delta(
      OMHI[Disease == "Healthy"],
      OMHI[Disease == "Disease"]
    ),
    
    balanced_accuracy = calc_balanced_accuracy(health01, OMHI),
    .groups = "drop"
  ) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    label = paste0(
      format_p(p_adj),
      "\n",
      "d = ", sprintf("%.2f", cliff_delta),
      "\n",
      "Balanced Accuracy = ", sprintf("%.1f%%", 100 * balanced_accuracy)
    )
  )

write.csv(
  stat_C,
  file.path(OUT_DIR, "Fig2C_OMHI_violin_species_stats.csv"),
  row.names = FALSE
)

y_min_all <- min(violin_dat$OMHI, na.rm = TRUE)
y_max_all <- max(violin_dat$OMHI, na.rm = TRUE)
y_rng_all <- y_max_all - y_min_all
if(!is.finite(y_rng_all) || y_rng_all == 0) y_rng_all <- 1

label_pos <- stat_C %>%
  mutate(
    y_text = y_max_all + 0.22 * y_rng_all
  )

p_C <- ggplot(violin_dat, aes(x = Disease_n, y = OMHI, fill = Disease)) +
  geom_violin(
    aes(fill = Disease),
    width = 0.85,
    alpha = 0.6,
    color = NA,
    trim = FALSE
  ) +
  geom_violin(
    aes(color = Disease),
    width = 0.85,
    fill = NA,
    linewidth = 0.6,
    trim = FALSE
  ) +
  geom_jitter(
    aes(color = Disease),
    width = 0.10,
    size = 0.35,
    alpha = 0.45,
    stroke = 0
  ) +
  geom_boxplot(
    width = 0.18,
    outlier.shape = NA,
    color = "grey35",
    fill = NA,
    linewidth = 0.45
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "grey65",
    linewidth = 0.55
  ) +
  geom_blank(
    data = label_pos,
    aes(x = 1.5, y = y_text),
    inherit.aes = FALSE
  ) +
  geom_text(
    data = label_pos,
    aes(x = 1.5, y = y_text, label = label),
    inherit.aes = FALSE,
    size = 2.9,
    lineheight = 0.92
  ) +
  scale_fill_manual(
    values = c("Healthy" = COL_HEALTH, "Disease" = COL_DISEASE)
  ) +
  scale_color_manual(
    values = c("Healthy" = COL_HEALTH, "Disease" = COL_DISEASE)
  ) +
  facet_wrap(~ set, nrow = 1, scales = "free_x") +
  labs(
    x = NULL,
    y = "OMHI"
  ) +
  coord_cartesian(
  ylim = c(y_min_all, y_max_all + 0.42 * y_rng_all),
  clip = "off"
) +
  theme_fig(10) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 8.5),
    strip.text = element_text(size = 9.5, face = "bold")
  )

print(p_C)

save_panel(
  p_C,
  "Fig2C_species_OMHI_violin",
  width = 4.4,
  height = 3.6
)

# =========================================================
# 8. Output prompts 
# =========================================================
cat("\nComplete. Output directory:\n")
cat(OUT_DIR, "\n")

cat("\nKey ROC metrics:\n")
cat("LOOCV OOF AUC =", round(auc_loocv, 3), "\n")
cat("10-fold CV OOF AUC =", round(auc_cv, 3), "\n")
cat("Mean 10-fold fold-wise AUC =", round(mean_cv_auc, 3), "\n")
cat("External validation AUC =", round(auc_test, 3), "\n")

cat("\nKey OMHI distribution indicators:\n")
print(stat_C)