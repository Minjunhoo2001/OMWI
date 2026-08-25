
# =========================================================
# Focused analysis of OMWI and host/lifestyle variables
# Train + test version.
# Final output in this script:
#   1) Overall OMWI violin plots
#   2) OR forest plots for OMWI < 0
#   3) Healthy vs disease comparisons stratified by smoking, age, and BMI
#   4) A 3 x 3 combined figure arranged as:
#      smoking row, age row, and BMI row
#   5) Exact Trend P values shown in the upper-right corner
# =========================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(forcats)
  library(broom)
  library(purrr)
  library(stringr)
  library(readr)
  library(ggpubr)
  library(rstatix)
  library(patchwork)
})

options(stringsAsFactors = FALSE)

# =========================================================
# 0. Parameters
# =========================================================
ROOT_DIR <- "/share/home/HeMinjun/metagenomic/Oral_proxy_v2"

SPLIT_DIR <- file.path(ROOT_DIR, "01.split64")
OMWI_DIR  <- file.path(ROOT_DIR, "04.lasso64", "omwi_scores")

OUT_DIR   <- "/share/home/HeMinjun/metagenomic/GitHub/results/Supplementary_Figure_3"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT_DIR, "plots_focus"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT_DIR, "plots_focus", "categorical"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT_DIR, "plots_focus", "OR"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT_DIR, "plots_focus", "stratified"), showWarnings = FALSE, recursive = TRUE)

LEVEL_USE <- "species"
OMWI_COL  <- "OMWI"

MIN_N_CONT    <- 30
MIN_N_CAT     <- 20
MIN_PER_GRP   <- 5
MAX_CAT_LEVEL <- 8
P_CUTOFF      <- 0.05

# ---------------------------------------------------------
# Global color palette
# ---------------------------------------------------------
COL_1 <- "#4DBBD5"
COL_2 <- "#FEE5D9"
COL_3 <- "#FCAE91"
COL_4 <- "#FB6A4A"
COL_5 <- "#CB181D"

COL_HEALTH  <- COL_1
COL_DISEASE <- COL_5
COL_FEMALE  <- "#3C5488"

PRIOR_MEANINGFUL_VARS <- c(
  "gender",
  "RBC_group",
  "WBC_group",
  "IgG_group",
  "IgM_group"
)

# Reference groups for OR models
OR_REF_MAP <- list(
  smoker = "Never",
  H_pylori_status = "Negative",
  BMI_group = "18.5-24.9",
  age_group = "18-44",
  gender = "Male",
  WBC_group = "4-10",
  IgG_group = "7-16",
  IgM_group = "0.4-2.3"
)

# =========================================================
# 1. Utility functions
# =========================================================
pick_first_existing_col <- function(df, candidates) {
  x <- intersect(candidates, colnames(df))
  if (length(x) == 0) return(NA_character_)
  x[1]
}

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

clean_empty_to_na <- function(df) {
  for (j in colnames(df)) {
    x <- df[[j]]
    if (is.character(x) || is.factor(x)) {
      x <- as.character(x)
      x <- iconv(x, from = "", to = "UTF-8", sub = "")
      x <- trimws(x)
      x[x %in% c("", "NA", "N/A", "Unknown", "unknown", "missing", "Missing", "NULL", "null")] <- NA
      df[[j]] <- x
    }
  }
  df
}

make_valid_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9_\\-]+", "_", x)
  x <- gsub("_+", "_", x)
  x
}

is_probably_numeric <- function(x) {
  if (is.numeric(x)) return(TRUE)
  x2 <- suppressWarnings(as.numeric(as.character(x)))
  mean(!is.na(x2)) > 0.8
}

theme_omwi_clean <- function(base_size = 10.5) {
  theme_classic(base_size = base_size) +
    theme(
      legend.position = "none",
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      axis.line = element_line(color = "black", linewidth = 0.2),
      axis.ticks = element_line(color = "black", linewidth = 0.2),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.2),
      strip.text = element_text(face = "bold", color = "black"),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.2),
      legend.title = element_text(face = "bold", color = "black"),
      legend.text = element_text(color = "black"),
      plot.title = element_text(face = "bold", color = "black", hjust = 0.5),
      plot.subtitle = element_text(face = "plain", color = "black", hjust = 0.5)
    )
}

fmt_trend_label <- function(p) {
  if (is.na(p)) return("Trend P = NA")
  if (p == 0) return(paste0("Trend P = ", formatC(.Machine$double.xmin, format = "E", digits = 2)))
  paste0("Trend P = ", formatC(p, format = "E", digits = 2))
}

make_palette <- function(n) {
  base_cols <- c(COL_1, COL_2, COL_3, COL_4, COL_5)
  if (n <= length(base_cols)) {
    return(base_cols[seq_len(n)])
  } else {
    return(colorRampPalette(base_cols)(n))
  }
}

get_manual_palette_for_var <- function(sub, var_name = NULL) {
  levs <- levels(droplevels(factor(sub$group)))
  n <- length(levs)
  
  if (!is.null(var_name) && var_name %in% c("gender", "Gender")) {
    pal <- c("Male" = COL_1, "Female" = COL_FEMALE)
    pal <- pal[names(pal) %in% levs]
    return(pal[levs])
  }
  
  if (all(c("Healthy", "Disease") %in% levs) && length(levs) == 2) {
    pal <- c("Healthy" = COL_HEALTH, "Disease" = COL_DISEASE)
    pal <- pal[levs]
    return(pal)
  }
  
  if (!is.null(var_name) && var_name %in% c("BMI_group", "BMI group", "Abnormal BMI")) {
    pal <- c(
      "18.5-24.9" = COL_1,
      "<18.5"     = COL_2,
      "25-29.9"   = COL_4,
      "≥30"       = COL_5
    )
    pal <- pal[levs]
    return(pal)
  }
  
  if (!is.null(var_name) && var_name %in% c("smoker", "Smoking status")) {
    pal <- c(
      "Never"   = COL_1,
      "Past"    = COL_4,
      "Current" = COL_5
    )
    pal <- pal[levs]
    return(pal)
  }
  
  if (!is.null(var_name) && var_name %in% c("age_group", "age group", "Age group", "Age category")) {
    pal <- c(
      "18-44" = COL_1,
      "45-59" = COL_3,
      "60-74" = COL_4,
      "75-90" = COL_5
    )
    pal <- pal[levs]
    return(pal)
  }
  
  if (!is.null(var_name) && var_name %in% c("H_pylori_status", "H.pylori", "H.pylori status")) {
    pal <- c(
      "Negative" = COL_1,
      "Positive" = COL_5
    )
    pal <- pal[levs]
    return(pal)
  }
  
  cols <- make_palette(n)
  names(cols) <- levs
  cols
}

get_palette_for_levels <- function(levels_use, var_name = NULL) {
  levels_use <- as.character(levels_use)
  levels_use <- levels_use[!is.na(levels_use)]
  sub <- tibble(group = factor(levels_use, levels = levels_use))
  
  pal <- get_manual_palette_for_var(sub, var_name = var_name)
  pal <- pal[levels_use]
  
  if (any(is.na(pal))) {
    miss <- names(pal)[is.na(pal)]
    fill_cols <- make_palette(length(miss))
    names(fill_cols) <- miss
    pal[miss] <- fill_cols
  }
  
  pal
}

format_p <- function(p) {
  ifelse(is.na(p), NA_character_,
         ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

format_es <- function(x, digits = 3) {
  ifelse(is.na(x), NA_character_, sprintf(paste0("%.", digits, "f"), x))
}

interpret_rbc_magnitude <- function(x) {
  case_when(
    is.na(x) ~ NA_character_,
    abs(x) < 0.10 ~ "negligible",
    abs(x) < 0.30 ~ "small",
    abs(x) < 0.50 ~ "moderate",
    TRUE ~ "large"
  )
}

interpret_eps2_magnitude <- function(x) {
  case_when(
    is.na(x) ~ NA_character_,
    x < 0.01 ~ "negligible",
    x < 0.08 ~ "small",
    x < 0.26 ~ "moderate",
    TRUE ~ "large"
  )
}

make_n_axis_labels <- function(sub, group_col = "group") {
  tb <- sub %>%
    dplyr::filter(!is.na(.data[[group_col]])) %>%
    dplyr::distinct(.data[[group_col]]) %>%
    dplyr::mutate(
      group_chr = as.character(.data[[group_col]]),
      label = group_chr
    )
  
  labs <- tb$label
  names(labs) <- tb$group_chr
  labs
}

# =========================================================
# 1.1 Statistical helper functions
# =========================================================
get_overall_p <- function(sub) {
  sub <- sub %>% filter(!is.na(group), !is.na(OMWI))
  if (nlevels(droplevels(sub$group)) < 2) return(NA_real_)
  
  if (nlevels(droplevels(sub$group)) == 2) {
    out <- tryCatch(wilcox.test(OMWI ~ group, data = sub)$p.value, error = function(e) NA_real_)
  } else {
    out <- tryCatch(kruskal.test(OMWI ~ group, data = sub)$p.value, error = function(e) NA_real_)
  }
  out
}

get_overall_effect <- function(sub) {
  sub <- sub %>% filter(!is.na(group), !is.na(OMWI))
  sub$group <- droplevels(factor(sub$group))
  k <- nlevels(sub$group)
  n <- nrow(sub)
  
  if (k < 2) return(NULL)
  
  if (k == 2) {
    eff <- tryCatch({
      rstatix::wilcox_effsize(sub, OMWI ~ group)
    }, error = function(e) NULL)
    
    if (is.null(eff) || nrow(eff) == 0) {
      return(tibble(
        effect_type = "rank_biserial",
        effect_size = NA_real_,
        effect_magnitude = NA_character_,
        n = n
      ))
    }
    
    return(tibble(
      effect_type = "rank_biserial",
      effect_size = eff$effsize[1],
      effect_magnitude = if ("magnitude" %in% colnames(eff)) eff$magnitude[1] else interpret_rbc_magnitude(eff$effsize[1]),
      n = n
    ))
  }
  
  kr <- tryCatch(kruskal.test(OMWI ~ group, data = sub), error = function(e) NULL)
  if (is.null(kr)) {
    return(tibble(
      effect_type = "epsilon_squared",
      effect_size = NA_real_,
      effect_magnitude = NA_character_,
      n = n
    ))
  }
  
  H <- as.numeric(kr$statistic)
  eps2 <- (H - k + 1) / (n - k)
  eps2 <- max(0, eps2)
  
  tibble(
    effect_type = "epsilon_squared",
    effect_size = eps2,
    effect_magnitude = interpret_eps2_magnitude(eps2),
    n = n
  )
}

get_pairwise_p <- function(sub, group_col = "group", value_col = "OMWI") {
  sub <- sub %>% filter(!is.na(.data[[group_col]]), !is.na(.data[[value_col]]))
  if (n_distinct(sub[[group_col]]) < 2) return(NULL)
  
  pw <- tryCatch({
    sub %>%
      rstatix::wilcox_test(
        as.formula(paste0(value_col, " ~ ", group_col)),
        p.adjust.method = "BH"
      )
  }, error = function(e) NULL)
  
  pw
}

get_pairwise_effect <- function(sub, group_col = "group", value_col = "OMWI") {
  sub <- sub %>%
    filter(!is.na(.data[[group_col]]), !is.na(.data[[value_col]])) %>%
    mutate(group_tmp = factor(.data[[group_col]]))
  
  lv <- levels(droplevels(sub$group_tmp))
  if (length(lv) < 2) return(NULL)
  
  combs <- combn(lv, 2, simplify = FALSE)
  
  res <- lapply(combs, function(cc) {
    tmp <- sub %>%
      filter(group_tmp %in% cc) %>%
      mutate(group_tmp = droplevels(group_tmp))
    
    med1 <- median(tmp[[value_col]][tmp$group_tmp == cc[1]], na.rm = TRUE)
    med2 <- median(tmp[[value_col]][tmp$group_tmp == cc[2]], na.rm = TRUE)
    mean1 <- mean(tmp[[value_col]][tmp$group_tmp == cc[1]], na.rm = TRUE)
    mean2 <- mean(tmp[[value_col]][tmp$group_tmp == cc[2]], na.rm = TRUE)
    
    eff <- tryCatch({
      rstatix::wilcox_effsize(tmp, as.formula(paste0(value_col, " ~ group_tmp")))
    }, error = function(e) NULL)
    
    tibble(
      group1 = cc[1],
      group2 = cc[2],
      n1 = sum(tmp$group_tmp == cc[1]),
      n2 = sum(tmp$group_tmp == cc[2]),
      median1 = med1,
      median2 = med2,
      mean1 = mean1,
      mean2 = mean2,
      delta_median = med2 - med1,
      delta_mean = mean2 - mean1,
      effect_type = "rank_biserial",
      effect_size = if (!is.null(eff) && nrow(eff) > 0) eff$effsize[1] else NA_real_,
      effect_magnitude = if (!is.null(eff) && nrow(eff) > 0 && "magnitude" %in% colnames(eff)) {
        eff$magnitude[1]
      } else if (!is.null(eff) && nrow(eff) > 0) {
        interpret_rbc_magnitude(eff$effsize[1])
      } else {
        NA_character_
      }
    )
  })
  
  bind_rows(res)
}

# =========================================================
# 1.2 OR functions: outcome is OMWI < 0
# =========================================================
fit_or_by_variable_binary <- function(dat, var, ref_level = NULL, outcome_type = c("lt0")) {
  outcome_type <- match.arg(outcome_type)
  if (!all(c(var, "OMWI") %in% colnames(dat))) return(NULL)
  
  sub <- dat %>%
    dplyr::select(sample_id, OMWI, all_of(var)) %>%
    dplyr::filter(!is.na(OMWI), !is.na(.data[[var]])) %>%
    dplyr::mutate(
      outcome = ifelse(OMWI < 0, 1, 0),
      group = as.character(.data[[var]])
    )
  
  grp_tab <- table(sub$group)
  keep_grp <- names(grp_tab[grp_tab >= MIN_PER_GRP])
  sub <- sub %>% dplyr::filter(group %in% keep_grp)
  
  if (nrow(sub) < MIN_N_CAT || dplyr::n_distinct(sub$group) < 2) return(NULL)
  
  if (!is.null(ref_level) && ref_level %in% unique(sub$group)) {
    levs <- c(ref_level, setdiff(unique(sub$group), ref_level))
  } else {
    levs <- unique(sub$group)
  }
  
  sub$group <- factor(sub$group, levels = levs)
  
  fit <- tryCatch(
    glm(outcome ~ group, data = sub, family = binomial()),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)
  
  tidy_res <- broom::tidy(fit, conf.int = TRUE, exponentiate = TRUE)
  
  res_ref <- tibble(
    variable = var,
    outcome_type = outcome_type,
    term = paste0("group", levs[1]),
    group = levs[1],
    OR = 1,
    conf.low = 1,
    conf.high = 1,
    p.value = NA_real_,
    is_ref = TRUE
  )
  
  res_nonref <- tidy_res %>%
    dplyr::filter(term != "(Intercept)") %>%
    dplyr::mutate(
      variable = var,
      outcome_type = outcome_type,
      group = stringr::str_remove(term, "^group"),
      OR = estimate,
      is_ref = FALSE
    ) %>%
    dplyr::select(variable, outcome_type, term, group, OR, conf.low, conf.high, p.value, is_ref)
  
  res <- bind_rows(res_ref, res_nonref) %>%
    left_join(
      sub %>% count(group, name = "n"),
      by = "group"
    ) %>%
    mutate(
      label = paste0(group, " (n=", n, ")"),
      sig = dplyr::case_when(
        is.na(p.value) ~ "",
        p.value < 0.001 ~ "***",
        p.value < 0.01  ~ "**",
        p.value < 0.05  ~ "*",
        TRUE ~ ""
      )
    )
  
  res
}

plot_or_forest_binary <- function(or_df, var_label = NULL, out_file = NULL, outcome_type = c("lt0")) {
  outcome_type <- match.arg(outcome_type)
  if (is.null(or_df) || nrow(or_df) == 0) return(NULL)
  if (is.null(var_label)) var_label <- unique(or_df$variable)[1]
  
  x_title <- "Odds ratio for OMWI < 0 (95% CI)"
  
  group_levels <- unique(or_df$group)
  group_levels <- group_levels[!is.na(group_levels)]
  pal <- get_palette_for_levels(group_levels, var_name = unique(or_df$variable)[1])
  
  plot_df <- or_df %>%
    mutate(
      group = factor(group, levels = group_levels),
      label = paste0(group, " (n=", n, ")"),
      label = factor(label, levels = rev(label))
    )
  
  x_star_right <- max(c(plot_df$conf.high, plot_df$OR, 16), na.rm = TRUE) * 1.12
  if (!is.finite(x_star_right)) x_star_right <- 16

  star_df <- plot_df %>%
    filter(!is_ref, sig != "") %>%
    mutate(x_star = x_star_right)
  
  p <- ggplot(plot_df, aes(x = OR, y = label)) +
    geom_vline(xintercept = 1, linetype = 2, color = "grey60", linewidth = 0.5) +
    geom_errorbarh(
      data = plot_df %>% filter(!is_ref),
      aes(xmin = conf.low, xmax = conf.high, color = group),
      height = 0.18,
      linewidth = 0.65,
      show.legend = FALSE
    ) +
    geom_point(
      data = plot_df %>% filter(!is_ref),
      aes(color = group),
      size = 2.8,
      show.legend = FALSE
    ) +
    geom_point(
      data = plot_df %>% filter(is_ref),
      fill = "white",
      color = "black",
      size = 2.8,
      shape = 21,
      stroke = 0.8,
      show.legend = FALSE
    )+
    geom_text(
      data = star_df,
      aes(x = x_star, y = label, label = sig),
      color = "black",
      size = 3.4,
      hjust = 0,
      na.rm = TRUE,
      inherit.aes = FALSE
    ) +
    scale_color_manual(values = pal, drop = FALSE) +
    scale_fill_manual(values = pal, drop = FALSE) +
    scale_x_log10(
      breaks = c(0.25, 0.5, 1, 2, 4, 8, 16),
      labels = c("0.25", "0.5", "1", "2", "4", "8", "16")
    ) +
    labs(
      title = NULL,
      subtitle = NULL,
      x = x_title,
      y = NULL
    ) +
    coord_cartesian(xlim = c(NA, x_star_right * 1.05), clip = "off") +
    theme_omwi_clean() +
    theme(
      axis.text.y = element_text(color = "black"),
      axis.text.x = element_text(color = "black"),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.margin = margin(5.5, 18, 5.5, 5.5)
    )
  
  if (!is.null(out_file)) {
    ggsave(out_file, p, width = 4, height = 4, units = "in", dpi = 600, bg = "white")
  }
  
  p
}

# =========================================================
# 2. Read metadata; keep train and test samples only
# =========================================================
meta_train <- fread(file.path(SPLIT_DIR, "meta_train.csv"), data.table = FALSE, check.names = FALSE)
meta_test  <- fread(file.path(SPLIT_DIR, "meta_test.csv"),  data.table = FALSE, check.names = FALSE)

meta_train$dataset_split <- "train"
meta_test$dataset_split  <- "test"

meta_all <- bind_rows(meta_train, meta_test)
meta_all <- clean_empty_to_na(meta_all)

if (!"sample_id" %in% colnames(meta_all)) {
  stop("Cannot find sample_id column in metadata.")
}

# =========================================================
# 3. Read OMWI scores; keep train and test samples only
# =========================================================
omwi_train_file <- file.path(OMWI_DIR, paste0("OMWI_train_", LEVEL_USE, "_alpha050.csv"))
omwi_test_file  <- file.path(OMWI_DIR, paste0("OMWI_test_", LEVEL_USE, "_alpha050.csv"))

if (!file.exists(omwi_train_file)) stop("Missing file: ", omwi_train_file)
if (!file.exists(omwi_test_file))  stop("Missing file: ", omwi_test_file)

omwi_train <- fread(omwi_train_file, data.table = FALSE, check.names = FALSE)
omwi_test  <- fread(omwi_test_file,  data.table = FALSE, check.names = FALSE)

omwi_train$dataset_split <- "train"
omwi_test$dataset_split  <- "test"

omwi_all <- bind_rows(omwi_train, omwi_test)

if (!"sample_id_use" %in% colnames(omwi_all)) {
  stop("Cannot find sample_id_use column in OMWI files.")
}
if (!(OMWI_COL %in% colnames(omwi_all))) {
  stop("Cannot find OMWI column: ", OMWI_COL)
}

omwi_all <- omwi_all %>%
  mutate(sample_id = as.character(sample_id_use)) %>%
  select(sample_id, dataset_split, all_of(OMWI_COL), everything())

# =========================================================
# 4. Merge metadata with OMWI scores
# =========================================================
dat <- meta_all %>%
  mutate(sample_id = as.character(sample_id)) %>%
  left_join(
    omwi_all %>% select(sample_id, OMWI = all_of(OMWI_COL), prob_healthy, level, prev_cutoff),
    by = "sample_id"
  )

if (all(is.na(dat$OMWI))) {
  stop("All OMWI values are NA after merging. Please check sample_id consistency.")
}

dat$OMWI <- safe_numeric(dat$OMWI)
dat$OMWI_low <- ifelse(!is.na(dat$OMWI) & dat$OMWI < 0, 1, ifelse(!is.na(dat$OMWI), 0, NA))

# =========================================================
# 5. Preprocess key variables
# =========================================================
if ("age" %in% colnames(dat)) {
  dat$age_raw <- trimws(as.character(dat$age))
}

if ("BMI" %in% colnames(dat)) {
  dat$BMI_raw <- trimws(as.character(dat$BMI))
  dat$BMI <- safe_numeric(dat$BMI_raw)
}

if ("RBC" %in% colnames(dat)) {
  dat$RBC_raw <- trimws(as.character(dat$RBC))
  dat$RBC <- safe_numeric(dat$RBC_raw)
}

if ("WBC" %in% colnames(dat)) {
  dat$WBC_raw <- trimws(as.character(dat$WBC))
  dat$WBC <- safe_numeric(dat$WBC_raw)
}

if ("IgG" %in% colnames(dat)) {
  dat$IgG_raw <- trimws(as.character(dat$IgG))
  dat$IgG <- safe_numeric(dat$IgG_raw)
}

if ("IgM" %in% colnames(dat)) {
  dat$IgM_raw <- trimws(as.character(dat$IgM))
  dat$IgM <- safe_numeric(dat$IgM_raw)
}

# -----------------------------
# 5.1 Harmonize the healthy/disease grouping variable
# -----------------------------
if ("disease" %in% colnames(dat)) {
  dat$disease_raw <- trimws(as.character(dat$disease))
  
  dat$health01 <- dplyr::case_when(
    tolower(dat$disease_raw) %in% c("healthy", "health", "control", "ctrl", "0") ~ "Healthy",
    tolower(dat$disease_raw) %in% c("disease", "case", "1") ~ "Disease",
    TRUE ~ NA_character_
  )
  
  dat$health01 <- factor(dat$health01, levels = c("Healthy", "Disease"))
} else {
  stop("Cannot find the disease column in metadata; healthy/disease stratified analysis cannot be performed.")
}

# -----------------------------
# 5.2 Harmonize smoking status
# -----------------------------
if ("smoker" %in% colnames(dat)) {
  dat$smoker <- trimws(as.character(dat$smoker))
  
  dat$smoker <- dplyr::case_when(
    tolower(dat$smoker) %in% c("never", "nerver", "no", "non-smoker", "nonsmoker") ~ "Never",
    tolower(dat$smoker) %in% c("past", "former", "ex-smoker", "previous") ~ "Past",
    tolower(dat$smoker) %in% c("current", "yes", "smoker") ~ "Current",
    TRUE ~ NA_character_
  )
  
  dat$smoker <- factor(dat$smoker, levels = c("Never", "Past", "Current"))
}

# -----------------------------
# 5.3 Harmonize H. pylori status
# -----------------------------
hp_col <- pick_first_existing_col(
  dat,
  c("H.pylori", "H_pylori", "hpylori", "Hpylori", "h_pylori",
    "Helicobacter_pylori", "helicobacter_pylori", "H.pylori_status", "Hp_status")
)

if (!is.na(hp_col)) {
  dat$H_pylori_raw <- trimws(as.character(dat[[hp_col]]))
  
  dat$H_pylori_status <- dplyr::case_when(
    tolower(dat$H_pylori_raw) %in% c("negative", "neg", "-", "0", "no", "n") ~ "Negative",
    tolower(dat$H_pylori_raw) %in% c("positive", "pos", "+", "1", "yes", "y") ~ "Positive",
    TRUE ~ NA_character_
  )
  
  dat$H_pylori_status <- factor(dat$H_pylori_status, levels = c("Negative", "Positive"))
}

# -----------------------------
# 5.4 Classify BMI using WHO-style categories
# -----------------------------
if ("BMI" %in% colnames(dat)) {
  dat$BMI_group <- dplyr::case_when(
    !is.na(dat$BMI) & dat$BMI >= 18.5 & dat$BMI < 25 ~ "18.5-24.9",
    !is.na(dat$BMI) & dat$BMI < 18.5 ~ "<18.5",
    !is.na(dat$BMI) & dat$BMI >= 25 & dat$BMI < 30 ~ "25-29.9",
    !is.na(dat$BMI) & dat$BMI >= 30 ~ "≥30",
    TRUE ~ NA_character_
  )
  
  dat$BMI_group <- factor(dat$BMI_group, levels = c("18.5-24.9", "<18.5", "25-29.9", "≥30"))
}

# -----------------------------
# 5.5 Age groups
# -----------------------------
if ("age" %in% colnames(dat)) {
  dat$age_num_strict <- ifelse(
    grepl("^[0-9]+(\\.[0-9]+)?$", dat$age_raw),
    dat$age_raw,
    NA_character_
  )
  
  dat$age_num_strict <- as.numeric(dat$age_num_strict)
  
  dat$age_group <- dplyr::case_when(
    !is.na(dat$age_num_strict) & dat$age_num_strict >= 18 & dat$age_num_strict <= 44 ~ "18-44",
    !is.na(dat$age_num_strict) & dat$age_num_strict >= 45 & dat$age_num_strict <= 59 ~ "45-59",
    !is.na(dat$age_num_strict) & dat$age_num_strict >= 60 & dat$age_num_strict <= 74 ~ "60-74",
    !is.na(dat$age_num_strict) & dat$age_num_strict >= 75 & dat$age_num_strict <= 90 ~ "75-90",
    TRUE ~ NA_character_
  )
  
  dat$age_group <- factor(dat$age_group, levels = c("18-44", "45-59", "60-74", "75-90"))
}

# -----------------------------
# 5.6 Harmonize gender
# -----------------------------
if ("gender" %in% colnames(dat)) {
  dat$gender <- as.character(dat$gender)
  dat$gender <- trimws(dat$gender)
  
  dat$gender <- dplyr::case_when(
    tolower(dat$gender) %in% c("male", "m", "man", "boy", "1", "男") ~ "Male",
    tolower(dat$gender) %in% c("female", "f", "woman", "girl", "2", "女") ~ "Female",
    TRUE ~ NA_character_
  )
  
  dat$gender <- factor(dat$gender, levels = c("Male", "Female"))
}

# -----------------------------
# 5.7 RBC groups based on sex-specific reference ranges
# -----------------------------
if (all(c("RBC", "gender") %in% colnames(dat))) {
  dat$RBC_group <- dplyr::case_when(
    dat$gender == "Male"   & !is.na(dat$RBC) & dat$RBC < 4.0 ~ "Male:<4.0",
    dat$gender == "Male"   & !is.na(dat$RBC) & dat$RBC >= 4.0 & dat$RBC <= 5.5 ~ "Male:4.0-5.5",
    dat$gender == "Male"   & !is.na(dat$RBC) & dat$RBC > 5.5 ~ "Male:>5.5",
    
    dat$gender == "Female" & !is.na(dat$RBC) & dat$RBC < 3.5 ~ "Female:<3.5",
    dat$gender == "Female" & !is.na(dat$RBC) & dat$RBC >= 3.5 & dat$RBC <= 5.0 ~ "Female:3.5-5.0",
    dat$gender == "Female" & !is.na(dat$RBC) & dat$RBC > 5.0 ~ "Female:>5.0",
    
    TRUE ~ NA_character_
  )
  
  dat$RBC_group <- factor(
    dat$RBC_group,
    levels = c(
      "Male:<4.0", "Male:4.0-5.5", "Male:>5.5",
      "Female:<3.5", "Female:3.5-5.0", "Female:>5.0"
    )
  )
}

# -----------------------------
# 5.8 WBC groups
# -----------------------------
if ("WBC" %in% colnames(dat)) {
  dat$WBC_group <- dplyr::case_when(
    !is.na(dat$WBC) & dat$WBC < 4 ~ "<4",
    !is.na(dat$WBC) & dat$WBC >= 4 & dat$WBC <= 10 ~ "4-10",
    !is.na(dat$WBC) & dat$WBC > 10 ~ ">10",
    TRUE ~ NA_character_
  )
  
  dat$WBC_group <- factor(dat$WBC_group, levels = c("<4", "4-10", ">10"))
}

# -----------------------------
# 5.9 IgG groups
# -----------------------------
if ("IgG" %in% colnames(dat)) {
  dat$IgG_group <- dplyr::case_when(
    !is.na(dat$IgG) & dat$IgG < 7 ~ "<7",
    !is.na(dat$IgG) & dat$IgG >= 7 & dat$IgG <= 16 ~ "7-16",
    !is.na(dat$IgG) & dat$IgG > 16 ~ ">16",
    TRUE ~ NA_character_
  )
  
  dat$IgG_group <- factor(dat$IgG_group, levels = c("<7", "7-16", ">16"))
}

# -----------------------------
# 5.10 IgM groups
# -----------------------------
if ("IgM" %in% colnames(dat)) {
  dat$IgM_group <- dplyr::case_when(
    !is.na(dat$IgM) & dat$IgM < 0.4 ~ "<0.4",
    !is.na(dat$IgM) & dat$IgM >= 0.4 & dat$IgM <= 2.3 ~ "0.4-2.3",
    !is.na(dat$IgM) & dat$IgM > 2.3 ~ ">2.3",
    TRUE ~ NA_character_
  )
  
  dat$IgM_group <- factor(dat$IgM_group, levels = c("<0.4", "0.4-2.3", ">2.3"))
}

# -----------------------------
# 5.11 Strata for Healthy vs Disease robustness analysis
# -----------------------------
if ("smoker" %in% colnames(dat)) {
  dat$smoking_strata_hd <- dplyr::case_when(
    as.character(dat$smoker) == "Never" ~ "Never",
    as.character(dat$smoker) == "Past" ~ "Past",
    as.character(dat$smoker) == "Current" ~ "Current",
    TRUE ~ NA_character_
  )
  dat$smoking_strata_hd <- factor(dat$smoking_strata_hd, levels = c("Never", "Past", "Current"))
}

if ("age_num_strict" %in% colnames(dat)) {
  dat$age_strata_hd <- dplyr::case_when(
    !is.na(dat$age_num_strict) & dat$age_num_strict >= 18 & dat$age_num_strict <= 44 ~ "18-44",
    !is.na(dat$age_num_strict) & dat$age_num_strict >= 45 & dat$age_num_strict <= 59 ~ "45-59",
    !is.na(dat$age_num_strict) & dat$age_num_strict >= 60 & dat$age_num_strict <= 74 ~ "60-74",
    !is.na(dat$age_num_strict) & dat$age_num_strict >= 75 & dat$age_num_strict <= 90 ~ "75-90",
    TRUE ~ NA_character_
  )
  dat$age_strata_hd <- factor(dat$age_strata_hd, levels = c("18-44", "45-59", "60-74", "75-90"))
}

if ("BMI" %in% colnames(dat)) {
  dat$BMI_strata_hd <- dplyr::case_when(
    !is.na(dat$BMI) & dat$BMI >= 18.5 & dat$BMI < 25 ~ "Normal",
    !is.na(dat$BMI) & dat$BMI >= 25 & dat$BMI < 30 ~ "Overweight",
    !is.na(dat$BMI) & dat$BMI >= 30 ~ "Obese",
    TRUE ~ NA_character_
  )
  dat$BMI_strata_hd <- factor(dat$BMI_strata_hd, levels = c("Normal", "Overweight", "Obese"))
}

# =========================================================
# 6. Define variables for focused analysis
# =========================================================
focus_categorical_vars <- c(
  "smoker",
  "H_pylori_status",
  "BMI_group",
  "age_group",
  "gender",
  "RBC_group",
  "WBC_group",
  "IgG_group",
  "IgM_group"
)
focus_categorical_vars <- intersect(focus_categorical_vars, colnames(dat))

prior_vars_exist <- intersect(PRIOR_MEANINGFUL_VARS, colnames(dat))
prior_cont_vars <- prior_vars_exist[sapply(prior_vars_exist, function(v) is_probably_numeric(dat[[v]]))]
prior_cat_vars  <- setdiff(prior_vars_exist, prior_cont_vars)

focus_categorical_vars <- unique(c(focus_categorical_vars, prior_cat_vars))

# =========================================================
# 7. Overall violin plot function
# =========================================================
plot_group_violin <- function(sub, var, xlab = var, out_file = NULL) {
  sub <- sub %>% filter(!is.na(group), !is.na(OMWI))
  sub$group <- droplevels(factor(sub$group))
  
  n_grp <- nlevels(sub$group)
  if (n_grp < 2) return(NULL)
  
  pal <- get_manual_palette_for_var(sub, var)
  pw <- get_pairwise_p(sub, group_col = "group", value_col = "OMWI")
  
  trend_p <- tryCatch(
    suppressWarnings(cor.test(as.numeric(sub$group), sub$OMWI, method = "spearman")$p.value),
    error = function(e) NA_real_
  )
  trend_lab <- fmt_trend_label(trend_p)
  
  axis_labs <- make_n_axis_labels(sub, "group")
  
  median_df <- sub %>%
    group_by(group) %>%
    summarise(median_OMWI = median(OMWI, na.rm = TRUE), .groups = "drop") %>%
    mutate(group = factor(group, levels = levels(sub$group)))
  
  y_max <- max(sub$OMWI, na.rm = TRUE)
  y_min <- min(sub$OMWI, na.rm = TRUE)
  y_rng <- y_max - y_min
  if (is.na(y_rng) || y_rng == 0) y_rng <- 0.2
  
  p <- ggplot(sub, aes(x = group, y = OMWI, fill = group, color = group)) +
    geom_violin(width = 0.88, trim = FALSE, color = NA, alpha = 0.7) +
    geom_jitter(width = 0.06, size = 0.9, alpha = 0.5, show.legend = FALSE) +
    geom_boxplot(
      width = 0.20,
      outlier.shape = NA,
      fill = "transparent",
      color = "black",
      linewidth = 0.45
    ) +
    geom_line(
      data = median_df,
      aes(x = group, y = median_OMWI, group = 1),
      inherit.aes = FALSE,
      color = "black",
      linewidth = 0.75
    ) +
    geom_point(
      data = median_df,
      aes(x = group, y = median_OMWI),
      inherit.aes = FALSE,
      shape = 21,
      fill = "white",
      color = "black",
      stroke = 0.8,
      size = 2.5
    ) +
    scale_fill_manual(values = pal) +
    scale_color_manual(values = pal) +
    scale_x_discrete(labels = axis_labs) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.28))) +
    labs(
      title = NULL,
      subtitle = NULL,
      x = xlab,
      y = "OMWI"
    ) +
    annotate(
      "text",
      x = Inf, y = Inf,
      label = trend_lab,
      hjust = 1.03, vjust = 1.2,
      size = 3.2
    ) +
    theme_omwi_clean() +
    theme(
      axis.text.x = element_text(
        color = "black",
        angle = 25,
        hjust = 1,
        vjust = 1
      ),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.margin = margin(5.5, 6, 5.5, 5.5)
    )
  
  if (!is.null(pw) && nrow(pw) > 0) {
    p_col <- c("p.adj", "p", "p.value", "pval", "p.adj.signif", "p.signif")
    p_col <- p_col[p_col %in% names(pw)][1]
    
    if (!is.na(p_col)) {
      if (p_col %in% c("p.adj.signif", "p.signif")) {
        pw$p_label <- as.character(pw[[p_col]])
      } else {
        pw$p_label <- dplyr::case_when(
          is.na(pw[[p_col]])  ~ NA_character_,
          pw[[p_col]] < 0.001 ~ "***",
          pw[[p_col]] < 0.01  ~ "**",
          pw[[p_col]] < 0.05  ~ "*",
          TRUE ~ "ns"
        )
      }
    } else {
      pw$p_label <- NA_character_
    }
    
    pw_plot <- pw %>%
      dplyr::filter(!is.na(p_label), p_label != "ns") %>%
      dplyr::mutate(
        group1 = as.character(group1),
        group2 = as.character(group2)
      )
    
    if (nrow(pw_plot) > 0) {
      levs <- levels(sub$group)
      
      pw_plot <- pw_plot %>%
        dplyr::mutate(
          x1 = match(group1, levs),
          x2 = match(group2, levs)
        ) %>%
        dplyr::arrange(x1, x2) %>%
        dplyr::mutate(
          y.bracket = y_max + (0.06 + 0.09 * (row_number() - 1)) * y_rng,
          y.label   = y.bracket + 0.03 * y_rng,
          xmid      = (x1 + x2) / 2
        )
      
      p <- p +
        geom_segment(
          data = pw_plot,
          aes(x = x1, xend = x2, y = y.bracket, yend = y.bracket),
          inherit.aes = FALSE,
          linewidth = 0.4,
          color = "black"
        ) +
        geom_segment(
          data = pw_plot,
          aes(x = x1, xend = x1, y = y.bracket, yend = y.bracket - 0.02 * y_rng),
          inherit.aes = FALSE,
          linewidth = 0.4,
          color = "black"
        ) +
        geom_segment(
          data = pw_plot,
          aes(x = x2, xend = x2, y = y.bracket, yend = y.bracket - 0.02 * y_rng),
          inherit.aes = FALSE,
          linewidth = 0.4,
          color = "black"
        ) +
        geom_text(
          data = pw_plot,
          aes(x = xmid, y = y.label, label = p_label),
          inherit.aes = FALSE,
          size = 3.4
        )
    }
  }
  
  if (!is.null(out_file)) {
    ggsave(out_file, p, width = 4, height = 4, dpi = 600, bg = "white")
  }
  
  return(list(
    plot = p,
    trend_p = trend_p,
    pairwise = pw
  ))
}

# =========================================================
# 8. Two-panel layout function: overall violin + OR
# =========================================================
plot_violin_or_pair <- function(violin_plot, or_plot, out_file = NULL, widths = c(1.15, 1.0)) {
  if (is.null(violin_plot) || is.null(or_plot)) return(NULL)
  
  p_left <- violin_plot +
    labs(title = NULL, subtitle = NULL) +
    theme(plot.title = element_blank(), plot.subtitle = element_blank())
  
  p_right <- or_plot +
    labs(title = NULL, subtitle = NULL) +
    theme(plot.title = element_blank(), plot.subtitle = element_blank())
  
  p_out <- p_left + p_right +
    patchwork::plot_layout(nrow = 1, widths = widths)
  
  if (!is.null(out_file)) {
    ggsave(
      filename = out_file,
      plot = p_out,
      width = 9,
      height = 3,
      units = "in",
      dpi = 600,
      bg = "white"
    )
  }
  
  return(p_out)
}

add_panel_title <- function(p, title) {
  if (is.null(p)) return(NULL)
  p +
    labs(title = title) +
    theme(
      plot.title = element_text(
        face = "bold",
        color = "black",
        hjust = 0,
        size = 10.5,
        margin = margin(b = 3)
      )
    )
}

# =========================================================
# 8.1 Stratified healthy vs disease OMWI comparison
# =========================================================
plot_stratified_health_comparison <- function(dat, strat_var, strat_label, out_prefix = NULL) {
  if (!all(c("sample_id", "OMWI", "health01", strat_var) %in% colnames(dat))) return(NULL)
  
  sub <- dat %>%
    dplyr::select(sample_id, OMWI, health01, stratum = all_of(strat_var)) %>%
    dplyr::filter(!is.na(OMWI), !is.na(health01), !is.na(stratum)) %>%
    dplyr::mutate(
      health01 = factor(as.character(health01), levels = c("Healthy", "Disease")),
      stratum = factor(as.character(stratum), levels = unique(as.character(stratum)))
    )
  
  keep_strata <- sub %>%
    dplyr::count(stratum, health01, name = "n") %>%
    tidyr::complete(stratum, health01 = factor(c("Healthy", "Disease"), levels = c("Healthy", "Disease")), fill = list(n = 0)) %>%
    dplyr::group_by(stratum) %>%
    dplyr::summarise(
      has_both = all(n >= MIN_PER_GRP),
      .groups = "drop"
    ) %>%
    dplyr::filter(has_both) %>%
    dplyr::pull(stratum) %>%
    as.character()
  
  sub <- sub %>%
    dplyr::filter(as.character(stratum) %in% keep_strata) %>%
    dplyr::mutate(stratum = droplevels(stratum))
  
  if (nrow(sub) < MIN_N_CAT || dplyr::n_distinct(sub$stratum) < 1) return(NULL)
  
  stat_df <- sub %>%
    dplyr::group_by(stratum) %>%
    dplyr::group_modify(~ {
      tmp <- .x %>% dplyr::mutate(health01 = droplevels(health01))
      p_val <- tryCatch(
        wilcox.test(OMWI ~ health01, data = tmp)$p.value,
        error = function(e) NA_real_
      )
      eff <- tryCatch(
        rstatix::wilcox_effsize(tmp, OMWI ~ health01),
        error = function(e) NULL
      )
      med_healthy <- median(tmp$OMWI[tmp$health01 == "Healthy"], na.rm = TRUE)
      med_disease <- median(tmp$OMWI[tmp$health01 == "Disease"], na.rm = TRUE)
      q1_healthy <- quantile(tmp$OMWI[tmp$health01 == "Healthy"], 0.25, na.rm = TRUE, names = FALSE)
      q3_healthy <- quantile(tmp$OMWI[tmp$health01 == "Healthy"], 0.75, na.rm = TRUE, names = FALSE)
      q1_disease <- quantile(tmp$OMWI[tmp$health01 == "Disease"], 0.25, na.rm = TRUE, names = FALSE)
      q3_disease <- quantile(tmp$OMWI[tmp$health01 == "Disease"], 0.75, na.rm = TRUE, names = FALSE)
      tibble(
        n_healthy = sum(tmp$health01 == "Healthy"),
        n_disease = sum(tmp$health01 == "Disease"),
        median_healthy = med_healthy,
        q1_healthy = q1_healthy,
        q3_healthy = q3_healthy,
        IQR_healthy = q3_healthy - q1_healthy,
        median_IQR_healthy = paste0(
          format_es(med_healthy), " [",
          format_es(q1_healthy), ", ",
          format_es(q3_healthy), "]"
        ),
        median_disease = med_disease,
        q1_disease = q1_disease,
        q3_disease = q3_disease,
        IQR_disease = q3_disease - q1_disease,
        median_IQR_disease = paste0(
          format_es(med_disease), " [",
          format_es(q1_disease), ", ",
          format_es(q3_disease), "]"
        ),
        delta_median_disease_minus_healthy = med_disease - med_healthy,
        direction = dplyr::case_when(
          is.na(med_healthy) | is.na(med_disease) ~ NA_character_,
          med_disease < med_healthy ~ "Disease lower",
          med_disease > med_healthy ~ "Disease higher",
          TRUE ~ "No median difference"
        ),
        p_value = p_val,
        p_label = dplyr::case_when(
          is.na(p_val) ~ "P = NA",
          p_val < 0.001 ~ "P < 0.001",
          TRUE ~ paste0("P = ", signif(p_val, 3))
        ),
        sig = dplyr::case_when(
          is.na(p_val) ~ "NA",
          p_val < 0.001 ~ "***",
          p_val < 0.01 ~ "**",
          p_val < 0.05 ~ "*",
          TRUE ~ "ns"
        ),
        effect_type = "rank_biserial",
        effect_size = if (!is.null(eff) && nrow(eff) > 0) eff$effsize[1] else NA_real_,
        effect_magnitude = if (!is.null(eff) && nrow(eff) > 0 && "magnitude" %in% colnames(eff)) {
          eff$magnitude[1]
        } else if (!is.null(eff) && nrow(eff) > 0) {
          interpret_rbc_magnitude(eff$effsize[1])
        } else {
          NA_character_
        }
      )
    }) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      variable = strat_var,
      variable_label = strat_label
    )
  
  pal <- c("Healthy" = COL_HEALTH, "Disease" = COL_DISEASE)
  n_strata <- dplyr::n_distinct(sub$stratum)
  plot_width <- max(2.8, 2.1 * n_strata)
  
  p <- ggplot(sub, aes(x = health01, y = OMWI, fill = health01, color = health01)) +
    geom_violin(width = 0.82, trim = FALSE, color = NA, alpha = 0.72) +
    geom_jitter(width = 0.06, size = 0.75, alpha = 0.45, show.legend = FALSE) +
    geom_boxplot(
      width = 0.18,
      outlier.shape = NA,
      fill = "transparent",
      color = "black",
      linewidth = 0.4
    ) +
    ggpubr::stat_compare_means(
      comparisons = list(c("Healthy", "Disease")),
      method = "wilcox.test",
      label = "p.signif",
      tip.length = 0.015,
      bracket.size = 0.35,
      size = 4,
      hide.ns = FALSE
    ) +
    facet_wrap(~ stratum, nrow = 1, scales = "free_x") +
    scale_fill_manual(values = pal, drop = FALSE) +
    scale_color_manual(values = pal, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.30))) +
    labs(
      title = NULL,
      subtitle = NULL,
      x = strat_label,
      y = "OMWI"
    ) +
    theme_omwi_clean(base_size = 10) +
    theme(
      axis.text.x = element_text(color = "black", angle = 25, hjust = 1, vjust = 1),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.margin = margin(5.5, 8, 5.5, 5.5)
    )
  
  if (!is.null(out_prefix)) {
    ggsave(paste0(out_prefix, ".png"), p, width = plot_width, height = 3.2, units = "in", dpi = 600, bg = "white")
    ggsave(paste0(out_prefix, ".pdf"), p, width = plot_width, height = 3.2, units = "in", dpi = 600, bg = "white")
  }
  
  list(plot = p, stats = stat_df)
}

# =========================================================
# 9. Result containers
# =========================================================
focus_results <- list()
focus_summary_list <- list()
pairwise_list <- list()
or_list <- list()
overall_plot_list <- list()
or_plot_list <- list()

# =========================================================
# 10. Overall analysis module for categorical variables
# =========================================================
run_categorical_block <- function(dat, var, plot_label = var, xlab = var, out_stub = var) {
  if (!var %in% colnames(dat)) return(NULL)
  
  sub <- dat %>%
    select(sample_id, OMWI, all_of(var)) %>%
    filter(!is.na(OMWI), !is.na(.data[[var]])) %>%
    mutate(group = .data[[var]])
  
  grp_tab <- table(sub$group)
  keep_grp <- names(grp_tab[grp_tab >= MIN_PER_GRP])
  
  sub <- sub %>% filter(group %in% keep_grp)
  sub$group <- droplevels(factor(sub$group))
  if (nrow(sub) < MIN_N_CAT || nlevels(sub$group) < 2) return(NULL)
  
  overall_p <- get_overall_p(sub)
  overall_eff <- get_overall_effect(sub)
  
  tmp_sum <- sub %>%
    group_by(group) %>%
    summarise(
      n = n(),
      median_OMWI = median(OMWI, na.rm = TRUE),
      IQR_OMWI = IQR(OMWI, na.rm = TRUE),
      mean_OMWI = mean(OMWI, na.rm = TRUE),
      sd_OMWI = sd(OMWI, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(variable = var)
  
  tmp_pw <- get_pairwise_p(sub)
  if (!is.null(tmp_pw)) tmp_pw <- tmp_pw %>% mutate(variable = var)
  
  tmp_pw_eff <- get_pairwise_effect(sub)
  if (!is.null(tmp_pw_eff)) {
    tmp_pw_eff <- tmp_pw_eff %>% mutate(variable = var)
  }
  
  if (!is.null(tmp_pw) && !is.null(tmp_pw_eff)) {
    tmp_pw <- tmp_pw %>%
      left_join(tmp_pw_eff, by = c("group1", "group2", "variable"))
  }
  
  focus_results[[var]] <<- tibble(
    variable = var,
    variable_type = "categorical",
    overall_test = ifelse(nlevels(droplevels(sub$group)) == 2, "Wilcoxon", "Kruskal-Wallis"),
    overall_p = overall_p
  ) %>%
    bind_cols(overall_eff)
  
  focus_summary_list[[var]] <<- tmp_sum
  pairwise_list[[var]] <<- tmp_pw
  
  invisible(list(summary = tmp_sum, pairwise = tmp_pw, overall_p = overall_p))
}

# =========================================================
# 11. Run overall statistics first
# =========================================================
run_categorical_block(
  dat = dat,
  var = "smoker",
  plot_label = "smoker",
  xlab = "Smoking status",
  out_stub = "smoker"
)

run_categorical_block(
  dat = dat,
  var = "H_pylori_status",
  plot_label = "H.pylori",
  xlab = "H.pylori status",
  out_stub = "H_pylori_status"
)

if ("BMI_group" %in% colnames(dat)) {
  sub <- dat %>%
    select(sample_id, OMWI, BMI, BMI_group) %>%
    filter(!is.na(OMWI), !is.na(BMI_group)) %>%
    mutate(group = BMI_group)
  
  grp_tab <- table(sub$group)
  keep_grp <- names(grp_tab[grp_tab >= MIN_PER_GRP])
  
  sub <- sub %>% filter(group %in% keep_grp)
  sub$group <- factor(sub$group, levels = c("18.5-24.9", "<18.5", "25-29.9", "≥30"))
  sub <- sub %>% filter(!is.na(group))
  
  if (nrow(sub) >= MIN_N_CAT && nlevels(droplevels(sub$group)) >= 2) {
    overall_p <- get_overall_p(sub)
    overall_eff <- get_overall_effect(sub)
    
    tmp_sum <- sub %>%
      group_by(group) %>%
      summarise(
        n = n(),
        median_OMWI = median(OMWI, na.rm = TRUE),
        IQR_OMWI = IQR(OMWI, na.rm = TRUE),
        mean_OMWI = mean(OMWI, na.rm = TRUE),
        sd_OMWI = sd(OMWI, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(variable = "BMI_group")
    
    tmp_pw <- get_pairwise_p(sub)
    if (!is.null(tmp_pw)) tmp_pw <- tmp_pw %>% mutate(variable = "BMI_group")
    
    tmp_pw_eff <- get_pairwise_effect(sub)
    if (!is.null(tmp_pw_eff)) tmp_pw_eff <- tmp_pw_eff %>% mutate(variable = "BMI_group")
    
    if (!is.null(tmp_pw) && !is.null(tmp_pw_eff)) {
      tmp_pw <- tmp_pw %>%
        left_join(tmp_pw_eff, by = c("group1", "group2", "variable"))
    }
    
    focus_results[["BMI_group"]] <- tibble(
      variable = "BMI_group",
      variable_type = "categorical",
      overall_test = ifelse(nlevels(droplevels(sub$group)) == 2, "Wilcoxon", "Kruskal-Wallis"),
      overall_p = overall_p
    ) %>%
      bind_cols(overall_eff)
    
    focus_summary_list[["BMI_group"]] <- tmp_sum
    pairwise_list[["BMI_group"]] <- tmp_pw
  }
}

if ("age_group" %in% colnames(dat)) {
  sub <- dat %>%
    select(sample_id, OMWI, age_num_strict, age_group) %>%
    filter(!is.na(OMWI), !is.na(age_num_strict), !is.na(age_group)) %>%
    mutate(group = age_group)
  
  grp_tab <- table(sub$group)
  keep_grp <- names(grp_tab[grp_tab >= MIN_PER_GRP])
  
  sub <- sub %>% filter(group %in% keep_grp)
  sub$group <- factor(sub$group, levels = c("18-44", "45-59", "60-74", "75-90"))
  sub <- sub %>% filter(!is.na(group))
  
  if (nrow(sub) >= MIN_N_CAT && nlevels(droplevels(sub$group)) >= 2) {
    overall_p <- get_overall_p(sub)
    overall_eff <- get_overall_effect(sub)
    
    tmp_sum <- sub %>%
      group_by(group) %>%
      summarise(
        n = n(),
        median_OMWI = median(OMWI, na.rm = TRUE),
        IQR_OMWI = IQR(OMWI, na.rm = TRUE),
        mean_OMWI = mean(OMWI, na.rm = TRUE),
        sd_OMWI = sd(OMWI, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(variable = "age_group")
    
    tmp_pw <- get_pairwise_p(sub)
    if (!is.null(tmp_pw)) tmp_pw <- tmp_pw %>% mutate(variable = "age_group")
    
    tmp_pw_eff <- get_pairwise_effect(sub)
    if (!is.null(tmp_pw_eff)) tmp_pw_eff <- tmp_pw_eff %>% mutate(variable = "age_group")
    
    if (!is.null(tmp_pw) && !is.null(tmp_pw_eff)) {
      tmp_pw <- tmp_pw %>%
        left_join(tmp_pw_eff, by = c("group1", "group2", "variable"))
    }
    
    focus_results[["age_group"]] <- tibble(
      variable = "age_group",
      variable_type = "categorical",
      overall_test = ifelse(nlevels(droplevels(sub$group)) == 2, "Wilcoxon", "Kruskal-Wallis"),
      overall_p = overall_p
    ) %>%
      bind_cols(overall_eff)
    
    focus_summary_list[["age_group"]] <- tmp_sum
    pairwise_list[["age_group"]] <- tmp_pw
  }
}

run_categorical_block(
  dat = dat,
  var = "gender",
  plot_label = "gender",
  xlab = "Gender",
  out_stub = "gender"
)

for (v in setdiff(prior_cat_vars, c("smoker", "H_pylori_status", "BMI_group", "age_group", "gender"))) {
  sub <- dat %>%
    select(sample_id, OMWI, all_of(v)) %>%
    mutate(group = as.character(.data[[v]])) %>%
    filter(!is.na(OMWI), !is.na(group), group != "")
  
  if (nrow(sub) < MIN_N_CAT) next
  
  grp_tab <- table(sub$group)
  keep_grp <- names(grp_tab[grp_tab >= MIN_PER_GRP])
  
  sub <- sub %>% filter(group %in% keep_grp)
  if (n_distinct(sub$group) < 2) next
  if (n_distinct(sub$group) > MAX_CAT_LEVEL) next
  
  ord <- names(sort(table(sub$group), decreasing = FALSE))
  sub$group <- factor(sub$group, levels = ord)
  
  overall_p <- get_overall_p(sub)
  overall_eff <- get_overall_effect(sub)
  
  tmp_sum <- sub %>%
    group_by(group) %>%
    summarise(
      n = n(),
      median_OMWI = median(OMWI, na.rm = TRUE),
      IQR_OMWI = IQR(OMWI, na.rm = TRUE),
      mean_OMWI = mean(OMWI, na.rm = TRUE),
      sd_OMWI = sd(OMWI, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(variable = v)
  
  tmp_pw <- get_pairwise_p(sub)
  if (!is.null(tmp_pw)) tmp_pw <- tmp_pw %>% mutate(variable = v)
  
  tmp_pw_eff <- get_pairwise_effect(sub)
  if (!is.null(tmp_pw_eff)) tmp_pw_eff <- tmp_pw_eff %>% mutate(variable = v)
  
  if (!is.null(tmp_pw) && !is.null(tmp_pw_eff)) {
    tmp_pw <- tmp_pw %>%
      left_join(tmp_pw_eff, by = c("group1", "group2", "variable"))
  }
  
  focus_results[[v]] <- tibble(
    variable = v,
    variable_type = "categorical",
    overall_test = ifelse(nlevels(droplevels(sub$group)) == 2, "Wilcoxon", "Kruskal-Wallis"),
    overall_p = overall_p
  ) %>%
    bind_cols(overall_eff)
  
  focus_summary_list[[v]] <- tmp_sum
  pairwise_list[[v]] <- tmp_pw
}

# =========================================================
# 12. Output overall violin plots, OMWI < 0 OR plots, and paired layouts
# =========================================================
or_vars <- unique(c(
  "smoker", "H_pylori_status", "BMI_group", "age_group", "gender",
  "RBC_group", "WBC_group", "IgG_group", "IgM_group"
))
or_vars <- intersect(or_vars, colnames(dat))

for (v in or_vars) {
  ref_level <- if (v %in% names(OR_REF_MAP)) OR_REF_MAP[[v]] else NULL
  
  plot_label <- dplyr::case_when(
    v == "BMI_group" ~ "Abnormal BMI",
    v == "age_group" ~ "Age group",
    v == "H_pylori_status" ~ "H.pylori status",
    v == "gender" ~ "Gender",
    TRUE ~ v
  )
  
  xlab_now <- dplyr::case_when(
    v == "smoker" ~ "Smoking status",
    v == "BMI_group" ~ "Abnormal BMI",
    v == "age_group" ~ "Age category",
    v == "H_pylori_status" ~ "H.pylori status",
    v == "gender" ~ "Gender",
    TRUE ~ plot_label
  )
  
  # -------- Overall violin plot --------
  sub_violin <- dat %>%
    dplyr::select(sample_id, OMWI, all_of(v)) %>%
    dplyr::filter(!is.na(OMWI), !is.na(.data[[v]])) %>%
    dplyr::mutate(group = .data[[v]])
  
  grp_tab <- table(sub_violin$group)
  keep_grp <- names(grp_tab[grp_tab >= MIN_PER_GRP])
  
  sub_violin <- sub_violin %>%
    dplyr::filter(group %in% keep_grp)
  
  if (nrow(sub_violin) < MIN_N_CAT || dplyr::n_distinct(sub_violin$group) < 2) next
  
  if (!is.null(ref_level) && ref_level %in% unique(as.character(sub_violin$group))) {
    levs <- c(ref_level, setdiff(unique(as.character(sub_violin$group)), ref_level))
    sub_violin$group <- factor(as.character(sub_violin$group), levels = levs)
  } else {
    sub_violin$group <- droplevels(factor(sub_violin$group))
  }
  
  violin_res <- plot_group_violin(
    sub = sub_violin,
    var = plot_label,
    xlab = xlab_now,
    out_file = file.path(
      OUT_DIR, "plots_focus", "categorical",
      paste0(make_valid_filename(v), "_overall_violin.png")
    )
  )
  if (!is.null(violin_res)) {
    overall_plot_list[[v]] <- violin_res$plot
  }
  plot_group_violin(
    sub = sub_violin,
    var = plot_label,
    xlab = xlab_now,
    out_file = file.path(
      OUT_DIR, "plots_focus", "categorical",
      paste0(make_valid_filename(v), "_overall_violin.pdf")
    )
  )
  
  # -------- OR plot: OMWI < 0 only --------
  or_df <- fit_or_by_variable_binary(
    dat = dat,
    var = v,
    ref_level = ref_level,
    outcome_type = "lt0"
  )
  if (is.null(or_df)) next
  
  or_list[[paste0(v, "_lt0")]] <- or_df
  
  or_plot <- plot_or_forest_binary(
    or_df = or_df,
    var_label = plot_label,
    outcome_type = "lt0",
    out_file = file.path(
      OUT_DIR, "plots_focus", "OR",
      paste0(make_valid_filename(v), "_OR_lt0_forest.png")
    )
  )
  if (!is.null(or_plot)) {
    or_plot_list[[v]] <- or_plot
  }
  plot_or_forest_binary(
    or_df = or_df,
    var_label = plot_label,
    outcome_type = "lt0",
    out_file = file.path(
      OUT_DIR, "plots_focus", "OR",
      paste0(make_valid_filename(v), "_OR_lt0_forest.pdf")
    )
  )
  
  # -------- Paired layout: overall violin + OR --------
  if (!is.null(violin_res) && !is.null(or_plot)) {
    plot_violin_or_pair(
      violin_plot = violin_res$plot,
      or_plot = or_plot,
      out_file = file.path(
        OUT_DIR, "plots_focus", "categorical",
        paste0(make_valid_filename(v), "_overall_violin_plus_OR_lt0.png")
      )
    )
    plot_violin_or_pair(
      violin_plot = violin_res$plot,
      or_plot = or_plot,
      out_file = file.path(
        OUT_DIR, "plots_focus", "categorical",
        paste0(make_valid_filename(v), "_overall_violin_plus_OR_lt0.pdf")
      )
    )
  }
}

# =========================================================
# 13. Stratified healthy vs disease robustness figures
# =========================================================
stratified_plot_list <- list()
stratified_stats_list <- list()

stratified_specs <- tibble::tribble(
  ~var, ~label, ~file_stub,
  "smoking_strata_hd", "Smoking status", "smoking",
  "age_strata_hd", "Age group", "age",
  "BMI_strata_hd", "BMI category", "BMI"
)

for (i in seq_len(nrow(stratified_specs))) {
  spec <- stratified_specs[i, ]
  if (!spec$var %in% colnames(dat)) next
  
  strat_res <- plot_stratified_health_comparison(
    dat = dat,
    strat_var = spec$var,
    strat_label = spec$label,
    out_prefix = file.path(
      OUT_DIR, "plots_focus", "stratified",
      paste0("healthy_vs_disease_by_", spec$file_stub)
    )
  )
  
  if (!is.null(strat_res)) {
    stratified_plot_list[[spec$file_stub]] <- strat_res$plot
    stratified_stats_list[[spec$file_stub]] <- strat_res$stats %>%
      dplyr::mutate(panel = spec$file_stub)
  }
}

stratified_stats_df <- bind_rows(stratified_stats_list)

if (length(stratified_plot_list) > 0) {
  stratified_combined <- patchwork::wrap_plots(stratified_plot_list, ncol = 1)
  stratified_combined <- stratified_combined +
    patchwork::plot_annotation(
      tag_levels = "A",
      theme = theme(plot.tag = element_text(face = "bold", size = 12))
    )
  
  ggsave(
    filename = file.path(OUT_DIR, "plots_focus", "stratified", "stratified_healthy_vs_disease_OMWI.png"),
    plot = stratified_combined,
    width = 5,
    height = 7,
    units = "in",
    dpi = 600,
    bg = "white"
  )
  ggsave(
    filename = file.path(OUT_DIR, "plots_focus", "stratified", "stratified_healthy_vs_disease_OMWI.pdf"),
    plot = stratified_combined,
    width = 5,
    height = 7,
    units = "in",
    dpi = 600,
    bg = "white"
  )
}

# =========================================================
# 14. Combined 3 x 3 host/lifestyle figure
# =========================================================
combined_9panel_specs <- list(
  list(key = "smoking", overall = "smoker", stratified = "smoking", or = "smoker",
       titles = c("Smoking: overall OMWI trend", "Smoking: disease comparison by smoking", "Smoking: OR")),
  list(key = "age", overall = "age_group", stratified = "age", or = "age_group",
       titles = c("Age: overall OMWI trend", "Age: disease comparison by age", "Age: OR")),
  list(key = "BMI", overall = "BMI_group", stratified = "BMI", or = "BMI_group",
       titles = c("BMI: overall OMWI trend", "BMI: disease comparison by BMI", "BMI: OR"))
)

combined_9panel_list <- list()

for (spec in combined_9panel_specs) {
  p_overall <- overall_plot_list[[spec$overall]]
  p_strata  <- stratified_plot_list[[spec$stratified]]
  p_or      <- or_plot_list[[spec$or]]

  if (!is.null(p_overall)) {
    combined_9panel_list[[length(combined_9panel_list) + 1]] <-
      add_panel_title(p_overall, spec$titles[1])
  }
  if (!is.null(p_strata)) {
    combined_9panel_list[[length(combined_9panel_list) + 1]] <-
      add_panel_title(p_strata, spec$titles[2])
  }
  if (!is.null(p_or)) {
    combined_9panel_list[[length(combined_9panel_list) + 1]] <-
      add_panel_title(p_or, spec$titles[3])
  }
}

if (length(combined_9panel_list) == 9) {
  combined_9panel <- patchwork::wrap_plots(combined_9panel_list, ncol = 3, byrow = TRUE)
  combined_9panel <- combined_9panel +
    patchwork::plot_layout(widths = c(1.05, 1.55, 1.10), heights = c(1, 1, 1))
  combined_9panel <- combined_9panel +
    patchwork::plot_annotation(
      tag_levels = "A",
      theme = theme(plot.tag = element_text(face = "bold", color = "black", size = 12))
    )

  ggsave(
    filename = file.path(OUT_DIR, "plots_focus", "host_lifestyle_OMWI_3x3.png"),
    plot = combined_9panel,
    width = 13.5,
    height = 10.5,
    units = "in",
    dpi = 600,
    bg = "white"
  )
  ggsave(
    filename = file.path(OUT_DIR, "plots_focus", "host_lifestyle_OMWI_3x3.pdf"),
    plot = combined_9panel,
    width = 13.5,
    height = 10.5,
    units = "in",
    dpi = 600,
    bg = "white"
  )
}

# =========================================================
# 15. Export result tables
# =========================================================
focus_result_df <- bind_rows(focus_results)
focus_summary_df <- bind_rows(focus_summary_list)
focus_pairwise_df <- bind_rows(pairwise_list)
or_df_all <- bind_rows(or_list)

if (nrow(focus_result_df) > 0) {
  focus_result_df <- focus_result_df %>%
    mutate(candidate_for_main_figure = ifelse(overall_p < P_CUTOFF, "YES", "NO")) %>%
    arrange(overall_p)
  
  write_csv(focus_result_df, file.path(OUT_DIR, "08_focus_variables_overall_tests_alphwith_effectsize.csv"))
}

if (nrow(focus_summary_df) > 0) {
  write_csv(focus_summary_df, file.path(OUT_DIR, "08b_focus_variables_group_summary.csv"))
}

if (nrow(focus_pairwise_df) > 0) {
  write_csv(focus_pairwise_df, file.path(OUT_DIR, "08c_focus_variables_pairwise_tests_with_effectsize.csv"))
}

if (nrow(or_df_all) > 0) {
  write_csv(or_df_all, file.path(OUT_DIR, "08g_OR_OMWI_lt0_by_variable.csv"))
}

if (nrow(stratified_stats_df) > 0) {
  write_csv(stratified_stats_df, file.path(OUT_DIR, "08h_stratified_healthy_vs_disease_OMWI_tests.csv"))
}

# =========================================================
# 16. Save the full analysis dataset
# =========================================================
write_csv(dat, file.path(OUT_DIR, "00_analysis_dataset_with_OMWI.csv"))

cat("=====================================================\n")
cat("OMWI focus analysis (train + test) completed.\n")
cat("Outputs: overall violin plots, OMWI < 0 OR plots, stratified plots, and the 3 x 3 combined figure.\n")
cat("Output directory: ", OUT_DIR, "\n")
cat("=====================================================\n")

