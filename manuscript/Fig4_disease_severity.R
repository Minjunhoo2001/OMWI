#!/usr/bin/env Rscript
# ============================================================================
# Fig. 4 | Oral v3.0 | six-panel manuscript figure (a-f)
# a/c/e: OMWI distribution across progression/activity stages
# b/d/f: odds ratios per 1 healthy-reference SD decrease in OMWI
#
# Adapted from: 08.fig4_revised-v2.0.R
# V3 DATA SOURCE: data/model_outputs/final/
#   species_train_in_sample_scores.csv        (final model, train in-sample)
#   external_independent_predictions.csv       (frozen model, independent test)
# CAUTION: these are not cross-validation out-of-fold predictions.
# All inputs and outputs are INSIDE ROOT_DIR. Never mix older Oral_proxy2.0 scores.
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(broom)
})
options(stringsAsFactors = FALSE)

# ---- 0. ONLY ROOT TO MODIFY -------------------------------------------------
ROOT_DIR <- Sys.getenv("OMWI_ROOT", getwd())

V3_MODEL_DIR <- file.path(ROOT_DIR, "data/model_outputs")
INPUT_DIR    <- file.path(V3_MODEL_DIR, "final")
TRAIN_FILE   <- file.path(INPUT_DIR, "species_train_in_sample_scores.csv")
TEST_FILE    <- file.path(INPUT_DIR, "external_independent_predictions.csv")
OUT_DIR      <- file.path(ROOT_DIR, "Fig4")

# Plot text follows the Fig.4 reference PDF. Underlying numeric column is OMHI.
INDEX_LABEL <- "OMWI"  # change to "OMHI" if this is the chosen manuscript name
SHOW_N_IN_FOREST <- FALSE  # reference PDF has plain group labels, no (n=...)

# Output size (inches). 2 columns x 3 rows, panels a/b/c/d/e/f.
FULL_WIDTH  <- 11.0
FULL_HEIGHT <- 13
DPI <- 600

for (f in c(TRAIN_FILE, TEST_FILE)) {
  if (!file.exists(f)) {
    stop("找不到新版数据文件:\n", f,
         "\n请确认 data/model_outputs 已解压，并包含 final 文件夹。",
         call. = FALSE)
  }
}
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- 1. Palette and ggplot theme --------------------------------------------
COL_1 <- "#4DBBD5"
COL_2 <- "#FEE5D9"
COL_3 <- "#FCAE91"
COL_4 <- "#FB6A4A"
COL_5 <- "#CB181D"

GROUP_COLORS <- c(
  "Healthy"                 = COL_1,
  "Colorectal polyps"       = COL_3,
  "Colorectal cancer"       = COL_5,
  "Sup./Atr. gastritis"     = COL_3,
  "Intestinal metaplasia"   = COL_5,
  "RA low/mod"              = COL_3,
  "RA high"                 = COL_5
)

PANEL_LEVELS <- list(
  Colorectal = c("Healthy", "Colorectal polyps", "Colorectal cancer"),
  Gastric    = c("Healthy", "Sup./Atr. gastritis", "Intestinal metaplasia"),
  RA         = c("Healthy", "RA low/mod", "RA high")
)

publication_theme <- function(base_size = 10.2) {
  theme_classic(base_size = base_size) +
    theme(
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.50),
      panel.grid = element_blank(),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      axis.line = element_blank(),
      legend.position = "none",
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.margin = margin(13, 11, 10, 10),
      plot.tag = element_text(face = "bold", size = 13, color = "black"),
      plot.tag.position = c(-0.08, 1.05)
    )
}

star_label <- function(p) {
  vapply(p, function(x) {
    if (!is.finite(x)) ""
    else if (x < 1e-4) "****"
    else if (x < 1e-3) "***"
    else if (x < 1e-2) "**"
    else if (x < 0.05) "*"
    else ""
  }, character(1))
}

trend_label <- function(p) {
  if (!is.finite(p)) return("Trend P = NA")
  paste0("Trend P = ", sprintf("%.2E", p))
}

# ---- 2. Load V3 score tables without guessing old model files ---------------
read_v3_score <- function(file, set_name) {
  # Original model output contains a few Latin-1 characters in unrelated
  # metadata fields. fread's encoding parameter is supported on Windows.
  x <- data.table::fread(file, data.table = FALSE,
                         check.names = FALSE, encoding = "Latin-1")
  must <- c("sample_id", "study_name", "group", "OMHI", "Disease activity")
  missing <- setdiff(must, colnames(x))
  if (length(missing)) {
    stop("数据列不匹配: ", basename(file), "; 缺少 ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  x %>%
    transmute(
      sample_id = trimws(as.character(.data$sample_id)),
      study_name = trimws(as.character(.data$study_name)),
      group_original = trimws(as.character(.data$group)),
      activity_original = trimws(as.character(.data[["Disease activity"]])),
      OMHI = suppressWarnings(as.numeric(.data$OMHI)),
      set = set_name
    )
}

raw_scores <- bind_rows(
  read_v3_score(TRAIN_FILE, "Training (in-sample)"),
  read_v3_score(TEST_FILE, "External validation")
)

# Do not silently drop duplicated samples or mix cohorts.
if (anyNA(raw_scores$sample_id) || any(raw_scores$sample_id == "")) {
  stop("sample_id 存在缺失值; 请核对输入数据。", call. = FALSE)
}
if (anyDuplicated(raw_scores$sample_id)) {
  duplicated_ids <- unique(raw_scores$sample_id[duplicated(raw_scores$sample_id)])
  stop("训练集/外部集存在重复 sample_id: ",
       paste(head(duplicated_ids, 10), collapse = ", "), call. = FALSE)
}

score_data <- raw_scores %>%
  mutate(
    group_lower = tolower(trimws(group_original)),
    group = case_when(
      group_lower %in% c("control", "healthy")          ~ "Healthy",
      group_lower == "colorectal polyps"                ~ "Colorectal polyps",
      group_lower == "colorectal cancer"                ~ "Colorectal cancer",
      group_lower == "superficial gastritis"            ~ "superficial gastritis",
      group_lower == "atrophic gastritis"               ~ "atrophic gastritis",
      group_lower == "intestinal metaplasia"           ~ "Intestinal metaplasia",
      group_lower == "rheumatoid arthritis"            ~ "Rheumatoid Arthritis",
      TRUE                                              ~ group_original
    ),
    activity = tolower(trimws(activity_original))
  ) %>%
  filter(is.finite(OMHI), !is.na(group), group != "",
         !is.na(study_name), study_name != "")

cat("\nV3 full train/test score table:\n")
print(table(score_data$set))


# ---- 3. Stage definitions ------------------------------------

make_crc_panel <- function(df) {
  groups <- c("Colorectal polyps", "Colorectal cancer")
  studies <- unique(df$study_name[df$group %in% groups])
  
  df %>%
    filter(
      study_name %in% studies,
      group %in% c("Healthy", groups)
    ) %>%
    mutate(
      panel = "Colorectal",
      group_axis = as.character(group),
      stage_num = match(group_axis, PANEL_LEVELS$Colorectal)
    )
}

make_gastric_panel <- function(df) {
  original_groups <- c(
    "superficial gastritis",
    "atrophic gastritis",
    "Intestinal metaplasia"
  )
  
  studies <- unique(
    df$study_name[df$group %in% original_groups]
  )
  
  df %>%
    filter(
      study_name %in% studies,
      group %in% c("Healthy", original_groups)
    ) %>%
    mutate(
      panel = "Gastric",
      group_axis = dplyr::if_else(
        group %in% c(
          "superficial gastritis",
          "atrophic gastritis"
        ),
        "Sup./Atr. gastritis",
        as.character(group)
      ),
      stage_num = match(group_axis, PANEL_LEVELS$Gastric)
    )
}

# Normalize explicit clinical labels without redefining disease activity.
normalize_label <- function(x) {
  y <- tolower(trimws(as.character(x)))
  y <- gsub("[_/\\.-]+", " ", y)
  y <- gsub("[[:space:]]+", " ", y)
  y[is.na(x)] <- ""
  trimws(y)
}

make_ra_panel <- function(df) {
  # Important: use explicit RA phenotype labels, not study name alone.
  # Mild, severe and remission are NOT automatically recoded as low/high.
  phenotype <- normalize_label(df$group_original)
  ra_pattern <- "\\brheumatoid[[:space:]]+arthritis\\b|\\bra\\b"
  ra_candidates <- df[grepl(ra_pattern, phenotype, perl = TRUE), , drop = FALSE]
  
  if (nrow(ra_candidates) == 0L) {
    stop(
      "RA 数据为 0：V3 的 group 列未找到明确的 rheumatoid arthritis / RA 标签。",
      "\n请查看 Fig4_all_original_group_activity_counts.csv，",
      "确认 RA 是否存在于这两份新版分数文件中。",
      call. = FALSE
    )
  }
  
  act <- normalize_label(ra_candidates$activity_original)
  grp <- normalize_label(ra_candidates$group_original)
  # Use the group label only when the Disease activity field is missing.
  act_missing <- act %in% c("", "na", "n a", "not available", "unknown")
  stage_text <- ifelse(act_missing, grp, act)
  
  ra_candidates$group_axis <- dplyr::case_when(
    grepl("\\bhigh\\b", stage_text, perl = TRUE) &
      !grepl("\\blow\\b|\\bmoderate\\b", stage_text, perl = TRUE) ~ "RA high",
    grepl("\\blow\\b|\\bmoderate\\b", stage_text, perl = TRUE) &
      !grepl("\\bhigh\\b", stage_text, perl = TRUE) ~ "RA low/mod",
    TRUE ~ NA_character_
  )
  
  # Export all candidates (including unmatched activity labels) for auditing.
  data.table::fwrite(
    ra_candidates %>%
      select(sample_id, study_name, set, group_original,
             activity_original, group_axis, OMHI),
    file.path(OUT_DIR, "Fig4_RA_original_activity_audit.csv")
  )
  cat("\nRA candidates and their original disease-activity labels:\n")
  print(dplyr::as_tibble(ra_candidates %>%
                           count(group_original, activity_original, group_axis, sort = TRUE)),
        n = Inf)
  
  ra <- ra_candidates %>% filter(!is.na(group_axis))
  if (nrow(ra) == 0L) {
    stop(
      "发现 RA 病例，但 Disease activity 无法按明确的 low/moderate/high 标签分组。",
      "\n请查看 Fig4_RA_original_activity_audit.csv 并核对原始临床编码；",
      "不要把 mild/severe/remission 自动当作 low/high。",
      call. = FALSE
    )
  }
  studies <- unique(ra$study_name)
  healthy <- df %>%
    filter(study_name %in% studies, group == "Healthy") %>%
    mutate(group_axis = "Healthy")
  
  bind_rows(healthy, ra) %>%
    mutate(panel = "RA", stage_num = match(group_axis, PANEL_LEVELS$RA))
}

# Always save the raw phenotype/activity labels before selecting stages.
data.table::fwrite(
  score_data %>%
    count(set, study_name, group_original, activity_original, name = "n", sort = TRUE),
  file.path(OUT_DIR, "Fig4_all_original_group_activity_counts.csv")
)

# Build each panel separately for easier debugging
crc_data <- make_crc_panel(score_data)
gastric_data <- make_gastric_panel(score_data)
ra_data <- make_ra_panel(score_data)

cat("\nCRC groups:\n")
print(table(crc_data$group_axis))

cat("\nGastric groups:\n")
print(table(gastric_data$group_axis))

cat("\nRA groups:\n")
print(table(ra_data$group_axis))

# Combine
panel_data <- bind_rows(
  crc_data,
  gastric_data,
  ra_data
) %>%
  filter(!is.na(stage_num)) %>%
  mutate(
    panel = factor(
      panel,
      levels = c("Colorectal", "Gastric", "RA")
    ),
    group_axis = as.character(group_axis)
  )

cat("\nFinal panel sample counts:\n")
print(table(panel_data$panel, panel_data$group_axis))


# Require all three prespecified categories per panel before statistical tests.
for (nm in names(PANEL_LEVELS)) {
  d <- panel_data %>% filter(panel == nm)
  counts <- table(factor(d$group_axis, levels = PANEL_LEVELS[[nm]]))
  cat("\n", nm, " sample counts: ",
      paste(names(counts), counts, sep = "=", collapse = "; "), "\n", sep = "")
  if (any(counts < 5)) {
    stop(
      "Fig4 ", nm, " 缺少预设的三组，或其中某组 n<5。不能计算有效的趋势/OR。",
      "\n请检查 Fig4_all_original_group_activity_counts.csv；若为 RA，",
      "另检查 Fig4_RA_original_activity_audit.csv。",
      call. = FALSE
    )
  }
}

# ---- 4. Trend Spearman and pairwise Wilcoxon --------------------------------
trend_stats <- bind_rows(lapply(names(PANEL_LEVELS), function(nm) {
  d <- panel_data %>% filter(panel == nm)
  t <- suppressWarnings(cor.test(d$stage_num, d$OMHI,
                                 method = "spearman", exact = FALSE))
  data.frame(panel = nm, n = nrow(d), rho = unname(t$estimate),
             p_trend = t$p.value, label = trend_label(t$p.value))
}))

pairwise_stats <- bind_rows(lapply(names(PANEL_LEVELS), function(nm) {
  d <- panel_data %>% filter(panel == nm)
  lev <- PANEL_LEVELS[[nm]]
  pair_inds <- list(c(1, 2), c(2, 3), c(1, 3))
  bind_rows(lapply(seq_along(pair_inds), function(i) {
    ij <- pair_inds[[i]]
    u <- d$OMHI[d$group_axis == lev[ij[1]]]
    v <- d$OMHI[d$group_axis == lev[ij[2]]]
    p <- suppressWarnings(wilcox.test(u, v, exact = FALSE)$p.value)
    data.frame(panel = nm, comparison = i,
               group1 = lev[ij[1]], group2 = lev[ij[2]],
               x1 = ij[1], x2 = ij[2], n1 = length(u), n2 = length(v),
               p_value = p, significance = star_label(p))
  }))
}))

fwrite(trend_stats, file.path(OUT_DIR, "Fig4_Spearman_trend_statistics.csv"))
fwrite(pairwise_stats, file.path(OUT_DIR, "Fig4_pairwise_Wilcoxon_statistics.csv"))

# ---- 5. Violin + boxplot + median trajectory + significance -----------------
plot_violin_panel <- function(panel_name) {
  d <- panel_data %>% filter(panel == panel_name)
  lev <- PANEL_LEVELS[[panel_name]]
  d$group_axis <- factor(d$group_axis, levels = lev)
  med <- d %>%
    group_by(group_axis) %>%
    summarise(median_OMHI = median(OMHI), .groups = "drop")
  ymax <- max(d$OMHI)
  ymin <- min(d$OMHI)
  spread <- max(ymax - ymin, 0.1)
  
  # Leave enough headroom for bracket levels and trend test.
  p <- ggplot(d, aes(x = group_axis, y = OMHI, fill = group_axis)) +
    geom_violin(trim = FALSE, width = 0.9, alpha = 0.68,
                color = NA) +
    geom_jitter(aes(color = group_axis),
                width = 0.055, height = 0, alpha = 0.28,
                size = 0.55, show.legend = FALSE) +
    geom_boxplot(width = 0.17, outlier.shape = NA,
                 fill = "transparent", linewidth = 0.48,
                 color = "black") +
    geom_line(data = med,
              aes(x = group_axis, y = median_OMHI, group = 1),
              inherit.aes = FALSE, linewidth = 0.65,
              color = "black") +
    geom_point(data = med,
               aes(x = group_axis, y = median_OMHI),
               inherit.aes = FALSE, shape = 21, size = 2.4,
               fill = "white", color = "black", stroke = 0.70) +
    scale_fill_manual(values = GROUP_COLORS, drop = FALSE) +
    scale_color_manual(values = GROUP_COLORS, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.02))) +
    labs(x = NULL, y = INDEX_LABEL) +
    coord_cartesian(ylim = c(ymin - 0.05 * spread, ymax + 0.63 * spread),
                    clip = "off") +
    publication_theme(10.2) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1,
                                     size = 9.5),
          plot.margin = margin(17, 10, 24, 13))
  
  stat <- pairwise_stats %>% filter(panel == panel_name, significance != "")
  if (nrow(stat) > 0) {
    for (i in seq_len(nrow(stat))) {
      st <- stat[i, ]
      y <- ymax + c(0.10, 0.22, 0.34)[st$comparison] * spread
      tick <- 0.024 * spread
      p <- p +
        annotate("segment", x = st$x1, xend = st$x2,
                 y = y, yend = y, linewidth = 0.48, color = "black") +
        annotate("segment", x = st$x1, xend = st$x1,
                 y = y - tick, yend = y, linewidth = 0.48) +
        annotate("segment", x = st$x2, xend = st$x2,
                 y = y - tick, yend = y, linewidth = 0.48) +
        annotate("text", x = (st$x1 + st$x2) / 2,
                 y = y + 0.032 * spread,
                 label = st$significance, size = 3.7, vjust = 0)
    }
  }
  tlabel <- trend_stats$label[trend_stats$panel == panel_name]
  p + annotate("text", x = 3.08, y = ymax + 0.535 * spread,
               label = tlabel, hjust = 1, vjust = 0,
               size = 2.95, color = "black")
}

# ---- 6. Logistic OR: disease group vs matching healthy controls -------------
# Per original Fig.4 code: -1 standardized unit of the index, where SD is
# calculated in the *healthy subjects within the respective panel*.
# Unadjusted binary logistic regression and profile-likelihood 95% CI.
run_case_vs_healthy_or <- function(panel_name, case_group) {
  d_all <- panel_data %>% filter(panel == panel_name)
  healthy <- d_all %>% filter(group_axis == "Healthy")
  ref_mean <- mean(healthy$OMHI)
  ref_sd <- sd(healthy$OMHI)
  if (!is.finite(ref_sd) || ref_sd <= 0) {
    stop(panel_name, " 的 Healthy 组 SD 为 0/NA，无法计算 OR")
  }
  d <- d_all %>%
    filter(group_axis %in% c("Healthy", case_group)) %>%
    mutate(case = as.integer(group_axis == case_group),
           index_drop_1sd = -(OMHI - ref_mean) / ref_sd)
  fit <- glm(case ~ index_drop_1sd, data = d, family = binomial())
  if (!isTRUE(fit$converged)) {
    warning(panel_name, " / ", case_group,
            " logistic 回归未收敛，请检查 OR 稳定性。")
  }
  x <- suppressMessages(suppressWarnings(
    broom::tidy(fit, conf.int = TRUE, exponentiate = TRUE)
  ))
  x <- x[x$term == "index_drop_1sd", , drop = FALSE]
  if (nrow(x) != 1 || any(!is.finite(c(x$estimate, x$conf.low, x$conf.high)))) {
    stop(panel_name, " / ", case_group,
         " 的 OR/CI 不可用。请检查模型分离或数值不稳定。")
  }
  data.frame(
    panel = panel_name, group = case_group,
    OR = x$estimate, CI_low = x$conf.low, CI_high = x$conf.high,
    p_value = x$p.value, significance = star_label(x$p.value),
    n_case = sum(d$case == 1), n_healthy = sum(d$case == 0),
    healthy_ref_mean = ref_mean, healthy_ref_sd = ref_sd,
    formula = "case ~ index_drop_1sd",
    ci_method = "profile likelihood 95% CI"
  )
}

or_stats <- bind_rows(lapply(names(PANEL_LEVELS), function(nm) {
  groups <- PANEL_LEVELS[[nm]][-1]
  bind_rows(lapply(groups, function(g) run_case_vs_healthy_or(nm, g)))
}))

or_stats <- or_stats %>%
  mutate(OR_CI = sprintf("%.2f (%.2f-%.2f)", OR, CI_low, CI_high))
fwrite(or_stats, file.path(OUT_DIR, "Fig4_unadjusted_OR_per_healthySD.csv"))

# ---- 7. Forest panels b, d, f -----------------------------------------------
plot_or_panel <- function(panel_name) {
  lev <- PANEL_LEVELS[[panel_name]]
  x <- or_stats %>% filter(panel == panel_name)
  healthy_n <- unique(x$n_healthy)[1]
  ref <- data.frame(panel = panel_name, group = "Healthy", OR = 1,
                    CI_low = 1, CI_high = 1, p_value = NA_real_,
                    significance = "", n_case = NA_integer_,
                    n_healthy = healthy_n, stringsAsFactors = FALSE)
  d <- bind_rows(ref, x) %>%
    mutate(
      label = if (SHOW_N_IN_FOREST) {
        ifelse(group == "Healthy",
               paste0("Healthy (n=", n_healthy, ")"),
               paste0(group, " (n=", n_case, ")"))
      } else group,
      is_healthy = group == "Healthy"
    )
  label_order <- if (SHOW_N_IN_FOREST) {
    vapply(lev, function(v) d$label[match(v, d$group)], character(1))
  } else lev
  d$label <- factor(d$label, levels = rev(label_order))
  d$label <- droplevels(d$label)
  
  # Expand in powers of two so no observed CI is clipped.
  xmax <- max(16, 2 ^ ceiling(log2(max(d$CI_high, na.rm = TRUE) * 1.32)))
  xmin <- min(0.25, 2 ^ floor(log2(min(d$CI_low, na.rm = TRUE) * 0.78)))
  tick_all <- 2 ^ (-5:10)
  tick <- tick_all[tick_all >= xmin & tick_all <= xmax]
  label_tick <- function(z) format(z, scientific = FALSE, trim = TRUE)
  
  starred <- d %>%
    filter(!is_healthy, significance != "") %>%
    mutate(star_x = pmin(CI_high * 1.09, xmax * 0.87))
  
  ggplot(d, aes(x = OR, y = label)) +
    geom_vline(xintercept = 1, linetype = "dashed",
               color = "grey45", linewidth = 0.55) +
    geom_segment(data = d %>% filter(!is_healthy),
                 aes(x = CI_low, xend = CI_high, y = label, yend = label,
                     color = group), linewidth = 0.78,
                 show.legend = FALSE) +
    geom_point(data = d %>% filter(!is_healthy),
               aes(color = group), size = 3.05, show.legend = FALSE) +
    geom_point(data = d %>% filter(is_healthy),
               shape = 21, size = 3, fill = "white",
               stroke = 0.75, color = "black") +
    geom_text(data = starred,
              aes(x = star_x, y = label, label = significance),
              inherit.aes = FALSE,
              hjust = 0, size = 3.9, color = "black") +
    scale_color_manual(values = GROUP_COLORS, drop = FALSE) +
    scale_x_log10(limits = c(xmin, xmax), breaks = tick,
                  labels = label_tick, expand = expansion(mult = c(0, 0))) +
    scale_y_discrete(drop = FALSE) +
    labs(x = paste0("Odds ratio per 1-SD decrease in ", INDEX_LABEL),
         y = NULL) +
    coord_cartesian(clip = "off") +
    publication_theme(10.2) +
    theme(axis.text.y = element_text(size = 9.0),
          axis.text.x = element_text(size = 9.1),
          axis.title.x = element_text(size = 9.6),
          plot.margin = margin(17, 12, 24, 11))
}

# ---- 8. Create all 6 panels and composite exactly in Fig4 order ------------
p_a <- plot_violin_panel("Colorectal")
p_b <- plot_or_panel("Colorectal")
p_c <- plot_violin_panel("Gastric")
p_d <- plot_or_panel("Gastric")
p_e <- plot_violin_panel("RA")
p_f <- plot_or_panel("RA")

fig4 <- patchwork::wrap_plots(p_a, p_b, p_c, p_d, p_e, p_f,
                              ncol = 2, byrow = TRUE) +
  patchwork::plot_layout(widths = c(1.23, 1.00), heights = c(1, 1, 1)) +
  patchwork::plot_annotation(tag_levels = "a")

# ---- 9. Save PDF/PNG, independent panels and source tables -------------------
save_pair <- function(p, stem, width, height) {
  pdf_file <- file.path(OUT_DIR, paste0(stem, ".pdf"))
  png_file <- file.path(OUT_DIR, paste0(stem, ".png"))
  if (capabilities("cairo")) {
    ggsave(pdf_file, p, width = width, height = height, units = "in",
           device = grDevices::cairo_pdf)
  } else {
    ggsave(pdf_file, p, width = width, height = height, units = "in",
           device = "pdf")
  }
  ggsave(png_file, p, width = width, height = height, units = "in",
         dpi = DPI, limitsize = FALSE)
}

save_pair(fig4, "Fig4_Oral_v3_combined", FULL_WIDTH, FULL_HEIGHT)
save_pair(p_a, "Fig4a_colorectal_violin", 5.65, 3.20)
save_pair(p_b, "Fig4b_colorectal_OR", 5.15, 3.20)
save_pair(p_c, "Fig4c_gastric_violin", 5.65, 3.20)
save_pair(p_d, "Fig4d_gastric_OR", 5.15, 3.20)
save_pair(p_e, "Fig4e_RA_violin", 5.65, 3.20)
save_pair(p_f, "Fig4f_RA_OR", 5.15, 3.20)

cat("\n==============================================\n")
cat("Fig.4 created successfully.\n")
cat("INPUT train : ", TRAIN_FILE, "\n", sep = "")
cat("INPUT test  : ", TEST_FILE, "\n", sep = "")
cat("OUTPUT      : ", OUT_DIR, "\n", sep = "")
cat("Fig4_Oral_v3_combined.pdf / png; individual a-f panels; CSV audit tables (including metadata and RA activity diagnostics).\n")
cat("NOTE: Training scores are in-sample; OR is unadjusted.\n")
cat("==============================================\n")
print(fig4)

