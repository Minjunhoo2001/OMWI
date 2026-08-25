#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)
  library(patchwork)
  library(scales)
  library(broom)
  library(forcats)
})

options(stringsAsFactors = FALSE)

# =========================================================
# 0. 路径
# =========================================================
OMWI_DIR <- "/share/home/HeMinjun/metagenomic/Oral_proxy_v2/04.lasso64/omwi_scores"
OUT_DIR  <- "/share/home/HeMinjun/metagenomic/GitHub/results/Figure4"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

LEVEL_TO_USE <- "species"

TRAIN_FILE <- file.path(OMWI_DIR, paste0("OMWI_train_", LEVEL_TO_USE, "_alpha050.csv"))
TEST_FILE  <- file.path(OMWI_DIR, paste0("OMWI_test_",  LEVEL_TO_USE, "_alpha050.csv"))

# =========================================================
# 1. 配色与主题
# =========================================================
COL_1 <- "#4DBBD5"
COL_2 <- "#FEE5D9"
COL_3 <- "#FCAE91"
COL_4 <- "#FB6A4A"
COL_5 <- "#CB181D"

group_fill_map <- c(
  "Healthy" = COL_1,
  "Colorectal polyps" = COL_3,
  "Colorectal cancer" = COL_5,
  "superficial gastritis" = COL_2,
  "atrophic gastritis" = COL_4,
  "Intestinal metaplasia" = COL_5,
  "Sup./Atr. gastritis" = COL_3,
  "Rheumatoid Arthritis_low" = COL_2,
  "Rheumatoid Arthritis_moderate" = COL_4,
  "RA high" = COL_5,
  "RA low/mod" = COL_3
)

panel_label_map <- c(
  Colorectal = "Colorectal progression",
  Gastric    = "Gastric progression",
  RA         = "RA activity"
)

group_order_list <- list(
  Colorectal = c("Healthy", "Colorectal polyps", "Colorectal cancer"),
  Gastric    = c("Healthy", "Sup./Atr. gastritis", "Intestinal metaplasia"),
  RA         = c("Healthy", "RA low/mod", "RA high")
)

# Global line style used for panel borders and neutral outlines
BORDER_COL <- "black"
BORDER_LWD <- 0.3

theme_cell <- function(base_size = 11){
  theme_classic(base_size = base_size) +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_blank(),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      axis.ticks = element_line(color = BORDER_COL, linewidth = BORDER_LWD),
      strip.background = element_rect(fill = "white", color = BORDER_COL, linewidth = BORDER_LWD),
      strip.text = element_text(face = "bold", color = "black"),
      panel.border = element_rect(fill = NA, color = BORDER_COL, linewidth = BORDER_LWD),
      legend.title = element_text(face = "bold", color = "black"),
      legend.text = element_text(color = "black"),
      plot.title = element_text(face = "bold", color = "black", hjust = 0.5)
    )
}

fmt_p <- function(p){
  dplyr::case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "<0.001",
    TRUE ~ sprintf("%.3f", p)
  )
}

p_to_star <- function(p){
  dplyr::case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE ~ ""
  )
}

# =========================================================
# 2. 读取 OMWI 数据（train + test）
# =========================================================
read_one_omwi <- function(path, set_name){
  if (!file.exists(path)) stop("File not found: ", path)
  df <- fread(path, data.table = FALSE, check.names = FALSE)
  df$set <- set_name
  df
}

score_df_all <- bind_rows(
  read_one_omwi(TRAIN_FILE, "train"),
  read_one_omwi(TEST_FILE,  "test")
)

need_cols <- c("OMWI", "group", "study_name")
miss_cols <- setdiff(need_cols, colnames(score_df_all))
if (length(miss_cols) > 0) stop("缺少必要列: ", paste(miss_cols, collapse = ", "))

score_df_all <- score_df_all %>%
  mutate(
    OMWI = as.numeric(OMWI),
    group = trimws(as.character(group)),
    study_name = trimws(as.character(study_name)),
    group = case_when(
      group %in% c("control", "Control", "healthy", "Healthy") ~ "Healthy",
      TRUE ~ group
    )
  ) %>%
  filter(
    !is.na(OMWI),
    !is.na(group), group != "",
    !is.na(study_name), study_name != ""
  )

if ("Disease activity" %in% colnames(score_df_all)) {
  score_df_all[["Disease activity"]] <- trimws(as.character(score_df_all[["Disease activity"]]))
}
if ("Stage" %in% colnames(score_df_all)) {
  score_df_all$Stage <- trimws(as.character(score_df_all$Stage))
}
if ("stage" %in% colnames(score_df_all)) {
  score_df_all$stage <- trimws(as.character(score_df_all$stage))
}

cat("\n========== unique groups ==========\n")
print(sort(unique(score_df_all$group)))

# =========================================================
# 3. 定义三个 panel 数据
# =========================================================
crc_groups <- c("Healthy", "Colorectal polyps", "Colorectal cancer")

make_crc_panel <- function(df){
  studies_keep <- df %>%
    filter(group %in% c("Colorectal polyps", "Colorectal cancer")) %>%
    distinct(study_name) %>%
    pull(study_name)
  
  df %>%
    filter(study_name %in% studies_keep, group %in% crc_groups) %>%
    mutate(
      panel_group = "Colorectal",
      panel_label = panel_label_map[["Colorectal"]],
      group_axis = factor(group, levels = group_order_list[["Colorectal"]]),
      stage_num = c(1, 2, 3)[match(group, crc_groups)]
    )
}

gastric_groups <- c("Healthy", "superficial gastritis", "atrophic gastritis", "Intestinal metaplasia")

make_gastric_panel <- function(df){
  studies_keep <- df %>%
    filter(group %in% c("superficial gastritis", "atrophic gastritis", "Intestinal metaplasia")) %>%
    distinct(study_name) %>%
    pull(study_name)
  
  df %>%
    filter(study_name %in% studies_keep, group %in% gastric_groups) %>%
    mutate(
      panel_group = "Gastric",
      panel_label = panel_label_map[["Gastric"]],
      group_axis = case_when(
        group %in% c("superficial gastritis", "atrophic gastritis") ~ "Sup./Atr. gastritis",
        TRUE ~ group
      ),
      group_axis = factor(group_axis, levels = group_order_list[["Gastric"]]),
      stage_num = c(1, 2, 3)[match(
        as.character(group_axis),
        c("Healthy", "Sup./Atr. gastritis", "Intestinal metaplasia")
      )]
    )
}

make_ra_panel <- function(df){
  if (!("Disease activity" %in% colnames(df))) {
    warning("没有 Disease activity 列，跳过 RA 分组")
    return(NULL)
  }
  
  ra_raw <- df %>%
    filter(group == "Rheumatoid Arthritis") %>%
    filter(!is.na(`Disease activity`), `Disease activity` != "") %>%
    filter(`Disease activity` %in% c("low", "moderate", "high")) %>%
    mutate(
      group_axis = case_when(
        `Disease activity` %in% c("low", "moderate") ~ "RA low/mod",
        `Disease activity` == "high" ~ "RA high",
        TRUE ~ NA_character_
      )
    )
  
  studies_keep <- ra_raw %>% distinct(study_name) %>% pull(study_name)
  
  healthy_match <- df %>%
    filter(group == "Healthy", study_name %in% studies_keep) %>%
    mutate(group_axis = "Healthy")
  
  bind_rows(healthy_match, ra_raw) %>%
    mutate(
      panel_group = "RA",
      panel_label = panel_label_map[["RA"]],
      group_axis = factor(group_axis, levels = group_order_list[["RA"]]),
      stage_num = c(1, 2, 3)[match(
        as.character(group_axis),
        c("Healthy", "RA low/mod", "RA high")
      )]
    )
}

panel_df <- bind_rows(
  make_crc_panel(score_df_all),
  make_gastric_panel(score_df_all),
  make_ra_panel(score_df_all)
) %>%
  filter(!is.na(group_axis))

cat("\n========== panel-wise sample counts ==========\n")
print(table(panel_df$panel_group, panel_df$group_axis))

cat("\n========== healthy reference ==========\n")
cat("Healthy reference is calculated within each progression panel.\n")
# =========================================================
# 4. 趋势检验
# =========================================================
trend_df <- panel_df %>%
  group_by(panel_group, panel_label) %>%
  summarise(
    n = n(),
    rho = suppressWarnings(cor(stage_num, OMWI, method = "spearman", use = "complete.obs")),
    p_trend = tryCatch(
      suppressWarnings(cor.test(stage_num, OMWI, method = "spearman")$p.value),
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    p_label = case_when(
      is.na(p_trend) ~ "Trend P = NA",
      p_trend < 0.001 ~ "Trend P < 0.001",
      TRUE ~ paste0("Trend P = ", signif(p_trend, 3))
    )
  )

# =========================================================
# 5. 左图：violin
# =========================================================
comparisons_crc <- list(
  c("Healthy", "Colorectal polyps"),
  c("Colorectal polyps", "Colorectal cancer"),
  c("Healthy", "Colorectal cancer")
)

comparisons_gastric <- list(
  c("Healthy", "Sup./Atr. gastritis"),
  c("Sup./Atr. gastritis", "Intestinal metaplasia"),
  c("Healthy", "Intestinal metaplasia")
)

comparisons_ra <- list(
  c("Healthy", "RA low/mod"),
  c("RA low/mod", "RA high"),
  c("Healthy", "RA high")
)

plot_violin_one <- function(df_sub, comparisons_list, trend_label){
  ymax <- max(df_sub$OMWI, na.rm = TRUE)
  ymin <- min(df_sub$OMWI, na.rm = TRUE)
  yrng <- ymax - ymin
  if (yrng <= 0) yrng <- 1
  
  y_positions <- ymax + c(0.08, 0.18, 0.28) * yrng
  
  median_df <- df_sub %>%
    group_by(group_axis) %>%
    summarise(
      median_OMWI = median(OMWI, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      group_axis = factor(group_axis, levels = levels(df_sub$group_axis))
    )
  
  ggplot(df_sub, aes(x = group_axis, y = OMWI, fill = group_axis)) +
    geom_violin(width = 0.88, trim = FALSE, color = NA, alpha = 0.65) +
    geom_violin(
    aes(color = group_axis),
    width = 0.85,
    fill = NA,
    linewidth = 0.6,
    trim = FALSE
  ) +
    geom_jitter(aes(color = group_axis),
                width = 0.06, size = 0.9, alpha = 0.35,
                show.legend = FALSE) +
    geom_boxplot(width = 0.20, outlier.shape = NA,
                 fill = "transparent", color = BORDER_COL,
                 linewidth = BORDER_LWD) +
    geom_line(
      data = median_df,
      aes(x = group_axis, y = median_OMWI, group = 1),
      inherit.aes = FALSE,
      color = BORDER_COL,
      linewidth = BORDER_LWD
    ) +
    geom_point(
      data = median_df,
      aes(x = group_axis, y = median_OMWI),
      inherit.aes = FALSE,
      shape = 21,
      fill = "white",
      color = BORDER_COL,
      stroke = BORDER_LWD,
      size = 2.5
    ) +
    
    ggpubr::stat_compare_means(
      comparisons = comparisons_list,
      method = "wilcox.test",
      label = "p.signif",
      hide.ns = TRUE,
      size = 3.4,
      color = BORDER_COL,
      bracket.size = BORDER_LWD,
      label.y = y_positions[seq_along(comparisons_list)]
    ) +
    annotate("text", x = Inf, y = Inf, label = trend_label,
             hjust = 1.03, vjust = 1.2, size = 3.2, color = "black") +
    scale_fill_manual(values = group_fill_map, drop = FALSE) +
    scale_color_manual(values = group_fill_map, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.03, 0.24))) +
    labs(x = NULL, y = "OMWI", title = NULL) +
    theme_cell(10.5) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 25, hjust = 1),
      plot.margin = margin(5.5, 6, 5.5, 5.5)
    )
}

# =========================================================
# 6. 连续型 OR：每下降 1 SD 的 OMWI
# =========================================================
run_or_vs_healthy <- function(data, case_group, healthy_group = "Healthy", panel_name = NA_character_) {
  
  # 当前比较的数据：Healthy + 当前疾病组
  df <- data %>%
    filter(as.character(group_axis) %in% c(healthy_group, case_group)) %>%
    mutate(
      case = ifelse(as.character(group_axis) == case_group, 1, 0)
    ) %>%
    filter(!is.na(OMWI), !is.na(case))
  
  if (nrow(df) < 10 || length(unique(df$case)) < 2) {
    return(NULL)
  }
  
  # 关键：panel 内的 Healthy reference
  # 注意这里不是 df 里的 Healthy，而是当前 panel data 里的全部 Healthy
  healthy_ref_df <- data %>%
    filter(as.character(group_axis) == healthy_group) %>%
    filter(!is.na(OMWI))
  
  if (nrow(healthy_ref_df) < 5) {
    return(NULL)
  }
  
  healthy_ref_mean <- mean(healthy_ref_df$OMWI, na.rm = TRUE)
  healthy_ref_sd   <- sd(healthy_ref_df$OMWI, na.rm = TRUE)
  
  if (is.na(healthy_ref_sd) || healthy_ref_sd == 0) {
    return(NULL)
  }
  
  # 用 panel-specific healthy mean/sd 标准化
  df <- df %>%
    mutate(
      OMWI_drop1sd = -(OMWI - healthy_ref_mean) / healthy_ref_sd
    )
  
  if (all(is.na(df$OMWI_drop1sd))) {
    return(NULL)
  }
  
  fit <- tryCatch(
    glm(case ~ OMWI_drop1sd, data = df, family = binomial()),
    error = function(e) NULL
  )
  
  if (is.null(fit)) return(NULL)
  
  res <- broom::tidy(fit, conf.int = TRUE, exponentiate = TRUE) %>%
    filter(term == "OMWI_drop1sd") %>%
    transmute(
      panel_group = panel_name,
      group = case_group,
      OR = estimate,
      CI_low = conf.low,
      CI_high = conf.high,
      p_value = p.value,
      n_case = sum(df$case == 1),
      n_healthy_in_model = sum(df$case == 0),
      healthy_ref_n = nrow(healthy_ref_df),
      healthy_ref_mean = healthy_ref_mean,
      healthy_ref_sd = healthy_ref_sd
    )
  
  return(res)
}

or_crc <- bind_rows(lapply(
  c("Colorectal polyps", "Colorectal cancer"),
  function(g) run_or_vs_healthy(
    data = panel_df %>% filter(panel_group == "Colorectal"),
    case_group = g,
    healthy_group = "Healthy",
    panel_name = "Colorectal"
  )
))

or_gastric <- bind_rows(lapply(
  c("Sup./Atr. gastritis", "Intestinal metaplasia"),
  function(g) run_or_vs_healthy(
    data = panel_df %>% filter(panel_group == "Gastric"),
    case_group = g,
    healthy_group = "Healthy",
    panel_name = "Gastric"
  )
))

or_ra <- bind_rows(lapply(
  c("RA low/mod", "RA high"),
  function(g) run_or_vs_healthy(
    data = panel_df %>% filter(panel_group == "RA"),
    case_group = g,
    healthy_group = "Healthy",
    panel_name = "RA"
  )
))

or_all <- bind_rows(or_crc, or_gastric, or_ra)

or_all_show <- or_all %>%
  mutate(
    OR_CI = paste0(sprintf("%.2f", OR), " (", sprintf("%.2f", CI_low), "–", sprintf("%.2f", CI_high), ")"),
    P = fmt_p(p_value),
    sig = p_to_star(p_value),
    healthy_ref_mean = round(healthy_ref_mean, 4),
    healthy_ref_sd = round(healthy_ref_sd, 4)
  )

# =========================================================
# 7. 右图：连续型 OR 森林图（星号版）
# =========================================================
plot_or_forest_one <- function(or_df, panel_name, row_order){
  
  df_sub0 <- or_df %>%
    filter(panel_group == panel_name) %>%
    mutate(sig = p_to_star(p_value))
  
  if (nrow(df_sub0) == 0) return(NULL)
  
  healthy_n <- unique(df_sub0$healthy_ref_n)[1]
  
  ref_row <- data.frame(
    panel_group = panel_name,
    group = "Healthy",
    OR = 1,
    CI_low = 1,
    CI_high = 1,
    p_value = NA_real_,
    n_case = healthy_n,
    n_healthy_in_model = healthy_n,
    healthy_ref_n = healthy_n,
    healthy_ref_mean = NA_real_,
    healthy_ref_sd = NA_real_,
    sig = "",
    stringsAsFactors = FALSE
  )
  
  df_sub <- bind_rows(ref_row, df_sub0) %>%
    mutate(
      # Keep the OR-panel y-axis clean; sample sizes are retained in the output tables.
      group_label = group,
      is_ref = group == "Healthy"
    )
  
  level_labels <- c("Healthy", row_order)
  
  # 关键：Healthy 放第一行，不要 rev()
  df_sub$group_label <- factor(df_sub$group_label, levels = rev(level_labels))
  
  xmax_data <- max(df_sub$CI_high, df_sub$OR, na.rm = TRUE)
  xmin_data <- min(df_sub$CI_low, df_sub$OR, na.rm = TRUE)
  
  xmax_plot <- max(1.5, xmax_data * 1.45)
  xmin_plot <- min(0.25, xmin_data * 0.85)
  
  df_star <- df_sub %>%
    filter(!is_ref, sig != "") %>%
    mutate(
      x_star = pmax(CI_high * 1.08, OR * 1.12)
    )
  
  p <- ggplot(df_sub, aes(x = OR, y = group_label)) +
    geom_vline(
      xintercept = 1,
      linetype = 2,
      color = BORDER_COL,
      linewidth = BORDER_LWD
    ) +
    
    geom_errorbarh(
      data = df_sub %>% filter(!is_ref),
      aes(xmin = CI_low, xmax = CI_high, color = group),
      height = 0.16,
      linewidth = BORDER_LWD,
      show.legend = FALSE
    ) +
    
    geom_point(
      data = df_sub %>% filter(!is_ref),
      aes(color = group),
      size = 3.2,
      show.legend = FALSE
    ) +
    
    geom_point(
      data = df_sub %>% filter(is_ref),
      size = 3.2,
      shape = 21,
      stroke = BORDER_LWD,
      fill = "white",
      color = BORDER_COL
    ) +
    
    geom_text(
      data = df_star,
      aes(x = x_star, y = group_label, label = sig),
      inherit.aes = FALSE,
      hjust = 0,
      size = 4.4,
      color = BORDER_COL
    ) +
    
    scale_color_manual(values = group_fill_map, drop = FALSE) +
    
    scale_x_log10(
      limits = c(xmin_plot, xmax_plot),
      breaks = c(0.25, 0.5, 1, 2, 4, 8, 16),
      labels = c("0.25", "0.5", "1", "2", "4", "8", "16")
    ) +
    
    labs(
      x = "Odds ratio per 1-SD decrease in OMWI",
      y = NULL,
      subtitle = NULL
    ) +
    
    coord_cartesian(clip = "off") +
    theme_cell(10.5) +
    theme(
      axis.text.y = element_text(color = "black"),
      panel.grid = element_blank(),
      plot.margin = margin(5.5, 28, 5.5, 8),
      plot.subtitle = element_text(hjust = 0.5)
    )
  
  return(p)
}

p_crc_or <- plot_or_forest_one(
  or_all, "Colorectal",
  row_order = c("Colorectal polyps", "Colorectal cancer")
)
p_crc_or
p_gas_or <- plot_or_forest_one(
  or_all, "Gastric",
  row_order = c("Sup./Atr. gastritis", "Intestinal metaplasia")
)

p_ra_or <- plot_or_forest_one(
  or_all, "RA",
  row_order = c("RA low/mod", "RA high")
)

# =========================================================
# 8. 生成三套图（左 violin + 右 OR）
# =========================================================
p_crc_left <- plot_violin_one(
  panel_df %>% filter(panel_group == "Colorectal"),
  comparisons_crc,
  trend_df %>% filter(panel_group == "Colorectal") %>% pull(p_label)
)

fig_crc <- p_crc_left + p_crc_or +
  plot_layout(widths = c(1.20, 0.92)) +
  plot_annotation(title = "Colorectal progression")

p_gas_left <- plot_violin_one(
  panel_df %>% filter(panel_group == "Gastric"),
  comparisons_gastric,
  trend_df %>% filter(panel_group == "Gastric") %>% pull(p_label)
)

fig_gas <- p_gas_left + p_gas_or +
  plot_layout(widths = c(1.20, 0.92)) +
  plot_annotation(title = "Gastric precancerous progression")

p_ra_left <- plot_violin_one(
  panel_df %>% filter(panel_group == "RA"),
  comparisons_ra,
  trend_df %>% filter(panel_group == "RA") %>% pull(p_label)
)

fig_ra <- p_ra_left + p_ra_or +
  plot_layout(widths = c(1.20, 0.92)) +
  plot_annotation(title = "Rheumatoid arthritis activity")

figure_progression_all <- fig_crc / fig_gas / fig_ra +
  plot_layout(heights = c(1, 1, 1))

# =========================================================
# 9. 输出图
# =========================================================
ggsave(
  file.path(OUT_DIR, "Figure_progression_crc_continuousOR_star.pdf"),
  fig_crc, width = 10.2, height = 4.3, units = "in", device = cairo_pdf
)
ggsave(
  file.path(OUT_DIR, "Figure_progression_crc_continuousOR_star.png"),
  fig_crc, width = 10.2, height = 4.3, units = "in", dpi = 600
)

ggsave(
  file.path(OUT_DIR, "Figure_progression_gastric_continuousOR_star.pdf"),
  fig_gas, width = 10.2, height = 4.3, units = "in", device = cairo_pdf
)
ggsave(
  file.path(OUT_DIR, "Figure_progression_gastric_continuousOR_star.png"),
  fig_gas, width = 10.2, height = 4.3, units = "in", dpi = 600
)

ggsave(
  file.path(OUT_DIR, "Figure_progression_ra_continuousOR_star.pdf"),
  fig_ra, width = 10.2, height = 4.3, units = "in", device = cairo_pdf
)
ggsave(
  file.path(OUT_DIR, "Figure_progression_ra_continuousOR_star.png"),
  fig_ra, width = 10.2, height = 4.3, units = "in", dpi = 600
)

ggsave(
  file.path(OUT_DIR, "Figure_progression_all_continuousOR_star.pdf"),
  figure_progression_all, width = 8, height = 8, units = "in", device = cairo_pdf
)
ggsave(
  file.path(OUT_DIR, "Figure_progression_all_continuousOR_star.png"),
  figure_progression_all, width = 8, height = 8, units = "in", dpi = 600
)

# =========================================================
# 10. 输出数据表
# =========================================================
fwrite(panel_df, file.path(OUT_DIR, "panel_data_used.csv"))
fwrite(trend_df, file.path(OUT_DIR, "panel_trend_statistics.csv"))
fwrite(or_all, file.path(OUT_DIR, "OR_per1SDdecrease_all_panels.csv"))
fwrite(or_all_show, file.path(OUT_DIR, "OR_per1SDdecrease_all_panels_show.csv"))

cat("\n[DONE] 输出目录: ", OUT_DIR, "\n")

