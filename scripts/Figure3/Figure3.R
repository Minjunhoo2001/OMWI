#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
  library(tibble)
})

options(stringsAsFactors = FALSE)

# =========================================================
# 0. Path and figure settings
# =========================================================
ROOT_DIR <- "/share/home/HeMinjun/metagenomic/Oral_proxy_v2/04.lasso64"

OMHI_DIR <- file.path(ROOT_DIR, "omhi_scores")
OUT_DIR  <- "/share/home/HeMinjun/metagenomic/GitHub/results/Figure3"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

LEVEL <- "species"
FIGURE_PREFIX <- "Fig3"

# Output sizes. These keep the final sizes used in the original script.
PANEL_WIDTH  <- 8
PANEL_HEIGHT <- 5
COMBINED_WIDTH  <- 8
COMBINED_HEIGHT <- 9

# =========================================================
# 1. Colors and theme
# =========================================================
COL_GREY <- "#7F7F7F"

system_color_map <- c(
  "Healthy" = "#2CA25F",
  "Oral Disease" = "#C7A6D8",
  "Digestive Disease" = "#D97C8A",
  "Respiratory Disease" = "#8FA8C9",
  "Immune Disease" = "#D8C97A",
  "Metabolic Disease" = "#4198AC",
  "Abnormal BMI" = "#FFC100"
)

theme_cell <- function(base_size = 10) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", linewidth = 0.1),
      axis.line = element_line(color = "black", linewidth = 0.1),
      axis.ticks = element_line(color = "black", linewidth = 0.1),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      plot.title = element_blank(),
      legend.title = element_blank(),
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(color = "black", face = "bold")
    )
}

# =========================================================
# 2. Helper functions
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
      "File not found in: ", dir,
      "\nCandidate patterns:\n",
      paste(patterns, collapse = "\n")
    )
  }

  NULL
}

pick_col <- function(df, candidates, required = TRUE) {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) == 0) {
    if (required) {
      stop("Column not found. Candidate columns: ", paste(candidates, collapse = ", "))
    }
    return(NULL)
  }
  hit[1]
}

check_required_cols <- function(df, cols, label = "data") {
  miss <- setdiff(cols, colnames(df))
  if (length(miss) > 0) {
    stop(label, " is missing required columns: ", paste(miss, collapse = ", "))
  }
}

standardize_binary_label <- function(x) {
  x0 <- trimws(tolower(as.character(x)))
  num <- suppressWarnings(as.numeric(x0))

  out <- rep(NA_real_, length(x0))
  out[!is.na(num) & num == 1] <- 1
  out[!is.na(num) & num == 0] <- 0

  out[x0 %in% c("healthy", "health", "control", "ctrl", "normal")] <- 1
  out[x0 %in% c("disease", "case", "patient", "unhealthy", "oral disease", "nonhealthy", "non-healthy")] <- 0

  out
}

standardize_set_label <- function(x) {
  x0 <- str_squish(tolower(as.character(x)))
  dplyr::case_when(
    x0 %in% c("training", "train") ~ "Training",
    x0 %in% c("external validation", "external_validation", "external", "test", "validation") ~ "External validation",
    TRUE ~ as.character(x)
  )
}

recode_group_or_system <- function(x) {
  x <- str_squish(as.character(x))
  x <- ifelse(tolower(x) %in% c("control", "healthy", "health"), "Healthy", x)
  x <- ifelse(x == "BMI Group", "Abnormal BMI", x)
  x
}

get_health01 <- function(df, set_name) {
  if ("y_healthy" %in% colnames(df)) return(standardize_binary_label(df$y_healthy))
  if ("health01" %in% colnames(df)) return(standardize_binary_label(df$health01))
  if ("disease" %in% colnames(df)) return(standardize_binary_label(df$disease))
  if ("label" %in% colnames(df)) return(standardize_binary_label(df$label))

  stop(set_name, ": no label column found. Expected one of y_healthy, health01, disease, or label.")
}

read_score_one <- function(file, set_name) {
  x <- fread(file, data.table = FALSE, check.names = FALSE)
  x$.source_file <- basename(file)
  x$set <- set_name

  check_required_cols(x, c("prob_healthy", "OMHI"), label = set_name)

  x$prob_healthy <- as.numeric(x$prob_healthy)
  x$OMHI <- as.numeric(x$OMHI)
  x$health01 <- get_health01(x, set_name)

  x %>%
    mutate(
      set = standardize_set_label(set),
      Disease = factor(
        ifelse(health01 == 1, "Healthy", "Disease"),
        levels = c("Healthy", "Disease")
      )
    ) %>%
    filter(
      is.finite(prob_healthy),
      is.finite(OMHI),
      !is.na(health01),
      !is.na(Disease)
    )
}

make_fill_map <- function(levels_vec, color_map, default = "#BDBDBD") {
  levels_vec <- as.character(levels_vec)
  out <- unname(color_map[levels_vec])
  out[is.na(out)] <- default
  setNames(out, levels_vec)
}

save_plot <- function(plot, prefix, width, height) {
  ggsave(
    filename = file.path(OUT_DIR, paste0(prefix, ".png")),
    plot = plot,
    width = width,
    height = height,
    dpi = 600,
    bg = "white"
  )

  ggsave(
    filename = file.path(OUT_DIR, paste0(prefix, ".pdf")),
    plot = plot,
    width = width,
    height = height,
    bg = "white"
  )
}

# =========================================================
# 3. Read OMHI score files
# =========================================================
best_file <- file.path(ROOT_DIR, "LASSO_best_model_by_level_full.csv")

if (file.exists(best_file)) {
  best_tab <- fread(best_file, data.table = FALSE)
  best_species <- best_tab %>% filter(level == LEVEL)

  if (nrow(best_species) != 1) {
    stop("The best-model table should contain exactly one row for level = ", LEVEL, ": ", best_file)
  }

  alpha_val <- best_species$alpha[1]
  prev_val  <- best_species$prev_cutoff[1]
} else {
  warning("Best-model table not found. Defaulting to alpha = 0.5.")
  alpha_val <- 0.5
  prev_val <- NA_real_
}

atag <- alpha_tag(alpha_val)

cat("\nBest species model:\n")
cat("alpha =", alpha_val, "\n")
cat("prev_cutoff =", prev_val, "\n")
cat("alpha tag =", atag, "\n")

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

cat("Training file:", basename(train_file), "\n")
cat("External validation file:", basename(test_file), "\n")

train_df <- read_score_one(train_file, "Training")
test_df  <- read_score_one(test_file, "External validation")

score_df_all <- bind_rows(train_df, test_df) %>%
  mutate(
    set = factor(
      standardize_set_label(set),
      levels = c("Training", "External validation")
    )
  )

cat("\nSample counts by set and disease status:\n")
print(table(score_df_all$set, score_df_all$Disease, useNA = "ifany"))

# =========================================================
# 4. Prepare group and system labels
# =========================================================
group_col <- pick_col(score_df_all, c("group", "Group", "group_label"), required = FALSE)
system_col <- pick_col(score_df_all, c("system", "System", "system_label"), required = FALSE)

if (is.null(group_col) || is.null(system_col)) {
  stop("Both group and system columns are required for Figure 3.")
}

plot_df <- score_df_all %>%
  filter(is.finite(OMHI)) %>%
  mutate(
    group2 = recode_group_or_system(.data[[group_col]]),
    system2 = recode_group_or_system(.data[[system_col]]),
    set = factor(standardize_set_label(set), levels = c("Training", "External validation"))
  ) %>%
  filter(
    !is.na(group2), group2 != "",
    !is.na(system2), system2 != ""
  )

# =========================================================
# 5. Ordering and statistics
# =========================================================
make_orders_for_one_set <- function(dat_set) {
  system_rank_df <- dat_set %>%
    group_by(system2) %>%
    summarise(med_system = median(OMHI, na.rm = TRUE), .groups = "drop")

  other_systems <- system_rank_df %>%
    filter(system2 != "Healthy") %>%
    arrange(desc(med_system)) %>%
    pull(system2)

  system_levels <- c("Healthy", other_systems) %>% unique()

  group_system_df <- dat_set %>%
    count(group2, system2, name = "n") %>%
    group_by(group2) %>%
    slice_max(order_by = n, n = 1, with_ties = FALSE) %>%
    ungroup()

  group_rank_df <- dat_set %>%
    group_by(group2) %>%
    summarise(med_group = median(OMHI, na.rm = TRUE), .groups = "drop") %>%
    left_join(group_system_df, by = "group2")

  group_levels <- "Healthy"
  for (ss in setdiff(system_levels, "Healthy")) {
    tmp_groups <- group_rank_df %>%
      filter(system2 == ss, group2 != "Healthy") %>%
      arrange(desc(med_group)) %>%
      pull(group2)

    group_levels <- c(group_levels, tmp_groups)
  }

  list(
    system_levels = unique(system_levels),
    group_levels = unique(group_levels),
    group_system_df = group_system_df
  )
}

make_p_table <- function(dat, xvar, ref_level = "Healthy", min_n = 10, star_by = "FDR") {
  empty_out <- tibble(
    label = character(),
    p_value = numeric(),
    FDR = numeric(),
    p_signif = character()
  )

  lvls <- unique(as.character(dat[[xvar]]))
  lvls <- lvls[!is.na(lvls) & lvls != ref_level]

  ref_dat <- dat %>% filter(.data[[xvar]] == ref_level)
  if (nrow(ref_dat) < min_n) return(empty_out)

  out <- lapply(lvls, function(xx) {
    sub_dat <- dat %>% filter(.data[[xvar]] == xx)
    if (nrow(sub_dat) < min_n) return(NULL)

    pval <- tryCatch(
      wilcox.test(ref_dat$OMHI, sub_dat$OMHI, exact = FALSE)$p.value,
      error = function(e) NA_real_
    )

    tibble(label = xx, p_value = pval)
  }) %>%
    bind_rows()

  if (nrow(out) == 0) return(empty_out)

  out %>%
    mutate(
      FDR = p.adjust(p_value, method = "fdr"),
      p_for_star = if (star_by == "FDR") FDR else p_value,
      p_signif = case_when(
        is.na(p_for_star) ~ "",
        p_for_star < 0.001 ~ "***",
        p_for_star < 0.01 ~ "**",
        p_for_star < 0.05 ~ "*",
        TRUE ~ ""
      )
    ) %>%
    select(label, p_value, FDR, p_signif)
}

# =========================================================
# 6. Plot one dataset
# =========================================================
plot_one_set_group_system <- function(dat_all, set_name, out_tag,
                                      x_limits = c(-4, 4.5), min_n = 10) {
  dat <- dat_all %>%
    filter(set == set_name, is.finite(OMHI))

  if (nrow(dat) == 0) {
    warning(set_name, ": no available data.")
    return(NULL)
  }

  ord <- make_orders_for_one_set(dat)

  # System-level data and order.
  system_df <- dat %>%
    group_by(system2) %>%
    filter(n() >= min_n) %>%
    ungroup()

  keep_systems <- unique(as.character(system_df$system2))
  system_levels <- c(
    "Healthy",
    ord$system_levels[ord$system_levels %in% keep_systems & ord$system_levels != "Healthy"]
  ) %>%
    unique()

  system_df <- system_df %>%
    mutate(system2 = factor(system2, levels = rev(system_levels))) %>%
    filter(!is.na(system2))

  # Group-level data and order.
  group_df <- dat %>%
    left_join(ord$group_system_df, by = "group2", suffix = c("", "_map")) %>%
    mutate(system_map = system2_map) %>%
    filter(system_map %in% keep_systems) %>%
    group_by(group2) %>%
    filter(n() >= min_n) %>%
    ungroup()

  keep_groups <- unique(as.character(group_df$group2))
  group_levels <- c(
    "Healthy",
    ord$group_levels[ord$group_levels %in% keep_groups & ord$group_levels != "Healthy"]
  ) %>%
    unique()

  group_df <- group_df %>%
    mutate(group2 = factor(group2, levels = rev(group_levels))) %>%
    filter(!is.na(group2))

  group_fill_map <- ord$group_system_df %>%
    filter(group2 %in% keep_groups, system2 %in% keep_systems) %>%
    mutate(fill_col = make_fill_map(system2, system_color_map)[system2]) %>%
    select(group2, fill_col) %>%
    tibble::deframe()

  missing_groups <- setdiff(levels(group_df$group2), names(group_fill_map))
  if (length(missing_groups) > 0) {
    group_fill_map <- c(
      group_fill_map,
      setNames(rep("#BDBDBD", length(missing_groups)), missing_groups)
    )
  }

  sys_fill_map <- make_fill_map(levels(system_df$system2), system_color_map)

  # Significance labels. Tests are two-sided Wilcoxon rank-sum tests,
  # with Benjamini-Hochberg FDR correction used for stars.
  star_x <- x_limits[2] - 0.12

  p_group_tbl <- make_p_table(
    group_df %>% mutate(group2_chr = as.character(group2)),
    xvar = "group2_chr",
    ref_level = "Healthy",
    min_n = min_n,
    star_by = "FDR"
  )

  p_system_tbl <- make_p_table(
    system_df %>% mutate(system2_chr = as.character(system2)),
    xvar = "system2_chr",
    ref_level = "Healthy",
    min_n = min_n,
    star_by = "FDR"
  )

  group_anno <- group_df %>%
    mutate(group2_chr = as.character(group2)) %>%
    distinct(group2_chr) %>%
    left_join(p_group_tbl, by = c("group2_chr" = "label")) %>%
    mutate(
      x = star_x,
      y = factor(group2_chr, levels = rev(group_levels)),
      p_signif = ifelse(is.na(p_signif) | group2_chr == "Healthy", "", p_signif)
    )

  system_anno <- system_df %>%
    mutate(system2_chr = as.character(system2)) %>%
    distinct(system2_chr) %>%
    left_join(p_system_tbl, by = c("system2_chr" = "label")) %>%
    mutate(
      x = star_x,
      y = factor(system2_chr, levels = rev(system_levels)),
      p_signif = ifelse(is.na(p_signif) | system2_chr == "Healthy", "", p_signif)
    )

  p_group <- ggplot(group_df, aes(y = group2, x = OMHI, fill = group2)) +
    geom_violin(scale = "width", trim = FALSE, color = NA, alpha = 0.78) +
    geom_boxplot(
      width = 0.16,
      outlier.shape = NA,
      fill = NA,
      color = "black",
      linewidth = 0.45
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = COL_GREY, linewidth = 0.7) +
    geom_text(
      data = group_anno %>% filter(p_signif != ""),
      aes(x = x, y = y, label = p_signif),
      inherit.aes = FALSE,
      size = 5,
      hjust = 1
    ) +
    scale_fill_manual(values = group_fill_map) +
    coord_cartesian(xlim = x_limits, clip = "off") +
    labs(x = "OMHI", y = NULL) +
    theme_cell(10) +
    theme(
      legend.position = "none",
      axis.text.y = element_text(size = 12),
      plot.margin = margin(5.5, 5.5, 5.5, 5.5)
    )

  p_system <- ggplot(system_df, aes(y = system2, x = OMHI, fill = system2)) +
    geom_violin(scale = "width", trim = FALSE, color = NA, alpha = 0.78) +
    geom_boxplot(
      width = 0.16,
      outlier.shape = NA,
      fill = NA,
      color = "black",
      linewidth = 0.45
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = COL_GREY, linewidth = 0.7) +
    geom_text(
      data = system_anno %>% filter(p_signif != ""),
      aes(x = x, y = y, label = p_signif),
      inherit.aes = FALSE,
      size = 5,
      hjust = 1
    ) +
    scale_fill_manual(values = sys_fill_map) +
    coord_cartesian(xlim = x_limits, clip = "off") +
    labs(x = "OMHI", y = NULL) +
    theme_cell(10) +
    theme(
      legend.position = "none",
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.margin = margin(5.5, 0, 5.5, 5.5)
    )

  sys_label_df <- tibble(
    system2 = factor(rev(system_levels), levels = rev(system_levels))
  ) %>%
    mutate(
      system_label = case_when(
        system2 == "Metabolic Disease" ~ "Metabolic\nDisease",
        system2 == "Digestive Disease" ~ "Digestive\nDisease",
        system2 == "Respiratory Disease" ~ "Respiratory\nDisease",
        system2 == "Immune Disease" ~ "Immune\nDisease",
        system2 == "Oral Disease" ~ "Oral\nDisease",
        system2 == "Abnormal BMI" ~ "Abnormal\nBMI",
        TRUE ~ as.character(system2)
      ),
      x = 1
    )

  p_sys_label <- ggplot(sys_label_df, aes(x = x, y = system2)) +
    geom_tile(aes(fill = system2), width = 2.2, height = 0.9, color = NA) +
    geom_text(aes(label = system_label), size = 4, fontface = "bold", lineheight = 0.9) +
    scale_fill_manual(values = sys_fill_map) +
    coord_cartesian(xlim = c(-0.15, 2.15), clip = "off") +
    labs(x = NULL, y = NULL) +
    theme_void() +
    theme(
      legend.position = "none",
      plot.margin = margin(5.5, 5.5, 5.5, 0)
    )

  p_out <- p_group + p_system + p_sys_label +
    plot_layout(widths = c(1.05, 1.05, 0.55))

  list(
    plot = p_out,
    group_plot = p_group,
    system_plot = p_system,
    system_label_plot = p_sys_label,
    group_stats = p_group_tbl,
    system_stats = p_system_tbl
  )
}

# =========================================================
# 7. Generate Figure 3
# =========================================================
res_A <- plot_one_set_group_system(
  dat_all = plot_df,
  set_name = "Training",
  out_tag = paste0(FIGURE_PREFIX, "A_OMHI_by_group_system_Training"),
  x_limits = c(-4, 4.5),
  min_n = 10
)

res_B <- plot_one_set_group_system(
  dat_all = plot_df,
  set_name = "External validation",
  out_tag = paste0(FIGURE_PREFIX, "B_OMHI_by_group_system_ExternalValidation"),
  x_limits = c(-4, 4.5),
  min_n = 10
)

if (is.null(res_A) || is.null(res_B)) {
  stop("Figure 3 could not be generated because one of the datasets is empty.")
}

save_plot(
  res_A$plot,
  paste0(FIGURE_PREFIX, "A_OMHI_by_group_system_Training"),
  width = PANEL_WIDTH,
  height = PANEL_HEIGHT
)

save_plot(
  res_B$plot,
  paste0(FIGURE_PREFIX, "B_OMHI_by_group_system_ExternalValidation"),
  width = PANEL_WIDTH,
  height = PANEL_HEIGHT
)

p_AB <- res_A$plot / res_B$plot +
  plot_layout(heights = c(1.4, 0.6))

save_plot(
  p_AB,
  paste0(FIGURE_PREFIX, "AB_OMHI_by_group_system_vertical"),
  width = COMBINED_WIDTH,
  height = COMBINED_HEIGHT
)

fwrite(res_A$group_stats, file.path(OUT_DIR, paste0(FIGURE_PREFIX, "A_group_wilcoxon_fdr.csv")))
fwrite(res_A$system_stats, file.path(OUT_DIR, paste0(FIGURE_PREFIX, "A_system_wilcoxon_fdr.csv")))
fwrite(res_B$group_stats, file.path(OUT_DIR, paste0(FIGURE_PREFIX, "B_group_wilcoxon_fdr.csv")))
fwrite(res_B$system_stats, file.path(OUT_DIR, paste0(FIGURE_PREFIX, "B_system_wilcoxon_fdr.csv")))

cat("\nFigure 3 outputs saved to:\n")
cat(OUT_DIR, "\n")
