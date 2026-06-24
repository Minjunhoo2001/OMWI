suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tidyr)
  library(scales)
  library(patchwork)
})

# =========================================================
# 0. File paths
# =========================================================
loocv_file <- "/share/home/HeMinjun/metagenomic/Oral_proxy_v2/04.lasso64/loocv_oof/OOF_LOOCV_species_alpha050_best_model.csv"
cv_file    <- "/share/home/HeMinjun/metagenomic/Oral_proxy_v2/04.lasso64/cv_oof/OOF_10fold_repeated_species_alpha050_best_model.csv"
train_file <- "/share/home/HeMinjun/metagenomic/Oral_proxy_v2/04.lasso64/omhi_scores/OMHI_train_species_alpha050.csv"
validation_file <- "/share/home/HeMinjun/metagenomic/Oral_proxy_v2/04.lasso64/omhi_scores/OMHI_test_species_alpha050.csv"

out_dir <- "/share/home/HeMinjun/metagenomic/GitHub/results/Figure2"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# =========================================================
# 1. Helper functions
# =========================================================
read_and_standardize <- function(file, dataset_name) {
  read_csv(file, show_col_types = FALSE) %>%
    transmute(
      sample_id = as.character(sample_id),
      y_healthy = as.numeric(y_healthy),
      OMHI = as.numeric(OMHI),
      prob_healthy = suppressWarnings(as.numeric(prob_healthy)),
      dataset = dataset_name
    ) %>%
    filter(!is.na(sample_id), !is.na(y_healthy), !is.na(OMHI))
}

collapse_duplicate_samples <- function(df) {
  if (anyDuplicated(df$sample_id) == 0) {
    return(df)
  }

  df %>%
    group_by(sample_id) %>%
    summarise(
      y_healthy = first(y_healthy),
      OMHI = mean(OMHI, na.rm = TRUE),
      prob_healthy = mean(prob_healthy, na.rm = TRUE),
      dataset = first(dataset),
      .groups = "drop"
    )
}

calc_cutoff_metrics <- function(df, cutoff, dataset_name) {
  sub <- df %>%
    filter(abs(OMHI) >= cutoff)

  if (nrow(sub) == 0) {
    return(tibble(
      dataset = dataset_name,
      cutoff = cutoff,
      balanced_acc = NA_real_
    ))
  }

  sub <- sub %>%
    mutate(pred_healthy = ifelse(OMHI > 0, 1, 0))

  healthy_dat <- sub %>% filter(y_healthy == 1)
  disease_dat <- sub %>% filter(y_healthy == 0)

  healthy_acc <- ifelse(
    nrow(healthy_dat) > 0,
    100 * mean(healthy_dat$pred_healthy == 1),
    NA_real_
  )

  nonhealthy_acc <- ifelse(
    nrow(disease_dat) > 0,
    100 * mean(disease_dat$pred_healthy == 0),
    NA_real_
  )

  balanced_acc <- ifelse(
    is.na(healthy_acc) | is.na(nonhealthy_acc),
    NA_real_,
    (healthy_acc + nonhealthy_acc) / 2
  )

  tibble(
    dataset = dataset_name,
    cutoff = cutoff,
    balanced_acc = balanced_acc
  )
}

make_panel_a_data <- function(df, breaks, labels) {
  df %>%
    mutate(
      dataset = "Overall",
      bin = cut(
        OMHI,
        breaks = breaks,
        labels = labels,
        right = FALSE,
        include.lowest = TRUE
      ),
      health_status = ifelse(y_healthy == 1, "Healthy", "Disease")
    ) %>%
    count(dataset, bin, health_status, .drop = FALSE) %>%
    mutate(
      dataset = factor(dataset, levels = "Overall"),
      bin = factor(bin, levels = labels),
      health_status = factor(health_status, levels = c("Disease", "Healthy"))
    ) %>%
    group_by(dataset, bin) %>%
    mutate(
      total_n = sum(n),
      prop = ifelse(total_n > 0, n / total_n, 0)
    ) %>%
    ungroup()
}

# =========================================================
# 2. Read and prepare data
# =========================================================
train_df <- read_and_standardize(train_file, "Training")
validation_df <- read_and_standardize(validation_file, "Validation") %>%
  collapse_duplicate_samples()
loocv_df <- read_and_standardize(loocv_file, "LOOCV")
cv_df <- read_and_standardize(cv_file, "10-fold CV") %>%
  collapse_duplicate_samples()

overall_df <- bind_rows(train_df, validation_df) %>%
  mutate(dataset = "Overall")

# =========================================================
# 3. Calculate cutoff-based metrics
# =========================================================
cutoffs <- c(0, 0.5, 1)

metrics_df <- bind_rows(
  lapply(cutoffs, function(x) calc_cutoff_metrics(train_df, x, "Training")),
  lapply(cutoffs, function(x) calc_cutoff_metrics(loocv_df, x, "LOOCV")),
  lapply(cutoffs, function(x) calc_cutoff_metrics(cv_df, x, "10-fold CV"))
) %>%
  mutate(dataset = factor(dataset, levels = c("Training", "LOOCV", "10-fold CV")))

# =========================================================
# 4. Panel A: sample proportion by OMHI bin
# =========================================================
breaks_train <- c(-Inf, -2, -1.5, -1, -0.5, 0, 0.5, 1, Inf)
labels_train <- c(
  "< -2", "[-2, -1.5)", "[-1.5, -1)", "[-1, -0.5)",
  "[-0.5, 0)", "[0, 0.5)", "[0.5, 1)", ">= 1"
)

panelA_train <- make_panel_a_data(
  df = overall_df,
  breaks = breaks_train,
  labels = labels_train
)

p_A <- ggplot(panelA_train, aes(x = bin, y = prop, fill = health_status)) +
  geom_col(position = "stack", width = 0.84, color = "#4D4D4D", linewidth = 0.25) +
  geom_hline(yintercept = 0.5, linetype = "dotted", color = "#999999", linewidth = 0.6) +
  geom_vline(xintercept = 5.5, linetype = "dashed", color = "#666666", linewidth = 0.8) +
  scale_fill_manual(
    values = c("Disease" = "#E64B35", "Healthy" = "#4DBBD5"),
    breaks = c("Healthy", "Disease")
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1.08),
    breaks = seq(0, 1, 0.2),
    expand = c(0, 0)
  ) +
  labs(
    title = "Overall",
    x = NULL,
    y = "Sample proportion by OMHI bin",
    fill = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "#DDDDDD", linewidth = 0.4),
    panel.border = element_rect(color = "black", linewidth = 0.1),
      axis.line = element_line(color = "black", linewidth = 0.1),
      axis.ticks = element_line(color = "black", linewidth = 0.1),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black"),
    axis.text.y = element_text(color = "black", size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )

print(p_A)

# =========================================================
# 5. Panel B: balanced accuracy by OMHI magnitude cutoff
# =========================================================
p_B <- ggplot(metrics_df, aes(x = cutoff, y = balanced_acc, group = 1)) +
  geom_line(color = "#00A087", linewidth = 1.0) +
  geom_point(color = "#2CA25F", size = 2.6) +
  geom_text(
    aes(label = ifelse(is.na(balanced_acc), "", sprintf("%.1f%%", balanced_acc))),
    vjust = -0.9, color = "black", size = 3.5
  ) +
  facet_wrap(~dataset, nrow = 1) +
  scale_x_continuous(
    breaks = cutoffs,
    expand = expansion(mult = c(0.25, 0.25))
  ) +
  scale_y_continuous(
    name = "Balanced accuracy (%)",
    limits = c(68, 92),
    breaks = seq(70, 90, 5),
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  labs(
    x = expression("|OMHI| cutoff"),
    title = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "#DDDDDD", linewidth = 0.4),
    panel.border = element_rect(color = "black", linewidth = 0.2),
      axis.line = element_line(color = "black", linewidth = 0.2),
      axis.ticks = element_line(color = "black", linewidth = 0.2),
    strip.background = element_rect(fill = "#F0F0F0", color = "#666666"),
    strip.text = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  )

print(p_B)

# =========================================================
# 6. Combine panels and save outputs
# =========================================================
fig3_final <- (p_B | p_A) + plot_layout(widths = c(1.5, 1.0))
print(fig3_final)

ggsave(
  file.path(out_dir, "Figure2DE.png"),
  fig3_final, width = 14.5, height = 4, dpi = 320, bg = "white"
)

ggsave(
  file.path(out_dir, "Figure2DE.pdf"),
  fig3_final, width = 14.5, height = 4, bg = "white"
)

write_csv(panelA_train, file.path(out_dir, "Figure2E_overall_bins.csv"))
write_csv(metrics_df, file.path(out_dir, "Figure2D_balanced_accuracy.csv"))

message("Done. Files saved in: ", out_dir)
