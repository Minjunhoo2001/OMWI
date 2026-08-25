#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(ragg)
  library(scales)
  library(tidyr)
})

options(stringsAsFactors = FALSE)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else "scripts/Figure2/Figure2F_LODO_stability.R"
repo_root <- normalizePath(file.path(dirname(script_file), "..", ".."), mustWork = FALSE)

get_arg <- function(flag, default) {
  args <- commandArgs(trailingOnly = TRUE)
  hit <- match(flag, args)
  if (!is.na(hit) && hit < length(args)) args[[hit + 1]] else default
}

input_file <- get_arg(
  "--input",
  file.path(repo_root, "results", "Figure2", "Fig2F_LODO_stability", "LODO_stability_heatmap_data.csv")
)
out_dir <- get_arg(
  "--out-dir",
  file.path(repo_root, "results", "Figure2", "Fig2F_LODO_stability")
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

dat <- read.csv(input_file, check.names = FALSE)
required <- c("model_label", "training_cv_auc", "pooled_external_auc", "spearman_vs_final", "model_type")
missing <- setdiff(required, names(dat))
if (length(missing)) stop("Missing columns: ", paste(missing, collapse = ", "))

dat$model_label <- factor(dat$model_label, levels = dat$model_label)

auc_df <- dat %>%
  select(model_label, model_type, training_cv_auc, pooled_external_auc) %>%
  pivot_longer(
    cols = c(training_cv_auc, pooled_external_auc),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = recode(
      metric,
      training_cv_auc = "Training CV AUC",
      pooled_external_auc = "Pooled external AUC"
    ),
    metric = factor(metric, levels = c("Pooled external AUC", "Training CV AUC")),
    value_label = sprintf("%.2f", value),
    text_colour = ifelse(value >= 0.76, "white", "black")
  )

rho_df <- dat %>%
  transmute(
    model_label,
    model_type,
    metric = factor("External-score\nSpearman r", levels = "External-score\nSpearman r"),
    value = spearman_vs_final,
    value_label = sprintf("%.3f", value),
    text_colour = ifelse(value >= 0.94, "white", "black")
  )

final_auc <- filter(auc_df, model_type == "Final model")
final_rho <- filter(rho_df, model_type == "Final model")

base_theme <- theme_classic(base_size = 8, base_family = "Arial") +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_text(colour = "black", face = "bold", size = 7.2),
    legend.title = element_text(face = "bold", size = 7.5),
    legend.text = element_text(size = 7),
    plot.margin = margin(2, 4, 1, 4)
  )

p_auc <- ggplot(auc_df, aes(model_label, metric, fill = value)) +
  geom_tile(colour = "white", linewidth = 0.65) +
  geom_tile(data = final_auc, fill = NA, colour = "black", linewidth = 0.85) +
  geom_text(aes(label = value_label, colour = text_colour), size = 2.25, fontface = "bold", show.legend = FALSE) +
  scale_colour_identity() +
  scale_fill_gradientn(
    colours = c("#FBEFEF", "#F5B7B1", "#EC7063", "#CB4335", "#7B241C"),
    values = rescale(c(0.60, 0.65, 0.70, 0.75, 0.84)),
    limits = c(0.60, 0.84),
    oob = squish,
    name = "AUC"
  ) +
  scale_x_discrete(drop = FALSE, position = "top") +
  labs(x = NULL, y = NULL, tag = "f") +
  base_theme +
  theme(
    axis.text.x = element_blank(),
    plot.tag = element_text(face = "bold", size = 9)
  )

p_rho <- ggplot(rho_df, aes(model_label, metric, fill = value)) +
  geom_tile(colour = "white", linewidth = 0.65) +
  geom_tile(data = final_rho, fill = NA, colour = "black", linewidth = 0.85) +
  geom_text(aes(label = value_label, colour = text_colour), size = 2.15, fontface = "bold", show.legend = FALSE) +
  scale_colour_identity() +
  scale_fill_gradientn(
    colours = c("#FFF5EB", "#FDD0A2", "#FDAE6B", "#F16913", "#D94801"),
    values = rescale(c(0.80, 0.85, 0.90, 0.95, 1.00)),
    limits = c(0.80, 1.00),
    oob = squish,
    name = "Spearman r"
  ) +
  scale_x_discrete(drop = FALSE) +
  labs(x = "Excluded training cohort (healthy/disease)", y = NULL) +
  base_theme +
  theme(
    axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, colour = "black", size = 6.4),
    axis.title.x = element_text(size = 7.2, margin = margin(t = 4))
  )

fig <- p_auc / p_rho + plot_layout(heights = c(2.0, 1.35))

pdf_file <- file.path(out_dir, "Fig2F_LODO_stability_heatmap.pdf")
png_file <- file.path(out_dir, "Fig2F_LODO_stability_heatmap.png")
ggsave(pdf_file, fig, width = 13.2, height = 3.9, device = grDevices::cairo_pdf, bg = "white")
ggsave(png_file, fig, width = 13.2, height = 3.9, dpi = 600, device = ragg::agg_png, bg = "white")
capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))

cat("Saved Figure 2f to:\n", pdf_file, "\n", png_file, "\n", sep = "")
