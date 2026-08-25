#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ragg)
  library(stringr)
})

options(stringsAsFactors = FALSE)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else "scripts/Figure3/Figure3F_system_AUC.R"
repo_root <- normalizePath(file.path(dirname(script_file), "..", ".."), mustWork = FALSE)

get_arg <- function(flag, default) {
  args <- commandArgs(trailingOnly = TRUE)
  hit <- match(flag, args)
  if (!is.na(hit) && hit < length(args)) args[[hit + 1]] else default
}

input_file <- get_arg(
  "--input",
  file.path(repo_root, "results", "Figure3", "Fig3F_system_AUC", "Fig3F_system_AUC_source_data.tsv")
)
out_dir <- get_arg(
  "--out-dir",
  file.path(repo_root, "results", "Figure3", "Fig3F_system_AUC")
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

plot_df <- read.delim(input_file, check.names = FALSE) %>%
  mutate(
    set = case_when(
      str_detect(str_to_lower(set), "train") ~ "Training",
      str_detect(str_to_lower(set), "external|test") ~ "External validation",
      TRUE ~ as.character(set)
    ),
    set = factor(set, levels = c("Training", "External validation")),
    auc = as.numeric(auc),
    ci_low = as.numeric(ci_low),
    ci_high = as.numeric(ci_high)
  ) %>%
  filter(!is.na(auc), !is.na(ci_low), !is.na(ci_high)) %>%
  group_by(set) %>%
  arrange(desc(auc), .by_group = TRUE) %>%
  mutate(rank_in_set = row_number()) %>%
  ungroup() %>%
  mutate(
    system_plot = paste(set, system, sep = "___"),
    auc_ci_label = paste0(sprintf("%.2f", auc), "\n(", sprintf("%.2f", ci_low), "-", sprintf("%.2f", ci_high), ")"),
    tile_y = 0.512
  )

level_order <- plot_df %>% arrange(set, rank_in_set) %>% pull(system_plot)
plot_df$system_plot <- factor(plot_df$system_plot, levels = level_order)
x_labels <- setNames(plot_df$system, plot_df$system_plot)

system_colours <- c(
  "Digestive disease" = "#D97C8A",
  "Oral disease" = "#C7A6D8",
  "Immune disease" = "#D8C97A",
  "Respiratory disease" = "#8FA8C9",
  "Metabolic disease" = "#4198AC"
)
missing_colours <- setdiff(unique(plot_df$system), names(system_colours))
if (length(missing_colours)) stop("Missing system colours for: ", paste(missing_colours, collapse = ", "))

p_auc <- ggplot(plot_df, aes(system_plot, auc)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey60", linewidth = 0.45) +
  geom_tile(aes(y = tile_y, fill = system), width = 0.72, height = 0.018, colour = NA) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high, colour = system), width = 0.10, linewidth = 0.70) +
  geom_point(aes(colour = system), size = 2.6) +
  geom_text(aes(label = auc_ci_label), vjust = -0.55, colour = "black", size = 2.35, lineheight = 0.88) +
  facet_wrap(~set, nrow = 1, scales = "free_x") +
  scale_colour_manual(values = system_colours, drop = FALSE) +
  scale_fill_manual(values = system_colours, drop = FALSE) +
  scale_x_discrete(labels = x_labels, drop = TRUE) +
  scale_y_continuous(
    limits = c(0.50, 1.00),
    breaks = seq(0.5, 1.0, by = 0.1),
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(x = NULL, y = "AUC", tag = "f") +
  theme_bw(base_size = 8, base_family = "Arial") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", linewidth = 0.45),
    axis.line = element_line(colour = "black", linewidth = 0.35),
    axis.ticks = element_line(colour = "black", linewidth = 0.35),
    axis.text = element_text(colour = "black"),
    axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1, size = 6.8),
    strip.background = element_blank(),
    strip.text = element_text(colour = "black", face = "bold", size = 8.5),
    legend.position = "none",
    plot.tag = element_text(face = "bold", size = 9),
    plot.margin = margin(5, 4, 4, 4)
  )

write.table(
  plot_df %>% select(set, system, auc, ci_low, ci_high, rank_in_set),
  file = file.path(out_dir, "Fig3F_system_AUC_plot_data.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

pdf_file <- file.path(out_dir, "Fig3F_system_AUC.pdf")
png_file <- file.path(out_dir, "Fig3F_system_AUC.png")
ggsave(pdf_file, p_auc, width = 7.2, height = 3.9, device = grDevices::cairo_pdf, bg = "white")
ggsave(png_file, p_auc, width = 7.2, height = 3.9, dpi = 600, device = ragg::agg_png, bg = "white")
capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))

cat("Saved Figure 3f AUC panel to:\n", pdf_file, "\n", png_file, "\n", sep = "")
