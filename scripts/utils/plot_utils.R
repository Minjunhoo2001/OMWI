# ============================================================
# plot_utils.R
# Plotting utilities for OMHI PCoA analyses
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
  library(scales)
  library(grid)
  library(dplyr)
})

theme_paper <- function(base_size = 13) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(
        size = base_size + 1,
        face = "bold",
        hjust = 0.5
      ),
      axis.title = element_text(size = base_size, face = "bold"),
      axis.text = element_text(size = base_size - 1, color = "black"),
      legend.title = element_text(size = base_size - 1, face = "bold"),
      legend.text = element_text(size = base_size - 2),
      legend.key.size = unit(0.9, "lines"),
      strip.text = element_text(size = base_size, face = "bold"),
      plot.margin = margin(8, 8, 8, 8)
    )
}

add_stat_text <- function(p, df, label_text, size = 4.0) {
  xr <- range(df$PCoA1, na.rm = TRUE)
  yr <- range(df$PCoA2, na.rm = TRUE)

  p +
    annotate(
      "text",
      x = xr[2] - 0.02 * diff(xr),
      y = yr[2] - 0.02 * diff(yr),
      label = label_text,
      hjust = 1,
      vjust = 1,
      size = size
    )
}

make_disease_colors <- function(levels_vec) {
  base_cols <- c(
    "healthy" = "#4DBBD5",
    "disease" = "#E64B35",
    "Healthy" = "#4DBBD5",
    "Disease" = "#E64B35"
  )

  cols <- base_cols[levels_vec]

  if (any(is.na(cols))) {
    cols <- setNames(
      hue_pal()(length(levels_vec)),
      levels_vec
    )
  }

  cols
}

make_system_colors <- function(levels_vec) {
  base_cols <- c(
    "Healthy" = "#4DBBD5",
    "Oral Disease" = "#C7A6D8",
    "Digestive Disease" = "#D97C8A",
    "Respiratory Disease" = "#8FA8C9",
    "Immune Disease" = "#D8C97A",
    "Metabolic Disease" = "#4198AC",
    "BMI Group" = "#FFC100"
  )

  missing_levels <- setdiff(levels_vec, names(base_cols))

  if (length(missing_levels) > 0) {
    extra_cols <- hue_pal()(length(missing_levels))
    names(extra_cols) <- missing_levels
    base_cols <- c(base_cols, extra_cols)
  }

  base_cols[levels_vec]
}

plot_pcoa_by_group <- function(
    df,
    group_col,
    color_values,
    eig_pct,
    title,
    stat_label = NULL,
    shape_col = NULL,
    base_size = 12
) {
  aes_base <- aes(
    x = PCoA1,
    y = PCoA2,
    color = .data[[group_col]]
  )

  if (!is.null(shape_col)) {
    aes_base$shape <- rlang::expr(.data[[shape_col]])
  }

  p <- ggplot(df, aes_base) +
    geom_point(size = 2.2, alpha = 0.85) +
    stat_ellipse(
      aes(color = .data[[group_col]]),
      level = 0.95,
      linetype = 2,
      linewidth = 0.55,
      show.legend = FALSE
    ) +
    scale_color_manual(values = color_values) +
    labs(
      title = title,
      x = paste0("PCoA1 (", sprintf("%.1f", eig_pct[1]), "%)"),
      y = paste0("PCoA2 (", sprintf("%.1f", eig_pct[2]), "%)")
    ) +
    guides(
      color = guide_legend(
        override.aes = list(size = 3.5, alpha = 1)
      )
    ) +
    theme_paper(base_size = base_size)

  if (!is.null(shape_col)) {
    n_shape <- length(unique(df[[shape_col]]))
    p <- p + scale_shape_manual(values = rep(0:25, length.out = n_shape))
  }

  if (!is.null(stat_label)) {
    p <- add_stat_text(p, df, stat_label, size = 3.8)
  }

  p
}

plot_pcoa_centroid <- function(
    df,
    center_df,
    group_col,
    color_values,
    eig_pct,
    title,
    stat_label = NULL,
    base_size = 12
) {
  seg_df <- df %>%
    left_join(center_df, by = group_col)

  p <- ggplot() +
    geom_segment(
      data = seg_df,
      aes(
        x = cx,
        y = cy,
        xend = PCoA1,
        yend = PCoA2,
        color = .data[[group_col]]
      ),
      alpha = 0.45,
      linewidth = 0.25,
      show.legend = FALSE
    ) +
    geom_point(
      data = df,
      aes(
        x = PCoA1,
        y = PCoA2,
        color = .data[[group_col]]
      ),
      size = 1.6,
      alpha = 0.18
    ) +
    geom_point(
      data = center_df,
      aes(
        x = cx,
        y = cy,
        color = .data[[group_col]]
      ),
      size = 4.3,
      alpha = 1
    ) +
    scale_color_manual(values = color_values) +
    labs(
      title = title,
      x = paste0("PCoA1 (", sprintf("%.1f", eig_pct[1]), "%)"),
      y = paste0("PCoA2 (", sprintf("%.1f", eig_pct[2]), "%)")
    ) +
    guides(
      color = guide_legend(
        ncol = 1,
        override.aes = list(size = 3.5, alpha = 1)
      )
    ) +
    theme_paper(base_size = base_size) +
    theme(legend.position = "right")

  if (!is.null(stat_label)) {
    p <- add_stat_text(p, df, stat_label, size = 3.8)
  }

  p
}

save_pdf_png <- function(plot, out_prefix, width, height) {
  ggsave(
    filename = paste0(out_prefix, ".pdf"),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    useDingbats = FALSE
  )

  ggsave(
    filename = paste0(out_prefix, ".png"),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 300
  )
}