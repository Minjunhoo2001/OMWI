# ============================================================================
# Supplementary Fig. X | OMWI transfer to dental samples (Nature-style figure)
# a | Dental ROC: all n=480 vs unseen-study n=419
# b | Paired saliva vs dental OMWI (n=129)
# Frozen predictions only: no retraining, no reoptimizing the cutoff.
# Compatibility fix: no patchwork '& theme()'; no deprecated label.size.
# ============================================================================

for (p in c("ggplot2", "pROC")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
library(ggplot2)

# ---- 1. Input/output paths --------------------------------------------------
base_dir <- file.path(Sys.getenv("OMWI_ROOT",getwd()), "data/dental_validation")
b_file   <- file.path(base_dir, "B", "frozen_predictions_Dental_only.csv")
c_file   <- file.path(base_dir, "C", "Saliva_Dental_subject_pairs.csv")
out_dir  <- file.path(base_dir, "figure")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(b_file), file.exists(c_file))

# ---- 2. Import and verify sample coverage ----------------------------------
dental <- read.csv(b_file, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(all(c("sample_id", "y_true", "score", "study_seen_in_training",
                "panel_all_zero") %in% names(dental)))
stopifnot(!anyNA(dental[, c("sample_id", "y_true", "score")]))
stopifnot(anyDuplicated(dental$sample_id) == 0)
stopifnot(all(dental$y_true %in% c(0, 1)), all(is.finite(dental$score)))

as_true <- function(x) tolower(trimws(as.character(x))) %in% c("true", "t", "1")
dental_unseen <- dental[!as_true(dental$study_seen_in_training), , drop = FALSE]
stopifnot(nrow(dental) == 480L, nrow(dental_unseen) == 419L)

cat("All dental n:", nrow(dental), "\n")
cat("Unseen-study dental n:", nrow(dental_unseen), "\n")
cat("All-zero panel samples retained:", sum(as_true(dental$panel_all_zero)), "\n")

# ---- 3. ROC calculated from each sample's frozen OMWI score ----------------
# 0: non-healthy, 1: healthy. Higher OMWI means more healthy.
r_all <- pROC::roc(dental$y_true, dental$score,
                   levels = c(0, 1), direction = "<", quiet = TRUE)
r_unseen <- pROC::roc(dental_unseen$y_true, dental_unseen$score,
                      levels = c(0, 1), direction = "<", quiet = TRUE)
a_all <- as.numeric(pROC::auc(r_all))
a_unseen <- as.numeric(pROC::auc(r_unseen))

lab_all <- sprintf("All dental (n = %d), AUC = %.3f", nrow(dental), a_all)
lab_unseen <- sprintf("Unseen studies (n = %d), AUC = %.3f",
                      nrow(dental_unseen), a_unseen)
roc_xy <- function(r, label) {
  data.frame(fpr = 1 - r$specificities,
             tpr = r$sensitivities,
             series = label)
}
roc_df <- rbind(roc_xy(r_all, lab_all), roc_xy(r_unseen, lab_unseen))
roc_df$series <- factor(roc_df$series, levels = c(lab_all, lab_unseen))

# ---- 4. Matched paired observations ----------------------------------------
pairs <- read.csv(c_file, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(all(c("study_name", "score_saliva", "score_dental") %in% names(pairs)))
stopifnot(nrow(pairs) == 129L)
stopifnot(all(is.finite(pairs$score_saliva)), all(is.finite(pairs$score_dental)))

rho <- cor(pairs$score_saliva, pairs$score_dental, method = "spearman")
mean_diff <- mean(pairs$score_dental - pairs$score_saliva)
median_diff <- median(pairs$score_dental - pairs$score_saliva)

study_counts <- table(pairs$study_name)
expected_studies <- c("BelstromD_2021", "ZhangDY_2015")
stopifnot(setequal(names(study_counts), expected_studies))
legend_labels <- c(
  BelstromD_2021 = sprintf("Belstrom et al. (n = %d)", study_counts[["BelstromD_2021"]]),
  ZhangDY_2015 = sprintf("Zhang et al. (n = %d)", study_counts[["ZhangDY_2015"]])
)
pairs$group <- factor(pairs$study_name, levels = expected_studies)

# Same x/y limits so gray diagonal is truly y = x.
limits <- range(c(pairs$score_saliva, pairs$score_dental))
limits <- c(floor(limits[1] - .20), ceiling(limits[2] + .20))

# ---- 5. Nature-inspired, publication-friendly graphical system ------------
# Restrained colors and very little decorative ink.
col_navy <- "#365C8D"
col_orange <- "#DA684B"
col_cyan <- "#4CA9BD"
col_gray <- "#AEB4B9"
col_text <- "#24292E"

nature_theme <- theme_classic(base_family = "sans", base_size = 8.5) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "#343A40", linewidth = .38),
    axis.ticks = element_line(color = "#343A40", linewidth = .33),
    axis.ticks.length = grid::unit(1.6, "mm"),
    axis.text = element_text(size = 7.6, color = col_text),
    axis.title = element_text(size = 8.6, color = col_text),
    plot.title = element_text(size = 9, face = "bold", color = col_text,
                              hjust = 0, margin = margin(b = 6)),
    legend.position = "bottom",
    legend.direction = "vertical",
    legend.justification = "left",
    legend.title = element_blank(),
    legend.text = element_text(size = 7.1, color = col_text),
    legend.key.height = grid::unit(3.6, "mm"),
    legend.key.width = grid::unit(10.5, "mm"),
    legend.spacing.y = grid::unit(.2, "mm"),
    plot.margin = margin(7, 9, 5, 8)
  )

# ---- 6. Panel a: ROC ---------------------------------------------------------
p_a <- ggplot(roc_df, aes(fpr, tpr, color = series, group = series)) +
  geom_abline(slope = 1, intercept = 0, color = "#CACED1",
              linewidth = .4, linetype = "dashed") +
  geom_path(linewidth = .8, lineend = "butt") +
  scale_color_manual(values = setNames(c(col_navy, col_orange),
                                       c(lab_all, lab_unseen)), drop = FALSE) +
  scale_x_continuous(breaks = seq(0, 1, .25), limits = c(0, 1),
                     expand = expansion(mult = 0)) +
  scale_y_continuous(breaks = seq(0, 1, .25), limits = c(0, 1),
                     expand = expansion(mult = 0)) +
  coord_equal() +
  labs(title = "a   Transfer to dental samples",
       x = "1 - Specificity", y = "Sensitivity") +
  nature_theme

# ---- 7. Panel b: paired samples --------------------------------------------
label_b <- sprintf("n = %d paired subjects\nSpearman's rho = %.3f\nMean dental - saliva = %.2f",
                   nrow(pairs), rho, mean_diff)

p_b <- ggplot(pairs, aes(score_saliva, score_dental, color = group)) +
  geom_hline(yintercept = 0, color = "#E7E9EA", linewidth = .33,
             linetype = "dotted") +
  geom_vline(xintercept = 0, color = "#E7E9EA", linewidth = .33,
             linetype = "dotted") +
  geom_abline(slope = 1, intercept = 0, color = col_gray,
              linewidth = .55, linetype = "dashed") +
  geom_point(size = 1.2, alpha = .87, stroke = 0) +
  annotate("text", x = limits[1] + .18, y = limits[2] - .16,
           label = label_b, hjust = 0, vjust = 1,
           size = 2.45, lineheight = 1.25, color = col_text) +
  scale_color_manual(values = c(BelstromD_2021 = col_cyan,
                                ZhangDY_2015 = col_orange),
                     breaks = expected_studies,
                     labels = unname(legend_labels[expected_studies]),
                     drop = FALSE) +
  scale_x_continuous(limits = limits, breaks = seq(-4, 4, 2),
                     expand = expansion(mult = 0)) +
  scale_y_continuous(limits = limits, breaks = seq(-4, 4, 2),
                     expand = expansion(mult = 0)) +
  coord_equal() +
  labs(title = "b   Matched saliva–dental samples",
       x = "Saliva OMWI", y = "Dental OMWI") +
  nature_theme

# ---- 8. Arrange using base grid (NO patchwork dependency) -----------------
# The previous error arose from patchwork's '& theme()' operator.
# Base grid bypasses that operator completely.
draw_fig <- function() {
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(
    layout = grid::grid.layout(1, 2,
                               widths = grid::unit(c(1, 1), "null"))))
  print(p_a, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1),
        newpage = FALSE)
  print(p_b, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2),
        newpage = FALSE)
  grid::popViewport()
}

# 180 mm wide (~Nature double-column width), 88 mm tall.
fig_width <- 180 / 25.4
fig_height <- 88 / 25.4

# Preview in RStudio Plots pane.
draw_fig()

pdf_file <- file.path(out_dir, "Supplementary_Fig_X_Nature.pdf")
png_file <- file.path(out_dir, "Supplementary_Fig_X_Nature_600dpi.png")

# PDF with editable/vector lines and texts.
grDevices::pdf(pdf_file, width = fig_width, height = fig_height,
               family = "sans", useDingbats = FALSE)
draw_fig()
grDevices::dev.off()

# PNG with 600 dpi, white background.
grDevices::png(png_file, width = fig_width, height = fig_height,
               units = "in", res = 600, bg = "white")
draw_fig()
grDevices::dev.off()

# ---- 9. Exact figure metrics ------------------------------------------------
roc_metrics <- data.frame(
  group = c("All dental", "Unseen studies"),
  n = c(nrow(dental), nrow(dental_unseen)),
  healthy_n = c(sum(dental$y_true == 1), sum(dental_unseen$y_true == 1)),
  nonhealthy_n = c(sum(dental$y_true == 0), sum(dental_unseen$y_true == 0)),
  n_all_zero_panel = c(sum(as_true(dental$panel_all_zero)),
                       sum(as_true(dental_unseen$panel_all_zero))),
  AUC = c(a_all, a_unseen)
)
write.csv(roc_metrics, file.path(out_dir, "Supplementary_Fig_X_Nature_ROC_metrics.csv"),
          row.names = FALSE)
write.csv(data.frame(n_pairs = nrow(pairs), spearman_rho = rho,
                     mean_dental_minus_saliva = mean_diff,
                     median_dental_minus_saliva = median_diff),
          file.path(out_dir, "Supplementary_Fig_X_Nature_paired_metrics.csv"),
          row.names = FALSE)

cat("\nROC AUCs:", sprintf("%.3f, %.3f", a_all, a_unseen), "\n")
cat("Paired rho:", sprintf("%.3f", rho), "\n")
cat("Mean dental - saliva:", sprintf("%.2f", mean_diff), "\n")
cat("Saved PDF: ", pdf_file, "\n", sep = "")
cat("Saved PNG: ", png_file, "\n", sep = "")
