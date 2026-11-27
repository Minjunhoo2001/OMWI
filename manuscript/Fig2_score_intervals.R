#!/usr/bin/env Rscript
# ================================================================
# Fig. 2f | Healthy / non-healthy proportions across OMWI bins
#          Training vs independent external validation (separate facets)
#
# Original Fig2g design retained: seven score bins (extreme left-tail merged), cyan Healthy,
# red Non-healthy, stacked proportions, 50% horizontal reference,
# and vertical dashed boundary at OMWI = 0.
#
# Input: frozen stability-model final predictions. NO model retraining.
# Training scores are IN-SAMPLE; external scores are independent.
# ================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# 0. Fixed Windows paths; independent of working directory / getwd().
ROOT_DIR <- Sys.getenv("OMWI_ROOT", getwd())
MODEL_DIR <- file.path(ROOT_DIR, "data/model_outputs")
OUT_DIR <- file.path(ROOT_DIR, "fig2")
TABLE_DIR <- file.path(OUT_DIR, "tables")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)

TRAIN_FILE <- file.path(MODEL_DIR, "final", "species_train_in_sample_scores.csv")
EXTERNAL_FILE <- file.path(MODEL_DIR, "final", "external_independent_predictions.csv")

# 1. Read exactly the final frozen model outputs (no older model files).
read_scores <- function(path, cohort) {
  if (!file.exists(path)) stop("Missing input CSV: ", path)
  tab <- data.table::fread(path, data.table = FALSE, check.names = FALSE,
                           select = c("sample_id", "disease", "OMHI"))
  required <- c("sample_id", "disease", "OMHI")
  missing <- setdiff(required, names(tab))
  if (length(missing)) stop(basename(path), " missing: ", paste(missing, collapse = ", "))
  
  tab$sample_id <- trimws(as.character(tab$sample_id))
  if (anyNA(tab$sample_id) || any(!nzchar(tab$sample_id)) || anyDuplicated(tab$sample_id))
    stop("Invalid or duplicate sample IDs in: ", path)
  
  tab$OMWI <- suppressWarnings(as.numeric(tab$OMHI))
  if (any(!is.finite(tab$OMWI))) stop("Non-finite OMWI values in: ", path)
  
  raw_label <- tolower(trimws(as.character(tab$disease)))
  tab$health <- ifelse(raw_label %in% c("healthy", "health", "control", "ctrl", "0"),
                       "Healthy",
                       ifelse(raw_label %in% c("disease", "diseased", "non-healthy", "nonhealthy", "case", "1"),
                              "Non-healthy", NA_character_))
  if (anyNA(tab$health)) {
    stop("Unexpected disease labels in ", basename(path), ": ",
         paste(utils::head(unique(tab$disease[is.na(tab$health)]), 10), collapse = ", "))
  }
  tab$dataset <- cohort
  tab[, c("sample_id", "OMWI", "health", "dataset"), drop = FALSE]
}

train <- read_scores(TRAIN_FILE, "Training")
external <- read_scores(EXTERNAL_FILE, "External validation")
if (length(intersect(train$sample_id, external$sample_id)))
  stop("Training and external validation contain overlapping sample IDs.")

# Preflight: ensure the expected version, fail loudly if an old model was used.
if (nrow(train) != 1077L || nrow(external) != 716L) {
  stop("Unexpected model sample counts: Training=", nrow(train),
       ", External=", nrow(external),
       ". Expected 1077 / 716 for the current stability model.")
}

dat <- rbind(train, external)
dat$dataset <- factor(dat$dataset, levels = c("Training", "External validation"))
dat$health <- factor(dat$health, levels = c("Non-healthy", "Healthy"))

# 2. SAME ORIGINAL FIG2f BINNING, except merge the two extreme left-tail bins.
# Previously: < -2 and [-2, -1.5). Now: <= -1.5.
# Both cohorts use exactly the same boundaries; the model scores are unchanged.
BIN_BREAKS <- c(-Inf, -1.5, -1, -0.5, 0, 0.5, 1, Inf)
BIN_LABELS <- c("<= -1.5", "(-1.5, -1)", "[-1, -0.5)",
                "[-0.5, 0)", "[0, 0.5)", "[0.5, 1)", ">= 1")
dat$omwi_bin <- cut(dat$OMWI, breaks = BIN_BREAKS, labels = BIN_LABELS,
                    include.lowest = TRUE, right = FALSE)
# Include any exact -1.5 score in the first group to match '<= -1.5'.
dat$omwi_bin[dat$OMWI == -1.5] <- BIN_LABELS[1]
if (anyNA(dat$omwi_bin)) stop("OMWI bin assignment produced missing values.")

# 3. Counts and proportions PER BIN, WITHIN EACH DATASET (not pooled).
counts <- as.data.frame(table(dataset = dat$dataset,
                              omwi_bin = dat$omwi_bin,
                              health = dat$health),
                        stringsAsFactors = FALSE)
names(counts)[4] <- "n"
counts$n <- as.integer(counts$n)
counts$dataset <- factor(counts$dataset, levels = levels(dat$dataset))
counts$omwi_bin <- factor(counts$omwi_bin, levels = BIN_LABELS)
counts$health <- factor(counts$health, levels = levels(dat$health))

# The table() call keeps empty cells, so zero-count groups remain represented.
total_per_bin <- aggregate(n ~ dataset + omwi_bin, data = counts, FUN = sum)
names(total_per_bin)[3] <- "n_bin_total"
counts <- merge(counts, total_per_bin,
                by = c("dataset", "omwi_bin"), all.x = TRUE, sort = FALSE)
counts$proportion <- ifelse(counts$n_bin_total > 0,
                            counts$n / counts$n_bin_total, NA_real_)
counts$percentage <- 100 * counts$proportion
counts <- counts[order(counts$dataset, counts$omwi_bin, counts$health), ]
rownames(counts) <- NULL

if (sum(counts$n[counts$dataset == "Training"]) != nrow(train) ||
    sum(counts$n[counts$dataset == "External validation"]) != nrow(external))
  stop("Bin counts do not sum to the corresponding dataset sizes.")

# Save a complete, auditable Supplementary Table 2 source.
write_table <- function(df, filename) {
  data.table::fwrite(df, file.path(TABLE_DIR, filename), bom = TRUE, na = "NA")
}
write_table(counts, "Fig2g_7bins_Training_External_bin_counts_and_proportions.csv")
write_table(subset(counts, dataset == "Training"),
            "Fig2g_7bins_Training_bin_counts_and_proportions.csv")
write_table(subset(counts, dataset == "External validation"),
            "Fig2g_7bins_External_bin_counts_and_proportions.csv")

# A supplementary wide table with both statuses on one row per bin.
wide <- reshape(as.data.frame(counts),
                idvar = c("dataset", "omwi_bin", "n_bin_total"),
                timevar = "health", direction = "wide",
                v.names = c("n", "percentage"))
wide <- wide[order(wide$dataset, wide$omwi_bin), ]
write_table(wide, "Fig2g_7bins_Training_External_bin_summary_wide.csv")

# 4. EXACT ORIGINAL FIG2f styling: cyan below, red above.
# This block deliberately reproduces the original split-figure theme.
# NO alternative NC/minimal theme, white segment outlines, sample-size labels,
# legend rearrangement or adjusted typography.
FILL_COLORS <- c("Non-healthy" = "#E64B35", "Healthy" = "#4DBBD5")

plot_original_fig2f <- function(layout = c("vertical", "horizontal")) {
  layout <- match.arg(layout)
  p <- ggplot(counts, aes(x = omwi_bin, y = proportion, fill = health)) +
    geom_col(width = 0.84, color = "#555555", linewidth = 0.45,
             position = position_stack(reverse = FALSE), na.rm = TRUE) +
    geom_hline(yintercept = 0.5, linetype = "dashed",
               color = "#B5B5B5", linewidth = 0.55) +
    # 4 bins below zero; zero boundary lies between bin 4 and bin 5.
    geom_vline(xintercept = 4.5, linetype = "dashed",
               color = "#858585", linewidth = 0.65) +
    facet_wrap(~dataset, ncol = if (layout == "vertical") 1 else 2) +
    scale_fill_manual(values = FILL_COLORS,
                      breaks = c("Non-healthy", "Healthy"), drop = FALSE) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2),
                       labels = function(x) paste0(round(x * 100), "%"),
                       limits = c(0, 1), expand = expansion(mult = c(0, 0))) +
    scale_x_discrete(drop = FALSE) +
    labs(x = "OMWI bin", y = "Sample proportion by OMWI bin", fill = NULL) +
    theme_classic(base_size = 11) +
    theme(
      strip.background = element_rect(fill = "white", color = "#555555", linewidth = 0.45),
      strip.text = element_text(face = "bold", color = "black", size = 11),
      panel.border = element_rect(fill = NA, color = "#555555", linewidth = 0.45),
      axis.line = element_blank(),
      axis.text.x = element_text(angle = 48, hjust = 1, vjust = 1,
                                 color = "black", size = 9.5),
      axis.text.y = element_text(color = "black", size = 10),
      axis.title = element_text(color = "black", size = 11),
      axis.ticks = element_line(color = "black", linewidth = 0.4),
      panel.spacing.x = grid::unit(0.5, "lines"),
      panel.spacing.y = grid::unit(0.7, "lines"),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.text = element_text(size = 11, color = "black"),
      legend.key.width = grid::unit(0.8, "cm"),
      plot.margin = margin(7, 7, 5, 7)
    )
  p
}

p_vertical <- plot_original_fig2f("vertical")
p_horizontal <- plot_original_fig2f("horizontal")
print(p_vertical)

save_fig2g <- function(p, filename, width, height) {
  ggplot2::ggsave(file.path(OUT_DIR, paste0(filename, ".pdf")), p,
                  width = width, height = height, units = "in", bg = "white",
                  device = if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf)
  ggplot2::ggsave(file.path(OUT_DIR, paste0(filename, ".png")), p,
                  width = width, height = height, units = "in", dpi = 600, bg = "white")
}

# Original appearance, two ways of arranging the SAME two cohorts.
# Vertical is better when the original Fig2g was one wide panel.
save_fig2g(p_vertical,   "Fig2g_7bins_original_style_VERTICAL",   9, 8.0)
save_fig2g(p_horizontal, "Fig2g_7bins_original_style_HORIZONTAL", 6, 4)

message("\n[DONE] Original Fig2g appearance retained; only binning and facets changed.")
message("Input: ", MODEL_DIR)
message("Output: ", OUT_DIR)
message("Training n = ", nrow(train), "; External n = ", nrow(external))
message("7 bins. First bin <= -1.5 in both cohorts. Training scores are in-sample.")
