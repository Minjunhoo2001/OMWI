# ============================================================
# Script: pcoa_batch_effect_train_species.R
# Project: Oral Microbiome Health Index (OMHI)
#
# Purpose:
#   Evaluate study-level batch effects before and after batch
#   correction using Bray-Curtis PCoA and PERMANOVA.
#
# Inputs:
#   1. meta_train.csv
#   2. train_species.tsv
#   3. train_species.batch.normal.tsv
#
# Outputs:
#   1. Four-panel PCoA figure before and after batch correction
#   2. PERMANOVA summary table
#   3. Session information
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(cowplot)
  library(vegan)
})

options(stringsAsFactors = FALSE)

options(error = function() {
  traceback(2)
  quit(save = "no", status = 1)
})

# ============================================================
# 0. Locate project scripts
# ============================================================

script_dir <- dirname(normalizePath(sys.frame(1)$ofile %||% commandArgs()[1]))
utils_dir <- file.path(script_dir, "utils")

source(file.path(utils_dir, "io_utils.R"))
source(file.path(utils_dir, "ordination_utils.R"))
source(file.path(utils_dir, "plot_utils.R"))

# ============================================================
# 1. Parameters
# ============================================================

defaults <- list(
  raw_dir = "data/02.by_level",
  batch_dir = "data/03.batch_train",
  meta_dir = "data/01.split",
  out_dir = "results/pcoa_species_train_only",
  level = "species",
  batch_var = "study_name",
  permutations = "999"
)

args <- parse_args(defaults)

RAW_DIR <- args$raw_dir
BATCH_DIR <- args$batch_dir
META_DIR <- args$meta_dir
OUT_DIR <- args$out_dir

LEVEL <- args$level
BATCH_VAR <- args$batch_var
DISEASE_VAR <- args$disease_var
SYSTEM_VAR <- args$system_var
PERM_N <- as.integer(args$permutations)

make_output_dirs(OUT_DIR)

# ============================================================
# 2. Input files
# ============================================================

meta_file <- file.path(META_DIR, "meta_train.csv")
raw_file <- file.path(RAW_DIR, paste0("train_", LEVEL, ".tsv"))
adj_file <- file.path(BATCH_DIR, paste0("train_", LEVEL, ".batch.normal.tsv"))

message("[OMHI] Metadata: ", meta_file)
message("[OMHI] Raw abundance table: ", raw_file)
message("[OMHI] Batch-corrected table: ", adj_file)
message("[OMHI] Output directory: ", OUT_DIR)

meta <- read_metadata(
  meta_file,
  required_cols = c("sample_id", BATCH_VAR, DISEASE_VAR, SYSTEM_VAR)
)

raw_dt <- read_abundance_table(raw_file)
adj_dt <- read_abundance_table(adj_file)

# ============================================================
# 3. Prepare matrices and metadata
# ============================================================

raw_obj <- prepare_abundance_matrix(raw_dt, meta)
adj_obj <- prepare_abundance_matrix(adj_dt, meta)

common_samples <- intersect(colnames(raw_obj$mat), colnames(adj_obj$mat))

if (length(common_samples) < 10) {
  stop("Too few common samples between raw and adjusted matrices: ",
       length(common_samples))
}

raw_mat <- raw_obj$mat[, common_samples, drop = FALSE]
adj_mat <- adj_obj$mat[, common_samples, drop = FALSE]
meta_use <- raw_obj$meta[common_samples, , drop = FALSE]

meta_use[[BATCH_VAR]] <- factor(meta_use[[BATCH_VAR]])
meta_use[[DISEASE_VAR]] <- order_binary_factor(meta_use[[DISEASE_VAR]])
meta_use[[SYSTEM_VAR]] <- order_binary_factor(
  meta_use[[SYSTEM_VAR]],
  first_levels = c("Healthy")
)

message("[OMHI] Samples used: ", ncol(raw_mat))
message("[OMHI] Raw features used: ", nrow(raw_mat))
message("[OMHI] Adjusted features used: ", nrow(adj_mat))

sample_summary <- data.frame(
  level = LEVEL,
  n_samples = ncol(raw_mat),
  n_raw_features = nrow(raw_mat),
  n_adjusted_features = nrow(adj_mat),
  n_studies = length(unique(meta_use[[BATCH_VAR]])),
  n_disease_groups = length(unique(meta_use[[DISEASE_VAR]])),
  n_system_groups = length(unique(meta_use[[SYSTEM_VAR]]))
)

fwrite(
  sample_summary,
  file.path(OUT_DIR, "tables", paste0("sample_feature_summary_", LEVEL, ".csv"))
)

# ============================================================
# 4. PCoA
# ============================================================

raw_pcoa <- run_pcoa(raw_mat, distance_method = "bray")
adj_pcoa <- run_pcoa(adj_mat, distance_method = "bray")

raw_df <- cbind(raw_pcoa$points, meta_use)
adj_df <- cbind(adj_pcoa$points, meta_use)

# ============================================================
# 5. PERMANOVA
# ============================================================

adonis_raw_batch <- run_permanova(
  raw_pcoa$dist,
  meta_use,
  paste0("`", BATCH_VAR, "`"),
  permutations = PERM_N
)

adonis_adj_batch <- run_permanova(
  adj_pcoa$dist,
  meta_use,
  paste0("`", BATCH_VAR, "`"),
  permutations = PERM_N
)


# Optional adjusted models.
adonis_adj_disease_plus_batch <- run_permanova(
  adj_pcoa$dist,
  meta_use,
  paste0("`", DISEASE_VAR, "` + `", BATCH_VAR, "`"),
  permutations = PERM_N
)

adonis_adj_system_plus_batch <- run_permanova(
  adj_pcoa$dist,
  meta_use,
  paste0("`", SYSTEM_VAR, "` + `", BATCH_VAR, "`"),
  permutations = PERM_N
)

permanova_table <- bind_rows(
  tidy_permanova(adonis_raw_batch, "batch_before_correction"),
  tidy_permanova(adonis_adj_batch, "batch_after_correction"),
  tidy_permanova(adonis_adj_disease, "disease_after_correction"),
  tidy_permanova(adonis_adj_system, "system_after_correction"),
  tidy_permanova(adonis_adj_disease_plus_batch, "disease_plus_batch_after_correction"),
  tidy_permanova(adonis_adj_system_plus_batch, "system_plus_batch_after_correction")
)

fwrite(
  permanova_table,
  file.path(OUT_DIR, "tables", paste0("PERMANOVA_", LEVEL, "_train_summary.csv"))
)

label_raw_batch <- make_permanova_label(adonis_raw_batch)
label_adj_batch <- make_permanova_label(adonis_adj_batch)

# ============================================================
# 6. Colors
# ============================================================

batch_levels <- levels(meta_use[[BATCH_VAR]])

batch_cols <- scales::hue_pal()(length(batch_levels))
names(batch_cols) <- batch_levels

# ============================================================
# 7. Four-panel batch-effect figure
# ============================================================

center_raw <- make_center_df(raw_df, BATCH_VAR)
center_adj <- make_center_df(adj_df, BATCH_VAR)

p_before_disease <- plot_pcoa_by_group(
  df = raw_df,
  group_col = DISEASE_VAR,
  color_values = disease_cols,
  eig_pct = raw_pcoa$eig_pct,
  title = "Before batch correction",
  stat_label = label_raw_batch,
  shape_col = BATCH_VAR,
  base_size = 12
)

p_after_disease <- plot_pcoa_by_group(
  df = adj_df,
  group_col = DISEASE_VAR,
  color_values = disease_cols,
  eig_pct = adj_pcoa$eig_pct,
  title = "After batch correction",
  stat_label = label_adj_batch,
  shape_col = BATCH_VAR,
  base_size = 12
)

p_before_center <- plot_pcoa_centroid(
  df = raw_df,
  center_df = center_raw,
  group_col = BATCH_VAR,
  color_values = batch_cols,
  eig_pct = raw_pcoa$eig_pct,
  title = "Study centroids before correction",
  stat_label = label_raw_batch,
  base_size = 12
)

p_after_center <- plot_pcoa_centroid(
  df = adj_df,
  center_df = center_adj,
  group_col = BATCH_VAR,
  color_values = batch_cols,
  eig_pct = adj_pcoa$eig_pct,
  title = "Study centroids after correction",
  stat_label = label_adj_batch,
  base_size = 12
)

fig_batch_4panel <- plot_grid(
  p_before_disease,
  p_after_disease,
  p_before_center,
  p_after_center,
  labels = c("A", "B", "C", "D"),
  ncol = 2,
  label_size = 14
)

save_pdf_png(
  fig_batch_4panel,
  file.path(
    OUT_DIR,
    "figures",
    paste0("Fig_", LEVEL, "_train_batch_before_after_4panel")
  ),
  width = 13.8,
  height = 12.0
)

# ============================================================
# 8. Session information
# ============================================================

sink(file.path(OUT_DIR, "logs", paste0("sessionInfo_", LEVEL, "_pcoa.txt")))
cat("OMHI PCoA batch-effect analysis\n")
cat("Date: ", as.character(Sys.time()), "\n\n")
cat("Parameters:\n")
print(args)
cat("\nSession info:\n")
print(sessionInfo())
sink()

message("[OMHI] Analysis completed.")
message("[OMHI] Results written to: ", OUT_DIR)