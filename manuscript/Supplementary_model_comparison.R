# ============================================================
# Original Elastic Net vs 4 conventional ML models — external ROC
# Read only external validation predictions; do NOT retrain models.
# All five models: differmodel/final/*_external_predictions.csv
# Original Elastic Net only; never use Compact K or the final OMWI score.
# Output: differmodel/figure/ (PDF, 600-dpi PNG, AUC summary CSV)
# ============================================================

for (pkg in c("pROC", "ggplot2")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
library(pROC)
library(ggplot2)

# Use / for Windows paths
base_dir <- Sys.getenv("OMWI_ROOT", getwd())
ml_dir   <- file.path(base_dir, "differmodel", "final")
out_dir  <- file.path(base_dir, "differmodel", "figure")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --------------------------------
# 1. Read five original models (no feature-panel compression)
# --------------------------------
ml_models <- c("ElasticNet", "RandomForest", "SVM_linear", "SVM_RBF", "XGBoost")
# Confirm this filename corresponds to the original Elastic Net in Supplementary Table 5.
model_files <- setNames(file.path(ml_dir, paste0(ml_models, "_external_predictions.csv")), ml_models)
model_files["ElasticNet"] <- file.path(ml_dir, "ElasticNet_external_predictions.csv")
if (!file.exists(model_files["ElasticNet"])) {
  available <- list.files(ml_dir, pattern = "external.*predictions.*[.]csv$", ignore.case = TRUE)
  stop("Original Elastic Net prediction file not found: ", model_files["ElasticNet"],
       "\nAvailable prediction files: ", paste(available, collapse = ", "),
       "\nSet model_files['ElasticNet'] to the verified ORIGINAL model file; do not use Compact K.")
}
if (grepl("compact|stability|subset|omwi|omhi", basename(model_files["ElasticNet"]), ignore.case = TRUE)) {
  stop("Elastic Net input appears to be a compact/stability/OMWI model. Verify the original file.")
}

model_names <- c(
  ElasticNet = "Elastic Net",
  XGBoost = "XGBoost",
  SVM_RBF = "SVM (RBF)",
  RandomForest = "Random Forest",
  SVM_linear = "SVM (linear)"
)
model_colors <- c(
  ElasticNet = "#3C5488",    # retain original figure styling
  XGBoost = "#E64B35",
  SVM_RBF = "#00A087",
  RandomForest = "#7E6148",
  SVM_linear = "#8491B4"
)

preds <- list()
for (m in ml_models) {
  f <- unname(model_files[m])
  if (!file.exists(f)) stop("Missing model prediction file: ", f)
  d <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE,
                colClasses = "character")
  req <- c("sample_id", "y_true", "prob_healthy")
  if (!all(req %in% names(d))) {
    stop("Missing columns in ", f, ": ", paste(setdiff(req, names(d)), collapse = ", "))
  }
  if (anyNA(d$y_true) || !all(trimws(d$y_true) %in% c("0", "1"))) {
    stop("Invalid labels (expected 0 = non-healthy, 1 = healthy) in ", f)
  }
  preds[[m]] <- data.frame(sample_id = trimws(d$sample_id),
                           y_true = as.integer(d$y_true),
                           score = suppressWarnings(as.numeric(d$prob_healthy)),
                           stringsAsFactors = FALSE)
  cat("Input for ", m, ": ", f, "\n", sep = "")
}

# --------------------------------
# 2. Strictly match external samples
# --------------------------------
reference <- preds$ElasticNet[, c("sample_id", "y_true")]
if (!setequal(unique(reference$y_true), c(0L, 1L))) stop("Both outcome classes are required")
for (m in names(preds)) {
  d <- preds[[m]]
  if (anyNA(d) || any(!nzchar(d$sample_id)) || anyDuplicated(d$sample_id) ||
      any(!is.finite(d$score)) || any(d$score < 0 | d$score > 1)) {
    stop("Invalid or duplicated values in ", m)
  }
  if (!setequal(reference$sample_id, d$sample_id)) {
    stop("External sample IDs differ between original Elastic Net and ", m)
  }
  d <- d[match(reference$sample_id, d$sample_id), , drop = FALSE]
  if (!identical(as.integer(d$y_true), as.integer(reference$y_true))) {
    stop("Health/disease labels differ for model ", m)
  }
  preds[[m]] <- d
}
cat("Matched external samples:", nrow(reference),
    "(healthy:", sum(reference$y_true == 1),
    ", non-healthy:", sum(reference$y_true == 0), ")\n")

# --------------------------------
# 3. Calculate five ROC curves/AUCs
# --------------------------------
roc_data <- list()
auc_values <- numeric(length(preds))
names(auc_values) <- names(preds)

for (m in names(preds)) {
  d <- preds[[m]]
  roc_obj <- pROC::roc(response = d$y_true, predictor = d$score,
                       levels = c(0, 1), direction = "<", quiet = TRUE)
  auc_values[m] <- as.numeric(pROC::auc(roc_obj))
  xy <- data.frame(model = m,
                   FPR = 1 - as.numeric(roc_obj$specificities),
                   TPR = as.numeric(roc_obj$sensitivities))
  roc_data[[m]] <- xy[order(xy$FPR, xy$TPR), , drop = FALSE]
}
plot_order <- c("ElasticNet", "XGBoost", "SVM_RBF", "RandomForest", "SVM_linear")
roc_df <- do.call(rbind, roc_data[plot_order])
roc_df$model <- factor(roc_df$model, levels = plot_order)
legend_labels <- setNames(
  sprintf("%s (AUC = %.3f)", model_names[plot_order], auc_values[plot_order]),
  plot_order
)

# --------------------------------
# 4. Plot (Nature-like style)
# --------------------------------
p <- ggplot() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              linewidth = 0.5, color = "grey70") +
  geom_path(data = subset(roc_df, model != "ElasticNet"),
            aes(x = FPR, y = TPR, color = model, group = model),
            linewidth = 0.9, lineend = "round") +
  geom_path(data = subset(roc_df, model == "ElasticNet"),
            aes(x = FPR, y = TPR, color = model, group = model),
            linewidth = 1.4, lineend = "round") +
  scale_color_manual(values = model_colors, breaks = plot_order,
                     labels = legend_labels) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = c(0, 0)
  ) +
  coord_equal() +
  labs(x = "1 - Specificity", y = "Sensitivity", color = NULL) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.line = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 0.6),
    axis.ticks.length = grid::unit(0.18, "cm"),
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(color = "black", size = 11),
    legend.position = "right",
    legend.text = element_text(size = 10, color = "black"),
    legend.key = element_blank(),
    legend.key.width = grid::unit(1.0, "cm"),
    legend.background = element_blank(),
    plot.margin = margin(8, 10, 8, 8)
  ) +
  guides(color = guide_legend(ncol = 1))

print(p)

# --------------------------------
# 5. Save files
# --------------------------------
file_prefix <- "Original_ElasticNet_vs_4models_external_ROC"
ggsave(file.path(out_dir, paste0(file_prefix, ".png")),
       p, width = 6.3, height = 5.6, dpi = 600, bg = "white")
# Cairo PDF keeps TrueType fonts more reliably; most recent R installations support it.
if (capabilities("cairo")) {
  ggsave(file.path(out_dir, paste0(file_prefix, ".pdf")),
         p, width = 6.3, height = 5.6, device = grDevices::cairo_pdf, bg = "white")
} else {
  message("Cairo not available; saving standard vector PDF with serif font.")
  p_pdf <- p + theme(text = element_text(family = "serif"))
  ggsave(file.path(out_dir, paste0(file_prefix, ".pdf")),
         p_pdf, width = 6.3, height = 5.6, device = "pdf", bg = "white")
}

auc_table <- data.frame(
  model = unname(model_names[plot_order]),
  AUC = round(unname(auc_values[plot_order]), 6),
  n = nrow(reference),
  prediction_file = basename(unname(model_files[plot_order]))
)
write.csv(auc_table, file.path(out_dir, paste0(file_prefix, "_AUC.csv")),
          row.names = FALSE)
print(auc_table)
cat("Files saved under:", out_dir, "\n")

