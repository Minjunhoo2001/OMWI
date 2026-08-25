#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx)) return(default)
  if (idx == length(args)) stop("Missing value for ", flag)
  args[idx + 1]
}

input_file <- get_arg("--input")
coef_file <- get_arg("--coef", "models/omwi_species_coefficients.tsv")
universe_file <- get_arg("--universe", "models/omwi_species_universe.tsv")
output_file <- get_arg("--output", "omwi_scores.tsv")
pseudocount <- as.numeric(get_arg("--pseudocount", "1e-6"))

if (is.null(input_file)) {
  stop("Usage: Rscript calc_omwi.R --input species_abundance.tsv --coef models/omwi_species_coefficients.tsv --universe models/omwi_species_universe.tsv --output omwi_scores.tsv")
}

cat("[OMWI] Reading input species table:", input_file, "\n")
dat <- fread(input_file, data.table = FALSE, check.names = FALSE)

if (ncol(dat) < 2) {
  stop("Input table must contain one feature column and at least one sample column.")
}

feature_col <- colnames(dat)[1]
features <- as.character(dat[[feature_col]])

X <- dat[, -1, drop = FALSE]
for (j in seq_len(ncol(X))) {
  X[[j]] <- suppressWarnings(as.numeric(X[[j]]))
}

if (any(is.na(as.matrix(X)))) {
  stop("Input abundance table contains non-numeric or NA values.")
}

rownames(X) <- features

cat("[OMWI] Reading universe:", universe_file, "\n")
universe <- fread(universe_file, data.table = FALSE)
universe_taxa <- as.character(universe[[1]])
universe_taxa <- universe_taxa[!is.na(universe_taxa) & universe_taxa != ""]

cat("[OMWI] Reading coefficients:", coef_file, "\n")
coef_df <- fread(coef_file, data.table = FALSE)

if (!all(c("feature", "coef") %in% colnames(coef_df))) {
  stop("Coefficient file must contain columns: feature, coef")
}

coef_df$feature <- as.character(coef_df$feature)
coef_df$coef <- as.numeric(coef_df$coef)

intercept <- coef_df$coef[coef_df$feature == "(Intercept)"]
if (length(intercept) == 0) {
  intercept <- 0
} else {
  intercept <- intercept[1]
}

coef_taxa <- coef_df$feature[coef_df$feature != "(Intercept)"]
coef_values <- coef_df$coef[coef_df$feature != "(Intercept)"]
names(coef_values) <- coef_taxa

# clean MetaPhlAn-style names if needed
clean_name <- function(x) {
  x <- gsub("^.*\\|s__", "", x)
  x <- gsub("^s__", "", x)
  x <- gsub("\\|t__.*$", "", x)
  x
}

rownames(X) <- clean_name(rownames(X))

# If duplicated species names exist after cleaning, sum them
X_dt <- data.table(feature = rownames(X), X, check.names = FALSE)
X_dt <- X_dt[, lapply(.SD, sum, na.rm = TRUE), by = feature]
features2 <- X_dt$feature
X2 <- as.data.frame(X_dt[, -1, drop = FALSE], check.names = FALSE)
rownames(X2) <- features2

# Convert MetaPhlAn percentages to proportions if needed
sample_sums <- colSums(X2, na.rm = TRUE)
if (median(sample_sums, na.rm = TRUE) > 2) {
  cat("[OMWI] Input appears to be percentage scale. Dividing by 100.\n")
  X2 <- X2 / 100
}

# Align to fixed OMWI universe
X_uni <- matrix(
  0,
  nrow = length(universe_taxa),
  ncol = ncol(X2),
  dimnames = list(universe_taxa, colnames(X2))
)

matched_taxa <- intersect(rownames(X2), universe_taxa)
X_uni[matched_taxa, ] <- as.matrix(X2[matched_taxa, , drop = FALSE])

cat("[OMWI] Matched universe taxa:", length(matched_taxa), "/", length(universe_taxa), "\n")
cat("[OMWI] Missing universe taxa were assigned zero abundance.\n")

# CLR over fixed universe
X_pc <- X_uni + pseudocount
log_X <- log(X_pc)
clr_X <- sweep(log_X, 2, colMeans(log_X), FUN = "-")

# Calculate OMWI
missing_coef_taxa <- setdiff(coef_taxa, rownames(clr_X))
if (length(missing_coef_taxa) > 0) {
  stop("Some coefficient taxa are not in the OMWI universe: ", paste(missing_coef_taxa, collapse = ", "))
}

score <- intercept + as.numeric(t(coef_values) %*% clr_X[coef_taxa, , drop = FALSE])

out <- data.frame(
  sample_id = colnames(clr_X),
  OMWI = score,
  predicted_state = ifelse(score > 0, "health_associated", "non_healthy_associated"),
  matched_universe_taxa = length(matched_taxa),
  total_universe_taxa = length(universe_taxa),
  stringsAsFactors = FALSE
)

fwrite(out, output_file, sep = "\t")

cat("[OMWI] Done.\n")
cat("[OMWI] Output:", output_file, "\n")