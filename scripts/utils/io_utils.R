# ============================================================
# io_utils.R
# Input/output and data preparation utilities for OMWI analyses
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

parse_args <- function(defaults = list()) {
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) == 0) {
    return(defaults)
  }

  out <- defaults

  for (arg in args) {
    if (!grepl("^--", arg)) next

    kv <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    if (length(kv) != 2) {
      stop("Invalid argument format: ", arg, "\nUse --key=value")
    }

    key <- kv[1]
    value <- kv[2]
    out[[key]] <- value
  }

  out
}

make_output_dirs <- function(out_dir) {
  dirs <- c(
    out_dir,
    file.path(out_dir, "figures"),
    file.path(out_dir, "tables"),
    file.path(out_dir, "logs")
  )

  invisible(lapply(dirs, dir.create, showWarnings = FALSE, recursive = TRUE))
}

check_file_exists <- function(file, label = "file") {
  if (!file.exists(file)) {
    stop("Cannot find ", label, ": ", file)
  }
}

read_abundance_table <- function(file) {
  check_file_exists(file, "abundance table")

  dt <- fread(file, data.table = FALSE, check.names = FALSE)

  if (ncol(dt) < 2) {
    stop("Abundance table must contain at least two columns: ", file)
  }

  colnames(dt)[1] <- "Taxon"
  dt
}

read_metadata <- function(file, required_cols) {
  check_file_exists(file, "metadata")

  meta <- fread(file, data.table = FALSE, check.names = FALSE)

  missing_cols <- setdiff(required_cols, colnames(meta))
  if (length(missing_cols) > 0) {
    stop("Metadata is missing required columns: ",
         paste(missing_cols, collapse = ", "))
  }

  meta$sample_id <- trimws(as.character(meta$sample_id))
  meta
}

collapse_same_taxon <- function(dt) {
  x <- as.data.table(dt)
  num_cols <- setdiff(colnames(x), "Taxon")

  x2 <- x[, lapply(.SD, function(z) {
    sum(as.numeric(z), na.rm = TRUE)
  }), by = Taxon, .SDcols = num_cols]

  as.data.frame(x2, check.names = FALSE)
}

make_num_matrix <- function(dt, id_col = "Taxon") {
  cols <- setdiff(colnames(dt), id_col)
  x <- dt[, cols, drop = FALSE]

  mat <- do.call(cbind, lapply(x, function(z) {
    as.numeric(as.character(z))
  }))

  mat <- as.matrix(mat)
  rownames(mat) <- dt[[id_col]]
  colnames(mat) <- cols

  mat
}

prepare_abundance_matrix <- function(dt, meta, sample_col = "sample_id") {
  dt <- collapse_same_taxon(dt)
  mat <- make_num_matrix(dt, "Taxon")

  colnames(mat) <- trimws(colnames(mat))

  common_samples <- intersect(colnames(mat), meta[[sample_col]])

  if (length(common_samples) < 10) {
    stop("Too few matched samples between abundance table and metadata: ",
         length(common_samples))
  }

  mat <- mat[, common_samples, drop = FALSE]
  meta_sub <- meta[match(common_samples, meta[[sample_col]]), , drop = FALSE]

  rownames(meta_sub) <- meta_sub[[sample_col]]
  meta_sub <- meta_sub[colnames(mat), , drop = FALSE]

  if (max(mat, na.rm = TRUE) > 1) {
    mat <- mat / 100
  }

  keep_samples <- colSums(mat, na.rm = TRUE) > 0
  mat <- mat[, keep_samples, drop = FALSE]
  meta_sub <- meta_sub[colnames(mat), , drop = FALSE]

  feature_var <- apply(mat, 1, var, na.rm = TRUE)
  keep_features <- is.finite(feature_var) & feature_var > 0
  mat <- mat[keep_features, , drop = FALSE]

  mat <- sweep(mat, 2, colSums(mat, na.rm = TRUE), "/")
  mat[!is.finite(mat)] <- 0

  list(
    mat = mat,
    meta = meta_sub
  )
}