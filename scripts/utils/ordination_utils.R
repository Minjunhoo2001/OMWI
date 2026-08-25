# ============================================================
# ordination_utils.R
# Ordination and PERMANOVA utilities for OMWI analyses
# ============================================================

suppressPackageStartupMessages({
  library(vegan)
  library(dplyr)
})

run_pcoa <- function(mat, distance_method = "bray") {
  dist_obj <- vegdist(t(mat), method = distance_method)

  ord <- cmdscale(
    dist_obj,
    k = 2,
    eig = TRUE,
    add = TRUE
  )

  eig_pct <- ord$eig / sum(ord$eig, na.rm = TRUE) * 100

  points <- as.data.frame(ord$points)
  colnames(points) <- c("PCoA1", "PCoA2")

  list(
    dist = dist_obj,
    points = points,
    eig_pct = eig_pct
  )
}

run_permanova <- function(dist_obj, meta_df, formula_text, permutations = 999) {
  fml <- as.formula(paste0("dist_obj ~ ", formula_text))

  adonis2(
    fml,
    data = meta_df,
    permutations = permutations
  )
}

tidy_permanova <- function(adonis_obj, analysis_name) {
  out <- as.data.frame(adonis_obj)
  out$term <- rownames(out)
  rownames(out) <- NULL
  out$analysis <- analysis_name

  out %>%
    select(analysis, term, everything())
}

make_permanova_label <- function(adonis_obj, term = NULL) {
  ad <- as.data.frame(adonis_obj)

  if (is.null(term)) {
    terms <- rownames(ad)
    terms <- terms[!grepl("^Residual|^Total", terms, ignore.case = TRUE)]
    term <- terms[1]
  }

  r2 <- ad[term, "R2"]
  p <- ad[term, "Pr(>F)"]

  r2_txt <- ifelse(is.na(r2), "NA", sprintf("%.3f", r2))
  p_txt <- ifelse(
    is.na(p),
    "NA",
    ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
  )

  paste0("PERMANOVA\nR² = ", r2_txt, "\nP = ", p_txt)
}

make_center_df <- function(df, group_col) {
  df %>%
    group_by(.data[[group_col]]) %>%
    summarise(
      cx = mean(PCoA1, na.rm = TRUE),
      cy = mean(PCoA2, na.rm = TRUE),
      .groups = "drop"
    )
}

order_binary_factor <- function(x, first_levels = c("healthy", "Healthy")) {
  x_chr <- as.character(x)
  lv <- unique(x_chr)

  first <- lv[tolower(lv) %in% tolower(first_levels)]
  rest <- setdiff(lv, first)

  factor(x_chr, levels = c(first, rest))
}