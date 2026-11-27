#!/usr/bin/env Rscript
# ============================================================================
# Fig. 3 (a-f): Oral v3.0, one-file figure-generation workflow
# Reference layout: a/b/c across top; d on lower left; e above f on lower right.
#
# a  OMWI for each phenotype within each study
# b  Shannon diversity for the same samples and rows
# c  Species richness for the same samples and rows
# d  Training: OMWI by phenotype / disease system
# e  External validation: OMWI by phenotype / disease system
# f  System-specific AUC (95% CI) in training and external validation
#
# Uses the model's EXISTING scores. Does NOT retrain or tune the model.
# Input files are ONLY read from ROOT_DIR; outputs are ONLY written to OUT_DIR.
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(pROC)
  library(scales)
})
options(stringsAsFactors = FALSE)

# ---------------------- 0. Windows paths and parameters ----------------------
ROOT_DIR <- Sys.getenv("OMWI_ROOT", getwd())
OUT_DIR <- file.path(ROOT_DIR, "fig3")
RAW_DIR <- file.path(ROOT_DIR, "data/model_outputs")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# "in_sample" reproduces the score definition in the supplied Fig. 3 layout;
# "oof_mean" takes the per-sample average from 10 repeats of nested OOF.
# Do not label in-sample training results as cross-validated performance.
TRAIN_SCORE_MODE <- "in_sample"  # alternative: "oof_mean"

# If auto-discovery cannot find species abundances, list explicit source files
# (absolute paths INSIDE ROOT_DIR) here, e.g.:
SPECIES_FILES <- c(file.path(ROOT_DIR,"data/abundance","all_species.tsv"))


MIN_STUDY_GROUP_N <- 5L   # a/b/c: include study-phenotype groups with >=5 samples
MIN_DE_GROUP_N <- 10L    # d/e: group and system minimum sample size
MIN_F_CASE_N <- 5L       # f: min cases and healthy controls in each dataset
OMWI_XLIMIT <- c(-4, 4.5)  # d/e only; remove clipping with NULL
EXPORT_WIDTH <- 12      # inches; portrait figure, editable PDF
EXPORT_HEIGHT <- 16.1

SYSTEM_COLORS <- c(
  "Healthy" = "#2CA25F", "Abnormal BMI" = "#FFC20A",
  "Metabolic Disease" = "#0072B2", "Oral Disease" = "#C7A6D8",
  "Immune Disease" = "#D8C97A", "Respiratory Disease" = "#8FA8C9",
  "Digestive Disease" = "#D97C8A"
)
SYSTEM_ORDER <- names(SYSTEM_COLORS)
DISEASE_SYSTEMS <- setdiff(SYSTEM_ORDER, c("Healthy", "Abnormal BMI"))

# ----------------------------- 1. Basic checks ------------------------------
check_input <- function(path) {
  if (!file.exists(path)) stop("找不到输入文件: ", path, call. = FALSE)
  norm <- function(x) tolower(gsub("\\\\", "/", normalizePath(x, mustWork = FALSE)))
  if (!startsWith(norm(path), paste0(norm(ROOT_DIR), "/"))) {
    stop("输入路径必须在 Oral v3.0 下: ", path, call. = FALSE)
  }
  path
}
read_csv <- function(path) {
  data.table::fread(check_input(path), data.table = FALSE, check.names = FALSE)
}
first_col <- function(df, opts, required = TRUE) {
  x <- intersect(opts, names(df))
  if (length(x)) return(x[[1]])
  if (required) stop("没有找到这些列中的任意一列: ", paste(opts, collapse = ", "))
  NULL
}
clean_group <- function(x) {
  x <- trimws(as.character(x))
  x[tolower(x) %in% c("control", "healthy", "health", "normal", "ctrl")] <- "Healthy"
  x[tolower(x) %in% c("type ii diabetes", "type 2 diabetes mellitus")] <- "Type 2 diabetes"
  x[tolower(x) == "covid"] <- "COVID"
  x[tolower(x) == "rheumatoid arthritis"] <- "Rheumatoid Arthritis"
  x[tolower(x) == "adenoid cystic carcinoma"] <- "Adenoid Cystic Carcinoma"
  x[tolower(x) == "bmi low"] <- "Underweight"
  x[x %in% c("", "NA", "N/A", "Unknown")] <- NA_character_
  x
}
clean_system <- function(x, group) {
  x <- trimws(as.character(x))
  x[tolower(x) %in% c("bmi group", "abnormal bmi", "abnormal body weight")] <- "Abnormal BMI"
  for (nm in SYSTEM_ORDER) x[tolower(x) == tolower(nm)] <- nm
  x[group == "Healthy"] <- "Healthy"
  x
}
pretty_study <- function(x) {
  x <- sub("_([0-9]{4})$", " (\\1)", as.character(x))
  gsub("_", " ", x, fixed = TRUE)
}
stars <- function(p) {
  ifelse(is.na(p), "", ifelse(p < 1e-4, "****", ifelse(p < .001, "***",
                                                       ifelse(p < .01, "**", ifelse(p < .05, "*", "")))))
}
base_theme <- function(sz = 10) {
  theme_bw(base_size = sz) + theme(
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "grey45", linewidth = .48, fill = NA),
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black"),
    legend.position = "none",
    plot.title = element_blank(),
    plot.margin = margin(4, 4, 4, 4)
  )
}
save_both <- function(plot, prefix, w, h) {
  p_png <- file.path(OUT_DIR, paste0(prefix, ".png"))
  p_pdf <- file.path(OUT_DIR, paste0(prefix, ".pdf"))
  ggsave(p_png, plot = plot, width = w, height = h, units = "in",
         dpi = 600, bg = "white", limitsize = FALSE)
  pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
  ggsave(p_pdf, plot = plot, width = w, height = h, units = "in",
         device = pdf_device, bg = "white", limitsize = FALSE)
  message("输出：", p_pdf)
}

# ------------------------ 2. Model score tables ------------------------------
train_file <- file.path(RAW_DIR, "final", "species_train_in_sample_scores.csv")
external_file <- file.path(RAW_DIR, "final", "external_independent_predictions.csv")

read_one_scores <- function(path, set_name) {
  dat <- read_csv(path)
  id_col <- first_col(dat, c("sample_id", "sample_id_use", "sample_name"))
  study_col <- first_col(dat, c("study_name", "study_use", "study"))
  group_col <- first_col(dat, c("group", "Group", "group_norm"))
  system_col <- first_col(dat, c("system_use", "system", "System"))
  score_col <- first_col(dat, c("OMHI", "OMWI"))
  out <- data.frame(
    sample_id = trimws(as.character(dat[[id_col]])),
    study = trimws(as.character(dat[[study_col]])),
    group = clean_group(dat[[group_col]]),
    OMHI = suppressWarnings(as.numeric(dat[[score_col]])),
    set = set_name, stringsAsFactors = FALSE
  )
  out$system <- clean_system(dat[[system_col]], out$group)
  if (anyDuplicated(out$sample_id)) {
    stop(set_name, " 评分表中 sample_id 重复，请先检查，不应直接去重")
  }
  out <- out[is.finite(out$OMHI) & !is.na(out$group) & !is.na(out$system), ]
  out
}

train <- read_one_scores(train_file, "Training")
external <- read_one_scores(external_file, "External validation")

if (TRAIN_SCORE_MODE == "oof_mean") {
  oof_file <- file.path(RAW_DIR, "summary", "nested_repeated10fold_OOF.csv")
  oof <- read_csv(oof_file)
  stopifnot(all(c("sample_id", "repeat_id", "OMHI") %in% names(oof)))
  if (anyDuplicated(oof[c("sample_id", "repeat_id")])) {
    stop("OOF 存在重复的 sample_id + repeat_id")
  }
  n_rep <- as.data.table(oof)[, .N, by = sample_id]
  if (any(n_rep$N != 10L)) stop("OOF: 并非每个训练样本都有10次预测")
  avg <- as.data.table(oof)[, .(OMHI = mean(as.numeric(OMHI))), by = sample_id]
  mt <- match(train$sample_id, avg$sample_id)
  if (anyNA(mt)) stop("部分训练样本没有 OOF 分数")
  train$OMHI <- avg$OMHI[mt]
} else if (TRAIN_SCORE_MODE != "in_sample") {
  stop("TRAIN_SCORE_MODE 必须为 in_sample 或 oof_mean")
}
score_df <- dplyr::bind_rows(train, external)
if (anyDuplicated(score_df$sample_id)) {
  stop("训练集和外部集 sample_id 重叠，不允许重复计算")
}
message("模型评分样本：Training=", nrow(train),
        "; External validation=", nrow(external),
        "; training mode=", TRAIN_SCORE_MODE)
write.csv(score_df, file.path(OUT_DIR, "Fig3_model_scores_used.csv"), row.names = FALSE)

# ---------- 3. Raw species abundances -> Shannon/richness (for b/c) ----------
# Important: biological alpha diversity must use RAW abundance, not MMUPHin-
# adjusted or CLR-transformed values. The function rejects negative values.
find_species_files <- function(ids) {
  if (length(SPECIES_FILES)) return(vapply(SPECIES_FILES, check_input, ""))
  scan_dirs <- file.path(ROOT_DIR, c(
    "data/abundance", "02.by_levelsaliva64",
    "data/metadata", "01.split_saliva64", "01.split_saliva_only"
  ))
  scan_dirs <- scan_dirs[dir.exists(scan_dirs)]
  candidates <- unlist(lapply(scan_dirs, function(d) {
    list.files(d, pattern = "(species|sp_level).*\\.(tsv|csv|txt)(\\.gz)?$",
               recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  }), use.names = FALSE)
  # Deliberately exclude batch-adjusted, CLR or log-scaled tables.
  candidates <- candidates[!grepl("(batch\\.|mmuphin|\\bclr\\b|log2|log10)",
                                  candidates, ignore.case = TRUE)]
  if (!length(candidates)) return(character(0))
  found <- character(0)
  for (f in candidates) {
    peek <- tryCatch(data.table::fread(f, nrows = 80, data.table = FALSE,
                                       check.names = FALSE, showProgress = FALSE),
                     error = function(e) NULL)
    if (is.null(peek) || ncol(peek) < 2) next
    header_hit <- sum(names(peek) %in% ids)
    row_hit <- sum(as.character(peek[[1]]) %in% ids)
    if (header_hit >= 2 || row_hit >= 2) found <- c(found, f)
  }
  found
}

calculate_diversity <- function(mat) {
  mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"
  if (any(mat < 0, na.rm = TRUE)) {
    stop("物种丰度包含负数，可能选中了 CLR/批次校正表，请指定原始丰度表")
  }
  mat[!is.finite(mat)] <- 0
  totals <- colSums(mat)
  valid <- is.finite(totals) & totals > 0
  rel <- mat[, valid, drop = FALSE]
  rel <- sweep(rel, 2L, totals[valid], "/")
  shannon <- -colSums(ifelse(rel > 0, rel * log(pmax(rel, .Machine$double.xmin)), 0))
  richness <- colSums(mat[, valid, drop = FALSE] > 0)
  data.frame(sample_id = colnames(mat)[valid], Shannon = as.numeric(shannon),
             Richness = as.numeric(richness), stringsAsFactors = FALSE)
}

read_diversity_file <- function(path, ids) {
  dat <- read_csv(path)
  if (ncol(dat) < 2) stop("不是有效的物种丰度矩阵: ", path)
  # Format 1: features x samples, first column is species ID
  id_columns <- intersect(ids, names(dat))
  if (length(id_columns) >= 2) {
    feat <- as.character(dat[[1]])
    # MetaPhlAn all-rank tables: restrict to species if taxonomic pipes exist
    if (any(grepl("\\|s__", feat))) {
      keep <- grepl("\\|s__", feat) & !grepl("\\|t__", feat)
    } else if (any(grepl("^s__", feat))) {
      keep <- grepl("^s__", feat)
    } else {
      keep <- rep(TRUE, length(feat))
    }
    if (!any(keep)) stop("没有 species 层级数据: ", path)
    mat <- as.matrix(dat[keep, id_columns, drop = FALSE])
    # taxa x samples, colnames are original sample IDs
  } else {
    # Format 2: samples x features, first column is sample ID
    row_ids <- as.character(dat[[1]])
    keep_rows <- row_ids %in% ids
    if (sum(keep_rows) < 2) stop("文件不含足够的匹配 sample_id: ", path)
    taxa_cols <- setdiff(names(dat), names(dat)[1])
    if (any(grepl("\\|s__", taxa_cols))) {
      taxa_cols <- taxa_cols[grepl("\\|s__", taxa_cols) & !grepl("\\|t__", taxa_cols)]
    } else if (any(grepl("^s__", taxa_cols))) {
      taxa_cols <- taxa_cols[grepl("^s__", taxa_cols)]
    }
    mat <- t(as.matrix(dat[keep_rows, taxa_cols, drop = FALSE]))
    colnames(mat) <- row_ids[keep_rows]
  }
  if (anyDuplicated(colnames(mat))) stop("丰度矩阵 sample_id 重复: ", path)
  calculate_diversity(mat)
}

species_files <- find_species_files(score_df$sample_id)
if (!length(species_files)) {
  stop(paste0(
    "找不到可与评分 sample_id 匹配的原始 species 丰度矩阵。\n",
    "请在本脚本 SPECIES_FILES 中填入 Oral v3.0 下的实际文件路径；\n",
    "不可用 batch.normal / CLR 矩阵代替原始物种丰度来计算 Shannon。"
  ), call. = FALSE)
}
message("检测到 species 文件：\n", paste(species_files, collapse = "\n"))
diversity_parts <- lapply(species_files, read_diversity_file, ids = score_df$sample_id)
diversity_df <- dplyr::bind_rows(diversity_parts)
if (anyDuplicated(diversity_df$sample_id)) {
  # Input files can overlap (e.g. train/test combined and their split tables).
  # Check actual diversity metrics before retaining one per sample.
  consistency <- diversity_df %>% group_by(sample_id) %>%
    summarise(n_shannon = n_distinct(round(Shannon, 7)),
              n_richness = n_distinct(Richness), .groups = "drop")
  if (any(consistency$n_shannon > 1 | consistency$n_richness > 1)) {
    stop("不同 species 文件对相同 sample_id 计算出不同的多样性，请只指定最终原始丰度表")
  }
  diversity_df <- diversity_df[!duplicated(diversity_df$sample_id), ]
}
write.csv(diversity_df, file.path(OUT_DIR, "Fig3_shannon_richness_raw.csv"), row.names = FALSE)
score_div <- dplyr::left_join(score_df, diversity_df, by = "sample_id")
message("可用于 a/b/c 联合比较的样本数：", sum(complete.cases(score_div[c("OMHI", "Shannon", "Richness")])) )

# ------------------------- 4. Panels a/b/c ---------------------------------
abc <- score_div %>%
  filter(is.finite(OMHI), is.finite(Shannon), is.finite(Richness)) %>%
  group_by(study, group) %>% filter(n() >= MIN_STUDY_GROUP_N) %>% ungroup()
studies_with_both <- abc %>% group_by(study) %>%
  filter(any(group == "Healthy") & any(group != "Healthy")) %>%
  distinct(study) %>% pull(study)
abc <- abc %>% filter(study %in% studies_with_both)
if (!nrow(abc)) stop("a/b/c 没有满足样本量要求且同时含健康组和疾病组的研究")

rows <- abc %>% count(study, group, name = "n") %>%
  arrange(study, desc(group == "Healthy"), group) %>%
  mutate(row_id = n():1L,
         row_label = paste0(pretty_study(study), " ", group, " (", n, ")"))
abc <- abc %>% left_join(rows %>% select(study, group, row_id),
                         by = c("study", "group"))
# study-wide gray stripes are identical across a/b/c
study_bg <- rows %>% group_by(study) %>%
  summarise(ymin = min(row_id) - .5, ymax = max(row_id) + .5,
            .groups = "drop") %>%
  arrange(study) %>% mutate(band = row_number() %% 2 == 0)
study_bg <- study_bg %>% filter(band)

# Intra-study Wilcoxon tests, no cross-study pooling of different control groups.
study_wilcox <- function(dat, metric) {
  test_rows <- rows[rows$group != "Healthy", ]
  ans <- lapply(seq_len(nrow(test_rows)), function(k) {
    x <- test_rows[k, ]
    sub <- dat[dat$study == x$study, ]
    ref <- sub[[metric]][sub$group == "Healthy"]
    alt <- sub[[metric]][sub$group == x$group]
    p <- if (length(ref) >= MIN_STUDY_GROUP_N && length(alt) >= MIN_STUDY_GROUP_N)
      tryCatch(wilcox.test(alt, ref, exact = FALSE)$p.value,
               error = function(e) NA_real_) else NA_real_
    data.frame(study = x$study, group = x$group, row_id = x$row_id,
               metric = metric, n_healthy = length(ref), n_group = length(alt),
               p = p, star = stars(p))
  })
  do.call(rbind, ans)
}
abc_stats <- dplyr::bind_rows(lapply(c("OMHI", "Shannon", "Richness"),
                                     function(x) study_wilcox(abc, x)))
abc_stats <- abc_stats %>% group_by(metric) %>%
  mutate(FDR_BH = p.adjust(p, method = "BH")) %>% ungroup()
write.csv(abc_stats, file.path(OUT_DIR, "Fig3_abc_within_study_wilcoxon.csv"), row.names = FALSE)
write.csv(abc, file.path(OUT_DIR, "Fig3_abc_plot_data.csv"), row.names = FALSE)

make_abc <- function(metric, panel_tag, xlabel, show_y = FALSE) {
  xmin <- min(abc[[metric]], na.rm = TRUE)
  xmax <- max(abc[[metric]], na.rm = TRUE)
  ran <- max(xmax - xmin, 1e-4)
  anno_x <- xmax + .045 * ran
  stats <- abc_stats %>% filter(metric == !!metric, star != "") %>%
    mutate(x = anno_x)
  p <- ggplot(abc, aes(x = .data[[metric]], y = row_id,
                       group = row_id, fill = system)) +
    geom_rect(data = study_bg,
              aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
              inherit.aes = FALSE, fill = "grey92", colour = NA) +
    geom_boxplot(width = .58, orientation = "y", outlier.shape = NA, alpha = .77,
                 colour = "grey25", linewidth = .38) +
    geom_point(position = position_jitter(width = 0, height = .12),
               size = .38, alpha = .30, colour = "grey25", show.legend = FALSE) +
    geom_text(data = stats, aes(x = x, y = row_id, label = star),
              inherit.aes = FALSE, size = 2.6, hjust = 0) +
    scale_fill_manual(values = SYSTEM_COLORS, na.value = "grey70") +
    scale_y_continuous(breaks = rows$row_id,
                       labels = if (show_y) rows$row_label else rep("", nrow(rows)),
                       expand = expansion(add = c(.5, .5))) +
    scale_x_continuous(breaks = scales::pretty_breaks(4),
                       expand = expansion(mult = c(.05, .20))) +
    labs(x = xlabel, y = NULL, tag = panel_tag) +
    base_theme(8.7) +
    theme(axis.text.y = if (show_y) element_text(size = 7.15)
          else element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.x = element_text(size = 9.4),
          plot.tag = element_text(size = 15, face = "bold"),
          plot.tag.position = c(0.01, 1.045))
  if (metric == "OMHI") p <- p +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50", linewidth = .42)
  p
}
p_a <- make_abc("OMHI", "a", "Oral microbiome wellness index (OMWI)", TRUE)
p_b <- make_abc("Shannon", "b", "Shannon diversity")
p_c <- make_abc("Richness", "c", "Species richness (no. of species)")
p_top <- wrap_plots(p_a, p_b, p_c, nrow = 1, widths = c(3.7, 1.34, 1.55))

# ------------------------ 5. Panels d/e -------------------------------------
GROUP_PRIORITY <- c(
  "Healthy", "Obesity", "Overweight", "Underweight", "Type 1 diabetes",
  "Type 2 diabetes", "Adenoid Cystic Carcinoma", "Periodontitis",
  "Rheumatoid Arthritis", "COVID", "Pancreatic cancer",
  "Atrophic gastritis", "Colorectal cancer", "Superficial gastritis",
  "Intestinal metaplasia", "Pancreatitis", "Colorectal polyps"
)
short_system_label <- function(x) {
  z <- x
  z[x == "Abnormal BMI"] <- "Abnormal\nBMI"
  z[x == "Metabolic Disease"] <- "Metabolic\nDisease"
  z[x == "Digestive Disease"] <- "Digestive\nDisease"
  z[x == "Respiratory Disease"] <- "Respiratory\nDisease"
  z[x == "Immune Disease"] <- "Immune\nDisease"
  z[x == "Oral Disease"] <- "Oral\nDisease"
  z
}

compare_healthy <- function(data, variable) {
  ref <- data$OMHI[as.character(data[[variable]]) == "Healthy"]
  groups <- setdiff(unique(as.character(data[[variable]])), "Healthy")
  if (!length(groups)) return(data.frame(category = character(), p = numeric(), star = character()))
  res <- lapply(groups, function(g) {
    alt <- data$OMHI[as.character(data[[variable]]) == g]
    p <- if (length(ref) >= MIN_DE_GROUP_N && length(alt) >= MIN_DE_GROUP_N)
      tryCatch(wilcox.test(ref, alt, exact = FALSE)$p.value,
               error = function(e) NA_real_) else NA_real_
    data.frame(category = g, n_healthy = length(ref), n_group = length(alt), p = p)
  })
  ret <- do.call(rbind, res)
  ret$FDR_BH <- p.adjust(ret$p, method = "BH")
  ret$star <- stars(ret$p)
  ret
}

make_de <- function(data, panel_tag, panel_title) {
  dat <- data %>% filter(is.finite(OMHI))
  keep_group <- dat %>% count(group) %>% filter(n >= MIN_DE_GROUP_N) %>% pull(group)
  dat <- dat %>% filter(group %in% keep_group)
  if (!"Healthy" %in% dat$group) stop(panel_title, " 缺少 Healthy 组")
  # Each phenotype's predominant system (there should be one unambiguous mapping).
  group_map <- dat %>% count(group, system, sort = TRUE) %>%
    arrange(group, desc(n), system) %>% group_by(group) %>% slice(1) %>% ungroup()
  dat <- dat %>% left_join(group_map %>% select(group, main_system = system), by = "group")
  # DISPLAY ORDER ONLY (d/e): Healthy first; then systems by median OMWI
  # descending, and phenotypes sorted within each system by median OMWI.
  # Leave the original sample selection, scores, statistics and a/b/c untouched.
  sys_counts <- dat %>% count(system) %>% filter(n >= MIN_DE_GROUP_N)
  sys_order <- intersect(SYSTEM_ORDER, sys_counts$system)
  dat <- dat %>% filter(main_system %in% sys_order)
  
  # Match the system summary used by p_sys: preserve original sample membership.
  sys_medians <- dat %>%
    filter(system %in% sys_order) %>%
    group_by(system) %>%
    summarise(med_OMWI = median(OMHI), .groups = "drop") %>%
    arrange(desc(med_OMWI), system)
  sys_order <- c(intersect("Healthy", sys_medians$system),
                 setdiff(sys_medians$system, "Healthy"))
  
  # Match the original phenotype-to-system mapping, changing factor order only.
  group_order <- dat %>%
    group_by(main_system, group) %>%
    summarise(med_OMWI = median(OMHI), .groups = "drop") %>%
    mutate(system_rank = match(main_system, sys_order)) %>%
    arrange(system_rank, desc(med_OMWI), group) %>%
    pull(group)
  dat$grp_y <- factor(dat$group, levels = rev(group_order))
  sysdat <- dat %>% filter(system %in% sys_order)
  sysdat$sys_y <- factor(sysdat$system, levels = rev(sys_order))
  group_cols <- SYSTEM_COLORS[group_map$system[match(group_order, group_map$group)]]
  names(group_cols) <- group_order
  sys_cols <- SYSTEM_COLORS[sys_order]
  
  st_g <- compare_healthy(dat, "group")
  st_s <- compare_healthy(sysdat, "system")
  st_g$dataset <- panel_title; st_g$type <- "phenotype"
  st_s$dataset <- panel_title; st_s$type <- "system"
  stats <- dplyr::bind_rows(st_g, st_s)
  write.csv(stats, file.path(OUT_DIR,
                             paste0("Fig3_", panel_tag, "_wilcoxon.csv")), row.names = FALSE)
  
  x_right <- if (is.null(OMWI_XLIMIT)) max(dat$OMHI, na.rm = TRUE) else OMWI_XLIMIT[2]
  gs <- st_g[st_g$star != "", , drop = FALSE]
  ss <- st_s[st_s$star != "", , drop = FALSE]
  gs$grp_y <- factor(gs$category, levels = rev(group_order))
  ss$sys_y <- factor(ss$category, levels = rev(sys_order))
  
  p_group <- ggplot(dat, aes(x = OMHI, y = grp_y, fill = group)) +
    geom_violin(trim = FALSE, scale = "width", colour = NA, alpha = .78) +
    geom_boxplot(width = .16, outlier.shape = NA, fill = NA,
                 colour = "black", linewidth = .44) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55", linewidth = .52) +
    geom_text(data = gs, aes(x = x_right - .1, y = grp_y, label = star),
              inherit.aes = FALSE, hjust = 1, size = 3.3) +
    scale_fill_manual(values = group_cols) +
    labs(x = "OMWI", y = NULL, tag = panel_tag) + base_theme(9.0) +
    theme(axis.text.y = element_text(size = if (panel_tag == "d") 8 else 8.5),
          axis.title.x = element_text(size = 9.4),
          plot.tag = element_text(size = 15, face = "bold"),
          plot.tag.position = c(-.02, 1.085))
  if (!is.null(OMWI_XLIMIT)) p_group <- p_group + coord_cartesian(xlim = OMWI_XLIMIT)
  
  p_sys <- ggplot(sysdat, aes(x = OMHI, y = sys_y, fill = system)) +
    geom_violin(trim = FALSE, scale = "width", colour = NA, alpha = .78) +
    geom_boxplot(width = .16, outlier.shape = NA, fill = NA,
                 colour = "black", linewidth = .44) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55", linewidth = .52) +
    geom_text(data = ss, aes(x = x_right - .1, y = sys_y, label = star),
              inherit.aes = FALSE, hjust = 1, size = 3.3) +
    scale_fill_manual(values = sys_cols) +
    labs(x = "OMWI", y = NULL) + base_theme(9.0) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.title.x = element_text(size = 9.4))
  if (!is.null(OMWI_XLIMIT)) p_sys <- p_sys + coord_cartesian(xlim = OMWI_XLIMIT)
  
  strip <- data.frame(system = sys_order, y = factor(sys_order, levels = rev(sys_order)),
                      label = short_system_label(sys_order), stringsAsFactors = FALSE)
  p_strip <- ggplot(strip, aes(x = 1, y = y)) +
    geom_tile(aes(fill = system), width = .98, height = .96, colour = "white", linewidth = .2) +
    geom_text(aes(label = label), size = if (panel_tag == "d") 2.7 else 2.35,
              fontface = "bold", lineheight = .93) +
    scale_fill_manual(values = sys_cols) +
    coord_cartesian(xlim = c(.5, 1.5), clip = "off") +
    theme_void() + theme(legend.position = "none", plot.margin = margin(4, 2, 5, 0))
  
  out <- wrap_plots(p_group, p_sys, p_strip, nrow = 1,
                    widths = c(1.28, 1.07, .57)) +
    plot_annotation(title = panel_title,
                    theme = theme(plot.title = element_text(hjust = .5, size = 12.5)))
  out
}

p_d <- make_de(train, "d", "Training set")
p_e <- make_de(external, "e", "External validation set")

# ------------------------------ 6. Panel f ----------------------------------
# Use same independent (frozen) model scores as panels d/e.
# Positive class = disease, predictor = -OMWI; one-versus-healthy per system.
f_stats <- list()
for (ss in c("Training", "External validation")) {
  dat <- score_df[score_df$set == ss, ]
  for (sy in DISEASE_SYSTEMS) {
    sub <- dat[dat$system %in% c("Healthy", sy), ]
    if (sum(sub$system == sy) < MIN_F_CASE_N ||
        sum(sub$system == "Healthy") < MIN_F_CASE_N) next
    y <- as.integer(sub$system == sy)
    roc_obj <- pROC::roc(response = y, predictor = -sub$OMHI,
                         levels = c(0, 1), direction = "<", quiet = TRUE)
    ci <- as.numeric(pROC::ci.auc(roc_obj, conf.level = .95))
    f_stats[[length(f_stats) + 1L]] <- data.frame(
      set = ss, system = sy, n_case = sum(y == 1),
      n_healthy = sum(y == 0), AUC = as.numeric(pROC::auc(roc_obj)),
      ci_low = ci[1], ci_high = ci[3]
    )
  }
}
if (!length(f_stats)) stop("f：找不到足够的病例/健康对照来计算系统级 AUC")
f_df <- dplyr::bind_rows(f_stats) %>%
  group_by(set) %>% arrange(desc(AUC), .by_group = TRUE) %>% ungroup()
write.csv(f_df, file.path(OUT_DIR, "Fig3_f_system_AUC_from_raw_scores.csv"), row.names = FALSE)
f_df <- f_df %>%
  mutate(set = factor(set, levels = c("Training", "External validation")),
         system_id = paste(as.character(set), system, sep = "___"),
         label = sprintf("%.2f\n(%.2f–%.2f)", AUC, ci_low, ci_high))
level_order <- f_df$system_id
f_df$system_id <- factor(f_df$system_id, levels = unique(level_order))
axis_labels <- setNames(gsub(" Disease", " disease", f_df$system, fixed = TRUE),
                        as.character(f_df$system_id))

p_f <- ggplot(f_df, aes(x = system_id, y = AUC, colour = system)) +
  geom_hline(yintercept = .5, linetype = "dashed", colour = "grey55", linewidth = .45) +
  geom_tile(aes(y = .511, fill = system), width = .73, height = .014,
            inherit.aes = TRUE, colour = NA, show.legend = FALSE) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = .10,
                linewidth = .72) +
  geom_point(size = 2.6) +
  geom_text(aes(label = label), vjust = -0.7, colour = "black",
            size = 2.0, lineheight = .90, show.legend = FALSE) +
  facet_wrap(~set, nrow = 1, scales = "free_x") +
  scale_color_manual(values = SYSTEM_COLORS) +
  scale_fill_manual(values = SYSTEM_COLORS) +
  scale_x_discrete(labels = axis_labels) +
  scale_y_continuous(breaks = seq(.5, 1, .1)) +
  coord_cartesian(ylim = c(.5, 1.06), clip = "off") +
  labs(x = NULL, y = "AUC", tag = "f") + base_theme(9) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 7.3),
        axis.text.y = element_text(size = 8),
        strip.background = element_blank(),
        strip.text = element_text(size = 10.3, face = "bold"),
        plot.tag = element_text(size = 15, face = "bold"),
        plot.tag.position = c(-.025, 1.085),
        plot.margin = margin(5, 5, 14, 5))

# ----------------------- 7. Assemble Fig. 3 a–f ------------------------------
# Matches the supplied fig3 PDF: top is a/b/c; bottom is d and (e above f).
p_right <- wrap_plots(p_e, p_f, ncol = 1, heights = c(1.0, 1.16))
p_bottom <- wrap_plots(p_d, p_right, nrow = 1, widths = c(1.13, .94))
fig3 <- patchwork::wrap_plots(
  p_top, p_bottom,
  ncol = 1,
  heights = c(1.10, 1.0)
)

# 给整张图设置白色背景
fig3 <- fig3 +
  patchwork::plot_annotation(
    theme = ggplot2::theme(
      plot.background = ggplot2::element_rect(
        fill = "white",
        colour = NA
      )
    )
  )

# 显示完整组图
print(fig3)

# 保存 PDF 和 PNG
save_both(
  fig3,
  "Fig3_OMWI_all_panels_abcdef",
  EXPORT_WIDTH,
  EXPORT_HEIGHT
)

save_both(fig3, "Fig3_OMWI_all_panels_abcdef", EXPORT_WIDTH, EXPORT_HEIGHT)
# Individual panels for future layout edits; CSVs already saved above.
save_both(p_top, "Fig3_abc_by_study", 17.4, 8.0)
save_both(p_d, "Fig3_d_training", 9.1, 8.1)
save_both(p_e, "Fig3_e_external", 7.6, 4.6)
save_both(p_f, "Fig3_f_system_AUC", 7.6, 4.6)
message("全部完成；所有输出在：", OUT_DIR)
