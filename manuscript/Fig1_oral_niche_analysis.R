#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(ggplot2)
  library(vegan)
  library(scales)
  library(grid)
})

options(stringsAsFactors = FALSE)
set.seed(704)

# =========================================================
# 0. 路径与参数
# =========================================================
INPUT_DIR <- file.path(Sys.getenv("OMWI_ROOT", getwd()), "data/all_oral")
LEVEL_DIR <- file.path(INPUT_DIR, "level")

# metadata 允许是 csv / tsv / txt 或无扩展名。
META_CANDIDATES <- c(
  file.path(INPUT_DIR, "all_oral_metadata.csv")
)
META_FILE <- META_CANDIDATES[file.exists(META_CANDIDATES)][1]

ABUND_FILE <- file.path(LEVEL_DIR, "taxonomy_species.tsv")
OUT_DIR    <- file.path(
  file.path(Sys.getenv("OMWI_ROOT", getwd()), "fig1"),
  "pcoa_oral_niche"
)

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

PERM_N <- 999
SEED   <- 704

SITE_ORDER <- c(
  "Dental",
  "Saliva",
  "Tongue",
  "Pharyngeal swab",
  "Oral rinse"
)

SITE_COLORS <- c(
  "Dental"          = "#E97863",
  "Saliva"          = "#73ADD3",
  "Tongue"          = "#B18BD5",
  "Pharyngeal swab" = "#72B77B",
  "Oral rinse"      = "#D9AD32"
)

SAMPLE_ID_CANDIDATES <- c(
  "sample_id", "SampleID", "sample", "Sample",
  "sample_name", "Sample_Name", "sample_id_use",
  "run", "Run", "run_id", "ID"
)

SITE_COL_CANDIDATES <- c("body_site")

if (length(META_FILE) == 0 || is.na(META_FILE)) {
  stop(
    "在以下目录中找不到 metadata：\n", INPUT_DIR, "\n",
    "已尝试：\n", paste(META_CANDIDATES, collapse = "\n")
  )
}

if (!file.exists(ABUND_FILE)) {
  stop("找不到物种丰度表：\n", ABUND_FILE)
}

cat("Metadata: ", META_FILE, "\n", sep = "")
cat("Abundance: ", ABUND_FILE, "\n", sep = "")
cat("Output: ", OUT_DIR, "\n\n", sep = "")

# =========================================================
# 1. 通用函数
# =========================================================
pick_first_col <- function(df, candidates, label) {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) == 0) {
    stop(
      label, "列未找到。\n",
      "当前文件列名：\n", paste(colnames(df), collapse = ", "), "\n\n",
      "候选列名：\n", paste(candidates, collapse = ", ")
    )
  }
  hit[1]
}

standardize_site <- function(x) {
  x0 <- trimws(as.character(x))
  x1 <- str_to_lower(x0)
  
  case_when(
    str_detect(x1, "oral[ _-]*rinse|mouth[ _-]*rinse|mouthwash|oral wash|rinse") ~ "Oral rinse",
    str_detect(x1, "pharyn|throat") ~ "Pharyngeal swab",
    str_detect(x1, "tongue") ~ "Tongue",
    str_detect(x1, "saliva|salivary") ~ "Saliva",
    str_detect(x1, "dental|tooth|teeth|plaque|gingiv|periodontal|subging|supraging") ~ "Dental",
    x1 == "oral rinse" ~ "Oral rinse",
    x1 == "pharyngeal swab" ~ "Pharyngeal swab",
    x1 == "tongue" ~ "Tongue",
    x1 == "saliva" ~ "Saliva",
    x1 == "dental" ~ "Dental",
    TRUE ~ NA_character_
  )
}

format_p_plot <- function(p) {
  if (is.na(p)) return("NA")
  # 999 permutations 的最小可报告 P 值为 0.001。
  if (p <= 0.001) return("0.001")
  sprintf("%.3f", p)
}

# =========================================================
# 2. 读取 metadata
# =========================================================
meta_raw <- fread(
  META_FILE,
  data.table = FALSE,
  check.names = FALSE,
  encoding = "UTF-8"
)

sample_col <- pick_first_col(meta_raw, SAMPLE_ID_CANDIDATES, "样本 ID")
site_col   <- pick_first_col(meta_raw, SITE_COL_CANDIDATES, "口腔位点")

meta <- meta_raw %>%
  transmute(
    sample_id = trimws(as.character(.data[[sample_col]])),
    Site_raw  = trimws(as.character(.data[[site_col]])),
    Site      = standardize_site(.data[[site_col]])
  ) %>%
  filter(!is.na(sample_id), sample_id != "") %>%
  distinct(sample_id, .keep_all = TRUE)

cat("识别到的样本 ID 列：", sample_col, "\n", sep = "")
cat("识别到的口腔位点列：", site_col, "\n\n", sep = "")

cat("metadata 原始位点取值及数量：\n")
print(sort(table(meta$Site_raw, useNA = "ifany"), decreasing = TRUE))

unknown_site <- meta %>% filter(is.na(Site))
if (nrow(unknown_site) > 0) {
  cat("\n未纳入七类口腔位点的原始取值：\n")
  print(sort(table(unknown_site$Site_raw, useNA = "ifany"), decreasing = TRUE))
}

meta <- meta %>%
  filter(!is.na(Site)) %>%
  mutate(Site = factor(Site, levels = SITE_ORDER))

# =========================================================
# 3. 读取 taxonomy_species.tsv，并自动判断矩阵方向
# =========================================================
abund_raw <- fread(
  ABUND_FILE,
  data.table = FALSE,
  check.names = FALSE,
  encoding = "UTF-8"
)

if (nrow(abund_raw) < 2 || ncol(abund_raw) < 2) {
  stop("taxonomy_species.tsv 的行列数异常。")
}

meta_ids <- meta$sample_id
col_overlap <- intersect(trimws(colnames(abund_raw)), meta_ids)
first_col_values <- trimws(as.character(abund_raw[[1]]))
row_overlap <- intersect(first_col_values, meta_ids)

if (length(col_overlap) >= 20) {
  # 常见格式：第一列为物种，后续列为样本。
  taxon_col <- colnames(abund_raw)[1]
  sample_cols <- col_overlap
  
  abund_sub <- abund_raw[, c(taxon_col, sample_cols), drop = FALSE]
  colnames(abund_sub)[1] <- "Taxon"
  
  for (j in sample_cols) {
    abund_sub[[j]] <- suppressWarnings(as.numeric(abund_sub[[j]]))
  }
  
  abund_dt <- as.data.table(abund_sub)
  abund_dt <- abund_dt[
    , lapply(.SD, function(z) sum(z, na.rm = TRUE)),
    by = Taxon,
    .SDcols = sample_cols
  ]
  
  taxa <- make.unique(as.character(abund_dt$Taxon))
  mat <- as.matrix(abund_dt[, sample_cols, with = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- taxa
  colnames(mat) <- sample_cols
  
  cat("丰度表方向：物种为行、样本为列。\n")
  
} else if (length(row_overlap) >= 20) {
  # 兼容格式：第一列为样本，后续列为物种。
  sample_id_col <- colnames(abund_raw)[1]
  use_rows <- first_col_values %in% meta_ids
  tmp <- abund_raw[use_rows, , drop = FALSE]
  
  sample_ids <- trimws(as.character(tmp[[sample_id_col]]))
  feature_df <- tmp[, -1, drop = FALSE]
  
  for (j in seq_len(ncol(feature_df))) {
    feature_df[[j]] <- suppressWarnings(as.numeric(feature_df[[j]]))
  }
  
  mat <- t(as.matrix(feature_df))
  storage.mode(mat) <- "numeric"
  rownames(mat) <- make.unique(colnames(feature_df))
  colnames(mat) <- sample_ids
  
  cat("丰度表方向：样本为行、物种为列，已自动转置。\n")
  
} else {
  stop(
    "taxonomy_species.tsv 与 metadata 的样本 ID 无法匹配。\n",
    "丰度表列名与 metadata 的交集数：", length(col_overlap), "\n",
    "丰度表第一列与 metadata 的交集数：", length(row_overlap), "\n",
    "请检查 metadata 样本 ID 与 taxonomy_species.tsv 样本名是否一致。"
  )
}

# =========================================================
# 4. 对齐样本并标准化相对丰度
# =========================================================
common_samples <- intersect(colnames(mat), meta$sample_id)

if (length(common_samples) < 20) {
  stop("丰度表与 metadata 的共同样本过少：", length(common_samples))
}

mat <- mat[, common_samples, drop = FALSE]
meta_use <- meta[match(common_samples, meta$sample_id), , drop = FALSE]

# 去除全零样本和全零物种。
keep_sample <- colSums(mat, na.rm = TRUE) > 0
mat <- mat[, keep_sample, drop = FALSE]
meta_use <- meta_use[keep_sample, , drop = FALSE]

keep_taxon <- rowSums(mat > 0, na.rm = TRUE) > 0
mat <- mat[keep_taxon, , drop = FALSE]

# 无论原表是百分数还是比例，最终都按每个样本归一化到 1。
mat[is.na(mat)] <- 0
mat <- sweep(mat, 2, colSums(mat, na.rm = TRUE), FUN = "/")
mat[!is.finite(mat)] <- 0

# 严格保持 metadata 与矩阵列顺序一致。
meta_use <- meta_use[match(colnames(mat), meta_use$sample_id), , drop = FALSE]
meta_use$Site <- factor(meta_use$Site, levels = SITE_ORDER)

site_counts <- meta_use %>%
  count(Site, name = "n", .drop = FALSE) %>%
  filter(n > 0)

cat("\n匹配丰度表后的七个位点样本量：\n")
print(site_counts)

# 两个新位点的样本量核对：以丰度表成功匹配、且非全零的样本为准。
new_site_counts <- site_counts %>%
  filter(as.character(Site) %in% c("Hard palate", "Mucosa"))
cat("\n新增位点实际纳入 PCoA 的样本量：\n")
print(new_site_counts)
expected_new_n <- c("Hard palate" = 19L, "Mucosa" = 4L)
for (site_name in names(expected_new_n)) {
  actual_n <- site_counts$n[as.character(site_counts$Site) == site_name]
  if (length(actual_n) == 1L && actual_n != expected_new_n[[site_name]]) {
    warning(site_name, ": 预期 ", expected_new_n[[site_name]],
            " 个，实际成功匹配 ", actual_n,
            " 个；请检查丰度表匹配及全零样本过滤。")
  }
}

missing_sites <- setdiff(SITE_ORDER, as.character(site_counts$Site))
if (length(missing_sites) > 0) {
  stop(
    "匹配后仍缺少位点：", paste(missing_sites, collapse = ", "), "\n",
    "请检查 metadata 中该位点样本是否存在于 taxonomy_species.tsv。"
  )
}

write.csv(
  site_counts,
  file.path(OUT_DIR, "oral_site_sample_counts.csv"),
  row.names = FALSE
)

# =========================================================
# 5. Bray-Curtis PCoA
# =========================================================
cat("\n开始计算 Bray-Curtis 距离和 PCoA……\n")
bray <- vegdist(t(mat), method = "bray")

# add = TRUE 用于处理 Bray-Curtis 可能出现的负特征值。
ord <- cmdscale(bray, k = 2, eig = TRUE, add = TRUE)

positive_eig <- ord$eig[ord$eig > 0]
eig_pct <- 100 * ord$eig[1:2] / sum(positive_eig)

pcoa_df <- as.data.frame(ord$points)
colnames(pcoa_df) <- c("PCoA1", "PCoA2")
pcoa_df$sample_id <- rownames(pcoa_df)
pcoa_df <- pcoa_df %>%
  left_join(meta_use %>% select(sample_id, Site), by = "sample_id")

# 统一 PCoA 方向，使 Dental 位于右侧。
dental_x <- mean(pcoa_df$PCoA1[pcoa_df$Site == "Dental"], na.rm = TRUE)
proxy_x  <- mean(pcoa_df$PCoA1[pcoa_df$Site != "Dental"], na.rm = TRUE)
if (is.finite(dental_x) && is.finite(proxy_x) && dental_x < proxy_x) {
  pcoa_df$PCoA1 <- -pcoa_df$PCoA1
}

# 进一步统一纵轴方向：Tongue/Oral rinse 倾向位于 Pharyngeal swab 上方。
upper_y <- mean(
  pcoa_df$PCoA2[pcoa_df$Site %in% c("Tongue", "Oral rinse")],
  na.rm = TRUE
)
lower_y <- mean(
  pcoa_df$PCoA2[pcoa_df$Site == "Pharyngeal swab"],
  na.rm = TRUE
)
if (is.finite(upper_y) && is.finite(lower_y) && upper_y < lower_y) {
  pcoa_df$PCoA2 <- -pcoa_df$PCoA2
}

# 总体 PERMANOVA。
set.seed(SEED)
ad_overall <- adonis2(bray ~ Site, data = meta_use, permutations = PERM_N)
overall_r2 <- as.numeric(ad_overall[1, "R2"])
overall_p  <- as.numeric(ad_overall[1, "Pr(>F)"])

write.csv(
  as.data.frame(ad_overall),
  file.path(OUT_DIR, "PERMANOVA_overall_oral_site.csv"),
  row.names = TRUE
)

# 各位点文字位置：基于每个位点的中位数
label_df <- pcoa_df %>%
  group_by(Site) %>%
  summarise(
    x = median(PCoA1, na.rm = TRUE),
    y = median(PCoA2, na.rm = TRUE),
    .groups = "drop"
  )

xr <- range(pcoa_df$PCoA1, na.rm = TRUE)
yr <- range(pcoa_df$PCoA2, na.rm = TRUE)
xd <- diff(xr)
yd <- diff(yr)

label_df <- label_df %>%
  mutate(
    x = x + case_when(
      Site == "Dental"          ~  0.00 * xd,
      Site == "Oral rinse"      ~ -0.03 * xd,
      Site == "Tongue"          ~ -0.04 * xd,
      Site == "Saliva"          ~  0.04 * xd,
      Site == "Pharyngeal swab" ~ -0.02 * xd,
      TRUE ~ 0
    ),
    y = y + case_when(
      Site == "Dental"          ~  0.01 * yd,
      Site == "Oral rinse"      ~  0.08 * yd,
      Site == "Tongue"          ~  0.02 * yd,
      Site == "Saliva"          ~  0.01 * yd,
      Site == "Pharyngeal swab" ~ -0.03 * yd,
      TRUE ~ 0
    )
  )

p_pcoa <- ggplot(
  pcoa_df,
  aes(x = PCoA1, y = PCoA2, color = Site)
) +
  geom_point(size = 3, alpha = 0.62, stroke = 0) +
  geom_label(
    data = label_df,
    aes(x = x, y = y, label = Site),
    inherit.aes = FALSE,
    fill = "white",
    color = "black",
    label.size = 0.38,
    label.padding = unit(0.12, "lines"),
    fontface = "bold",
    size = 5
  ) +
  annotate(
    "text",
    x = xr[2] - 0.02 * xd,
    y = yr[2] - 0.02 * yd,
    label = paste0(
      "R² = ", sprintf("%.3f", overall_r2), "\n",
      "P = ", format_p_plot(overall_p)
    ),
    hjust = 1,
    vjust = 1,
    size = 8,
    fontface = "bold",
    color = "black"
  ) +
  scale_color_manual(
    values = SITE_COLORS,
    breaks = SITE_ORDER,
    drop = FALSE
  ) +
  labs(
    x = paste0("PCoA1 (", sprintf("%.2f", eig_pct[1]), "%)"),
    y = paste0("PCoA2 (", sprintf("%.2f", eig_pct[2]), "%)"),
    color = NULL
  ) +
  theme_bw(base_size = 18) +
  theme(
    panel.grid.major = element_line(color = "grey90", linewidth = 0.35),
    panel.grid.minor = element_line(color = "grey95", linewidth = 0.25),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.15),
    axis.title = element_text(face = "bold", color = "black"),
    axis.text = element_text(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.position = "right",
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.text = element_text(size = 8.5, color = "black"),
    plot.margin = margin(7, 7, 7, 7)
  )
p_pcoa
# =========================================================
# 6. 两两 PERMANOVA R² 下三角热图
# =========================================================
pair_list <- combn(SITE_ORDER, 2, simplify = FALSE)
bray_mat <- as.matrix(bray)

pairwise_res <- lapply(pair_list, function(pair_i) {
  keep <- meta_use$Site %in% pair_i
  
  # 用完整距离矩阵按样本下标截取，避免重新计算物种距离。
  sub_dist <- as.dist(bray_mat[keep, keep, drop = FALSE])
  sub_meta <- droplevels(meta_use[keep, , drop = FALSE])
  
  set.seed(SEED)
  ad <- adonis2(sub_dist ~ Site, data = sub_meta, permutations = PERM_N)
  
  data.frame(
    Site1 = pair_i[1],
    Site2 = pair_i[2],
    R2 = as.numeric(ad[1, "R2"]),
    P  = as.numeric(ad[1, "Pr(>F)"]),
    n1 = sum(sub_meta$Site == pair_i[1]),
    n2 = sum(sub_meta$Site == pair_i[2]),
    stringsAsFactors = FALSE
  )
}) %>%
  bind_rows()

write.csv(
  pairwise_res,
  file.path(OUT_DIR, "PERMANOVA_pairwise_oral_site.csv"),
  row.names = FALSE
)

# 删除下三角中完全空白的 Dental 行及最后一列（Mucosa）。
# 行：Saliva → Mucosa；列：Dental → Hard palate。
heat_rows <- SITE_ORDER[-1]
heat_cols <- SITE_ORDER[-length(SITE_ORDER)]

heat_df <- pairwise_res %>%
  transmute(
    row_site = Site2,
    col_site = Site1,
    R2 = R2,
    P = P
  ) %>%
  mutate(
    row_site = factor(row_site, levels = rev(heat_rows)),
    col_site = factor(col_site, levels = heat_cols),
    label = sprintf("%.3f", R2)
  )

count_map <- setNames(site_counts$n, as.character(site_counts$Site))
row_label_map <- setNames(
  paste0(heat_rows, "\n(n=", count_map[heat_rows], ")"),
  heat_rows
)

# 删除两端的空白行/列，下三角热图紧凑显示。

heat_max <- max(heat_df$R2, na.rm = TRUE)
heat_limit <- max(0.20, ceiling(heat_max * 100) / 100)

p_heat <- ggplot(heat_df, aes(x = col_site, y = row_site)) +
  geom_tile(
    aes(fill = R2),
    color = "white",
    linewidth = 1.0,
    width = 0.96,
    height = 0.96
  ) +
  geom_text(
    aes(label = label),
    size = 3.5,
    fontface = "bold",
    color = "black"
  ) +
  scale_fill_gradient(
    low = "#FFF5F0",
    high = "#EF3B2C",
    limits = c(0, heat_limit),
    oob = squish,
    name = expression(R^2)
  ) +
  scale_x_discrete(
    limits = heat_cols,
    drop = FALSE,
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_discrete(
    limits = rev(heat_rows),
    labels = row_label_map,
    drop = FALSE,
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  coord_fixed(clip = "off") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 10.5) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      face = "bold",
      color = "black"
    ),
    axis.text.y = element_text(
      face = "bold",
      color = "black",
      lineheight = 0.95
    ),
    axis.ticks = element_blank(),
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(color = "black"),
    legend.key.height = unit(2.0, "cm"),
    plot.margin = margin(8, 8, 8, 8)
  )
p_heat

# =========================================================
# 6.1 完整 7 × 7 下三角热图
# =========================================================

# 保留全部七个位点作为横轴和纵轴
heat_rows <- SITE_ORDER
heat_cols <- SITE_ORDER

# pairwise_res 中 Site1 总是排序靠前的位点，
# Site2 总是排序靠后的位点，因此映射后正好形成下三角。
heat_df <- pairwise_res %>%
  transmute(
    row_site = Site2,
    col_site = Site1,
    R2 = R2,
    P = P
  ) %>%
  mutate(
    col_site = factor(col_site, levels = heat_cols),
    
    # ggplot 的第一个 y 轴水平位于最下方，
    # 所以反转 levels，使 Dental 位于最上方。
    row_site = factor(row_site, levels = rev(heat_rows)),
    
    label = sprintf("%.3f", R2)
  )

# 检查：七个位点应该有 21 个唯一的两两比较
expected_pair_n <- choose(length(SITE_ORDER), 2)

stopifnot(nrow(pairwise_res) == expected_pair_n)
stopifnot(nrow(heat_df) == expected_pair_n)

message(
  "Pairwise comparisons retained: ",
  nrow(heat_df), "/", expected_pair_n
)

# =========================================================
# 样本量标签
# =========================================================
count_map <- setNames(
  site_counts$n,
  as.character(site_counts$Site)
)

row_label_map <- setNames(
  paste0(
    SITE_ORDER,
    "\n(n=", count_map[SITE_ORDER], ")"
  ),
  SITE_ORDER
)

# =========================================================
# 构建完整的 7 × 7 空白网格
# 作用：即使 Dental 行和 Mucosa 列没有色块，
# 也仍然完整保留在坐标轴中。
# =========================================================
full_grid <- tidyr::expand_grid(
  col_site = factor(SITE_ORDER, levels = SITE_ORDER),
  row_site = factor(SITE_ORDER, levels = rev(SITE_ORDER))
)

# =========================================================
# 对角线虚线
#
# 离散坐标内部的位置为 1、2、3、4、5、6、7。
# 这条线从矩阵左上角延伸到右下角，
# 穿过五个对角单元格的中心。
# =========================================================
n_site <- length(SITE_ORDER)

diag_df <- data.frame(
  x     = 0.5,
  y     = n_site + 0.5,
  xend  = n_site + 0.5,
  yend  = 0.5
)

# =========================================================
# 色阶范围
# =========================================================
heat_max <- max(heat_df$R2, na.rm = TRUE)
heat_limit <- max(0.20, ceiling(heat_max * 100) / 100)

# =========================================================
# 作图
# =========================================================
p_heat <- ggplot() +
  
  # 建立完整 7 × 7 绘图区，但不绘制任何背景方块
  geom_blank(
    data = full_grid,
    aes(x = col_site, y = row_site)
  ) +
  
  # 只绘制下三角的两两比较结果
  geom_tile(
    data = heat_df,
    aes(
      x = col_site,
      y = row_site,
      fill = R2
    ),
    color = "white",
    linewidth = 1.0,
    width = 0.96,
    height = 0.96
  ) +
  
  # 对角线虚线
  geom_segment(
    data = diag_df,
    aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend
    ),
    inherit.aes = FALSE,
    color = "grey55",
    linetype = "dashed",
    linewidth = 0.55,
    lineend = "butt"
  ) +
  
  # 下三角色块中的 R² 数值
  geom_text(
    data = heat_df,
    aes(
      x = col_site,
      y = row_site,
      label = label
    ),
    size = 3.5,
    fontface = "bold",
    color = "black"
  ) +
  
  scale_fill_gradient(
    low = "#FFF5F0",
    high = "#EF3B2C",
    limits = c(0, heat_limit),
    oob = scales::squish,
    name = expression(R^2)
  ) +
  
  # 横轴完整保留七个位点
  scale_x_discrete(
    limits = SITE_ORDER,
    drop = FALSE,
    expand = expansion(add = 0.5)
  ) +
  
  # 纵轴完整保留七个位点，并显示样本量
  scale_y_discrete(
    limits = rev(SITE_ORDER),
    labels = row_label_map,
    drop = FALSE,
    expand = expansion(add = 0.5)
  ) +
  
  coord_fixed(
    xlim = c(0.5, n_site + 0.5),
    ylim = c(0.5, n_site + 0.5),
    clip = "off",
    expand = FALSE
  ) +
  
  labs(
    x = NULL,
    y = NULL
  ) +
  
  theme_minimal(base_size = 10.5) +
  
  theme(
    panel.grid = element_blank(),
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      face = "bold",
      color = "black"
    ),
    
    axis.text.y = element_text(
      face = "bold",
      color = "black",
      lineheight = 0.95
    ),
    
    axis.ticks = element_blank(),
    
    legend.title = element_text(
      face = "bold",
      size = 11
    ),
    
    legend.text = element_text(
      color = "black"
    ),
    
    legend.key.height = unit(2.0, "cm"),
    
    plot.margin = margin(
      t = 8,
      r = 8,
      b = 8,
      l = 8
    )
  )

p_heat
# =========================================================
# 7. 导出结果
# =========================================================
ggsave(
  file.path(OUT_DIR, "Fig_oral_site_PCoA.pdf"),
  p_pcoa,
  width = 8,
  height = 5.1,
  units = "in"
)

ggsave(
  file.path(OUT_DIR, "Fig_oral_site_PCoA.png"),
  p_pcoa,
  width = 8,
  height = 5.1,
  units = "in",
  dpi = 600
)

ggsave(
  file.path(OUT_DIR, "Fig_oral_site_pairwise_R2_triangle.pdf"),
  p_heat,
  width = 6.8,
  height = 6.0,
  units = "in"
)

ggsave(
  file.path(OUT_DIR, "Fig_oral_site_pairwise_R2_triangle.png"),
  p_heat,
  width = 6.8,
  height = 6.0,
  units = "in",
  dpi = 600
)
ggsave(
  file.path(OUT_DIR, "PERMANOVA_pairwise_heatmap_full_triangle.png"),
  p_heat,
  width = 7.2,
  height = 6.5,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(OUT_DIR, "PERMANOVA_pairwise_heatmap_full_triangle.pdf"),
  p_heat,
  width = 7.2,
  height = 6.5,
  bg = "white"
)
write.csv(
  pcoa_df,
  file.path(OUT_DIR, "oral_site_PCoA_coordinates.csv"),
  row.names = FALSE
)

cat("\n==================== DONE ====================\n")
cat("总体 PERMANOVA: R² = ", sprintf("%.3f", overall_r2),
    ", P = ", format_p_plot(overall_p), "\n", sep = "")
cat("输出目录：", OUT_DIR, "\n", sep = "")

