suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(vegan)
  library(patchwork)
})

# =========================================================
# 0. 路径
# =========================================================
TAX_FILE <- "/share/home/HeMinjun/metagenomic/Oral_proxy_v2/00.input/level/taxonomy_species.tsv"
META_FILE <- "/share/home/HeMinjun/metagenomic/Oral_proxy_v2/00.input/metadata6.0.csv"

OUT_DIR <- "/share/home/HeMinjun/metagenomic/GitHub/results/PCoA_body_site_species"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

set.seed(704)

# =========================================================
# 1. 工具函数
# =========================================================
pick_col <- function(df, candidates, label) {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) == 0) {
    stop("找不到列：", label, "\n候选列名：", paste(candidates, collapse = ", "))
  }
  hit[1]
}

format_p <- function(p) {
  ifelse(is.na(p), "NA", ifelse(p < 0.001, "< 0.001", sprintf("%.3f", p)))
}

sig_label <- function(p) {
  case_when(
    is.na(p) ~ "",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

standardize_site <- function(x) {
  x0 <- stringr::str_squish(as.character(x))
  case_when(
    str_detect(x0, regex("dental|plaque|supragingival|subgingival", ignore_case = TRUE)) ~ "Dental",
    str_detect(x0, regex("saliva", ignore_case = TRUE)) ~ "Saliva",
    str_detect(x0, regex("tongue", ignore_case = TRUE)) ~ "Tongue",
    str_detect(x0, regex("pharyngeal|throat|oropharyngeal|swab", ignore_case = TRUE)) ~ "Pharyngeal swab",
    str_detect(x0, regex("rinse|mouthwash|oral rinse", ignore_case = TRUE)) ~ "Oral rinse",
    TRUE ~ x0
  )
}

read_species_table <- function(tax_file, meta_sample_ids) {
  tax <- fread(tax_file, data.table = FALSE, check.names = FALSE)
  
  sample_cols <- intersect(colnames(tax), meta_sample_ids)
  
  if (length(sample_cols) == 0) {
    stop("taxonomy_species.tsv 的列名和 metadata 的 sample ID 没有交集。请检查样本名是否一致。")
  }
  
  feature_col <- colnames(tax)[1]
  feature_names <- make.unique(as.character(tax[[feature_col]]))
  
  x_df <- tax[, sample_cols, drop = FALSE]
  for (j in seq_len(ncol(x_df))) {
    x_df[[j]] <- suppressWarnings(as.numeric(x_df[[j]]))
  }
  
  x_mat <- t(as.matrix(x_df))
  colnames(x_mat) <- feature_names
  rownames(x_mat) <- sample_cols
  x_mat[is.na(x_mat)] <- 0
  
  x_mat
}

run_pcoa_permanova <- function(X0, meta0, site_order, tag) {
  
  message("\n================ ", tag, " ================")
  
  meta0 <- meta0 %>%
    filter(body_site %in% site_order) %>%
    mutate(body_site = factor(body_site, levels = site_order))
  
  X <- X0[meta0$SampleID, , drop = FALSE]
  
  # 去掉全 0 样本和全 0 物种
  keep_sample <- rowSums(X, na.rm = TRUE) > 0
  X <- X[keep_sample, , drop = FALSE]
  meta0 <- meta0[keep_sample, , drop = FALSE]
  
  keep_taxa <- colSums(X > 0, na.rm = TRUE) > 0
  X <- X[, keep_taxa, drop = FALSE]
  
  message("[样本数] ", nrow(X))
  message("[物种数] ", ncol(X))
  print(table(meta0$body_site))
  
  # Bray-Curtis 推荐先转相对丰度
  X_rel <- decostand(X, method = "total")
  
  # =======================================================
  # PCoA
  # =======================================================
  bray_dist <- vegdist(X_rel, method = "bray")
  
  pcoa <- cmdscale(bray_dist, eig = TRUE, k = 2)
  eig_pos <- pcoa$eig[pcoa$eig > 0]
  
  pcoa1_pct <- pcoa$eig[1] / sum(eig_pos) * 100
  pcoa2_pct <- pcoa$eig[2] / sum(eig_pos) * 100
  
  coord_df <- data.frame(
    SampleID = rownames(X),
    PCoA1 = pcoa$points[, 1],
    PCoA2 = pcoa$points[, 2],
    body_site = as.character(meta0$body_site),
    stringsAsFactors = FALSE
  )
  
  center_df <- coord_df %>%
    group_by(body_site) %>%
    summarise(
      PCoA1 = mean(PCoA1),
      PCoA2 = mean(PCoA2),
      n = n(),
      .groups = "drop"
    ) %>%
    mutate(body_site = factor(body_site, levels = site_order)) %>%
    arrange(body_site)
  
  fwrite(coord_df, file.path(OUT_DIR, paste0("PCoA_", tag, "_coordinates.csv")))
  fwrite(center_df, file.path(OUT_DIR, paste0("PCoA_", tag, "_centers.csv")))
  
  # =======================================================
  # global PERMANOVA
  # =======================================================
  adonis_meta <- data.frame(body_site = meta0$body_site)
  rownames(adonis_meta) <- rownames(X)
  
  global_adonis <- adonis2(
    bray_dist ~ body_site,
    data = adonis_meta,
    permutations = 999
  )
  
  global_R2 <- global_adonis$R2[1]
  global_P  <- global_adonis$`Pr(>F)`[1]
  
  message("[Global PERMANOVA] R2 = ", sprintf("%.3f", global_R2),
          ", P = ", format_p(global_P))
  
  # =======================================================
  # pairwise PERMANOVA
  # =======================================================
  dist_mat <- as.matrix(bray_dist)
  
  pair_list <- combn(site_order, 2, simplify = FALSE)
  
  pair_df <- map_dfr(pair_list, function(pair) {
    g1 <- pair[1]
    g2 <- pair[2]
    
    idx <- meta0$body_site %in% c(g1, g2)
    
    sub_meta <- data.frame(
      body_site = droplevels(meta0$body_site[idx])
    )
    rownames(sub_meta) <- rownames(X)[idx]
    
    sub_dist <- as.dist(dist_mat[idx, idx])
    
    fit <- adonis2(
      sub_dist ~ body_site,
      data = sub_meta,
      permutations = 999
    )
    
    data.frame(
      group1 = g1,
      group2 = g2,
      n1 = sum(sub_meta$body_site == g1),
      n2 = sum(sub_meta$body_site == g2),
      F_model = fit$F[1],
      R2 = fit$R2[1],
      P = fit$`Pr(>F)`[1],
      stringsAsFactors = FALSE
    )
  }) %>%
    mutate(
      P_adj_BH = p.adjust(P, method = "BH"),
      sig = sig_label(P_adj_BH)
    )
  
  fwrite(pair_df, file.path(OUT_DIR, paste0("PERMANOVA_pairwise_", tag, ".csv")))
  
# =======================================================
# 画图
# =======================================================
site_color <- c(
  "Dental"          = "#E86F61",
  "Saliva"          = "#48A6D1",
  "Pharyngeal swab" = "#2EA66F",
  "Tongue"          = "#B792D6",
  "Oral rinse"      = "#C2B84B"
)

coord_df <- coord_df %>%
  mutate(body_site = factor(body_site, levels = site_order))

center_df <- center_df %>%
  mutate(
    body_site = factor(body_site, levels = site_order),
    label = as.character(body_site)
  )

global_lab <- paste0(
  "R\u00B2 = ", sprintf("%.3f", global_R2),
  "\nP = ", format_p(global_P)
)

x_rng <- range(coord_df$PCoA1, na.rm = TRUE)
y_rng <- range(coord_df$PCoA2, na.rm = TRUE)

# -------------------------------------------------------
# PCoA plot
# -------------------------------------------------------
p_pcoa <- ggplot(coord_df, aes(PCoA1, PCoA2, color = body_site)) +
  geom_point(
    shape = 16,
    size = 2.0,
    alpha = 0.55
  ) +
  geom_label(
    data = center_df,
    aes(PCoA1, PCoA2, label = label),
    fill = scales::alpha("white", 0),
    color = "black",
    fontface = "bold",
    size = 4.0,
    label.size = 0.7,
    label.padding = unit(0.1, "lines"),
    label.r = unit(0, "lines"),
    show.legend = FALSE
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = global_lab,
    hjust = 1.1,
    vjust = 1.8,
    fontface = "bold",
    size = 5
  ) +
  scale_color_manual(values = site_color, drop = FALSE) +
  labs(
    x = paste0("PCoA1 (", sprintf("%.2f", pcoa1_pct), "%)"),
    y = paste0("PCoA2 (", sprintf("%.2f", pcoa2_pct), "%)"),
    color = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.border = element_rect(color = "black", linewidth = 1),
    axis.title = element_text(face = "bold", color = "black"),
    axis.text = element_text(color = "black"),
    legend.position = "right",
    legend.justification = c(0, 0),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.text = element_text(size = 10),
    plot.margin = margin(8, 8, 8, 8)
  )

# -------------------------------------------------------
# heatmap 数据整理
# -------------------------------------------------------
n_df <- center_df %>%
  mutate(
    body_site = as.character(body_site),
    body_site_label = case_when(
      body_site == "Pharyngeal swab" ~ paste0("Pharyngeal\nswab\n(n=", n, ")"),
      TRUE ~ paste0(body_site, "\n(n=", n, ")")
    )
  ) %>%
  select(body_site, body_site_label)

label_levels <- n_df$body_site_label[match(site_order, n_df$body_site)]

pair_plot_df <- pair_df %>%
  mutate(
    i = match(group1, site_order),
    j = match(group2, site_order),
    row_site = ifelse(i > j, group1, group2),
    col_site = ifelse(i > j, group2, group1)
  ) %>%
  left_join(n_df, by = c("row_site" = "body_site")) %>%
  mutate(
    row_site_label = factor(body_site_label, levels = label_levels),
    col_site = factor(col_site, levels = site_order),
    R2_label = sprintf("%.3f", R2)
  )

# -------------------------------------------------------
# heatmap plot
# -------------------------------------------------------
p_heat <- ggplot(pair_plot_df, aes(x = col_site, y = row_site_label, fill = R2)) +
  geom_tile(color = "white", linewidth = 1.2) +
  geom_text(aes(label = R2_label), size = 3.6, fontface = "bold") +
  scale_fill_gradient(
    low = "#F7E6E2",
    high = "#F11E14",
    limits = c(0, max(0.20, max(pair_plot_df$R2, na.rm = TRUE))),
    breaks = c(0.05, 0.10, 0.15),
    name = expression(R^2)
  ) +
  scale_x_discrete(limits = site_order, drop = FALSE) +
  scale_y_discrete(limits = rev(label_levels), drop = FALSE) +
  coord_fixed() +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      color = "black",
      face = "bold"
    ),
    axis.text.y = element_text(
      color = "black",
      face = "bold",
      hjust = 1
    ),
    legend.position = "right",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    plot.margin = margin(8, 8, 8, 8)
  )

# =======================================================
# 分开导出：PCoA 单独，heatmap 单独
# =======================================================

# ---- PCoA ----
ggsave(
  file.path(OUT_DIR, paste0("PCoA_", tag, ".png")),
  p_pcoa,
  width = 7,
  height = 4.8,
  dpi = 600
)

ggsave(
  file.path(OUT_DIR, paste0("PCoA_", tag, ".pdf")),
  p_pcoa,
  width = 7,
  height = 4.8
)

# ---- heatmap ----
ggsave(
  file.path(OUT_DIR, paste0("PERMANOVA_heatmap_", tag, ".png")),
  p_heat,
  width = 4.8,
  height = 4.6,
  dpi = 600
)

ggsave(
  file.path(OUT_DIR, paste0("PERMANOVA_heatmap_", tag, ".pdf")),
  p_heat,
  width = 4.8,
  height = 4.6
)

return(list(
  p_pcoa = p_pcoa,
  p_heat = p_heat,
  coord_df = coord_df,
  center_df = center_df,
  pair_df = pair_df,
  global_adonis = global_adonis
))
}
# =========================================================
# 2. 读取 metadata
# =========================================================
meta <- fread(META_FILE, data.table = FALSE, check.names = FALSE)

sample_col <- pick_col(
  meta,
  c("SampleID", "sample_id", "sample", "Sample", "sample_name", "Sample_Name", "run", "Run"),
  "sample ID"
)

site_col <- pick_col(
  meta,
  c("body_site", "Body_site", "BodySite", "body_site_use", "oral_site", "niche", "sample_type", "SampleType"),
  "body site"
)

meta <- meta %>%
  mutate(
    SampleID = as.character(.data[[sample_col]]),
    body_site_raw = as.character(.data[[site_col]]),
    body_site = standardize_site(body_site_raw)
  )

# =========================================================
# 3. 读取 species abundance table
# =========================================================
X <- read_species_table(TAX_FILE, meta$SampleID)

common_samples <- intersect(rownames(X), meta$SampleID)

if (length(common_samples) == 0) {
  stop("taxonomy_species.tsv 和 metadata6.0.csv 没有共同样本。")
}

X <- X[common_samples, , drop = FALSE]
meta <- meta[match(common_samples, meta$SampleID), , drop = FALSE]

stopifnot(identical(rownames(X), meta$SampleID))

# =========================================================
# 4. 跑 all body sites
# =========================================================
site_order_all <- c(
  "Dental",
  "Saliva",
  "Tongue",
  "Pharyngeal swab",
  "Oral rinse"
)

res_all <- run_pcoa_permanova(
  X0 = X,
  meta0 = meta,
  site_order = site_order_all,
  tag = "body_site_species"
)

res_all$p_pcoa
res_all$p_heat


