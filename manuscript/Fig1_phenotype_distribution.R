#!/usr/bin/env Rscript
# =======================================================================
# Fig 1: Number of oral metagenomic samples and studies by phenotype
# This version computes all numbers from metadata (no hard-coded counts).
# Preserves the original square-point + lollipop + Studies column design.
# =======================================================================

# ---- 0. Packages ----
packages <- c("data.table", "dplyr", "ggplot2", "scales")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(scales)
})

# ---- 1. Paths and filter ----
INPUT_CSV <- file.path(Sys.getenv("OMWI_ROOT", getwd()), "data/metadata", "mix_sited_metadata.csv")
OUT_DIR <- file.path(Sys.getenv("OMWI_ROOT", getwd()), "fig1", "condition_dotplot")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# "Oral_proxy" = only target oral sites; "ALL" = all sites in CSV
SET_FILTER <- "Oral_proxy"
if (!file.exists(INPUT_CSV)) stop("输入CSV不存在，请检查 INPUT_CSV：", INPUT_CSV)

# ---- 2. Read and clean metadata ----
meta <- data.table::fread(
  INPUT_CSV, data.table = FALSE, check.names = FALSE,
  na.strings = c("", "NA", "N/A"), encoding = "unknown"
)
names(meta) <- sub("^\ufeff", "", names(meta))
required <- c("sample_id", "study_name", "Set", "system", "group")
if (!all(required %in% names(meta))) {
  stop("CSV缺少列：", paste(setdiff(required, names(meta)), collapse = ", "))
}

meta <- meta %>%
  mutate(across(all_of(required), ~ trimws(as.character(.x)))) %>%
  filter(!is.na(sample_id), sample_id != "",
         !is.na(study_name), study_name != "",
         !is.na(system), system != "",
         !is.na(group), group != "")

if (SET_FILTER != "ALL") meta <- meta %>% filter(Set == SET_FILTER)
if (nrow(meta) == 0) stop("筛选后没有样本：", SET_FILTER)

# Only standardize spelling / capitalization. Do not silently merge diagnoses.
meta <- meta %>%
  mutate(
    group = case_when(
      system == "Healthy" & tolower(group) %in% c("control", "healthy") ~ "Healthy",
      tolower(group) == "rheumatoid arthritis" ~ "Rheumatoid arthritis",
      tolower(group) == "adenoid cystic carcinoma" ~ "Adenoid cystic carcinoma",
      tolower(group) == "superficial gastritis" ~ "Superficial gastritis",
      tolower(group) == "atrophic gastritis" ~ "Atrophic gastritis",
      tolower(group) == "intestinal metaplasia" ~ "Intestinal metaplasia",
      toupper(group) == "COVID" ~ "COVID-19",
      TRUE ~ group
    )
  )

# Each row of the plot is a phenotype within a system.
# A study can count in multiple phenotype rows, so DO NOT sum Studies_Num.
plot_dat <- meta %>%
  group_by(system, group) %>%
  summarise(
    Samples_Num = n_distinct(sample_id),
    Studies_Num = n_distinct(study_name),
    .groups = "drop"
  ) %>%
  rename(System = system, Condition = group)

# Retain the original paper's categorical palette, including the new system.
system_color_map <- c(
  "Healthy" = "#2CA25F",
  "Oral Disease" = "#C7A6D8",
  "Respiratory Disease" = "#8FA8C9",
  "Digestive Disease" = "#D97C8A",
  "Immune Disease" = "#D8C97A",
  "Metabolic Disease" = "#4198AC",
  "BMI Group" = "#FFC100",
  "Neuropsychiatric Disease" = "#9F86C0"
)
system_order <- c(
  "Healthy", "Oral Disease", "Respiratory Disease",
  "Digestive Disease", "Immune Disease", "Metabolic Disease",
  "BMI Group", "Neuropsychiatric Disease"
)

unknown_systems <- setdiff(unique(plot_dat$System), names(system_color_map))
if (length(unknown_systems)) {
  stop("尚未给这些system配置颜色：", paste(unknown_systems, collapse = ", "))
}

plot_dat <- plot_dat %>%
  mutate(System = factor(System, levels = system_order)) %>%
  arrange(System, desc(Samples_Num), Condition) %>%
  mutate(
    # factor levels must be unique; use an internal unique key per row
    condition_key = paste(as.character(System), Condition, sep = " | "),
    condition_key = factor(condition_key, levels = rev(condition_key))
  )

# Export inspectable statistics.
write.csv(
  plot_dat %>% select(System, Condition, Studies_Num, Samples_Num),
  file.path(OUT_DIR, paste0("Fig1_condition_counts_", SET_FILTER, ".csv")),
  row.names = FALSE, fileEncoding = "UTF-8"
)
cat("\n------ 数据核查 ------\n")
cat("Set：", SET_FILTER, "\n")
cat("不同样本：", n_distinct(meta$sample_id), "\n")
cat("不同研究：", n_distinct(meta$study_name), "\n")
cat("疾病/条件：", nrow(plot_dat), "\n")
print(plot_dat %>% select(System, Condition, Studies_Num, Samples_Num), n = Inf)

# ---- 3. Design: square lollipop + Studies column ----
max_n <- max(plot_dat$Samples_Num)
divider_x <- max_n * 1.16
study_x <- max_n * 1.34
x_max <- max_n * 1.48

# Dynamic tick marks on the sqrt scale; keep the axis interpretable.
base_breaks <- c(0, 10, 25, 50, 100, 250, 500, 750, 1000, 1250, 1500, 2000, 3000)
x_breaks <- base_breaks[base_breaks <= max_n]
if (length(x_breaks) < 2) x_breaks <- c(0, max_n)

# Workaround for duplicated condition names across Systems:
# store unique keys internally and strip the System prefix for axis labels.

p <- ggplot(plot_dat, aes(x = Samples_Num, y = condition_key)) +
  geom_segment(
    aes(x = 0, xend = Samples_Num, yend = condition_key),
    color = "grey55", linewidth = 0.35
  ) +
  geom_point(
    aes(fill = System), shape = 22, size = 5,
    color = "white", stroke = 0.35
  ) +
  geom_vline(xintercept = divider_x, color = "black", linewidth = 0.3) +
  geom_text(
    aes(x = study_x, label = Studies_Num),
    size = 5, color = "grey20"
  ) +
  annotate(
    "text", x = study_x, y = nrow(plot_dat) + 0.9,
    label = "Studies", size = 5, fontface = "bold", color = "grey20"
  ) +
  scale_fill_manual(values = system_color_map, drop = TRUE) +
  scale_y_discrete(labels = function(x) sub("^.*[|] ", "", x),
                   expand = expansion(add = c(0.65, 1.6))) +
  scale_x_continuous(
    trans = "sqrt", breaks = x_breaks, labels = scales::comma,
    limits = c(0, x_max), expand = expansion(mult = c(0.005, 0.02))
  ) +
  labs(x = "Number of samples", y = NULL, fill = "System") +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 10) +
  theme(
    axis.text.y = element_text(size = 15, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 12, color = "black"),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.x = element_line(linewidth = 0.35, color = "black"),
    axis.ticks.x = element_line(linewidth = 0.3, color = "black"),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 15),
    legend.key.size = grid::unit(0.35, "cm"),
    legend.position = "bottom",
    plot.margin = margin(4, 12, 4, 4)
  )

print(p)

# ---- 4. Outputs ----
# More phenotypes than in old figure: make the panel taller automatically.
plot_h <- max(4.1, 1.5 + 0.245 * nrow(plot_dat))
base <- file.path(OUT_DIR, paste0("Fig1_condition_dotplot_", SET_FILTER))

ggsave(paste0(base, ".pdf"), p, width = 7, height = 6, bg = "white")
ggsave(paste0(base, ".png"), p,
       width = 7, height = 6, dpi = 600, bg = "white")
ggsave(paste0(base, ".tif"), p,
       width = 7, height = 6, dpi = 600,
       device = "tiff", compression = "lzw", bg = "white")
message("Saved: ", normalizePath(OUT_DIR, mustWork = FALSE))
