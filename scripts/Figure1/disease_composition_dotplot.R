suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
})

system_color_map <- c(
  "Healthy" = "#2CA25F",
  "Oral Disease" = "#C7A6D8",
  "Digestive Disease" = "#D97C8A",
  "Respiratory Disease" = "#8FA8C9",
  "Immune Disease" = "#D8C97A",
  "Metabolic Disease" = "#4198AC",
  "BMI Group" = "#FFC100"
)

plot_dat <- tribble(
  ~System, ~Condition, ~Studies_Num, ~Samples_Num,
  "Oral Disease",        "Periodontitis",              7, 130,
  "Oral Disease",        "Adenoid cystic carcinoma",   1, 21,
  "Respiratory Disease", "COVID",                      1, 224,
  "Metabolic Disease",   "Diabetes",                   4, 42,
  "Metabolic Disease",   "Hypertension",               1, 10,
  "BMI Group",           "Overweight",                 6, 50,
  "BMI Group",           "Obesity",                    6, 26,
  "BMI Group",           "BMI low",                    2, 7,
  "Digestive Disease",   "Pancreatic cancer",          2, 73,
  "Digestive Disease",   "Gastritis",                  1, 55,
  "Digestive Disease",   "Colorectal cancer",          2, 52,
  "Digestive Disease",   "Colorectal polyps",          1, 30,
  "Digestive Disease",   "Pancreatitis",               2, 24,
  "Digestive Disease",   "Intestinal metaplasia",      1, 23,
  "Immune Disease",      "Rheumatoid arthritis",       1, 75,
  "Healthy",             "Healthy",                   23, 1024
)

system_order <- c(
  "Healthy",
  "Oral Disease",
  "Respiratory Disease",
  "Digestive Disease",
  "Immune Disease",
  "Metabolic Disease",
  "BMI Group"
)

plot_dat <- plot_dat %>%
  mutate(System = factor(System, levels = system_order)) %>%
  arrange(System, desc(Samples_Num)) %>%
  mutate(
    Condition = factor(Condition, levels = rev(Condition))
  )

study_x <- 1350
divider_x <- 1150

p <- ggplot(plot_dat, aes(x = Samples_Num, y = Condition)) +
  geom_segment(
    aes(x = 0, xend = Samples_Num, yend = Condition),
    color = "grey50",
    linewidth = 0.35
  ) +
  geom_point(
    aes(fill = System),
    shape = 22,
    size = 3.0,
    color = "white",
    stroke = 0.35
  ) +
  geom_vline(
    xintercept = divider_x,
    color = "black",
    linewidth = 0.3
  ) +
  geom_text(
    aes(x = study_x, label = Studies_Num),
    size = 3.5,
    color = "grey20"
  ) +
  annotate(
    "text",
    x = study_x,
    y = length(unique(plot_dat$Condition)) + 0.8,
    label = "Studies",
    size = 3.5,
    fontface = "bold",
    color = "grey20"
  ) +
  scale_fill_manual(values = system_color_map) +
  scale_x_continuous(
    trans = "sqrt",
    breaks = c(0, 25, 50, 100, 250, 500, 1000),
    labels = comma,
    limits = c(0, 1450),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  labs(
    x = "Num of samples",
    y = NULL,
    fill = "System"
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 10) +
  theme(
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 10, color = "black"),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.x = element_line(linewidth = 0.35, color = "black"),
    axis.ticks.x = element_line(linewidth = 0.3, color = "black"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.35, "cm"),
    legend.position = "right",
    plot.margin = margin(4, 12, 4, 4)
  )

p

ggsave(
  "/share/home/HeMinjun/metagenomic/GitHub/results/Figure1/Disease_composition_dot/panel_condition_dotplot_study_column.pdf",
  p,
  width = 7,
  height = 4,
  device = cairo_pdf
)

ggsave(
  "/share/home/HeMinjun/metagenomic/GitHub/results/Figure1/Disease_composition_dot/panel_condition_dotplot_study_column.png",
  p,
  width = 7,
  height = 4,
  dpi = 600
)
