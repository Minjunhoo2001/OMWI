#!/usr/bin/env Rscript
# ================================================================
# Global distribution of oral metagenomic samples
# Bubble size: number of unique samples per country
# Bubble fill: number of distinct studies represented per country
# Source: all_oral_metadata.csv
# ================================================================

# ---- 0. Packages ----
packages <- c("data.table", "dplyr", "ggplot2", "sf", "rnaturalearth", "rnaturalearthdata", "scales")
missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(sf)
  library(rnaturalearth)
  library(scales)
})

# ---- 1. User settings ----
# If the specified path does not exist, a file-selection dialog will open.
INPUT_CSV <- file.path(Sys.getenv("OMWI_ROOT", getwd()), "data/metadata", "mix_sited_metadata.csv")
if (!file.exists(INPUT_CSV)) stop("Missing input metadata: ", INPUT_CSV)

# "Oral_proxy": use the oral proxy samples only (exclude tooth biofilm etc.)
# "ALL"       : use every sample in the supplied CSV (all oral sites)
SET_FILTER <- "Oral_proxy"

OUT_DIR <- file.path(file.path(Sys.getenv("OMWI_ROOT", getwd()), "fig1", "map_figure"))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- 2. Read metadata and derive counts ----
meta <- data.table::fread(INPUT_CSV, data.table = FALSE, check.names = FALSE,
                          encoding = "unknown", na.strings = c("", "NA", "N/A"))
names(meta) <- sub("^\ufeff", "", names(meta))
needed <- c("country", "study_name", "sample_id", "Set")
if (!all(needed %in% names(meta))) {
  stop("Missing required metadata columns: ",
       paste(setdiff(needed, names(meta)), collapse = ", "))
}

meta <- meta %>%
  mutate(across(all_of(needed), ~trimws(as.character(.x)))) %>%
  filter(!is.na(country), country != "",
         !is.na(study_name), study_name != "",
         !is.na(sample_id), sample_id != "")

if (SET_FILTER != "ALL") {
  meta <- meta %>% filter(Set == SET_FILTER)
}
if (nrow(meta) == 0) stop("No samples remained after filtering Set = ", SET_FILTER)

# Normalize country spellings that may occur in other versions of the file.
meta <- meta %>% mutate(country = dplyr::recode(
  country,
  "USA" = "United States",
  "United States of America" = "United States",
  "UK" = "United Kingdom",
  "FJI" = "Fiji",
  "Japanese" = "Japan"
))

country_counts <- meta %>%
  group_by(country) %>%
  summarise(
    SampleCount = n_distinct(sample_id),
    StudyCount  = n_distinct(study_name),
    .groups = "drop"
  ) %>%
  arrange(desc(SampleCount))

message("Selected dataset: ", SET_FILTER)
message("Unique samples: ", n_distinct(meta$sample_id),
        "; unique studies: ", n_distinct(meta$study_name),
        "; countries: ", n_distinct(meta$country))
print(country_counts, n = Inf)

write.csv(country_counts,
          file.path(OUT_DIR, paste0("country_sample_study_counts_", SET_FILTER, ".csv")),
          row.names = FALSE, fileEncoding = "UTF-8")

# ---- 3. Country marker and label coordinates ----
# Representative point locations and separately arranged labels:
# using fixed label positions prevents crowding in northwestern Europe.
# Coordinates are only for visualization, NOT sample-level locations.
positions <- tibble::tribble(
  ~country,         ~lon,    ~lat,  ~label_lon, ~label_lat,
  "United States",  -98.0,   39.0,    -165,       39,
  "China",          104.0,   35.0,      81,       37,
  "Japan",          138.0,   36.0,     151,       39,
  "Germany",         10.5,   51.0,      27,       55,
  "Denmark",          9.9,   56.0,      26,       70,
  "Italy",           12.5,   42.5,      29,       42,
  "Fiji",           178.0,  -17.7,     158,      -28,
  "Spain",           -3.7,   40.4,     -24,       36,
  "Singapore",      103.8,    1.3,      74,       -8,
  "France",           2.2,   46.5,     -30,       49,
  "Philippines",    122.5,   12.8,     142,       12,
  "Luxembourg",       6.1,   49.7,     -50,       58,
  "Netherlands",      5.3,   52.1,     -48,       68
)

plot_df <- country_counts %>% left_join(positions, by = "country")
if (anyNA(plot_df$lon) || anyNA(plot_df$lat)) {
  stop("Unknown country coordinates. Please add these countries to 'positions': ",
       paste(plot_df$country[is.na(plot_df$lon) | is.na(plot_df$lat)], collapse = ", "))
}

# ---- 4. World geometry ----
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
# This plot uses the same plate-carree-like geographic orientation as the original.

# Breakpoints include 600 for Oral_proxy; expand to 800 if ALL has >600.
size_top <- max(600, ceiling(max(plot_df$SampleCount) / 200) * 200)
size_ticks <- c(50, 100, 200, 400, 600, 800, 1000, 1200)
size_ticks <- size_ticks[size_ticks <= size_top]

# Use fine connectors mainly in Europe (where markers are clustered).
connector_countries <- c("Netherlands", "Luxembourg", "Denmark", "Germany",
                         "France", "Italy", "Spain")

# ---- 5. Make publication-ready bubble map ----
p <- ggplot() +
  geom_sf(data = world, fill = "#ECECEC", color = "white", linewidth = 0.20) +
  geom_segment(
    data = filter(plot_df, country %in% connector_countries),
    aes(x = lon, y = lat, xend = label_lon, yend = label_lat),
    linewidth = 0.28, color = "#888888"
  ) +
  geom_point(
    data = plot_df,
    aes(x = lon, y = lat, size = SampleCount, fill = StudyCount),
    shape = 21, color = "#333333", stroke = 0.38, alpha = 0.94
  ) +
  geom_text(
    data = plot_df,
    aes(x = label_lon, y = label_lat, label = country),
    family = "sans", size = 3.4, color = "#222222", hjust = 0
  ) +
  scale_size_area(
    name = "Sample Count", max_size = 9,
    limits = c(0, size_top), breaks = size_ticks
  ) +
  scale_fill_gradient(
    name = "Study Count", low = "#FEE5D9", high = "#CB181D",
    breaks = scales::breaks_pretty(n = 4)
  ) +
  scale_x_continuous(
    breaks = c(-120, -60, 0, 60, 120),
    labels = c("120°W", "60°W", "0°", "60°E", "120°E")
  ) +
  scale_y_continuous(
    breaks = c(-50, 0, 50),
    labels = c("50°S", "0°", "50°N")
  ) +
  coord_sf(xlim = c(-174, 190), ylim = c(-60, 80),
           default_crs = sf::st_crs(4326), expand = FALSE) +
  guides(
    fill = guide_colorbar(order = 1, barheight = grid::unit(3.0, "cm"),
                          barwidth = grid::unit(0.43, "cm")),
    size = guide_legend(order = 2,
                        override.aes = list(fill = "white", color = "#333333"))
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 10.5, base_family = "sans") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "#303030", fill = NA, linewidth = 0.6),
    axis.text = element_text(color = "#444444", size = 9.5),
    axis.ticks = element_blank(),
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(color = "black", size = 10.5, face = "bold"),
    legend.text = element_text(color = "#333333", size = 9.5),
    legend.spacing.y = grid::unit(0.55, "cm"),
    plot.margin = margin(8, 9, 8, 9)
  )

print(p)

# ---- 6. Export both raster and vector formats ----
output_base <- file.path(OUT_DIR, paste0("global_oral_map_", SET_FILTER))
ggsave(paste0(output_base, ".png"), p,
       width = 11.8, height = 5.6, units = "in", dpi = 600, bg = "white")
ggsave(paste0(output_base, ".tif"), p,
       width = 11.8, height = 5.6, units = "in", dpi = 600,
       device = "tiff", compression = "lzw", bg = "white")
ggsave(paste0(output_base, ".pdf"), p,
       width = 7, height = 7, units = "in", bg = "white")
message("Files saved under: ", normalizePath(OUT_DIR, mustWork = FALSE))
