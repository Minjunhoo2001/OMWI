#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(ragg)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(sf)
})

options(stringsAsFactors = FALSE)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else "scripts/Figure1/Figure1D_map.R"
repo_root <- normalizePath(file.path(dirname(script_file), "..", ".."), mustWork = FALSE)

get_arg <- function(flag, default) {
  args <- commandArgs(trailingOnly = TRUE)
  hit <- match(flag, args)
  if (!is.na(hit) && hit < length(args)) args[[hit + 1]] else default
}

input_file <- get_arg(
  "--input",
  file.path(repo_root, "results", "Figure1", "Fig1D_map", "map_source_data.csv")
)
out_dir <- get_arg(
  "--out-dir",
  file.path(repo_root, "results", "Figure1", "Fig1D_map")
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

map_data <- read.csv(input_file, check.names = FALSE) %>%
  select(Region, StudyCount, SampleCount) %>%
  mutate(
    Region = recode(
      Region,
      "USA" = "United States of America",
      "United State" = "United States",
      "Japanese" = "Japan",
      "FJI" = "Fiji"
    ),
    StudyCount = as.numeric(StudyCount),
    SampleCount = as.numeric(SampleCount)
  )

stopifnot(all(c("Region", "StudyCount", "SampleCount") %in% names(map_data)))

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
coords <- suppressWarnings(sf::st_centroid(world)) %>%
  sf::st_coordinates() %>%
  as.data.frame() %>%
  cbind(Region = world$name_long)
colnames(coords)[1:2] <- c("lon", "lat")

plot_data <- left_join(map_data, coords, by = "Region")
if (anyNA(plot_data$lon) || anyNA(plot_data$lat)) {
  stop("Map coordinates were not resolved for: ", paste(plot_data$Region[is.na(plot_data$lon) | is.na(plot_data$lat)], collapse = ", "))
}

write.csv(
  plot_data,
  file.path(out_dir, "Fig1D_map_plot_data.csv"),
  row.names = FALSE
)

p_map <- ggplot(world) +
  geom_sf(fill = "grey78", colour = "white", linewidth = 0.18) +
  geom_point(
    data = plot_data,
    aes(x = lon, y = lat, size = SampleCount, fill = StudyCount),
    shape = 21,
    colour = "grey20",
    stroke = 0.35,
    alpha = 0.90
  ) +
  geom_text_repel(
    data = plot_data,
    aes(x = lon, y = lat, label = Region),
    size = 3.0,
    family = "Arial",
    min.segment.length = 0,
    max.overlaps = Inf,
    seed = 101,
    box.padding = 0.18,
    point.padding = 0.12
  ) +
  scale_size_continuous(
    name = "Sample count",
    limits = c(0, 600),
    breaks = c(100, 200, 400, 600),
    range = c(1.2, 7.2)
  ) +
  scale_fill_gradient(
    name = "Study count",
    low = "#FDE0DD",
    high = "#DE2D26",
    limits = c(1, max(plot_data$StudyCount, na.rm = TRUE))
  ) +
  coord_sf(xlim = c(-180, 195), ylim = c(-65, 85), expand = FALSE) +
  labs(x = NULL, y = NULL, tag = "d") +
  theme_minimal(base_size = 8, base_family = "Arial") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.45),
    axis.text = element_text(colour = "grey35", size = 7),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_text(size = 7.5),
    legend.text = element_text(size = 7),
    plot.tag = element_text(face = "bold", size = 9),
    plot.margin = margin(4, 5, 2, 4)
  ) +
  guides(
    fill = guide_colourbar(order = 1, title.position = "left", barwidth = grid::unit(23, "mm"), barheight = grid::unit(2.5, "mm")),
    size = guide_legend(order = 2, title.position = "left", nrow = 1)
  )

pdf_file <- file.path(out_dir, "Fig1D_map_bubble.pdf")
png_file <- file.path(out_dir, "Fig1D_map_bubble.png")
ggsave(pdf_file, p_map, width = 10, height = 4, device = grDevices::cairo_pdf, bg = "white")
ggsave(png_file, p_map, width = 10, height = 4, dpi = 600, device = ragg::agg_png, bg = "white")
capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))

cat("Saved Figure 1d map to:\n", pdf_file, "\n", png_file, "\n", sep = "")
