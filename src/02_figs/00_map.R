# ------------------------------------------------------------------------------
# Map
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)
library(maps)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(geodata)
source("src/_theme_colors.R")

# Load data
world <- ne_countries(scale = "medium", returnclass = "sf")
lat_long <- read_tsv("data/raw/survival/lat_long_data.txt") |> 
  st_as_sf(coords = c("long", "lat"), crs = 4326, remove = FALSE) |> 
  mutate(location = factor(location, levels = c("Vermont, USA",
                                                "Montpellier, France",
                                                "Shiojiri, Japan",
                                                "Chiapas, Mexico",
                                                "Accra, Ghana",
                                                "Mumbai, India")))

# Create sf objects for the equator and tropics (spanning longitudes -180 to 180)
equator <- st_linestring(matrix(c(-180, 0, 180, 0), ncol = 2, byrow = TRUE))
tropic_n <- st_linestring(matrix(c(-180, 23.43602, 180, 23.43602), ncol = 2, byrow = TRUE))
tropic_s <- st_linestring(matrix(c(-180, -23.43602, 180, -23.43602), ncol = 2, byrow = TRUE))

# Combine the lines into an sf object
lines_sf <- st_sfc(equator, tropic_n, tropic_s, crs = 4326)
line_types <- c("solid", "dashed", "dashed") 
lines_sf <- st_sf(geometry = lines_sf, line_type = line_types)

# Create plot with locations plotted here
map <- ggplot() +
  geom_sf(data = world, 
          fill = "lightgreen", 
          alpha = 0.11, color = "grey80", size = 0.2) +
  geom_sf(data = lines_sf, aes(linetype = line_types), color = "grey75") +
  geom_sf(data = lat_long, 
             aes(fill = location),
          shape = 21,
          color = "grey99",
          size = 4) +
  geom_sf_label(data = lat_long, 
          aes(label = location,
              fill = location),
          color = "grey99",
          size = 4,
          nudge_y = -1000000) +
  coord_sf(crs = "+proj=robin") +
  scale_linetype_identity() +
  scale_fill_manual(values = colors_genotypes) +
  theme_minimal() +
  theme(legend.position = "none")

ggsave(
  plot = map,
  here::here("output/figs/metabolomics/map.pdf")
)
ggsave(
  plot = map,
  here::here("output/figs/metabolomics/map.svg")
)

saveRDS(map, here::here("output/figs/metabolomics/map.rds"))