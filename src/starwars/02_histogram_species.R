# ------------------------------------------------------------------------------
# Histogram of species counts
# March 21, 2023
# TS O'Leary
# ------------------------------------------------------------------------------

# Description ----
# Create a histogram of the counted species.

# Load libraries
library(tidyverse)

# Load data
dat <- readRDS(here::here("data/processed/starwars_species_count.rds"))

# Histogram plot with cowplot theme
dat %>%
  ggplot() +
  geom_histogram(aes(x = n),
    color = "grey40",
    fill = "grey80"
  ) +
  labs(x = "Number of observations", y = "Count") +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05))
  ) +
  cowplot::theme_minimal_hgrid(font_family = "Myriad Pro")

# Save plot
cowplot::ggsave2(here::here("output/figs/starwars_species_hist.tiff"),
  height = 5,
  width = 7
)
