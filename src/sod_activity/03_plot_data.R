# ------------------------------------------------------------------------------
# Plot SOD Activity
# April 02, 2023
# TS O'Leary
# ------------------------------------------------------------------------------

# Description -----
# Plot SOD Activity data

# Load libraries
library(tidyverse)

# Load data
dat <- readRDS(here::here("data/processed/sod_activity_02.rds"))

# Plot
dat %>%
  mutate(region = ifelse(Sample %in% c("VT10", "VT12"), "VT", "Trop")) %>%
  filter(Sample != "Standard") %>%
  group_by(Temp, region) %>%
  ggplot(aes(
    x = as.factor(Temp),
    y = sod_activity_U_mL,
    fill = region
  )) +
  geom_boxplot(
    width = 0.5,
    alpha = 0.5,
    color = "grey50",
    outlier.shape = NA
  ) +
  ggbeeswarm::geom_beeswarm(
    aes(color = region),
    alpha = 0.6,              
    size = 3,
    dodge.width = 0.5,
    cex = 2,
    color = "grey50",
    shape = 21,
    alpha = 0.9
  ) +
  scale_y_continuous(
    limits = c(0,  0.8),
    name = "Total SOD Activity (U/mL)",
    expand = expansion(mult = c(0, 0.05)),
  ) +
  scale_x_discrete(
    name = "Acute heat shock",
    breaks = c(25, 36),
    labels = c("25°C", "36°C")
  ) +
  cowplot::theme_minimal_hgrid(font_family = "Myriad Pro") +
  scale_fill_manual(
    name = "Region",
    values = c("coral", "deepskyblue"),
    labels = c("Tropical", "Vermont")
  )

# Save the plot
ggsave("output/figs/sod_activity_2022.tiff",
       width = 6,
       height = 4
)