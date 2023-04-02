# ------------------------------------------------------------------------------
# Group iris data by species
# March 21, 2023
# TS O'Leary
# ------------------------------------------------------------------------------

# Description -----
# Load iris data and group by species and save as a rds file

# Load libraries
library(tidyverse)

# Load data
dat <- read_csv(here::here("data/raw/iris.csv"))

# Analyze data
dat <- dat %>%
  group_by(species)

# Save data
saveRDS(dat, here::here("data/processed/species_grouped.rds"))