# ------------------------------------------------------------------------------
# Count species in the Star Wars data set
# March 21, 2023
# TS O'Leary
# ------------------------------------------------------------------------------

# Description -----
# Count the number of observations of each species and save as a .rds file

# Load libraries
library(tidyverse)

# Load data
dat <- read_csv(here::here("data/raw/starwars.csv"))

# Analyze data
dat <- dat %>%
  group_by(species) %>%
  count()

# Save data
saveRDS(dat, here::here("data/processed/starwars_species_count.rds"))
