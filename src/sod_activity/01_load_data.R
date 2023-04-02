# ------------------------------------------------------------------------------
# SOD Activity Load data
# April 02, 2023
# TS O'Leary
# ------------------------------------------------------------------------------

# Description -----
# Load the SOD activity data

# Load libraries
library(tidyverse)

# Load data
dat <- bind_rows(
  "2022-02-15" = read_csv("data/raw/sod_activity_2022-02-15.csv"),
  "2022-03-09" = read_csv("data/raw/sod_activity_2022-03-09.csv"),
  .id = "Date"
)

# Save data
saveRDS(dat, here::here("data/processed/sod_activity_01.rds"))
