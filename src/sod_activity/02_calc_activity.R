# ------------------------------------------------------------------------------
# Calculated SOD Activity
# April 02, 2023
# TS O'Leary
# ------------------------------------------------------------------------------

# Description -----
# Count the number of observations of each species.

# Load libraries
library(tidyverse)
library(kableExtra)

# Load data
dat <- readRDS(here::here("data/processed/sod_activity.rds"))

# Standard Curve
df_std <- dat %>%
  filter(Sample == "Standard" & Date == "2022-02-15" & Activity < 0.1) %>%
  group_by(Sample, Activity) %>%
  summarise(Abs = mean(Abs)) %>%
  mutate(Ratio = Abs[Activity == 0] / Abs)

# Run Linear Regression and save information
fit <- lm(df_std, formula = Ratio ~ Activity)

y_int <- fit$coefficients[1]
slope <- fit$coefficients[2]

# Calculate Activity
df_calc <- dat %>%
  filter(!is.na(Bio_rep) | Sample == "Standard" & Activity == 0) %>%
  mutate(Ratio = Abs[Sample == "Standard"] / Abs) %>%
  group_by(Sample, Temp, Bio_rep, Date) %>%
  summarise(Ratio = mean(Ratio)) %>%
  mutate(sod_activity_U_mL = ((Ratio - y_int) / slope) * (0.230 / 0.01))

# Print out table
df_calc %>%
  filter(Sample != "Standard") %>%
  kable(digits = 4) %>%
  kable_styling("striped") %>%
  scroll_box(height = "250px")

# Save activity
saveRDS(df_calc, here::here("data/processed/sod_activity_02.rds"))
