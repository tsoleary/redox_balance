# ------------------------------------------------------------------------------
# Differential abundant metabolites
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)
library(purrr)
library(lmerTest)
library(emmeans)
library(broom)
library(glue)

# Load data
dat <- readRDS(here::here("data/processed/metabolomics/dat.rds")) |>
  unnest(data) |>
  filter(Region != "QC") |>
  nest()

# Model formula with genotype as random effect ---------------------------------
model_formula <- NormIntensity ~ Region * Temperature + (1 | Genotype)

# Adjust p-values and find metabolites with significant terms ------------------
results <- dat |>
  mutate(model = map(data, ~lmer(model_formula, data = .))) |>  
  dplyr::select(Metabolite, model) |> 
  mutate(model = map(model, ~broom.mixed::tidy(.))) |> 
  unnest(model) |> 
  group_by(term) |> 
  mutate(p.value.adj = p.adjust(p.value, method = "fdr"))

# Differentially abundant metabolites ------------------------------------------
DAMs <- results |> 
  filter(term != "(Intercept)") |> 
  filter(p.value.adj < 0.05)


DAMs |> 
  group_by(term, statistic > 0) |> 
  tally()

DAMs |> 
  filter(term == "Temperature32°C",
         statistic > 0) |> 
  pull(Metabolite)

DAMs |> 
  filter(term == "RegionTropical:Temperature32°C") |> 
  pull(Metabolite)



# Save -------------------------------------------------------------------------
saveRDS(results, here::here("output/metabolomics/results.rds"))
saveRDS(DAMs, here::here("output/metabolomics/DAMs.rds"))

# Save to CSV for MetaboAnalyst ---
DAMs |> 
  ungroup() |> 
  filter(term == "Temperature32°C",
         !str_detect(Metabolite, "Ratio")) |> 
  select(Metabolite) |>
  distinct() |> 
  write_csv("output/metabolomics/DAMs_temp.csv")

DAMs |> 
  ungroup() |> 
  filter(term == "RegionTropical:Temperature32°C",
         !str_detect(Metabolite, "Ratio")) |> 
  select(Metabolite) |> 
  distinct() |> 
  write_csv("output/metabolomics/DAMs_regionXtemp.csv")


dat |> 
  filter(!str_detect(Metabolite, "Ratio")) |> 
  select(Metabolite) |> 
  write_csv("output/metabolomics/all_metabolites.csv")


# Post-hoc of Heat Shock Change within a region --------------------------------

DAMs_int <- DAMs |>
  filter(term == "RegionTropical:Temperature32°C") |>
  pull(Metabolite)

hs_effect_fct <- function(dat, Metabolite) {
  # Set levels for contrast of 32°C vs. 25°C
  dat <- dat %>%
    mutate(Temperature = factor(Temperature, levels = c("32°C", "25°C")))
  
  # Run model
  mod <- lmer(model_formula, data = dat)
  
  # Get emeans pairs
  hs_dat <- mod |>
    emmeans(~ Temperature | Region) |>
    pairs(infer = c(TRUE, TRUE)) |> 
    summary() |>
    as_tibble() |> 
    mutate(Metabolite = Metabolite)
  
  # Get effect sizes
  means <- mod |> 
    emmeans(~ Region*Temperature) |> 
    as_tibble() |> 
    select(Region, Temperature, emmean) |> 
    pivot_wider(names_from = Temperature, values_from = emmean) %>%
    rename(mean_25 = `25°C`, mean_32 = `32°C`) %>%
    mutate(log2_fold_change = log2(mean_32 / mean_25))
    
  # Join data
  hs_dat |> 
    left_join(means, by = "Region") |> 
    mutate(Metabolite = Metabolite)
  
}

# Get the 
dat |> 
  filter(Metabolite %in% DAMs_int) |> 
  mutate(out = map2(data, Metabolite, hs_effect_fct)) |> 
  pull(out) |> 
  list_rbind() |> 
  mutate(padj = p.adjust(p.value, method = "BH")) |> 
  filter(padj < 0.05) |> 
  select(Metabolite, everything())

