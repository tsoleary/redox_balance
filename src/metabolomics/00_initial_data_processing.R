# ------------------------------------------------------------------------------
# Initial data processing for downstream analysis
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)

# Load data
dat_raw <- read_csv(
  here::here("data/raw/metabolomics/redox_metabolomics_data_norm.csv"))  |> 
  separate(Sample, into = c("Genotype", "Temperature", "Replicate")) |> 
  mutate(Replicate = ifelse(is.na(Replicate), Temperature, Replicate)) |> 
  mutate(Temperature = ifelse(Temperature %in% c("25", "32"), 
                              paste0(Temperature, "°C"), "QC")) |> 
  mutate(Replicate = paste0("Rep_", Replicate)) |> 
  mutate(Region = case_when(Genotype %in% c("BO", "CH", "GA") ~ "Tropical",
                            Genotype %in% c("VT", "FR", "JP") ~ "Temperate",
                            Genotype == "QC" ~ "QC"))

# Tidy & nest data for analysis and plotting ------------
dat <- dat_raw |> 
  mutate(Locale = case_when(Genotype == "VT" ~ "Vermont, USA",
                            Genotype == "FR" ~ "Montpellier, France",
                            Genotype == "JP" ~ "Shiojiri, Japan",
                            Genotype == "CH" ~ "Chiapas, Mexico",
                            Genotype == "GA" ~ "Accra, Ghana",
                            Genotype == "BO" ~ "Mumbai, India")) |> 
  mutate(Locale = factor(Locale, levels = c("Vermont, USA", 
                                            "Montpellier, France", 
                                            "Shiojiri, Japan", 
                                            "Chiapas, Mexico", 
                                            "Accra, Ghana", 
                                            "Mumbai, India"))) |>
  select(Region, Locale, everything()) |> 
  pivot_longer(cols = Alanine:`GSSG/GSH-NEM Ratio`,
               names_to = "Metabolite",
               values_to = "NormIntensity") |> 
  group_by(Metabolite) |> 
  nest()

# Calculate log-transformed and pareto scaled data for pca ---------
dat_norm <- dat |> 
  unnest(cols = data) |> 
  group_by(Metabolite) |> 
  mutate(log_transformed = log10(NormIntensity),
         pareto_scaled = 
           (log_transformed - mean(log_transformed)) / 
           sqrt(sd(log_transformed))) 

# Save data
saveRDS(dat_raw, here::here("data/processed/metabolomics/dat_raw.rds"))
saveRDS(dat, here::here("data/processed/metabolomics/dat.rds"))
saveRDS(dat_norm, here::here("data/processed/metabolomics/dat_norm.rds"))

# Create a boxplot and density plot of NormIntensity and pareto scaled data --

# Box plot of NormIntensity data
p_box <- dat_norm |> 
  ggplot() +
  geom_boxplot(aes(y = Metabolite,
                   x = NormIntensity),
               fill = "grey90") +
  theme_minimal() +
  theme(axis.title.y = element_blank())

# Denisty plot of NormIntensity data
p_density <- dat_norm |> 
  ggplot() +
  geom_density(aes(x = NormIntensity),
               fill = "grey90") +
  theme_void()

# NormIntensity data
p_norm <- cowplot::plot_grid(
  p_density,
  p_box,
  align = "v",
  ncol = 1,
  rel_heights = c(1, 6)
)


# Box plot of pareto scaled data
p_box <- dat_norm |> 
  ggplot() +
  geom_boxplot(aes(y = Metabolite,
                   x = pareto_scaled),
               fill = "grey90") +
  theme_minimal() +
  theme(axis.title.y = element_blank())

# Density plot of Pareto scaled data
p_density <- dat_norm |> 
  ggplot() +
  geom_density(aes(x = pareto_scaled),
               fill = "grey90") +
  theme_void()

# Pareto scaled data
p_pareto <- cowplot::plot_grid(
  p_density,
  p_box,
  align = "v",
  ncol = 1,
  rel_heights = c(1, 6)
)

# Put both plots together
p_norm_all <- cowplot::plot_grid(
  p_norm,
  p_pareto,
  ncol = 2
)

# Save normalization plot
ggsave(plot = p_norm_all,
       here::here("output/figs/metabolomics/normalization_box_density.pdf"),
       height = 8,
       width = 10,
       units = "in")

