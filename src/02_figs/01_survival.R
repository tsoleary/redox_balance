# ------------------------------------------------------------------------------
# Plots of survival of the metabolomics genotypes
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)
library(drc)
source("src/_theme_colors.R")

# Set directory for output of figures
fig_dir <- "output/figs/metabolomics"

# Load in survival data - and clean it up
dat_surv_2018 <- read_tsv("data/raw/survival/lockwood_et_al_2018_survival.txt") |> 
  rename("n_eggs" = "eggs",
         "n_hatched" = "hatched") |> 
  filter(genotype %in% c("VTECK_10", "CH", "MU", "GH"))
dat_surv_fr_jp <- read_csv("data/raw/survival/fr_jp_hatch_2023-04-27.csv") |> 
  janitor::clean_names() |> 
  mutate(survival = n_hatched/n_eggs) |> 
  dplyr::select(genotype, temperature, n_eggs, n_hatched, survival) |> 
  mutate(region = "temperate") |> 
  filter(genotype %in% c("FRMO-01", "JPSC"))

# Combine data
dat_surv <- bind_rows(dat_surv_2018, dat_surv_fr_jp) |> 
  mutate(genotype = case_when(genotype == "VTECK_10" ~ "Vermont, USA",
                              genotype == "FRMO-01" ~ "Montpellier, France",
                              genotype == "JPSC" ~ "Shiojiri, Japan",
                              genotype == "CH" ~ "Chiapas, Mexico",
                              genotype == "GH" ~ "Accra, Ghana",
                              genotype == "MU" ~ "Mumbai, India")) |> 
  mutate(genotype = factor(genotype, levels = c("Vermont, USA",
                                                "Montpellier, France",
                                                "Shiojiri, Japan",
                                                "Chiapas, Mexico",
                                                "Accra, Ghana",
                                                "Mumbai, India")))

# Survival curves
dat_lt50 <- dat_surv |> 
  group_by(genotype) |> 
  nest() |> 
  mutate(fit = map(data, ~ drm(n_hatched/n_eggs ~ temperature, 
                               data = .x, 
                               weights = n_eggs,
                               fct = LL.3(names = c("slope", "upper limit", "LT50")),
                               type = "binomial"))) |> 
  mutate(tidy_fit = map(fit, broom::tidy)) |> 
  unnest(tidy_fit) |> 
  filter(term == "LT50")

saveRDS(
  dat_lt50,
  here::here("data/processed/survival/lt50.rds")
)

# Print out LT50s
lt_se <- dat_lt50 |> 
  dplyr::select(genotype, term, estimate, std.error) |> 
  arrange(estimate) |> 
  dplyr::select(genotype, std.error)

# Create plot
p_surv <- dat_surv |> 
  ggplot(aes(x = temperature,
             y = survival,
             color = genotype,
             linetype = region)) +
  geom_smooth(method = drm,
              method.args = list(fct = LL.3()),
              se = FALSE,
              fullrange = TRUE) +
  scale_color_manual(values = colors_genotypes) +
  scale_x_continuous(limits = c(25, 40),
                     labels = function(x) paste0(x, "°C"),
                     name = "Acute heat shock") +
  scale_y_continuous(limits = c(0, 1),
                     labels = scales::percent,
                     name = "Hatching success") +
  scale_linetype_manual(values = c(1, 2)) +
  guides(linetype = "none") +
  theme_minimal(base_size = 16) +
  theme(legend.title = element_blank())

# Save plot
ggsave(
  plot = p_surv,
  here::here(fig_dir, "survival.pdf"),
  width = 18,
  height = 9, 
  units = "cm"
)


# BioClim variables ----

# Load data
bio <- geodata::worldclim_global(var = "bio", res = 10,
                        path = here::here("data/processed/worldclim"))
lat_long <- read_tsv("data/raw/survival/lat_long_data.txt") |> 
  sf::st_as_sf(coords = c("long", "lat"), crs = 4326, remove = FALSE)
bio_values <- terra::extract(bio, terra::vect(lat_long))
lat_long <- bind_cols(lat_long, bio_values) |> 
  left_join(lt_se, by = c("location" = "genotype"))

# lat_long |> 
#   mutate(genotype = factor(genotype, levels = c("VT", "FR", "JP",
#                                                 "CH", "GA", "MU"))) |> 
#   ggplot() +
#   geom_point(aes(x = lat,
#                  y = LT50,
#                  fill = genotype),
#              shape = 21,
#              color = "grey90",
#              size = 6) +
#   scale_fill_manual(values = c("lightblue", "lightblue3", "lightblue4",
#                                 "orchid4", "orchid", "orchid2" )) +
#   scale_x_continuous(labels = function(x) paste0(x, "°N"),
#                      name = "Lattitude") +
#   scale_y_continuous(labels = function(x) paste0(x, "°C"),
#                      name = expression("LT"[50])) +
#   theme_minimal() +
#   theme()

p_lt_temp <- lat_long |> 
  mutate(genotype = factor(genotype, levels = c("VT", "FR", "JP",
                                                "CH", "GA", "MU"))) |> 
  ggplot(aes(x = wc2.1_10m_bio_5,
             y = LT50)) +
  stat_smooth(method = "lm",
              color = "grey50",
              linetype = 2,
              se = FALSE,
              fullrange = TRUE) +
  geom_errorbar(aes(ymin = LT50 - std.error,
                    ymax = LT50 + std.error,
                    color = genotype)) +
  geom_point(aes(fill = genotype),
             shape = 21,
             color = "grey80",
             size = 4) +
  scale_fill_manual(values = colors_genotypes) +
  scale_color_manual(values = colors_genotypes) +
  scale_x_continuous(labels = scales::label_number(accuracy = 1, suffix = "°C"),
                     limits = c(24, 35),
                     breaks = c(25,  30, 35),
                     name = "Maximum temperature of the warmest month") +
  scale_y_continuous(labels = scales::label_number(accuracy = 1, suffix = "°C"),
                     limits = c(33.5, 36.5),
                     breaks = c(34, 35, 36),
                     name = expression("Acute heat tolerance (LT"[50]*")")) +
  theme_minimal(base_size = 16) +
  theme()


lat_long %>% 
  lm(LT50 ~ wc2.1_10m_bio_5,
     data = .) |> 
  broom::tidy()
lat_long %>% 
  lm(LT50 ~ wc2.1_10m_bio_5,,
     data = .) |> 
  broom::glance()

cowplot::plot_grid(
  p_surv,
  p_lt_temp,
  align = "v",
  nrow = 2
)

ggsave(
  here::here("output/figs/metabolomics/survival_lt50.pdf"),
  width = 7, 
  height = 7,
  units = "in"
)
ggsave(
  here::here("output/figs/metabolomics/survival_lt50.svg"),
  width = 7, 
  height = 7,
  units = "in"
)
