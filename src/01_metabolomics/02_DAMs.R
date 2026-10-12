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
dat <- readRDS(here::here("data/processed/metabolomics/dat_norm.rds")) |> 
  filter(Region != "QC") |>
  nest()

# Model formula with genotype as random effect ---------------------------------
model_formula <- log_transformed ~ Region * Temperature +
  (1 + Temperature | Genotype)

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

results |> 
  #filter(term == "RegionTropical:Temperature32°C") |> 
  filter(p.value.adj < 0.05) |> 
  arrange(term, p.value.adj) |> 
  select(-group) |> 
  print(n = 26) 

DAMs |> 
  group_by(term, statistic > 0) |> 
  tally()

DAMs |> 
  filter(term == "Temperature32°C",
         statistic > 0) |> 
  pull(Metabolite)

# Save -------------------------------------------------------------------------
saveRDS(results, here::here("output/metabolomics/results_log.rds"))
write_csv(results |> 
            arrange(effect, desc(term), p.value.adj), 
          here::here("output/metabolomics/results_log.csv"))
saveRDS(DAMs, here::here("output/metabolomics/DAMs_log.rds"))

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
  dat <- dat |>
    mutate(Temperature = factor(Temperature, levels = c("32°C", "25°C")))

  # Run model
  mod <- lmer(model_formula, data = dat)

  # Get emmeans pairs
  mod |>
    emmeans(~ Temperature | Region) |>
    pairs(infer = c(TRUE, TRUE)) |>
    summary() |>
    as_tibble() |>
    mutate(fold_change      = 10^estimate,
           pct_change       = 100 * (10^estimate - 1),
           log2_fold_change = estimate * log2(10),
           Metabolite       = Metabolite)
}

posthoc <- dat |>
  filter(Metabolite %in% DAMs_int) |>
  mutate(out = map2(data, Metabolite, hs_effect_fct)) |>
  pull(out) |>
  list_rbind() |>
  # BH across the 2 regions x n metabolites tested here
  mutate(padj = p.adjust(p.value, method = "BH")) |>
  select(Metabolite, everything())

posthoc |>
  arrange(padj) |>
  print(n = Inf)

saveRDS(posthoc, here::here("output/metabolomics/posthoc_heat_within_region.rds"))
write_csv(posthoc, here::here("output/metabolomics/posthoc_heat_within_region.csv"))

# Genotype level heat shock response of the redox ratios -----------------------
line_effect_fct <- function(dat) {
  dat <- dat |>
    mutate(Temperature = factor(Temperature, levels = c("32°C", "25°C")))

  lm(update(model_formula, . ~ Genotype * Temperature), data = dat) |>
    emmeans(~ Temperature | Genotype) |>
    pairs(infer = c(TRUE, TRUE)) |>
    summary() |>
    as_tibble()
}

line_effects <- dat |>
  mutate(out = map(data, line_effect_fct)) |>
  dplyr::select(Metabolite, out) |>
  unnest(out) |>
  group_by(Metabolite) |>
  mutate(padj = p.adjust(p.value, method = "BH")) |>
  ungroup()

saveRDS(line_effects, here::here("output/metabolomics/line_effects.rds"))

res_geno_redox_ratio <- line_effects |>
  filter(str_detect(Metabolite, "Ratio")) |>
  dplyr::select(Metabolite, Genotype, estimate, 
                lower.CL, upper.CL, p.value, padj) |>
  arrange(Metabolite, Genotype)


write_csv(res_geno_redox_ratio,
          here::here("output/metabolomics/figs6_genotype_heat_redox_ratio.csv"))

res_geno_redox_ratio |> 
  filter(padj < 0.05)



# Significance letters for Figs. 4, 5, S5 ------------------------------------
group_letters_fct <- function(dat) {
  mod <- lmer(model_formula, data = dat)

  emmeans(mod, ~ Region * Temperature) |>
    multcomp::cld(Letters = letters, adjust = "BH", sort = FALSE) |>
    as_tibble() |>
    mutate(letter = trimws(.group)) |>
    dplyr::select(Region, Temperature, emmean, SE, df, letter)
}

group_letters <- dat |>
  mutate(out = map(data, group_letters_fct)) |>
  dplyr::select(Metabolite, out) |>
  unnest(out)

saveRDS(group_letters, here::here("output/metabolomics/group_letters.rds"))
write_csv(group_letters, here::here("output/metabolomics/group_letters.csv"))

# The metabolites plotted in Figs. 4, 5 and S5
group_letters |>
  filter(Metabolite %in% c("NADH/NAD+ Ratio", "NADPH/NADP+ Ratio",
                           "GSH/GSSG Ratio", "AMP", "CMP", "GMP",
                           "NADH", "NAD+", "NADPH", "NADP+",
                           "GSSG", "GSH")) |>
  arrange(Metabolite, Region, Temperature) |>
  print(n = Inf)


# Variance partitioning across metabolites ------------------------------------
# Variance attributable to each fixed effect and variance component

# Stacking order: fixed effects first, then the variance components
comp_levels <- c("Region × Temperature",
                 "Temperature",
                 "Region",
                 "Genotype (baseline)",
                 "Genotype × Temperature",
                 "Residual")

# Extract variance components ---------------------------------------------------
var_comp_fct <- function(dat) {

  # Sum-to-zero contrasts, so the fixed variance splits additively
  old <- options(contrasts = c("contr.sum", "contr.poly"))
  on.exit(options(old), add = TRUE)

  dat <- dat |>
    mutate(Region      = factor(Region, levels = c("Temperate", "Tropical")),
           Temperature = factor(Temperature, levels = c("25°C", "32°C")),
           Genotype    = factor(Genotype))

  mod <- lmer(model_formula, data = dat)

  # Fixed effects, split by term
  X <- getME(mod, "X")
  b <- fixef(mod)
  cn    <- colnames(X)
  has_r <- str_detect(cn, "^Region")      | str_detect(cn, ":Region")
  has_t <- str_detect(cn, "^Temperature") | str_detect(cn, ":Temperature")

  grp <- case_when(
    cn == "(Intercept)" ~ NA_character_,
    has_r & has_t       ~ "Region × Temperature",
    has_r               ~ "Region",
    has_t               ~ "Temperature")

  fx <- tibble(
    component = comp_levels[1:3],
    variance  = map_dbl(comp_levels[1:3], function(g) {
      cols <- which(!is.na(grp) & grp == g)
      var(drop(X[, cols, drop = FALSE] %*% b[cols]))
    }))

  # Random effects and residual
  rn <- as.data.frame(VarCorr(mod)) |>
    filter(is.na(var2)) |>
    transmute(
      component = case_when(
        grp  == "Residual"    ~ "Residual",
        var1 == "(Intercept)" ~ "Genotype (baseline)",
        TRUE                  ~ "Genotype × Temperature"),
      variance = vcov)

  bind_rows(fx, rn)
}

var_comp <- dat |>
  mutate(out = map(data, var_comp_fct)) |>
  dplyr::select(Metabolite, out) |>
  unnest(out) |>
  group_by(Metabolite) |>
  mutate(prop_var = variance / sum(variance)) |>
  ungroup() |>
  mutate(component = factor(component, levels = comp_levels))

write_csv(var_comp, here::here("output/metabolomics/variance_partition.csv"))

# Median share of each component
var_comp |>
  group_by(component) |>
  summarise(median_prop = median(prop_var), .groups = "drop") |>
  print()

# Plot --------------------------------------------------------------------------

# Order by the region x temperature share
ord <- var_comp |>
  filter(component == "Region × Temperature") |>
  arrange(prop_var) |>
  pull(Metabolite)

colors_var <- c("Region × Temperature"    = "#9B6AAE",
                "Temperature"                 = "#C1452E",
                "Region"                      = "#12617F",
                "Genotype (baseline)"         = "grey70",
                "Genotype × Temperature"  = "grey80",
                "Residual"                    = "grey90")

p_var <- var_comp |>
  mutate(Metabolite = factor(Metabolite, levels = ord)) |>
  ggplot(aes(y = Metabolite, x = prop_var, fill = component)) +
  # reverse = TRUE anchors the fixed effects at zero
  geom_col(width = 0.75,
           position = position_stack(reverse = TRUE),
           colour = "white",
           linewidth = 0.25) +
  scale_fill_manual(values = colors_var, name = NULL) +
  scale_x_continuous(name = "Proportion of variance",
                     labels = scales::percent,
                     expand = c(0, 0)) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_minimal(base_size = 11) +
  theme(axis.title.y       = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        panel.grid.major.x = element_line(linewidth = 0.25, colour = "grey88"),
        legend.position    = "top",
        legend.key.size    = unit(10, "pt"),
        plot.margin        = margin(5.5, 14, 5.5, 5.5, "pt"))

ggsave(plot = p_var,
       here::here("output/figs/final/fig_variance_partition.pdf"),
       height = 7,
       width = 7,
       units = "in")
ggsave(plot = p_var,
       here::here("output/figs/final/fig_variance_partition.png"),
       height = 7,
       width = 7,
       units = "in")
ggsave(plot = p_var,
       here::here("output/figs/final/fig_variance_partition.svg"),
       height = 7,
       width = 7,
       units = "in")
