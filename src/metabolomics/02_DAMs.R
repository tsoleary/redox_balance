# ------------------------------------------------------------------------------
# Differential abundant metabolites
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)
library(purrr)
library(lmerTest)

# Load data
dat <- readRDS(here::here("data/processed/metabolomics/dat.rds")) |> 
  unnest(data) |> 
  filter(Region != "QC") |> 
  nest()

# Model formulas to compare ----------------------------------------------------
model_formula_nested <- NormIntensity ~ 
  Region * Temperature + (1 | Region:Genotype)
model_formula_random <- NormIntensity ~ 
  Region * Temperature + (1 | Genotype)
model_formula_fixed <- NormIntensity ~ 
  Region * Temperature

# Run the three models ---------------------------------------------------------
dat_models <- dat |>
  mutate(model_nested = map(data, ~lmer(model_formula_nested, data = .)),
         model_random = map(data, ~lmer(model_formula_random, data = .)),
         model_fixed = map(data, ~lm(model_formula_fixed, data = .))) |> 
  mutate(model_nested_AIC = map_dbl(model_nested, ~AIC(.)),
         model_random_AIC = map_dbl(model_random, ~AIC(.)),
         model_fixed_AIC = map_dbl(model_fixed, ~AIC(.))) |> 
  mutate(best_model = pmap_chr(list(model_fixed_AIC, 
                                    model_random_AIC, 
                                    model_nested_AIC),
      ~ {AICs <- c(fixed = ..1, random = ..2, nested = ..3)
        names(AICs)[which.min(rank(AICs))]}))

# # Count up which model preforms the best across all metabolites ----------------
# dat_models |> 
#   group_by(best_model) |> 
#   tally()

# Adjust p-values and find metabolites with significant terms ------------------
results <- dat_models |> 
  select(Metabolite, model_random) |> 
  mutate(model_random_tidy = map(model_random, ~broom.mixed::tidy(.))) |> 
  unnest(model_random_tidy) |> 
  group_by(term) |> 
  mutate(p.value.adj = p.adjust(p.value, method = "fdr"))

# Differentially abundant metabolites ------------------------------------------
DAMs <- results |> 
  filter(term != "(Intercept)") |> 
  filter(p.value.adj < 0.05)
# 
# ################################################################################
# ### FIXED model results --------------------------------------------------------
# ################################################################################
# 
# # Adjust p-values and find metabolites with significant terms ------------------
# results_fixed <- dat_models |> 
#   select(Metabolite, model_fixed) |> 
#   mutate(model_fixed_tidy = map(model_fixed, ~broom.mixed::tidy(.))) |> 
#   unnest(model_fixed_tidy) |> 
#   group_by(term) |> 
#   mutate(p.value.adj = p.adjust(p.value, method = "fdr"))
# 
# # Differentially abundant metabolites ------------------------------------------
# DAMs_fixed <- results_fixed |> 
#   filter(term != "(Intercept)") |> 
#   filter(p.value.adj < 0.05)






# Save -------------------------------------------------------------------------
saveRDS(results, here::here("output/metabolomics/results.rds"))
saveRDS(DAMs, here::here("output/metabolomics/DAMs.rds"))


# Maybe let's try to create some plots to compare the models -------------------
results_all <- full_join(
  results_fixed, 
  results_random, 
  by = c("Metabolite", "term"), 
  suffix = c(".fixed", ".random")
)


results_all |> 
  filter(term != "(Intercept)") |> 
  ggplot() +
  geom_point(aes(x = statistic.fixed,
                 y = statistic.random)) +
  scale_y_log10() +
  scale_x_log10() +
  theme_minimal()


results_all |> 
  filter(term != "(Intercept)") |> 
  ggplot(aes(x = p.value.adj.fixed,
             y = p.value.adj.random)) +
  geom_point() +
  geom_vline(xintercept = 0.05) +
  geom_hline(yintercept = 0.05) +
  theme_minimal()

results_all |> 
  filter(term != "(Intercept)") |> 
  filter(p.value.adj.fixed < 0.05 | p.value.adj.random < 0.05) |>   
  ggplot(aes(x = p.value.adj.fixed,
             y = p.value.adj.random)) +
  ggrepel::geom_label_repel(aes(label = paste(Metabolite, term))) +
  geom_point() +
  geom_vline(xintercept = 0.05) +
  geom_hline(yintercept = 0.05) +
  theme_minimal()

results_all |> 
  filter(term != "(Intercept)") |> 
  ggplot() +
  geom_point(aes(x = p.value.adj.fixed,
                 y = p.value.adj.random)) +
  geom_vline(xintercept = 0.05) +
  geom_hline(yintercept = 0.05) +
  scale_y_log10() +
  scale_x_log10() +
  theme_minimal()

results_random |> 
  filter(effect == "ran_pars")


r2s <- dat_models |> 
  select(Metabolite, model_fixed, model_random) |> 
  mutate(
    r2_marginal_fixed = map_dbl(model_fixed, ~pluck(performance::r2(.), "R2"))) |> 
  mutate(
    r2_values_random = map(model_random, ~performance::r2(.)),
    r2_marginal_random = map_dbl(r2_values_random, "R2_marginal"),
    r2_conditional_random = map_dbl(r2_values_random, "R2_conditional")) |> 
  select(Metabolite, r2_marginal_fixed, r2_marginal_random, r2_conditional_random)


# Compare marginal r2
r2s |> 
ggplot(aes(x = r2_marginal_fixed, 
           y = r2_marginal_random)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  geom_point(size = 3, shape = 21, fill = "grey80") +
  theme_minimal()

r2s |> 
  ggplot(aes(x = r2_marginal_fixed, 
             y = r2_conditional_random)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  geom_point(size = 3, shape = 21, fill = "grey80") +
  ggrepel::geom_label_repel(aes(label = Metabolite)) +
  theme_minimal()



