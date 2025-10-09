# ------------------------------------------------------------------------------
# script_description
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)
source("src/_theme_colors.R")


# Load data
dat <- readRDS(here::here("data/processed/metabolomics/dat.rds"))
dat_lt50 <- readRDS(here::here("data/processed/survival/lt50.rds"))


################################################################################
################################################################################
################################################################################
################################################################################
# dat_redox <- readRDS(here::here("data/processed/metabolomics/dat_norm.rds")) |>
#   filter(Region  != "QC",
#          str_detect(Metabolite, "Ratio")) |>
#   mutate(log_int = ifelse(
#     Metabolite == "GSSG/GSH-NEM Ratio",
#     1 / pareto_scaled,
#     pareto_scaled),
#   log_int = log10(NormIntensity))
# 
# # nest by metabolite & genotype ---------------------------------------------
# geno_tests <- dat_redox |>
#   group_by(Metabolite, Genotype) |>
#   nest() |>
#   ungroup()
# 
# # paired t‐test of 32°C vs 25°C per genotype × metabolite --------------------
# results_geno <- geno_tests |>
#   mutate(
#     ttest = map(data, ~ {
#       df <- .x
#       x <- df$log_int[df$Temperature == "32°C"]
#       y <- df$log_int[df$Temperature == "25°C"]
#       t.test(x, y, paired = TRUE)
#     }),
#     tidy = map(ttest, broom::tidy)
#   ) |>
#   dplyr::select(Metabolite, Genotype, tidy) |>
#   unnest(tidy) |>
#   dplyr::select(Metabolite, Genotype,
#          estimate,
#          statistic,
#          p.value) |>
#   group_by(Metabolite) |>
#   dplyr::mutate(padj = p.adjust(p.value, method = "BH")) |>
#   ungroup() |>
#   filter(p.value < 0.05)
################################################################################
################################################################################
################################################################################
# 
# # Calculate redox shift
# dat_redox <- dat |> 
#   filter(str_detect(Metabolite, "Ratio")) |> 
#   unnest(cols = data) |> 
#   filter(Genotype != "QC") |> 
#   mutate(NormIntensity = ifelse(Metabolite == "GSSG/GSH-NEM Ratio", 
#                                 NormIntensity^(-1),
#                                 NormIntensity)) |> 
#   mutate(Metabolite = ifelse(Metabolite == "GSSG/GSH-NEM Ratio",
#                             "GSH/GSSG Ratio", Metabolite)) |> 
#   group_by(Locale, Metabolite)  |> 
#   mutate(log_int = ifelse(Metabolite == "GSH/GSSG Ratio", 
#                           log(NormIntensity), NormIntensity)) %>%
#   group_by(Locale, Metabolite, Temperature) %>%
#   summarise(
#     mean_log = mean(log_int),
#     se_log = sd(log_int)/sqrt(n()),
#     .groups = "drop_last") |> 
#   pivot_wider(
#     names_from = Temperature,
#     values_from = c(mean_log, se_log)) |> 
#   transmute(
#     Locale,
#     Metabolite,
#     redox_shift = `mean_log_32°C` - `mean_log_25°C`,
#     se   = sqrt(`se_log_32°C`^2 + `se_log_25°C`^2)) |> 
#   left_join(dat_lt50, by = c("Locale" = "genotype"))

# Pareto scaled data
dat_redox <- dat_norm |> 
  filter(str_detect(Metabolite, "Ratio")) |> 
  filter(Genotype != "QC") |> 
  mutate(pareto_scaled = ifelse(Metabolite == "GSSG/GSH-NEM Ratio", 
                                pareto_scaled*(-1),
                                pareto_scaled)) |> 
  mutate(Metabolite = ifelse(Metabolite == "GSSG/GSH-NEM Ratio",
                             "GSH/GSSG Ratio", Metabolite)) |> 
  group_by(Locale, Metabolite)  |> 
  group_by(Locale, Metabolite, Temperature) %>%
  summarise(
    mean = mean(pareto_scaled),
    se = sd(pareto_scaled)/sqrt(n()),
    .groups = "drop_last") |> 
  pivot_wider(
    names_from = Temperature,
    values_from = c(mean, se)) |> 
  transmute(
    Locale,
    Metabolite,
    redox_shift = `mean_32°C` - `mean_25°C`,
    se = sqrt(`se_32°C`^2 + `se_25°C`^2)) |> 
  left_join(dat_lt50, by = c("Locale" = "genotype"))


p_gsh <- dat_redox |> 
  filter(Metabolite == "GSH/GSSG Ratio") |> 
  ggplot(aes(x = redox_shift,
             y = estimate,
             color = Locale)) +
  geom_vline(xintercept = 0,
             color = "grey20") +
  geom_errorbar(aes(xmin = redox_shift - 1.96*se,
                    xmax = redox_shift + 1.96*se),
                alpha = c(0.5, 1, 1, 1, 0.5, 1),
                width = 0.05) +
  geom_point(aes(fill = Locale),
             shape = 21,
             alpha = c(0.5, 1, 1, 1, 0.5, 1),
             color = "grey20",
             size = 4) +
  # geom_smooth(method = "lm",
  #             se = FALSE,
  #             color = "grey20",
  #             fullrange = TRUE,
  #             linetype = 2) +
  # ggpmisc::stat_poly_eq(
  #   aes(label = after_stat(rr.label)),
  #   formula = y ~ x,
  #   geom = "label",
  #   method = "lm",
  #   parse = TRUE,
  #   size = 3,
  #   label.x = 1.85,
  #   label.y = 36.8,
  #   color = "grey20") +
  scale_fill_manual(values = colors_genotypes) +
  scale_color_manual(values = colors_genotypes) +
  scale_y_continuous(labels = scales::label_number(accuracy = 1, suffix = "°C"),
                     limits = c(34, 37),
                     breaks = 34:37,
                     name = expression("LT"[50])) +
  scale_x_continuous(limits = c(-2.25, 2.25),
                     name = "Heat-induced shift in redox ratio") +
  theme_minimal() +
  labs(title = "GSH/GSSG") +
  guides(color = guide_legend(reverse = TRUE),
         fill = guide_legend(reverse = TRUE)) +
  theme(legend.title = element_blank())


legend <- cowplot::get_legend(p_gsh) |> 
  cowplot::ggdraw()

p_nadh <- dat_redox |> 
  filter(Metabolite == "NADH/NAD+ Ratio") |> 
  ggplot(aes(x = redox_shift,
             y = estimate,
             color = Locale)) +
  geom_vline(xintercept = 0,
             color = "grey20") +
  geom_errorbar(aes(xmin = redox_shift - 1.96*se,
                    xmax = redox_shift + 1.96*se),
                alpha = c(0.5, 1, 1, 0.5, 0.5, 1),
                width = 0.05) +
  geom_point(aes(fill = Locale),
             shape = 21,
             alpha = c(0.5, 1, 1, 0.5, 0.5, 1),
             color = "grey20",
             size = 4) +
  # geom_smooth(method = "lm",
  #             se = FALSE,
  #             color = "grey20",
  #             fullrange = TRUE,
  #             linetype = 2) +
  # ggpmisc::stat_poly_eq(
  #   aes(label = after_stat(rr.label)),
  #   formula = y ~ x,
  #   geom = "label",
  #   method = "lm",
  #   parse = TRUE,
  #   size = 3,
  #   label.x = 0.185,
  #   label.y = 37,
  #   color = "grey20") +
  scale_fill_manual(values = colors_genotypes) +
  scale_color_manual(values = colors_genotypes) +
  scale_y_continuous(labels = scales::label_number(accuracy = 1, suffix = "°C"),
                     limits = c(34, 37),
                     breaks = 34:37,
                     name = expression("LT"[50])) +
  scale_x_continuous(limits = c(-2.25, 2.25),
                     name = "Heat-induced shift in redox ratio") +
  theme_minimal() +
  labs(title = "NADH/NAD+") +
  theme(legend.position = "none")


p_nadph <- dat_redox |> 
  filter(Metabolite == "NADPH/NADP+ Ratio") |> 
  ggplot(aes(x = redox_shift,
             y = estimate,
             color = Locale)) +
  geom_vline(xintercept = 0,
             color = "grey20") +
  geom_errorbar(aes(xmin = redox_shift - 1.96*se,
                    xmax = redox_shift + 1.96*se),
                alpha = c(0.5, 0.5, 0.5, 0.5, 0.5, 1),
                width = 0.05) +
  geom_point(aes(fill = Locale),
             shape = 21,
             alpha = c(0.5, 0.5, 0.5, 0.5, 0.5, 1),
             color = "grey20",
             size = 4) +
  # geom_smooth(method = "lm",
  #             se = FALSE,
  #             color = "grey20",
  #             fullrange = TRUE,
  #             linetype = 2) +
  # ggpmisc::stat_poly_eq(
  #   aes(label = after_stat(rr.label)),
  #   formula = y ~ x,
  #   geom = "label",
  #   method = "lm",
  #   parse = TRUE,
  #   stroke = 0,
  #   size = 3,
  #   label.x = 0.0175,
  #   label.y = 36,
  #   color = "grey20") +
  scale_fill_manual(values = colors_genotypes) +
  scale_color_manual(values = colors_genotypes) +
  scale_y_continuous(labels = scales::label_number(accuracy = 1, suffix = "°C"),
                     limits = c(34, 37),
                     breaks = 34:37,
                     name = expression("LT"[50])) +
  scale_x_continuous(limits = c(-2.25, 2.25),
                     name = "Heat-induced shift in redox ratio") +
  theme_minimal() +
  labs(title = "NADPH/NADP+") +
  theme(legend.position = "none")


plots <- cowplot::plot_grid(
  p_nadh,
  p_nadph,
  p_gsh +
    theme(legend.position = "none"),
  nrow = 3,
  labels = c("A", "B", "C")
)

cowplot::plot_grid(
  plots,
  legend,
  nrow = 1,
  rel_widths = c(1, 0.35)
)

ggsave(here::here("output/figs/final/redox_shift.png"),
       width = 5.5,
       height = 8, 
       units = "in"
)

ggsave(here::here("output/figs/final/redox_shift.svg"),
  width = 5.5,
  height = 8, 
  units = "in"
)



dat_redox |> 
  group_by(Metabolite) |> 
  nest() |> 
  mutate(fit = map(data, ~ lm(redox_shift ~ estimate, data = .x)),
         tidied = map(fit, broom::tidy)) |> 
  unnest(tidied) |> 
  filter(term == "estimate") |> 
  dplyr::select(Metabolite, estimate, std.error, statistic, p.value) 


