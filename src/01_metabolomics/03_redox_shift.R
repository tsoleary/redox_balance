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

# Per-line heat-shock effects, from 02_DAMs.R ----------------------------------
dat_redox <- read_csv(
  here::here("output/metabolomics/figs6_genotype_heat_redox_ratio.csv")) |>
  mutate(Locale = case_when(Genotype == "VT" ~ "Vermont, USA",
                            Genotype == "FR" ~ "Montpellier, France",
                            Genotype == "JP" ~ "Shiojiri, Japan",
                            Genotype == "CH" ~ "Chiapas, Mexico",
                            Genotype == "GA" ~ "Accra, Ghana",
                            Genotype == "BO" ~ "Mumbai, India"),
         Locale = factor(Locale, levels = c("Vermont, USA",
                                            "Montpellier, France",
                                            "Shiojiri, Japan",
                                            "Chiapas, Mexico",
                                            "Accra, Ghana",
                                            "Mumbai, India")),
         redox_shift = estimate,
         sig = padj < 0.05) |>
  dplyr::select(-estimate) |>
  left_join(dat_lt50 |> dplyr::select(genotype, estimate),
            by = c("Locale" = "genotype"))


p_gsh <- dat_redox |> 
  filter(Metabolite == "GSH/GSSG Ratio") |> 
  ggplot(aes(x = redox_shift,
             y = estimate,
             color = Locale)) +
  geom_vline(xintercept = 0,
             color = "grey20") +
  geom_errorbar(aes(xmin = lower.CL,
                    xmax = upper.CL,
                    alpha = sig),
                width = 0.05) +
  geom_point(aes(fill = Locale, alpha = sig),
             shape = 21,
             color = "grey20",
             size = 4) +
  scale_fill_manual(values = colors_genotypes) +
  scale_color_manual(values = colors_genotypes) +
  scale_alpha_manual(values = c(`FALSE` = 0.45, `TRUE` = 1), guide = "none") +
  scale_y_continuous(labels = scales::label_number(accuracy = 1, suffix = "°C"),
                     limits = c(34, 37),
                     breaks = 34:37,
                     name = expression("LT"[50])) +
  scale_x_continuous(limits = c(-1.5, 1.5),
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
  geom_errorbar(aes(xmin = lower.CL,
                    xmax = upper.CL,
                    alpha = sig),
                width = 0.05) +
  geom_point(aes(fill = Locale, alpha = sig),
             shape = 21,
             color = "grey20",
             size = 4) +
  scale_fill_manual(values = colors_genotypes) +
  scale_color_manual(values = colors_genotypes) +
  scale_alpha_manual(values = c(`FALSE` = 0.45, `TRUE` = 1), guide = "none") +
  scale_y_continuous(labels = scales::label_number(accuracy = 1, suffix = "°C"),
                     limits = c(34, 37),
                     breaks = 34:37,
                     name = expression("LT"[50])) +
  scale_x_continuous(limits = c(-1.5, 1.5),
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
  geom_errorbar(aes(xmin = lower.CL,
                    xmax = upper.CL,
                    alpha = sig),
                width = 0.05) +
  geom_point(aes(fill = Locale, alpha = sig),
             shape = 21,
             color = "grey20",
             size = 4) +
  scale_fill_manual(values = colors_genotypes) +
  scale_color_manual(values = colors_genotypes) +
  scale_alpha_manual(values = c(`FALSE` = 0.45, `TRUE` = 1), guide = "none") +
  scale_y_continuous(labels = scales::label_number(accuracy = 1, suffix = "°C"),
                     limits = c(34, 37),
                     breaks = 34:37,
                     name = expression("LT"[50])) +
  scale_x_continuous(limits = c(-1.5, 1.5),
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

ggsave(
  here::here("output/figs/final/fig_genotype_redox_shift.pdf"),
  height = 8,
  width = 5.5,
  units = "in"
)
ggsave(
  here::here("output/figs/final/fig_genotype_redox_shift.png"),
  height = 8,
  width = 5.5,
  units = "in"
)
ggsave(
  here::here("output/figs/final/fig_genotype_redox_shift.svg"),
  height = 8,
  width = 5.5,
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


