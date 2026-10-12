# ------------------------------------------------------------------------------
# Summary plots for metabolomics
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)
source("src/_theme_colors.R")

# Load data
dat <- readRDS(here::here("data/processed/metabolomics/dat.rds"))
dat_norm <- readRDS(here::here("data/processed/metabolomics/dat_norm.rds"))
results <- readRDS(here::here("output/metabolomics/results.rds"))
DAMs <- readRDS(here::here("output/metabolomics/DAMs.rds"))

# Set directory for output of figures
fig_dir <- "output/figs/metabolomics"

# Define plotting colors for this script
color_sig <- "red"

# Redox ratios -----------------------------------------------------------------

p_nadh <- dat$data[[which(dat$Metabolite == "NADH/NAD+ Ratio")]] |>
  filter(Region != "QC") |>
  ggplot(aes(y = log10(NormIntensity),
             x = Region,
             fill = Temperature)) +
  geom_boxplot() +
  ggbeeswarm::geom_beeswarm(shape = 21,
                            size = 3,
                            cex = 2.5,
                            alpha = 0.8,
                            dodge.width = 0.75) +
  scale_fill_manual(values = c("grey90", "firebrick"),
                    name = "Temperature") +
  scale_y_continuous(
    name = expression("log"[10]*"(NADH / NAD"^"+"*")")) +
  theme_minimal(base_size = 14) +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "top")

legend <- cowplot::get_legend(p_nadh) 

p_nadh <- p_nadh +
  theme(legend.position = "none")


p_nadph <- dat$data[[which(dat$Metabolite == "NADPH/NADP+ Ratio")]] |>
  filter(Region != "QC") |>
  ggplot(aes(y = log10(NormIntensity),
             x = Region,
             fill = Temperature)) +
  geom_boxplot() +
  ggbeeswarm::geom_beeswarm(shape = 21,
                            size = 3,
                            cex = 2.5,
                            alpha = 0.8,
                            dodge.width = 0.75) +
  scale_fill_manual(values = c("grey90", "firebrick")) +
  scale_y_continuous(
    name = expression("log"[10]*"(NADPH / NADP"^"+"*")")) +
  theme_minimal(base_size = 14) +
  theme(axis.title.x = element_blank(),
        legend.position = "none")


p_gsh <- dat$data[[which(dat$Metabolite == "GSH/GSSG Ratio")]] |>
  filter(Region != "QC") |>
  ggplot(aes(y = log10(NormIntensity),
             x = Region,
             fill = Temperature)) +
  geom_boxplot() +
  ggbeeswarm::geom_beeswarm(shape = 21,
                            size = 3,
                            cex = 2.5,
                            alpha = 0.8,
                            dodge.width = 0.75) +
  scale_fill_manual(values = c("grey90", "firebrick")) +
  scale_y_continuous(name = expression("log"[10]*"(GSH / GSSG)")) +
  theme_minimal(base_size = 14) +
  theme(axis.title.x = element_blank(),
        legend.position = "none")


# # Vertical
cowplot::plot_grid(
  cowplot::ggdraw(legend),
  p_nadh,
  p_nadph,
  p_gsh,
  labels = c("", "A", "B", "C"),
  rel_heights = c(0.1, 1, 1, 1),
  label_size = 20,
  nrow = 4
)


ggsave(
  here::here("output/figs/final/fig_redox_ratios.pdf"),
  height = 9,
  width = 4,
  units = "in"
)
ggsave(
  here::here("output/figs/final/fig_redox_ratios.png"),
  height = 9,
  width = 4,
  units = "in"
)
ggsave(
  here::here("output/figs/final/fig_redox_ratios.svg"),
  height = 9,
  width = 4,
  units = "in"
)


# Nucleotide monophosphates --------------------------------------------------------------

p_amp <- dat$data[[which(dat$Metabolite == "AMP")]] |>
  filter(Region != "QC") |>
  ggplot(aes(y = log10(NormIntensity),
             x = Region,
             fill = Temperature)) +
  geom_boxplot() +
  ggbeeswarm::geom_beeswarm(shape = 21,
                            size = 3,
                            cex = 2.5,
                            alpha = 0.8,
                            dodge.width = 0.75) +
  scale_fill_manual(values = c("grey90", "firebrick"),
                    name = "Temperature") +
  scale_y_continuous(name = expression("log"[10]*"(AMP)")) +
  theme_minimal(base_size = 16) +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "top")

legend <- cowplot::get_legend(p_amp) 

p_amp <- p_amp +
  theme(legend.position = "none")


p_gmp <- dat$data[[which(dat$Metabolite == "GMP")]] |>
  filter(Region != "QC") |>
  ggplot(aes(y = log10(NormIntensity),
             x = Region,
             fill = Temperature)) +
  geom_boxplot() +
  ggbeeswarm::geom_beeswarm(shape = 21,
                            size = 3,
                            cex = 2.5,
                            alpha = 0.8,
                            dodge.width = 0.75) +
  scale_fill_manual(values = c("grey90", "firebrick")) +
  scale_y_continuous(name = expression("log"[10]*"(GMP)")) +
  theme_minimal(base_size = 16) +
  theme(axis.title.x = element_blank(),
        legend.position = "none")

p_cmp <- dat$data[[which(dat$Metabolite == "CMP")]] |>
  filter(Region != "QC") |>
  ggplot(aes(y = log10(NormIntensity),
             x = Region,
             fill = Temperature)) +
  geom_boxplot() +
  ggbeeswarm::geom_beeswarm(shape = 21,
                            size = 3,
                            cex = 2.5,
                            alpha = 0.8,
                            dodge.width = 0.75) +
  scale_fill_manual(values = c("grey90", "firebrick"),
                    name = "Temperature") +
  scale_y_continuous(name = expression("log"[10]*"(CMP)")) +
  theme_minimal(base_size = 16) +
  theme(axis.title.x = element_blank(),
        legend.position = "none")


# # Vertical
cowplot::plot_grid(
  cowplot::ggdraw(legend),
  p_amp,
  p_cmp,
  p_gmp,
  labels = c("", "A", "B", "C"),
  rel_heights = c(0.1, 1, 1, 1),
  label_size = 20,
  nrow = 4
)

ggsave(
  here::here("output/figs/final/fig_nucleotide_monophosphates.pdf"),
  height = 9,
  width = 4,
  units = "in"
)
ggsave(
  here::here("output/figs/final/fig_nucleotide_monophosphates.png"),
  height = 9,
  width = 4,
  units = "in"
)
ggsave(
  here::here("output/figs/final/fig_nucleotide_monophosphates.svg"),
  height = 9,
  width = 4,
  units = "in"
)


# Individual plot examples -----------------------------------------------------

# Create plot just for the legend later...
plot_legend <- dat$data[[1]] |>
  filter(Region != "QC") |>
  ggplot(aes(y = log10(NormIntensity),
             x = Region,
             fill = Temperature)) +
  geom_boxplot() +
  ggbeeswarm::geom_beeswarm(shape = 21,
                            size = 3,
                            cex = 2,
                            dodge.width = 0.75) +
  scale_fill_manual(values = c("grey90", "firebrick")) +
  scale_y_continuous(
    name = expression("log"[10]*"(Normalized intensity)")) +
  theme_minimal() +
  theme(legend.key.size = unit(1, 'cm'), 
        legend.key.height = unit(1, 'cm'),
        legend.key.width = unit(1, 'cm'),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 12))

# Extract legend
legend <- cowplot::get_legend(plot_legend)



redox_couples <- c(
  "NADH", "NAD+",
  "NADPH", "NADP+",
  "GSH","GSSG"
)


dat |> 
  filter(Metabolite %in% redox_couples) |> 
  unnest(cols = data) |> 
  filter(Region != "QC") |> 
  mutate(Metabolite = factor(Metabolite, levels = redox_couples)) |> 
  arrange(Metabolite) |> 
  ggplot(aes(y = log10(NormIntensity),
             x = Region,
             fill = Temperature)) +
  geom_boxplot() +
  ggbeeswarm::geom_beeswarm(shape = 21,
                            size = 3,
                            cex = 3,
                            alpha = 0.8,
                            dodge.width = 0.75) +
  scale_fill_manual(values = c("grey90", "firebrick")) +
  scale_y_continuous(name = expression("log"[10]*"(Normalized intensity)")) +
  theme_minimal() +
  facet_wrap(~ Metabolite,
             scales = "free_y",
             nrow = 3) +
  theme(axis.title.x = element_blank(),
        legend.position = "top",
        legend.title = element_blank(),
        strip.background = element_rect(color = "grey95", fill = "grey95"))


ggsave(
  here::here("output/figs/final/fig_redox_metabolites.pdf"),
  height = 7.5,
  width = 6.5,
  units = "in"
)
ggsave(
  here::here("output/figs/final/fig_redox_metabolites.png"),
  height = 7.5,
  width = 6.5,
  units = "in"
)
ggsave(
  here::here("output/figs/final/fig_redox_metabolites.svg"),
  height = 7.5,
  width = 6.5,
  units = "in"
)


# Categorize the metabolites and reorder the data ------------------------------

redox_metabolites <- c(
  "GSH","GSSG", "GSH/GSSG Ratio",
  "NADPH", "NADP+", "NADPH/NADP+ Ratio", 
  "NADH", "NAD+", "NADH/NAD+ Ratio"
)

amino_acids <- c( 
  "Alanine", "Leucine", "Isoleucine", "Methionine",
  "Arginine", "Glutamic acid", 
  "Serine", "Tyrosine", "Asparagine", "Glutamine"
)

nucleo_metabolites <- c( 
  "Uracil", "Thymine", "Cytidine", "Adenosine", "Guanosine", 
  "UMP", "CMP", "AMP", "GMP",
  "CDP", "GDP", 
  "CTP", "ATP"
)

misc_intermediate_metabolites <- c(
  "Acetyl-CoA",
  "Acetylcarnitine",
  "Citicoline",
  "Argininosuccinic acid"
)

ordered_metabolites <- c(
  redox_metabolites,
  amino_acids,
  nucleo_metabolites,
  misc_intermediate_metabolites
)

# Reorder metabolites
dat <- dat |> 
  mutate(Metabolite = factor(Metabolite, levels = ordered_metabolites)) |> 
  arrange(Metabolite)


# All metabolites --- Region & Temperature plot --------------------------------

# Empty plot object to populate
plots <- NULL

# Loop through all plots
for (i in 1:nrow(dat)) {
  p <- dat$data[[i]] |> 
    ggplot(aes(y = log10(NormIntensity),
               x = Region,
               fill = Temperature)) +
    geom_boxplot() +
    ggbeeswarm::geom_beeswarm(shape = 21,
                              size = 3,
                              cex = 2,
                              dodge.width = 0.75) +
    scale_fill_manual(values = c("grey90", "firebrick")) +
    theme_minimal() +
    labs(title = dat$Metabolite[[i]], x = element_blank(),
         y = expression("log"[10]*"(Normalized intensity)")) +
    theme(axis.text.x = element_text(face = "bold", size = 12),
          axis.title.y = element_text(size = 10, face = "plain"),
          legend.position = "none",
          plot.margin = unit(c(1, 1, 1, 1), "cm")) 
  
  if (dat$Metabolite[[i]] %in% unique(DAMs$Metabolite)) {
    plots[[i]] <- p +
      theme(plot.title = element_text(face = "bold", size = 20, color = "forestgreen"))
  } else {
    plots[[i]] <- p +
      theme(plot.title = element_text(face = "bold", size = 20, color = "grey80"))
  }
}

# Combine all plots together
plots_all <- cowplot::plot_grid(plotlist = plots)

# Add legend
plots_legend <- cowplot::plot_grid(
  plots_all, 
  legend, 
  rel_widths = c(1, 0.1))

# Save
ggsave(plot = plots_legend,
       here::here(fig_dir, "gssg", "all.pdf"),
       height = 30,
       width = 30,
       units = "in")


# Redox metabolites --- Region & Temperature plot ------------------------------

# Empty plot object to populate
plots <- NULL

# Loop through all plots
for (i in which(dat$Metabolite %in% redox_metabolites)) {
  p <- dat$data[[i]] |> 
    ggplot(aes(y = log10(NormIntensity),
               x = Region,
               fill = Temperature)) +
    geom_boxplot() +
    ggbeeswarm::geom_beeswarm(shape = 21,
                              size = 3,
                              cex = 2,
                              dodge.width = 0.75) +
    scale_fill_manual(values = c("grey90", "firebrick")) +
    theme_minimal() +
    labs(title = dat$Metabolite[[i]], x = element_blank(),
         y = expression("log"[10]*"(Normalized intensity)")) +
    theme(axis.text.x = element_text(face = "bold", size = 12),
          axis.title.y = element_text(size = 8, face = "plain"),
          legend.position = "none",
          plot.margin = unit(c(1, 1, 1, 1), "cm")) 
  
  if (dat$Metabolite[[i]] %in% unique(DAMs$Metabolite)) {
    plots[[i]] <- p +
      theme(plot.title = element_text(face = "bold", size = 20, color = "forestgreen"))
  } else {
    plots[[i]] <- p +
      theme(plot.title = element_text(face = "bold", size = 20, color = "grey80"))
  }
}

# Combine all plots together
plots_all <- cowplot::plot_grid(plotlist = plots)

# Add legend
plots_legend <- cowplot::plot_grid(
  plots_all, 
  legend, 
  rel_widths = c(1, 0.1))

# Save
ggsave(plot = plots_legend,
       here::here(fig_dir, "gssg", "redox_metabolites.pdf"),
       height = 15,
       width = 15,
       units = "in")


# Redox metabolites --- Region & Temperature plot ------------------------------

# Empty plot object to populate
plots <- NULL

# Loop through all plots
for (i in which(dat$Metabolite %in% amino_acids)) {
  p <- dat$data[[i]] |> 
    ggplot(aes(y = log10(NormIntensity),
               x = Region,
               fill = Temperature)) +
    geom_boxplot() +
    ggbeeswarm::geom_beeswarm(shape = 21,
                              size = 3,
                              cex = 2,
                              dodge.width = 0.75) +
    scale_fill_manual(values = c("grey90", "firebrick")) +
    theme_minimal() +
    labs(title = dat$Metabolite[[i]], x = element_blank(),
         y = expression("log"[10]*"(Normalized intensity)")) +
    theme(axis.text.x = element_text(face = "bold", size = 12),
          axis.title.y = element_text(size = 10, face = "plain"),
          legend.position = "none",
          plot.margin = unit(c(1, 1, 1, 1), "cm")) 
  
  if (dat$Metabolite[[i]] %in% unique(DAMs$Metabolite)) {
    plots[[i]] <- p +
      theme(plot.title = element_text(face = "bold", size = 20, color = "forestgreen"))
  } else {
    plots[[i]] <- p +
      theme(plot.title = element_text(face = "bold", size = 20, color = "grey80"))
  }
}


# Combine all plots together
plots_all <- cowplot::plot_grid(
  plotlist = plots[which(dat$Metabolite %in% amino_acids)])

# Add legend
plots_legend <- cowplot::plot_grid(
  plots_all, 
  legend, 
  rel_widths = c(1, 0.1))

# Save
ggsave(plot = plots_legend,
       here::here(fig_dir, "amino_acids.pdf"),
       height = 15,
       width = 20,
       units = "in")


# Nucleotide related metabolites --- Region & Temperature plot -----------------

# Empty plot object to populate
plots <- NULL

# Loop through all plots
for (i in which(dat$Metabolite %in% nucleo_metabolites)) {
  p <- dat$data[[i]] |> 
    ggplot(aes(y = log10(NormIntensity),
               x = Region,
               fill = Temperature)) +
    geom_boxplot() +
    ggbeeswarm::geom_beeswarm(shape = 21,
                              size = 3,
                              cex = 2,
                              dodge.width = 0.75) +
    scale_fill_manual(values = c("grey90", "firebrick")) +
    theme_minimal() +
    labs(title = dat$Metabolite[[i]], x = element_blank(),
         y = expression("log"[10]*"(Normalized intensity)")) +
    theme(axis.text.x = element_text(face = "bold", size = 12),
          axis.title.y = element_text(size = 10, face = "plain"),
          legend.position = "none",
          plot.margin = unit(c(1, 1, 1, 1), "cm")) 
  
  if (dat$Metabolite[[i]] %in% unique(DAMs$Metabolite)) {
    plots[[i]] <- p +
      theme(plot.title = element_text(face = "bold", size = 20, color = "forestgreen"))
  } else {
    plots[[i]] <- p +
      theme(plot.title = element_text(face = "bold", size = 20, color = "grey80"))
  }
}

# Combine all plots together
plots_all <- cowplot::plot_grid(
  plotlist = plots[which(dat$Metabolite %in% nucleo_metabolites)],
  ncol = 5)

# Add legend
plots_legend <- cowplot::plot_grid(
  plots_all, 
  legend, 
  rel_widths = c(1, 0.1))

# Save
ggsave(plot = plots_legend,
       here::here(fig_dir, "nucleo_metabolites.pdf"),
       height = 15,
       width = 25,
       units = "in")


# Misc intermediate metabolites --- Region & Temperature plot ------------------

# Empty plot object to populate
plots <- NULL

# Loop through all plots
for (i in which(dat$Metabolite %in% misc_intermediate_metabolites)) {
  p <- dat$data[[i]] |> 
    ggplot(aes(y = log10(NormIntensity),
               x = Region,
               fill = Temperature)) +
    geom_boxplot() +
    ggbeeswarm::geom_beeswarm(shape = 21,
                              size = 3,
                              cex = 2,
                              dodge.width = 0.75) +
    scale_fill_manual(values = c("grey90", "firebrick")) +
    theme_minimal() +
    labs(title = dat$Metabolite[[i]], x = element_blank(),
         y = expression("log"[10]*"(Normalized intensity)")) +
    theme(axis.text.x = element_text(face = "bold", size = 12),
          axis.title.y = element_text(size = 10, face = "plain"),
          legend.position = "none",
          plot.margin = unit(c(1, 1, 1, 1), "cm")) 
  
  if (dat$Metabolite[[i]] %in% unique(DAMs$Metabolite)) {
    plots[[i]] <- p +
      theme(plot.title = element_text(face = "bold", size = 20, color = "forestgreen"))
  } else {
    plots[[i]] <- p +
      theme(plot.title = element_text(face = "bold", size = 20, color = "grey80"))
  }
}


# Combine all plots together
plots_all <- cowplot::plot_grid(
  plotlist = plots[which(dat$Metabolite %in% misc_intermediate_metabolites)],
  ncol = 2)

# Add legend
plots_legend <- cowplot::plot_grid(
  plots_all, 
  legend, 
  rel_widths = c(1, 0.1))

# Save
ggsave(plot = plots_legend,
       here::here(fig_dir, "misc_intermediate_metabolites.pdf"),
       height = 12,
       width = 12,
       units = "in")


# All other metabolites

redox_metabolites <- c(
  "GSH","GSSG", "GSH/GSSG Ratio",
  "NADPH", "NADP+", "NADPH/NADP+ Ratio", 
  "NADH", "NAD+", "NADH/NAD+ Ratio"
)


dat |> 
  filter(!(Metabolite %in% c(redox_metabolites))) |> 
  filter(!(Metabolite %in% c("AMP", "CMP", "GMP"))) |> 
  unnest(data) |>
  filter(Region != "QC") |>
  ggplot(aes(y = log10(NormIntensity),
             x = Region,
             fill = Temperature)) +
  geom_boxplot() +
  ggbeeswarm::geom_beeswarm(shape = 21,
                            size = 3,
                            cex = 2.5,
                            alpha = 0.8,
                            dodge.width = 0.75) +
  scale_fill_manual(values = c("grey90", "firebrick"),
                    name = "Temperature") +
  scale_y_continuous(
    name = expression("log"[10]*"(Normalized intensity)")) +
  theme_minimal(base_size = 16) +
  facet_wrap(~Metabolite, scale = "free_y") +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "top")
  
dat_norm |> 
  filter(!(Metabolite %in% c(redox_metabolites))) |> 
  filter(!(Metabolite %in% c("AMP", "CMP", "GMP"))) |> 
  filter(Region != "QC") |> 
  ggplot(aes(y = pareto_scaled,
             x = Region,
             fill = Temperature)) +
  geom_boxplot() +
  ggbeeswarm::geom_beeswarm(shape = 21,
                            size = 3,
                            cex = 2.5,
                            alpha = 0.8,
                            dodge.width = 0.75) +
  scale_fill_manual(values = c("grey90", "firebrick"),
                    name = "Temperature") +
  expand_limits(y = 0) +
  theme_minimal(base_size = 16) +
  facet_wrap(~Metabolite) +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "top")

# ALL METABOLITES ------------

redox_metabolites <- c(
  "NADH", "NAD+", "NADH/NAD+",
  "NADPH", "NADP+", "NADPH/NADP+", 
  "GSH", "GSSG", "GSH/GSSG"
)

amino_acids <- c( 
  "Alanine", "Leucine", "Isoleucine", "Methionine",
  "Arginine", "Glutamic acid", 
  "Serine", "Tyrosine", "Asparagine", "Glutamine"
)

nucleo_metabolites <- c( 
  "Uracil", "Thymine", "Cytidine", "Adenosine", "Guanosine", 
  "UMP", "CMP", "AMP", "GMP",
  "CDP", "GDP", 
  "CTP", "ATP"
)

misc_intermediate_metabolites <- c(
  "Acetyl-CoA",
  "Acetylcarnitine",
  "Citicoline",
  "ASA"
)


dat_norm |>
  # Short panel labels
  mutate(Metabolite = str_remove_all(Metabolite, " Ratio")) |>
  mutate(Metabolite = ifelse(Metabolite == "Argininosuccinic acid",
                             "ASA", Metabolite)) |>
  mutate(Metabolite = factor(Metabolite, levels = c(
    redox_metabolites, nucleo_metabolites,
    amino_acids, misc_intermediate_metabolites))) |>
  filter(Region != "QC") |>
  group_by(Metabolite, Region, Temperature) |> 
  summarise(avg = mean(pareto_scaled),
            se = sd(pareto_scaled)/sqrt(n())) |> 
  ggplot(aes(y = avg,
             x = Temperature,
             fill = Region,
             color = Region,
             group = Region)) +
  geom_errorbar(aes(ymin = avg - se,
                    ymax = avg + se),
                width = 0.1) +
  geom_hline(yintercept = 0, 
             color = "grey70",
             linetype = 2) +
  geom_line() +
  geom_point(shape = 21,
             alpha = 0.9,
             color = "grey50",
             size = 2) +
  scale_fill_manual(values = colors_genotypes[c(1, 4)]) +
  scale_color_manual(values = colors_genotypes[c(1, 4)]) +
  scale_y_continuous(name = "Pareto-scaled abundance") +
  scale_x_discrete(expand = c(0.1, 0.1)) +
  theme_minimal(base_size = 16) +
  facet_wrap(~Metabolite) +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "top")


ggsave(
  here::here("output/figs/final/fig_all_metabolites.pdf"),
  height = 9,
  width = 9,
  units = "in",
  bg = "white"
)
ggsave(
  here::here("output/figs/final/fig_all_metabolites.png"),
  height = 9,
  width = 9,
  units = "in",
  bg = "white"
)
ggsave(
  here::here("output/figs/final/fig_all_metabolites.svg"),
  height = 9,
  width = 9,
  units = "in",
  bg = "white"
)


# Temp only results
temp_only <- setdiff(
  DAMs |> filter(term == "Temperature32°C") |> pull(Metabolite),
  DAMs |> filter(term == "RegionTropical:Temperature32°C") |> pull(Metabolite))

int <- DAMs |> filter(term == "RegionTropical:Temperature32°C") |> pull(Metabolite)

dat_norm |> 
  filter(Metabolite %in% int) |>
  mutate(Metabolite = str_remove_all(Metabolite, " Ratio")) |>
  mutate(Metabolite = ifelse(Metabolite == "Argininosuccinic acid",
                             "ASA", Metabolite)) |>
  mutate(Metabolite = factor(Metabolite, levels = c(
    nucleo_metabolites, 
    amino_acids, misc_intermediate_metabolites,
    redox_metabolites))) |> 
  filter(Region != "QC") |> 

  group_by(Metabolite, Region, Temperature) |> 
  summarise(avg = mean(pareto_scaled),
            se = sd(pareto_scaled)/sqrt(n())) |> 
  ggplot(aes(y = avg,
             x = Temperature,
             fill = Region,
             color = Region,
             group = Region)) +
  geom_errorbar(aes(ymin = avg - se,
                    ymax = avg + se),
                width = 0.1) +
  geom_hline(yintercept = 0, 
             color = "grey70",
             linetype = 2) +
  geom_line() +
  geom_point(shape = 21,
             alpha = 0.9,
             color = "grey50",
             size = 2) +
  scale_fill_manual(values = colors_genotypes[c(1, 4)]) +
  scale_color_manual(values = colors_genotypes[c(1, 4)]) +
  scale_y_continuous(name = "Pareto-scaled abundance") +
  scale_x_discrete(expand = c(0.15, 0.15)) +
  theme_minimal(base_size = 16) +
  facet_wrap(~Metabolite, nrow = 2) +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "top")

ggsave(
  here::here("output/figs/metabolomics/interaction.pdf"),
  height = 4.5,
  width = 9,
  units = "in",
  bg = "white"
)
