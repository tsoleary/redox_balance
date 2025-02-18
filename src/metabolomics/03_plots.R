# ------------------------------------------------------------------------------
# Summary plots for metabolomics
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)

# Load data
dat <- readRDS(here::here("data/processed/metabolomics/dat.rds"))
results <- readRDS(here::here("output/metabolomics/results.rds"))
DAMs <- readRDS(here::here("output/metabolomics/DAMs.rds"))

# Set directory for output of figures
fig_dir <- "output/figs/metabolomics"

# Define plotting colors for this script
color_sig <- "red"

# Summary plots of the statistics ----------------------------------------------

# Volcano plot of region effect
p_volcano_reg <- results |> 
  filter(term == "RegionTropical") |> 
  ggplot(aes(x = statistic,
             y = -log10(p.value.adj))) +
  geom_point(aes(fill = p.value.adj < 0.05,
                 size = p.value.adj < 0.05),
             shape = 21,
             color = "grey20") +
  ggrepel::geom_label_repel(aes(label = 
                                  ifelse(p.value.adj < 0.05, Metabolite, NA)),
                            size = 3,
                            alpha = 0.9,
                            color = color_sig) +
  scale_fill_manual(values = c("grey80", color_sig)) +
  scale_size_manual(values = c(2, 3)) +
  geom_hline(yintercept = -log10(0.05),
             linetype = 2) +
  geom_vline(xintercept = 0) +
  scale_x_continuous(limits = c(-5, 5)) +
  labs(title = "Region") +
  theme_minimal() +
  theme(legend.position = "none")

# Volcano plot of temperature effect
p_volcano_temp <- results |> 
  filter(term == "Temperature32°C") |> 
  ggplot(aes(x = statistic,
             y = -log10(p.value.adj))) +
  geom_point(aes(fill = p.value.adj < 0.05,
                 size = p.value.adj < 0.05),
             shape = 21,
             color = "grey20") +
  scale_fill_manual(values = c("grey80", color_sig)) +
  scale_size_manual(values = c(2, 3)) +
  geom_hline(yintercept = -log10(0.05),
             linetype = 2) +
  geom_vline(xintercept = 0) +
  ggrepel::geom_label_repel(aes(label = 
                                  ifelse(p.value.adj < 0.05, Metabolite, NA)),
                            size = 3,
                            alpha = 0.8,
                            fill = "white",
                            color = color_sig) +
  scale_x_continuous(limits = c(-12, 12)) +
  labs(title = "Temperature") +
  theme_minimal() +
  theme(legend.position = "none")

# Volcano plot of interaction
p_volcano_int <- results |> 
  filter(term == "RegionTropical:Temperature32°C") |> 
  ggplot(aes(x = statistic,
             y = -log10(p.value.adj))) +
  geom_point(aes(fill = p.value.adj < 0.05,
                 size = p.value.adj < 0.05),
             shape = 21,
             color = "grey20") +
  scale_fill_manual(values = c("grey80", color_sig)) +
  scale_size_manual(values = c(2, 3)) +
  geom_hline(yintercept = -log10(0.05),
             linetype = 2) +
  geom_vline(xintercept = 0) +
  ggrepel::geom_label_repel(aes(label = 
                                  ifelse(p.value.adj < 0.05, Metabolite, NA)),
                            size = 3,
                            alpha = 0.8,
                            fill = "white",
                            color = color_sig) +
  scale_x_continuous(limits = c(-6.25, 6.25)) +
  theme_minimal() +
  labs(title = "Interaction (Region x Temperature)") +
  theme(legend.position = "none")


volcanos <- cowplot::plot_grid(
  p_volcano_reg,
  p_volcano_temp,
  p_volcano_int,
  nrow = 1
)

# Save
ggsave(plot = volcanos,
       here::here(fig_dir, "gssg", "volcano_plots.pdf"),
       height = 6,
       width = 15,
       units = "in")

# Individual plot examples -----------------------------------------------------

plot_legend <- dat$data[[1]] |>
  ggplot(aes(y = NormIntensity,
             x = Region,
             fill = Temperature)) +
  geom_boxplot() +
  ggbeeswarm::geom_beeswarm(shape = 21,
                            size = 3,
                            cex = 2,
                            dodge.width = 0.75) +
  scale_fill_manual(values = c("grey90", "firebrick")) +
  expand_limits(y = 0) +
  theme_minimal() +
  theme(legend.key.size = unit(1, 'cm'), 
        legend.key.height = unit(1, 'cm'),
        legend.key.width = unit(1, 'cm'),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 12))

# Extract legend
legend <- cowplot::get_legend(plot_legend)

# Categorize the metabolites and reorder the data ------------------------------

# redox_metabolites <- c(
#   "Glutathione-NEM","GSSG", "GSH-NEM/GSSG Ratio", 
#   "NADPH", "NADP+", "NADPH/NADP+ Ratio", 
#   "NADH", "NAD+", "NADH/NAD+ Ratio"
# )

redox_metabolites <- c(
  "Glutathione-NEM","GSSG", "GSSG/GSH-NEM Ratio", 
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
    ggplot(aes(y = NormIntensity,
                        x = Region,
                        fill = Temperature)) +
    geom_boxplot() +
    ggbeeswarm::geom_beeswarm(shape = 21,
                              size = 3,
                              cex = 2,
                              dodge.width = 0.75) +
    scale_fill_manual(values = c("grey90", "firebrick")) +
    expand_limits(y = 0) +
    theme_minimal() +
    labs(title = dat$Metabolite[[i]], x = element_blank(), y = "Normalized intensity") +
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
    ggplot(aes(y = NormIntensity,
               x = Region,
               fill = Temperature)) +
    geom_boxplot() +
    ggbeeswarm::geom_beeswarm(shape = 21,
                              size = 3,
                              cex = 2,
                              dodge.width = 0.75) +
    scale_fill_manual(values = c("grey90", "firebrick")) +
    expand_limits(y = 0) +
    theme_minimal() +
    labs(title = dat$Metabolite[[i]], x = element_blank(), y = "Normalized intensity") +
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

# Correct Glutathione ratio

breaks <- 10^(-10:10)
minor_breaks <- rep(1:9, 21)*(10^rep(-10:10, each = 9))

plots[[which(dat$Metabolite == "GSSG/GSH-NEM Ratio")]] <- 
  dat$data[[which(dat$Metabolite == "GSSG/GSH-NEM Ratio")]] |>
  ggplot(aes(y = NormIntensity^(-1),
             x = Region,
             fill = Temperature)) +
  geom_boxplot() +
  ggbeeswarm::geom_beeswarm(shape = 21,
                            size = 3,
                            cex = 2,
                            dodge.width = 0.75) +
  scale_fill_manual(values = c("grey90", "firebrick")) +
  labs(title = "GSH-NEM/GSSG Ratio", x = element_blank(), y = "Normalized intensity") +
  scale_y_continuous(trans = "log2",
                     breaks = breaks,
                     minor_breaks = minor_breaks) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 20, color = "forestgreen"),
        axis.text.x = element_text(face = "bold", size = 12),
        axis.title.y = element_text(size = 10, face = "plain"),
        legend.position = "none",
        plot.margin = unit(c(1, 1, 1, 1), "cm"))


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
    ggplot(aes(y = NormIntensity,
               x = Region,
               fill = Temperature)) +
    geom_boxplot() +
    ggbeeswarm::geom_beeswarm(shape = 21,
                              size = 3,
                              cex = 2,
                              dodge.width = 0.75) +
    scale_fill_manual(values = c("grey90", "firebrick")) +
    expand_limits(y = 0) +
    theme_minimal() +
    labs(title = dat$Metabolite[[i]], x = element_blank(), y = "Normalized intensity") +
    theme(axis.text.x = element_text(face = "bold", size = 12),
          axis.title.y = element_text(size = 10, face = "plain"),
          legend.position = "none",
          plot.margin = unit(c(1, 1, 1, 1), "cm")) 
  
  if (dat_ordered$Metabolite[[i]] %in% unique(DAMs$Metabolite)) {
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
    ggplot(aes(y = NormIntensity,
               x = Region,
               fill = Temperature)) +
    geom_boxplot() +
    ggbeeswarm::geom_beeswarm(shape = 21,
                              size = 3,
                              cex = 2,
                              dodge.width = 0.75) +
    scale_fill_manual(values = c("grey90", "firebrick")) +
    expand_limits(y = 0) +
    theme_minimal() +
    labs(title = dat$Metabolite[[i]], x = element_blank(), y = "Normalized intensity") +
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
    ggplot(aes(y = NormIntensity,
               x = Region,
               fill = Temperature)) +
    geom_boxplot() +
    ggbeeswarm::geom_beeswarm(shape = 21,
                              size = 3,
                              cex = 2,
                              dodge.width = 0.75) +
    scale_fill_manual(values = c("grey90", "firebrick")) +
    expand_limits(y = 0) +
    theme_minimal() +
    labs(title = dat$Metabolite[[i]], x = element_blank(), y = "Normalized intensity") +
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


