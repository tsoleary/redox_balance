# ------------------------------------------------------------------------------
# Principle component analysis
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)
library(ggpattern)

# Load data
dat <- readRDS(here::here("data/processed/metabolomics/dat.rds"))
dat_norm <- readRDS(here::here("data/processed/metabolomics/dat_norm.rds"))

# PCA on data ------------------------------------------------------------------

# Pivot data wider for prcomp
dat_wide <- dat_norm |> 
  select(-c(log_transformed, NormIntensity)) |> 
  pivot_wider(values_from = pareto_scaled,
              names_from = Metabolite)

# PCA on pareto scaled data
pca_fit <- dat_wide |> 
  select(where(is.numeric)) |> 
  prcomp(center = TRUE)

# Colors for plotting
region_colors <- c("forestgreen", "pink3")

# Percent variance explained by each PC
pca_scree <- pca_fit |> 
  broom::tidy(matrix = "eigenvalues") |> 
  filter(PC <= 10) |> 
  ggplot() +
  geom_col(aes(x = PC,
               y = percent),
           fill = "grey80",
           color = "grey20") +
  scale_y_continuous(name = "Percent variance explained",
                     labels = scales::percent,
                     expand = c(0, 0)) +
  scale_x_continuous(breaks = 1:10,
                     name = "Principle component") +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank())

# Extract explained variance for PC1 & PC2
pc1_p_var <- round(broom::tidy(pca_fit, matrix = "eigenvalues")$percent[1]*100,
                   digits = 1)
pc2_p_var <- round(broom::tidy(pca_fit, matrix = "eigenvalues")$percent[2]*100,
                   digits = 1)


# First two PCs plotted with Temperature and Region
pca_temp <- pca_fit |> 
  broom::augment(dat_wide) |> 
  ggplot(aes(x = .fittedPC1,
             y = .fittedPC2)) +
  geom_point(aes(fill = Temperature, shape = Region), 
             size = 4, stroke = 1, alpha = 0.8, color = "grey20") +
  scale_fill_manual(values = c("grey80", "firebrick", "grey20")) +
  scale_shape_manual(values = c(15, 21, 24)) +
  scale_x_continuous(name = paste0("PC1 (", pc1_p_var, "%)")) +
  scale_y_continuous(name = paste0("PC2 (", pc2_p_var, "%)")) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  theme_minimal() +
  theme(legend.title = element_blank())


# First two PCs plotted with Temperature and Region
pca_region <- pca_fit |> 
  broom::augment(dat_wide) |> 
  mutate(Region = factor(Region, levels = c("Tropical", "Temperate", "QC"))) |> 
  arrange(desc(Region)) |> 
  ggplot(aes(x = .fittedPC1,
             y = .fittedPC2,
             shape = Temperature, 
             fill = Region)) +
  # Controlling the order of the ellipse plotting
  stat_ellipse(data = function(x) subset(x, Region == "QC"),
               geom = "polygon", alpha = 0.5, color = "grey50") +
  stat_ellipse(data = function(x) subset(x, Region == "Tropical"),
               geom = "polygon", alpha = 0.5, color = "grey50") +
  stat_ellipse(data = function(x) subset(x, Region == "Temperate"),
               geom = "polygon", alpha = 0.5, color = "grey50") +
  geom_point(size = 4, stroke = 1, alpha = 0.8, color = "grey20") +
  scale_fill_manual(values = c("grey10", rev(region_colors))) +
  scale_shape_manual(values = c(21, 24, 15)) +
  scale_x_continuous(name = paste0("PC1 (", pc1_p_var, "%)")) +
  scale_y_continuous(name = paste0("PC2 (", pc2_p_var, "%)")) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  theme_minimal() +
  theme(legend.title = element_blank())

# Rotation plot -------

# Define arrow style for plotting
arrow_style <- arrow(
  angle = 20, 
  ends = "first", 
  type = "closed", 
  length = grid::unit(8, "pt")
)

# plot rotation matrix
lim_dims <- 0.4
pca_rotation <- pca_fit |> 
  broom::tidy(matrix = "rotation") |> 
  pivot_wider(names_from = "PC", names_prefix = "PC", values_from = "value") |> 
  ggplot(aes(PC1, PC2)) +
  geom_segment(xend = 0, yend = 0, arrow = arrow_style) +
  ggrepel::geom_label_repel(aes(label = ifelse(abs(PC1) > 0.21 | abs(PC2) > 0.21, column, NA)),
                            alpha = 0.8,
                            size = 2,
                            hjust = 1, 
                            nudge_x = -0.02,
                            color = "#904C2F") +
  scale_x_continuous(limits = c(-lim_dims-0.2, lim_dims + 0.2)) +
  scale_y_continuous(limits = c(-lim_dims, lim_dims)) +
  theme_void()

# Combine pca scatter plots
pca_plots <- cowplot::plot_grid(
  pca_temp,
  pca_region,
  nrow = 1
)

# Combine scree plot & rotation plot
pca_extra <- cowplot::plot_grid(
  pca_scree,
  pca_rotation,
  rel_widths = c(0.7, 1),
  nrow = 1
)

# Combine all pca plots
pca_all <- cowplot::plot_grid(
  pca_plots,
  pca_extra,
  nrow = 2,
  rel_heights = c(2, 1)
)

# Save PCA plots together
ggsave(plot = pca_all,
       here::here("output/figs/metabolomics/pca.pdf"),
       height = 8,
       width = 12,
       units = "in")

