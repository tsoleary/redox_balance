# ------------------------------------------------------------------------------
# Principle component analysis
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)

# Load data
dat <- readRDS(here::here("data/processed/metabolomics/dat.rds"))
dat_norm <- readRDS(here::here("data/processed/metabolomics/dat_norm.rds"))

# PCA on data ------------------------------------------------------------------

# Pivot data wider for prcomp
dat_wide <- dat_norm |>
  filter(Genotype != "QC") |> 
  select(-c(log_transformed, NormIntensity)) |> 
  pivot_wider(values_from = pareto_scaled,
              names_from = Metabolite)

# PCA on pareto scaled data
pca_fit <- dat_wide |> 
  select(where(is.numeric)) |> 
  prcomp(center = TRUE)

# # Colors for plotting
# region_colors <- c("forestgreen", "pink3")

# Percent variance explained by each PC
pca_scree <- pca_fit |> 
  broom::tidy(matrix = "eigenvalues") |> 
  filter(PC <= 10) |> 
  mutate(PC = factor(PC, levels = 10:1)) |> 
  ggplot() +
  geom_col(aes(y = PC,
               x = percent),
           fill = "grey80",
           color = "grey20") +
  scale_x_continuous(name = "Percent variance explained",
                     labels = scales::percent,
                     expand = c(0, 0),
                     position = "top") +
  scale_y_discrete(breaks = 1:10,
                   name = "Principle component",
                   label = function(x) paste0("PC ", x)) +
  theme_minimal(base_size = 14) +
  theme(axis.title.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())

# Extract explained variance for PC1 & PC2
pc1_p_var <- round(broom::tidy(pca_fit, matrix = "eigenvalues")$percent[1]*100,
                   digits = 1)
pc2_p_var <- round(broom::tidy(pca_fit, matrix = "eigenvalues")$percent[2]*100,
                   digits = 1)

# First two PCs plotted with Temperature and Region
p_pca <- pca_fit |> 
  broom::augment(dat_wide) |> 
  ggplot(aes(x = .fittedPC1,
             y = .fittedPC2,
             fill = Temperature, 
             shape = Region)) +
  stat_ellipse(data = function(x) subset(x, Region == "QC"),
               geom = "polygon", alpha = 0.5, color = "grey50") +
  stat_ellipse(data = function(x) subset(x, Region == "Temperate"),
               geom = "polygon", alpha = 0.5, color = "grey50") +
  stat_ellipse(data = function(x) subset(x, Region == "Tropical"),
               geom = "polygon", alpha = 0.5, color = "grey50") +
  geom_point(size = 4, stroke = 1, alpha = 0.8, color = "grey20") +
  scale_fill_manual(values = c("grey80", "firebrick", "grey20")) +
  scale_shape_manual(values = c( 21, 24)) +
  scale_x_continuous(name = paste0("PC1 (", pc1_p_var, ".0%)")) +
  scale_y_continuous(name = paste0("PC2 (", pc2_p_var, "%)")) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank())

ggsave(
  here::here("output/figs/metabolomics/presentation/pca_all.pdf"),
  width = 6,
  height = 5,
  units = "in"
)

# Pivot data wider for prcomp
dat_wide_qc <- dat_norm |>
  select(-c(log_transformed, NormIntensity)) |> 
  pivot_wider(values_from = pareto_scaled,
              names_from = Metabolite)

# PCA on pareto scaled data
pca_fit_qc <- dat_wide_qc |> 
  select(where(is.numeric)) |> 
  prcomp(center = TRUE)

# Plot with QC
p_pca_qc <- pca_fit_qc |> 
  broom::augment(dat_wide_qc) |> 
  mutate(Temperature = factor(Temperature, levels = c("25°C", "32°C", "QC"))) |> 
  arrange(Temperature) |> 
  ggplot(aes(x = .fittedPC1,
             y = .fittedPC2,
             fill = Temperature, 
             shape = Region)) +
  stat_ellipse(data = function(x) subset(x, Region == "Temperate"),
               geom = "polygon", alpha = 0.25, color = "grey50") +
  stat_ellipse(data = function(x) subset(x, Region == "Tropical"),
               geom = "polygon", alpha = 0.25, color = "grey50") +
  stat_ellipse(data = function(x) subset(x, Region == "QC"),
               geom = "polygon", alpha = 0.75, color = "grey50") +
  geom_point(aes(alpha = Temperature),
             size = 4, stroke = 1, color = "grey20") +
  scale_fill_manual(values = c("grey80", "firebrick", "grey20")) +
  scale_shape_manual(values = c(15, 21, 24)) +
  scale_alpha_manual(values = c(0.5, 0.5, 0.8)) +
  scale_x_continuous(name = "PC1") +
  scale_y_continuous(name = "PC2") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  theme_minimal(base_size = 14) +
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
  geom_segment(xend = 0, 
               yend = 0, 
               arrow = arrow_style,
               color = "grey20",
               alpha = 0.8) +
  ggrepel::geom_label_repel(aes(label = 
                                  ifelse(abs(PC1) > 0.21 | 
                                           abs(PC2) > 0.21, column, NA)),
                            alpha = 0.95,
                            size = 2,
                            hjust = 1, 
                            nudge_x = -0.02,
                            color = "grey20") +
  scale_x_continuous(limits = c(-lim_dims - 0.2, lim_dims + 0.2)) +
  scale_y_continuous(limits = c(-lim_dims, lim_dims)) +
  theme_minimal(base_size = 14)

# Combine scree plot & rotation plot
pca_extra <- cowplot::plot_grid(
  pca_scree,
  pca_rotation + theme(plot.margin = margin(t = 30,r = 20, l = 20)),
  rel_widths = c(0.7, 1),
  labels = c("B", "C"),
  nrow = 1
)

# Combine all pca plots
pca_all <- cowplot::plot_grid(
  p_pca_qc,
  pca_extra,
  nrow = 2,
  labels = c("A", ""),
  rel_heights = c(2, 1)
)

# cowplot::plot_grid(
#   pca_scree,
#   pca_rotation + theme(plot.margin = margin(t = 30,r = 20, l = 20)),
#   rel_widths = c(0.85, 1),
#   labels = c("A", "B"),
#   nrow = 1
# )


ggsave(plot = pca_all,
       here::here("output/figs/metabolomics/pca_sf.pdf"),
       height = 8,
       width = 7,
       units = "in")
ggsave(plot = pca_all,
       here::here("output/figs/metabolomics/pca_sf.png"),
       height = 8,
       width = 7,
       units = "in")

ggsave(plot = p_pca,
       here::here("output/figs/metabolomics/pca_only.pdf"),
       height = 5,
       width = 6,
       units = "in")
ggsave(plot = p_pca,
       here::here("output/figs/metabolomics/pca_only.svg"),
       height = 5,
       width = 6,
       units = "in")

 # Save PCA plots together
ggsave(plot = pca_all,
       here::here("output/figs/metabolomics/pca.pdf"),
       height = 8,
       width = 7,
       units = "in")

ggsave(plot = pca_all,
       here::here("output/figs/metabolomics/pca.png"),
       height = 8,
       width = 7,
       units = "in")


# Heatmap ----------------------------------------------------------------------

# Create matrix
dat_filt <- dat_wide |> 
  filter(Genotype != "QC")
num_mat <- dat_filt |> 
  mutate(SampleID = paste(Genotype, Temperature, Replicate, sep = "_")) |> 
  column_to_rownames("SampleID") |> 
  select(Alanine:`GSSG/GSH-NEM Ratio`) |> 
  mutate(across(everything(), as.numeric)) |> 
  as.matrix()

# Annotation for the samples
annotation_df <- data.frame(
  Region = factor(dat_filt$Region),
  Temperature = factor(dat_filt$Temperature),
  Genotype = factor(dat_filt$Genotype),
  row.names = rownames(num_mat),
  stringsAsFactors = TRUE
)

# Define colors for plotting
annotation_colors <- list(
  Region = c(Tropical = "#C5A4CC", Temperate = "#9DC2D2"),
  Temperature = c(`25°C` = "grey80", `32°C` = "firebrick"),
  Genotype = c(VT = "green4", FR = "green4", JP = "green4",
               CH = "lightblue", GA = "lightblue", BO = "yellow")
)

# Redraw with gaps and cluster-bar
p_hmap <- pheatmap::pheatmap(
  t(num_mat),
  annotation_col = annotation_df |> select(-Genotype),
  annotation_colors = annotation_colors,
  color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(100),
  scale = "none",
  fontsize = 8,
  fontsize_row = 8,
  clustering_distance_cols = "correlation",
  clustering_distance_rows = "correlation",
  treeheight_col = 20,
  treeheight_row = 20,
  cuttree_cols = 3,
  border_color = NA,
  gaps_col = gaps,
  show_colnames = FALSE) |>
  grid::grid.grabExpr()


ggsave(plot = p_hmap,
       here::here("output/figs/metabolomics/heatmap.svg"),
       height = 6,
       width = 8,
       units = "in"
)


# PC2 vs LT50
dat_lt50 <- readRDS(here::here("data/processed/survival/lt50.rds"))

p_pc2 <- pca_fit |> 
  broom::augment(dat_wide) |> 
  mutate(Region = case_when(Genotype %in% c("BO", "CH", "GA") ~ "Tropical",
                            Genotype %in% c("VT", "FR", "JP") ~ "Temperate")) |> 
  dplyr::select(Region, Locale:Label, .fittedPC1:.fittedPC36) |> 
  left_join(dat_lt50, by = c("Locale" = "genotype")) |> 
  group_by(Genotype, Temperature, Region, estimate) |> 
  #summarise(across(.fittedPC1:.fittedPC36, mean)) |> 
  ggplot(aes(y = .fittedPC2,
             x = estimate,
             fill = Temperature, 
             color = Temperature)) +
  geom_smooth(method = "lm",
              se = FALSE, 
              fullrange = TRUE,
              linetype = 2)+
  geom_point(aes(shape = Region),
             size = 4,
             stroke = 1,
             alpha = 0.8,
             color = "grey20") +
  ggpmisc::stat_poly_eq(aes(label = after_stat(rr.label)),
                        formula = y ~ x,
                        color = c("grey20", "grey99"),
                        geom = "label",
                        label.y = c(2.2, 1),
                        label.x = 36.85,
                        alpha = 0.8,
                        parse = TRUE,
                        method = "lm") +
  scale_fill_manual(values = c("grey80", "firebrick")) +
  scale_color_manual(values = c("grey80", "firebrick")) +
  scale_shape_manual(values = c(21, 24)) +
  scale_x_continuous(name = expression("LT"[50]),
                     limits = c(34, 37),
                     breaks = 34:37,
                     labels = function(x) paste0(x, "°C")) +
  scale_y_continuous(name = "PC2") +
  theme_minimal(base_size = 16) +
  theme(legend.position = "none")


p_pc4 <- pca_fit |> 
  broom::augment(dat_wide) |> 
  mutate(Region = case_when(Genotype %in% c("BO", "CH", "GA") ~ "Tropical",
                            Genotype %in% c("VT", "FR", "JP") ~ "Temperate")) |> 
  dplyr::select(Region, Locale:Label, .fittedPC1:.fittedPC36) |> 
  left_join(dat_lt50, by = c("Locale" = "genotype")) |> 
  group_by(Genotype, Temperature, Region, estimate) |> 
  #summarise(across(.fittedPC1:.fittedPC36, mean)) |> 
  ggplot(aes(y = .fittedPC4,
             x = estimate,
             fill = Temperature, 
             color = Temperature)) +
  geom_smooth(method = "lm",
              se = FALSE, 
              fullrange = TRUE,
              linetype = 2)+
  geom_point(aes(shape = Region),
             size = 4,
             stroke = 1,
             alpha = 0.8,
             color = "grey20") +
  ggpmisc::stat_poly_eq(aes(label = after_stat(rr.label)),
                        formula = y ~ x,
                        color = c("grey20", "grey99"),
                        geom = "label",
                        label.y = c(1, -1.1),
                        label.x = 36.85,
                        alpha = 0.8,
                        parse = TRUE,
                        method = "lm") +
  scale_fill_manual(values = c("grey80", "firebrick")) +
  scale_color_manual(values = c("grey80", "firebrick")) +
  scale_shape_manual(values = c(21, 24)) +
  scale_x_continuous(name = expression("LT"[50]),
                     limits = c(34, 37),
                     breaks = 34:37,
                     labels = function(x) paste0(x, "°C")) +
  scale_y_continuous(name = "PC4") +
  theme_minimal(base_size = 16) +
  theme(legend.position = "none")

cowplot::plot_grid(
  p_pc2,
  p_pc4,
  labels = c("A", "B"),
  nrow = 2
)

ggsave(here::here("output/figs/metabolomics/pc_lt_50.png"),
       height = 8,
       width = 5,
       units = "in"
)

p_pc2_loadings <- pca_fit |> 
  broom::tidy(matrix = "rotation") |> 
  filter(PC == 2) |> 
  slice_max(abs(value), n = 5) |> 
  mutate(column = fct_reorder(column, desc(abs(value)))) |> 
  ggplot() +
  geom_col(aes(x = column,
               y = value,
               fill = value < 0),
           width = 0.75,
           alpha = 0.75,
           color = "grey20") +
  geom_hline(yintercept = 0,
             color = "grey20",
             size = 1) +
  scale_x_discrete(position = "top") +
  scale_y_continuous(name = "PC2 loadings",
                     limits = c(-0.5, 0.5)) +
  #scale_fill_manual(values = c("#8F707A", "#708F85")) +
  scale_fill_manual(values = c("#F18F0E", "#0E70F1") |> 
                      colorspace::desaturate(amount = 0.5)) +
  theme_minimal(base_size = 16) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 8),
        plot.margin = margin(t = 0.5, r = 0, b = 0.5, l = 0.25, unit = "in"))

pca_fit |>
  broom::tidy(matrix = "rotation") |>
  filter(PC == 2) |>
  mutate(contrib = 100 * value^2 / sum(value^2)) |>
  arrange(desc(value)) 

p_pc4_loadings <- pca_fit |> 
  broom::tidy(matrix = "rotation") |> 
  filter(PC == 4) |> 
  slice_max(abs(value), n = 5) |> 
  mutate(column = fct_reorder(column, desc(abs(value)))) |> 
  ggplot() +
  geom_col(aes(x = column,
               y = value, 
               fill = value < 0),
           width = 0.75,
           alpha = 0.75,
           color = "grey20") +
  geom_hline(yintercept = 0,
             color = "grey20",
             size = 1) +
  scale_x_discrete(position = "top") +
  scale_y_continuous(name = "PC4 loadings",
                     limits = c(-0.5, 0.5)) +
  #scale_fill_manual(values = c("#8F707A", "#708F85")) +
  scale_fill_manual(values = c("#0E70F1") |> 
                      colorspace::desaturate(amount = 0.5)) +
  theme_minimal(base_size = 16) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 8),
        plot.margin = margin(t = 0.5, r = 0, b = 0.5, l = 0.25, unit = "in"))

pca_fit |>
  broom::tidy(matrix = "rotation") |>
  filter(PC == 4) |>
  mutate(contrib = 100 * value^2 / sum(value^2)) |>
  arrange(desc(contrib))

cowplot::plot_grid(
  p_pc2, p_pc2_loadings,
  p_pc4, p_pc4_loadings,
  labels = c("A", "B",
             "C", "D"),
  nrow = 2
)

ggsave(here::here("output/figs/metabolomics/pc_lt_50.png"),
       height = 8,
       width = 10,
       units = "in"
)
ggsave(here::here("output/figs/metabolomics/pc_lt_50.pdf"),
       height = 8,
       width = 10,
       units = "in"
)


# Correlation of principal component scores on LT50 for both treatments --------
pca_fit |> 
  broom::augment(dat_wide) |> 
  dplyr::select(Locale:Label, .fittedPC1:.fittedPC36) |>
  left_join(dat_lt50, by = c("Locale" = "genotype")) |>
  pivot_longer(cols = starts_with(".fittedPC"),
               names_to  = "PC",
               values_to = "score") |> 
  group_by(Temperature, PC) |> 
  nest() |> 
  mutate(model = map(data, ~ lm(score ~ estimate, data = .x)),
         tidied = map(model, broom::tidy)) |> 
  unnest(tidied) |> 
  filter(term == "estimate") |> 
  dplyr::select(Temperature, PC, estimate, std.error, statistic, p.value) |> 
  ungroup() |> 
  mutate(padj = p.adjust(p.value, method = "BH")) |> 
  arrange(padj) |> 
  filter(padj < 0.05) |> 
  group_by(PC) |> 
  add_tally() |> 
  filter(n == 2) |> 
  arrange(PC)
