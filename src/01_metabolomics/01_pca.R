# ------------------------------------------------------------------------------
# Principle component analysis
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)
library(lmerTest)

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

# Per-PC and cumulative variance
broom::tidy(pca_fit, matrix = "eigenvalues") |>
  write_csv(here::here("output/metabolomics/pca_eigenvalues.csv"))

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


ggsave(plot = pca_all,
       here::here("output/figs/final/fig_pca_qc_scree_loadings.pdf"),
       height = 8,
       width = 7,
       units = "in",
       bg = "white")
ggsave(plot = pca_all,
       here::here("output/figs/final/fig_pca_qc_scree_loadings.png"),
       height = 8,
       width = 7,
       units = "in",
       bg = "white")
ggsave(plot = pca_all,
       here::here("output/figs/final/fig_pca_qc_scree_loadings.svg"),
       height = 8,
       width = 7,
       units = "in",
       bg = "white")

# p_pca is panel A of fig_pca_drivers

 # pca_all saved above as fig_pca_qc_scree_loadings


# Heatmap ----------------------------------------------------------------------

# Create matrix
dat_filt <- dat_wide |> 
  filter(Genotype != "QC")
num_mat <- dat_filt |> 
  mutate(SampleID = paste(Genotype, Temperature, Replicate, sep = "_")) |> 
  column_to_rownames("SampleID") |> 
  # Metabolite columns are the numeric ones
  select(where(is.numeric)) |>
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
       here::here("output/figs/final/fig_metabolite_heatmap.pdf"),
       height = 6,
       width = 8,
       units = "in")
ggsave(plot = p_hmap,
       here::here("output/figs/final/fig_metabolite_heatmap.png"),
       height = 6,
       width = 8,
       units = "in")
ggsave(plot = p_hmap,
       here::here("output/figs/final/fig_metabolite_heatmap.svg"),
       height = 6,
       width = 8,
       units = "in")


# PC2 vs LT50
dat_lt50 <- readRDS(here::here("data/processed/survival/lt50.rds"))

p_pc2 <- pca_fit |> 
  broom::augment(dat_wide) |> 
  mutate(Region = case_when(Genotype %in% c("BO", "CH", "GA") ~ "Tropical",
                            Genotype %in% c("VT", "FR", "JP") ~ "Temperate")) |> 
  dplyr::select(Region, Locale:Label, .fittedPC1:.fittedPC36) |> 
  left_join(dat_lt50, by = c("Locale" = "genotype")) |> 
  group_by(Genotype, Temperature, Region, estimate) |> 
  summarise(across(.fittedPC1:.fittedPC36, mean), .groups = "drop") |>
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
  summarise(across(.fittedPC1:.fittedPC36, mean), .groups = "drop") |>
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


# PC scores averaged within genotype x temperature ----------------------------
n_pcs <- 6   # ~92% of the variance

pc_line_means <- pca_fit |>
  broom::augment(dat_wide) |>
  dplyr::select(Region, Genotype, Locale, Temperature, .fittedPC1:.fittedPC36) |>
  pivot_longer(cols = starts_with(".fittedPC"),
               names_to  = "PC",
               values_to = "score") |>
  mutate(PC = as.integer(str_remove(PC, "\\.fittedPC"))) |>
  filter(PC <= n_pcs) |>
  group_by(Region, Genotype, Locale, Temperature, PC) |>
  summarise(score = mean(score), .groups = "drop") |>
  left_join(dat_lt50 |> dplyr::select(genotype, LT50 = estimate),
            by = c("Locale" = "genotype"))

# Correlation of principal component scores on LT50 for both treatments --------
pc_lt50 <- pc_line_means |>
  nest(data = -c(Temperature, PC)) |>
  mutate(fit    = map(data, ~ lm(score ~ LT50, data = .x)),
         tidied = map(fit, broom::tidy),
         r2     = map_dbl(fit, ~ summary(.x)$r.squared)) |>
  unnest(tidied) |>
  filter(term == "LT50") |>
  dplyr::select(Temperature, PC, estimate, std.error, statistic, p.value, r2) |>
  mutate(padj = p.adjust(p.value, method = "BH")) |>
  arrange(PC, Temperature)

pc_lt50

# PCs associated with LT50 at both temperatures
pc_lt50 |>
  filter(padj < 0.05) |>
  add_count(PC) |>
  filter(n == 2) |>
  arrange(PC)

# Test for region and temperature differences in PC space ---------------------
# Same model as the DAM tests, applied to the PC scores.
pc_tests <- pca_fit |>
  broom::augment(dat_wide) |>
  dplyr::select(Region, Genotype, Temperature, .fittedPC1:.fittedPC36) |>
  pivot_longer(cols = starts_with(".fittedPC"),
               names_to  = "PC",
               values_to = "score") |>
  mutate(PC = as.integer(str_remove(PC, "\\.fittedPC"))) |>
  filter(PC <= n_pcs) |>
  nest(data = -PC) |>
  mutate(model  = map(data, ~ lmer(score ~ Region * Temperature +
                                     (1 + Temperature | Genotype), data = .x)),
         tidied = map(model, broom.mixed::tidy)) |>
  unnest(tidied) |>
  filter(effect == "fixed", term != "(Intercept)") |>
  group_by(term) |>
  mutate(p.value.adj = p.adjust(p.value, method = "fdr")) |>
  ungroup() |>
  dplyr::select(PC, term, estimate, std.error, statistic, df,
                p.value, p.value.adj) |>
  arrange(term, PC)

pc_tests

# Which metabolites drive the group differences on each PC? -------------------
# contribution = loading x group mean difference; these sum to the observed gap

pc_drivers <- function(pc, by) {
  lv <- if (by == "Region") c("Temperate", "Tropical") else c("25°C", "32°C")

  loadings <- pca_fit |>
    broom::tidy(matrix = "rotation") |>
    filter(PC == pc) |>
    dplyr::select(Metabolite = column, loading = value)

  dat_wide |>
    ungroup() |>
    dplyr::select(all_of(by), where(is.numeric)) |>
    pivot_longer(-all_of(by), names_to = "Metabolite", values_to = "value") |>
    group_by(Metabolite, .data[[by]]) |>
    summarise(mean = mean(value), .groups = "drop") |>
    pivot_wider(names_from = all_of(by), values_from = mean) |>
    mutate(diff = .data[[lv[2]]] - .data[[lv[1]]]) |>
    left_join(loadings, by = "Metabolite") |>
    mutate(contribution = loading * diff,
           pct_of_gap   = 100 * contribution / sum(contribution),
           higher_in    = factor(if_else(diff > 0, lv[2], lv[1]), levels = lv)) |>
    arrange(desc(abs(contribution)))
}

pal_region <- c(Temperate = "#4B88A2", Tropical = "#92589D")
pal_temp   <- c(`25°C` = "grey70", `32°C` = "firebrick")

n_show <- 10   # top N by absolute value

p_drivers_panel <- function(d, xlab, pal) {
  d <- d |>
    slice_max(abs(contribution), n = n_show) |>
    mutate(Metabolite = fct_reorder(Metabolite, pct_of_gap))

  # Pad the left only when a bar runs that way, to hold its label
  pad_left <- if (any(d$pct_of_gap < 0)) 0.12 else 0.02

  d |>
    ggplot(aes(x = pct_of_gap, y = Metabolite, fill = higher_in)) +
    geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.3) +
    geom_col(width = 0.7, alpha = 0.8) +
    geom_text(aes(label = sprintf("%.1f", pct_of_gap),
                  hjust = if_else(pct_of_gap >= 0, -0.25, 1.25)),
              size = 3.2, colour = "grey30") +
    scale_fill_manual(values = pal, name = NULL, drop = FALSE) +
    # Fill both keys from the palette, so a level missing from the top n still
    # shows its colour instead of a blank swatch
    guides(fill = guide_legend(
      override.aes = list(fill = unname(pal[levels(d$higher_in)]),
                          alpha = 0.8))) +
    scale_x_continuous(name = xlab,
                       expand = expansion(mult = c(pad_left, 0.12))) +
    theme_minimal(base_size = 11) +
    theme(axis.title.y       = element_blank(),
          axis.text.y        = element_text(size = 11),
          panel.grid.major.y = element_blank(),
          panel.grid.minor   = element_blank(),
          legend.position    = "top",
          legend.key.size    = unit(9, "pt"))
}

pc1_temp   <- pc_drivers(1, "Temperature")
pc2_region <- pc_drivers(2, "Region")
pc2_temp   <- pc_drivers(2, "Temperature")
pc3_temp   <- pc_drivers(3, "Temperature")

write_csv(pc2_region, here::here("output/metabolomics/pc2_region_drivers.csv"))
write_csv(pc1_temp,   here::here("output/metabolomics/pc1_temp_drivers.csv"))
write_csv(pc2_temp,   here::here("output/metabolomics/pc2_temp_drivers.csv"))
write_csv(pc3_temp,   here::here("output/metabolomics/pc3_temp_drivers.csv"))

# Scores on any pair of PCs, for the supplement --------------------------------
pc_var <- function(pc) {
  sprintf("%.1f", 100 * broom::tidy(pca_fit, matrix = "eigenvalues")$percent[pc])
}

p_scores_panel <- function(pc_x, pc_y) {
  pca_fit |>
    broom::augment(dat_wide) |>
    mutate(.pcx = .data[[paste0(".fittedPC", pc_x)]],
           .pcy = .data[[paste0(".fittedPC", pc_y)]]) |>
    ggplot(aes(x = .pcx, y = .pcy, fill = Temperature, shape = Region)) +
    stat_ellipse(data = function(x) subset(x, Region == "Temperate"),
                 geom = "polygon", alpha = 0.5, color = "grey50") +
    stat_ellipse(data = function(x) subset(x, Region == "Tropical"),
                 geom = "polygon", alpha = 0.5, color = "grey50") +
    geom_point(size = 3, stroke = 1, alpha = 0.8, color = "grey20") +
    scale_fill_manual(values = c("grey80", "firebrick")) +
    scale_shape_manual(values = c(21, 24)) +
    scale_x_continuous(name = paste0("PC", pc_x, " (", pc_var(pc_x), "%)")) +
    scale_y_continuous(name = paste0("PC", pc_y, " (", pc_var(pc_y), "%)")) +
    guides(fill = guide_legend(override.aes = list(shape = 21))) +
    theme_minimal(base_size = 11) +
    theme(legend.title     = element_blank(),
          legend.position  = "top",
          legend.key.size  = unit(9, "pt"),
          panel.grid.minor = element_blank())
}

# Spread of scores by region on each PC ----------------------------------------
pc_spread <- pca_fit |>
  broom::augment(dat_wide) |>
  dplyr::select(Region, Temperature, .fittedPC1:.fittedPC3) |>
  pivot_longer(starts_with(".fittedPC"), names_to = "PC", values_to = "score") |>
  mutate(PC = as.integer(str_remove(PC, "\\.fittedPC"))) |>
  group_by(PC, Region) |>
  summarise(n = n(), mean = mean(score), sd = sd(score), .groups = "drop") |>
  pivot_wider(names_from = Region, values_from = c(n, mean, sd)) |>
  mutate(sd_ratio = sd_Temperate / sd_Tropical)

pc_spread

write_csv(pc_spread, here::here("output/metabolomics/pc_spread_by_region.csv"))

# Top loadings on a PC, for the supplement -------------------------------------
p_loadings_panel <- function(pc) {
  pca_fit |>
    broom::tidy(matrix = "rotation") |>
    filter(PC == pc) |>
    dplyr::select(Metabolite = column, value) |>
    slice_max(abs(value), n = n_show) |>
    mutate(Metabolite = fct_reorder(Metabolite, value)) |>
    ggplot(aes(x = value, y = Metabolite, fill = value < 0)) +
    geom_col(width = 0.7, alpha = 0.8) +
    geom_vline(xintercept = 0, colour = "grey30", linewidth = 0.4) +
    scale_fill_manual(values = c(`FALSE` = "#F18F0E", `TRUE` = "#0E70F1"),
                      guide = "none") +
    scale_x_continuous(name = paste0("PC", pc, " loading")) +
    theme_minimal(base_size = 11) +
    theme(axis.title.y       = element_blank(),
          axis.text.y        = element_text(size = 11),
          panel.grid.major.y = element_blank(),
          panel.grid.minor   = element_blank())
}

# Main figure: PCA scores and the drivers of each axis ------------------------
p_fig_pca <- cowplot::plot_grid(
  p_pca,
  p_drivers_panel(pc1_temp,
                  "Contribution to temperature difference in PC1 (%)", pal_temp),
  p_drivers_panel(pc2_region,
                  "Contribution to region difference in PC2 (%)", pal_region),
  p_drivers_panel(pc2_temp,
                  "Contribution to temperature difference in PC2 (%)", pal_temp),
  labels = c("A", "B", "C", "D"),
  ncol = 2
)

ggsave(plot = p_fig_pca,
       here::here("output/figs/final/fig_pca_drivers.pdf"),
       height = 8,
       width = 10,
       units = "in",
       bg = "white")
ggsave(plot = p_fig_pca,
       here::here("output/figs/final/fig_pca_drivers.png"),
       height = 8,
       width = 10,
       units = "in",
       bg = "white")
ggsave(plot = p_fig_pca,
       here::here("output/figs/final/fig_pca_drivers.svg"),
       height = 8,
       width = 10,
       units = "in",
       bg = "white")

# Supplement: PC3 scores, loadings for PC1-3, PC3 temperature drivers ---------
p_fig_loadings <- cowplot::plot_grid(
  p_scores_panel(2, 3),
  p_loadings_panel(1),
  p_loadings_panel(2),
  p_loadings_panel(3),
  p_drivers_panel(pc3_temp,
                  "Contribution to temperature difference in PC3 (%)", pal_temp),
  labels = c("A", "B", "C", "D", "E"),
  ncol = 2
)

ggsave(plot = p_fig_loadings,
       here::here("output/figs/final/fig_pc_loadings.pdf"),
       height = 12,
       width = 10,
       units = "in",
       bg = "white")
ggsave(plot = p_fig_loadings,
       here::here("output/figs/final/fig_pc_loadings.png"),
       height = 12,
       width = 10,
       units = "in",
       bg = "white")
ggsave(plot = p_fig_loadings,
       here::here("output/figs/final/fig_pc_loadings.svg"),
       height = 12,
       width = 10,
       units = "in",
       bg = "white")


saveRDS(pc_lt50,  here::here("output/metabolomics/pc_lt50_line_means.rds"))
saveRDS(pc_tests, here::here("output/metabolomics/pc_tests.rds"))
write_csv(pc_lt50,  here::here("output/metabolomics/pc_lt50_line_means.csv"))
write_csv(pc_tests, here::here("output/metabolomics/pc_tests.csv"))
