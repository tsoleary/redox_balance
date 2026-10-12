# ------------------------------------------------------------------------------
# Initial data processing for downstream analysis
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)

# Rename the glutathione couple and flip the ratio to GSH/GSSG -----------------
rename_met <- function(x) {
  case_when(x == "GSSG/GSH-NEM Ratio" ~ "GSH/GSSG Ratio",
            x == "Glutathione-NEM"    ~ "GSH",
            TRUE                      ~ x)
}

# Load data
dat_raw <- read_csv(
  here::here("data/raw/metabolomics/redox_metabolomics_data_norm.csv"))  |> 
  separate(Sample, into = c("Genotype", "Temperature", "Replicate")) |> 
  mutate(Replicate = ifelse(is.na(Replicate), Temperature, Replicate)) |> 
  mutate(Temperature = ifelse(Temperature %in% c("25", "32"), 
                              paste0(Temperature, "°C"), "QC")) |> 
  mutate(Replicate = paste0("Rep_", Replicate)) |> 
  mutate(Region = case_when(Genotype %in% c("BO", "CH", "GA") ~ "Tropical",
                            Genotype %in% c("VT", "FR", "JP") ~ "Temperate",
                            Genotype == "QC" ~ "QC"))

# Tidy & nest data for analysis and plotting ------------
dat <- dat_raw |> 
  mutate(Locale = case_when(Genotype == "VT" ~ "Vermont, USA",
                            Genotype == "FR" ~ "Montpellier, France",
                            Genotype == "JP" ~ "Shiojiri, Japan",
                            Genotype == "CH" ~ "Chiapas, Mexico",
                            Genotype == "GA" ~ "Accra, Ghana",
                            Genotype == "BO" ~ "Mumbai, India")) |> 
  mutate(Locale = factor(Locale, levels = c("Vermont, USA", 
                                            "Montpellier, France", 
                                            "Shiojiri, Japan", 
                                            "Chiapas, Mexico", 
                                            "Accra, Ghana", 
                                            "Mumbai, India"))) |>
  select(Region, Locale, everything()) |> 
  pivot_longer(cols = Alanine:`GSSG/GSH-NEM Ratio`,
               names_to = "Metabolite",
               values_to = "NormIntensity") |>
  # Flip, then rename
  mutate(NormIntensity = if_else(Metabolite == "GSSG/GSH-NEM Ratio",
                                 1 / NormIntensity,
                                 NormIntensity),
         Metabolite = rename_met(Metabolite)) |>
  group_by(Metabolite) |>
  nest()

# Calculate log-transformed and pareto scaled data for pca ---------
dat_norm <- dat |> 
  unnest(cols = data) |> 
  group_by(Metabolite) |> 
  mutate(log_transformed = log10(NormIntensity),
         pareto_scaled = 
           (log_transformed - mean(log_transformed)) / 
           sqrt(sd(log_transformed))) 

# Save data -- dat_raw keeps the original metabolite names
saveRDS(dat_raw, here::here("data/processed/metabolomics/dat_raw.rds"))
saveRDS(dat, here::here("data/processed/metabolomics/dat.rds"))
saveRDS(dat_norm, here::here("data/processed/metabolomics/dat_norm.rds"))

# Fig. S1 -- per-sample distributions before and after normalization ----------
# Ratios and process blanks excluded.

dat_core_raw <- read_csv(
  here::here("data/raw/metabolomics/redox_metabolomics_data_raw.csv")) |>
  filter(!str_detect(Sample, "^Pb")) |>
  pivot_longer(-c(Sample, Label), names_to = "Metabolite",
               values_to = "RawIntensity") |>
  # Match dat_norm's names and orientation for the join
  mutate(RawIntensity = if_else(Metabolite == "GSSG/GSH-NEM Ratio",
                                1 / RawIntensity,
                                RawIntensity),
         Metabolite = rename_met(Metabolite))

norm_factors <- read_csv(
  here::here("data/raw/metabolomics/redox_metabolomics_norm_factors.csv")) |>
  dplyr::select(Sample, Factor)

stage_levels <- c("log10 raw intensity", "log10 normalized intensity")

dat_stages <- dat_norm |>
  ungroup() |>
  filter(!str_detect(Metabolite, "Ratio"),
         Region != "QC") |>
  # rebuild the core's sample id: BO_25_1, QC_1
  mutate(rep_n  = str_remove(Replicate, "^Rep_"),
         Sample = if_else(Region == "QC",
                          paste0("QC_", rep_n),
                          paste(Genotype, str_remove(Temperature, "°C"), rep_n,
                                sep = "_")),
         # readable version for the axis: VT 25°C Rep 1
         SampleLabel = if_else(Region == "QC",
                               paste("QC Rep", rep_n),
                               paste(Genotype, Temperature, "Rep", rep_n))) |>
  dplyr::select(Sample, SampleLabel, Region, Temperature, Metabolite,
                `log10 normalized intensity` = log_transformed) |>
  left_join(dat_core_raw |>
              mutate(`log10 raw intensity` = log10(RawIntensity)) |>
              dplyr::select(Sample, Metabolite, `log10 raw intensity`),
            by = c("Sample", "Metabolite")) |>
  pivot_longer(all_of(stage_levels), names_to = "stage", values_to = "value") |>
  mutate(stage = factor(stage, levels = stage_levels),
         group = if_else(Region == "QC", "QC", Temperature))

# Order samples by group
samp_order <- dat_stages |>
  distinct(SampleLabel, group) |>
  arrange(group, SampleLabel) |>
  pull(SampleLabel)

dat_stages <- dat_stages |>
  mutate(SampleLabel = factor(SampleLabel, levels = samp_order))

fill_group <- c("25°C" = "grey80", "32°C" = "firebrick")

# Zoom the x axis to the bulk of the data
x_rng <- quantile(dat_stages$value, c(0.01, 0.99), na.rm = TRUE)

x_labs <- list(
  "log10 raw intensity"        = expression("log"[10]*"(Raw intensity)"),
  "log10 normalized intensity" = expression("log"[10]*"(Normalized intensity)")
)

box_panel <- function(dat, which_stage, show_labels) {
  p <- dat |>
    filter(stage == which_stage) |>
    ggplot(aes(x = value, y = SampleLabel, fill = group)) +
    geom_boxplot(outlier.shape = 21, outlier.size = 0.5,
                 outlier.alpha = 0.4, linewidth = 0.2) +
    scale_fill_manual(values = fill_group, name = NULL) +
    scale_x_continuous(name = x_labs[[which_stage]]) +
    scale_y_discrete(name = NULL, limits = rev(samp_order)) +
    coord_cartesian(xlim = x_rng) +
    theme_minimal(base_size = 9) +
    theme(panel.grid.major.y = element_line(linewidth = 0.15),
          panel.grid.minor   = element_blank(),
          legend.position    = "none")

  if (show_labels) {
    p + theme(axis.text.y = element_text(size = 5))
  } else {
    p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  }
}

# One density per sample, 25°C drawn over 32°C
dens_panel <- function(dat, which_stage) {
  d <- filter(dat, stage == which_stage)
  ggplot(mapping = aes(x = value, group = SampleLabel, colour = group)) +
    geom_line(data = filter(d, group == "32°C"), stat = "density",
              linewidth = 0.25, alpha = 0.7) +
    geom_line(data = filter(d, group == "25°C"), stat = "density",
              linewidth = 0.25, alpha = 0.7) +
    scale_colour_manual(values = fill_group, guide = "none") +
    coord_cartesian(xlim = x_rng) +
    theme_void()
}

stage_col <- function(st, show_labels) {
  cowplot::plot_grid(dens_panel(dat_stages, st),
                     box_panel(dat_stages, st, show_labels),
                     ncol = 1, rel_heights = c(1, 9),
                     align = "v", axis = "lr")
}

fig_legend <- cowplot::get_legend(
  box_panel(dat_stages, stage_levels[1], FALSE) +
    theme(legend.position = "top")
)

p_norm_all <- cowplot::plot_grid(
  fig_legend,
  cowplot::plot_grid(stage_col(stage_levels[1], TRUE),
                     stage_col(stage_levels[2], FALSE),
                     labels = c("A", "B"), ncol = 2,
                     rel_widths = c(1.25, 1)),
  ncol = 1, rel_heights = c(0.04, 1)
)

ggsave(plot = p_norm_all,
       here::here("output/figs/final/fig_normalization.pdf"),
       height = 9,
       width = 8,
       units = "in")
ggsave(plot = p_norm_all,
       here::here("output/figs/final/fig_normalization.png"),
       height = 9,
       width = 8,
       units = "in")
ggsave(plot = p_norm_all,
       here::here("output/figs/final/fig_normalization.svg"),
       height = 9,
       width = 8,
       units = "in")

# Numbers for the legend
dat_stages |>
  group_by(stage, group, Sample) |>
  summarise(med = median(value), .groups = "drop_last") |>
  summarise(n = n(), sd_of_medians = sd(med), .groups = "drop") |>
  print()

norm_factors |>
  summarise(min = min(Factor), max = max(Factor),
            median = median(Factor), sd = sd(Factor),
            fold_range = max(Factor) / min(Factor)) |>
  print()

# Fig. S2 -- per-metabolite distributions, log10 and Pareto scaled ------------

met_order <- dat_norm |>
  ungroup() |>
  distinct(Metabolite) |>
  arrange(desc(Metabolite)) |>
  pull(Metabolite)

dat_mets <- dat_norm |>
  ungroup() |>
  filter(Region != "QC") |>
  mutate(Metabolite = factor(Metabolite, levels = met_order))

met_panel <- function(dat, var, xlab) {
  box <- dat |>
    ggplot(aes(x = {{ var }}, y = Metabolite)) +
    geom_boxplot(outlier.shape = 21, outlier.size = 0.6,
                 outlier.alpha = 0.4, linewidth = 0.25, fill = "grey85") +
    scale_x_continuous(name = xlab) +
    scale_y_discrete(name = NULL) +
    theme_minimal(base_size = 9) +
    theme(panel.grid.major.y = element_line(linewidth = 0.2))

  dens <- dat |>
    ggplot(aes(x = {{ var }})) +
    geom_density(fill = "grey85", linewidth = 0.3) +
    scale_x_continuous(name = NULL) +
    scale_y_continuous(name = NULL) +
    theme_void()

  cowplot::plot_grid(dens, box, ncol = 1, rel_heights = c(1, 6),
                     align = "v", axis = "lr")
}

p_mets <- cowplot::plot_grid(
  met_panel(dat_mets, log_transformed,
            expression("log"[10]*"(Normalized intensity)")),
  met_panel(dat_mets, pareto_scaled, "Pareto scaled"),
  labels = c("A", "B"),
  ncol = 2
)

ggsave(plot = p_mets,
       here::here("output/figs/final/fig_metabolite_distributions.pdf"),
       height = 7,
       width = 9,
       units = "in")
ggsave(plot = p_mets,
       here::here("output/figs/final/fig_metabolite_distributions.png"),
       height = 7,
       width = 9,
       units = "in")
ggsave(plot = p_mets,
       here::here("output/figs/final/fig_metabolite_distributions.svg"),
       height = 7,
       width = 9,
       units = "in")
