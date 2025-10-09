# ------------------------------------------------------------------------------
# Results of the Differential Metabolite analysis
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

# Counts of DAMS
DAMs |> 
  group_by(term, statistic > 0) |> 
  filter(p.value.adj < 0.05) |> 
  tally()

# Summary plots of the statistics ----------------------------------------------

# Volcano plot of region effect
p_volcano_reg <- results |> 
  filter(term == "RegionTropical") |> 
  ggplot(aes(x = statistic,
             y = -log10(p.value.adj))) +
  geom_hline(yintercept = -log10(0.05),
             linetype = 2) +
  geom_vline(xintercept = 0) +
  geom_point(aes(fill = p.value.adj < 0.05,
                 size = p.value.adj < 0.05),
             shape = 21,
             alpha = 0.7,
             color = "grey20") +
  ggrepel::geom_label_repel(aes(label = 
                                  ifelse(p.value.adj < 0.05, Metabolite, NA)),
                            size = 2,
                            alpha = 0.9,
                            color = color_sig) +
  scale_fill_manual(values = c("grey80", color_sig)) +
  scale_size_manual(values = c(2, 3)) +
  scale_x_continuous(limits = c(-5, 5),
                     name = "Test statistic") +
  labs(title = "Region effect") +
  scale_y_continuous(name = expression("-log"[10]* "(adj. p-value)")) +
  theme_minimal() +
  theme(legend.position = "none")

# Volcano plot of temperature effect
p_volcano_temp <- results |> 
  filter(term == "Temperature32°C") |> 
  ggplot(aes(x = statistic,
             y = -log10(p.value.adj))) +
  geom_hline(yintercept = -log10(0.05),
             linetype = 2) +
  geom_vline(xintercept = 0) +
  geom_point(aes(fill = p.value.adj < 0.05,
                 size = p.value.adj < 0.05),
             shape = 21,
             alpha = 0.7,
             color = "grey20") +
  scale_fill_manual(values = c("grey80", color_sig)) +
  scale_size_manual(values = c(2, 3)) +
  ggrepel::geom_label_repel(aes(label = 
                                  ifelse(p.value.adj < 0.05, Metabolite, NA)),
                            size = 2,
                            alpha = 0.7,
                            force = 10, force_pull = 0.1,
                            fill = "white",
                            color = color_sig) +
  scale_x_continuous(limits = c(-12, 12),
                     name = "Test statistic") +
  scale_y_continuous(name = expression("-log"[10]* "(adj. p-value)")) +
  labs(title = "Temperature effect") +
  theme_minimal() +
  theme(legend.position = "none")

# Volcano plot of interaction
p_volcano_int <- results |> 
  filter(term == "RegionTropical:Temperature32°C") |> 
  ggplot(aes(x = statistic,
             y = -log10(p.value.adj))) +
  geom_hline(yintercept = -log10(0.05),
             linetype = 2) +
  geom_vline(xintercept = 0) +
  geom_point(aes(fill = p.value.adj < 0.05,
                 size = p.value.adj < 0.05),
             shape = 21,
             alpha = 0.7,
             color = "grey20") +
  scale_fill_manual(values = c("grey80", color_sig)) +
  scale_size_manual(values = c(2, 3)) +
  ggrepel::geom_label_repel(aes(label = 
                                  ifelse(p.value.adj < 0.05, Metabolite, NA)),
                            size = 2,
                            alpha = 0.7,
                            force = 40, force_pull = 0.075,
                            fill = "white",
                            max.time = 10,
                            max.iter = 100000,
                            color = color_sig) +
  scale_x_continuous(limits = c(-6.25, 6.25),
                     name = "Test statistic") +
  scale_y_continuous(name = expression("-log"[10]* "(adj. p-value)")) +
  theme_minimal() +
  labs(title = "Interaction (Region x Temperature)") +
  theme(legend.position = "none")



p_top <- cowplot::plot_grid(
  p_volcano_reg,
  p_volcano_temp,
  labels = c("A", "B"),
  nrow = 1
)

p_bottom <- cowplot::plot_grid(
  NULL, p_volcano_int, NULL,
  rel_widths = c(1, 2, 1),
  nrow = 1,
  labels = c("", "C", ""))

volcanos <- cowplot::plot_grid(
  p_top,
  p_bottom,
  nrow = 2
)


# Save
ggsave(plot = volcanos,
       here::here(fig_dir, "volcano_plots.pdf"),
       height = 7,
       width = 7,
       units = "in")
ggsave(plot = volcanos,
       here::here(fig_dir, "volcano_plots.png"),
       height = 7,
       width = 7,
       units = "in")
 
