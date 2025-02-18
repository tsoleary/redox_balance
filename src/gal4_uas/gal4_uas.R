# ------------------------------------------------------------------------------
# Gal4/UAS overexpression phenotype results
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)

# Load data
dat <- read_csv("data/raw/gal4_uas.csv") |> 
  mutate(gene = factor(gene, levels = c("Atg6", "Sod3", "Ucp4A", "Ctrl"))) |> 
  mutate(ctrl = ifelse(str_detect(gene, "Ctrl"), "ctrl", "overexpression")) |> 
  mutate(gal4_sex = ifelse(gal4_sex == "female", 
                           "♀ Hsp70-Gal4 x  ♂ UAS-Gene",
                           "♀ UAS-Gene x  ♂ Hsp70-Gal4"))


# Count total number of eggs
dat |> 
  summarize(total_eggs = sum(eggs))

# Number of eggs per gene
dat |> 
  group_by(genotype_id, temperature) |> 
  summarise(total_eggs = sum(eggs))

# Number of eggs per gene and direction of cross
dat |> 
  group_by(gene, temperature, gal4_sex) |> 
  summarise(total_eggs = sum(eggs)) |> 
  arrange(total_eggs)


# Analyze data
dat_avg <- dat |> 
  group_by(gene, temperature, ctrl, gal4_sex) |> 
  summarise(eggs = sum(eggs),
            hatched = sum(hatched)) |> 
  mutate(survival = hatched/eggs,
         se = sqrt(survival * (1 - survival) / eggs))

dat_avg_all <- dat_avg |> 
  group_by(gene, temperature, ctrl) |> 
  summarise(eggs = sum(eggs),
            hatched = sum(hatched)) |> 
  mutate(survival = hatched/eggs,
         se = sqrt(survival * (1 - survival) / eggs))

# Plot data
gene_colors <- c("orchid", "forestgreen", "lightblue", "grey50")

dat_avg |> 
  ggplot(aes(x = temperature,
             y = survival,
             fill = gene,
             color = gene)) +
  geom_line(aes(linetype = ctrl),
            size = 1) +
  geom_errorbar(aes(ymin = survival - se,
                    ymax = survival + se),
                width = 0.15) +
  geom_point(shape = 21,
             size = 4,
             alpha = 0.8,
             color = "grey80") +
  ggrepel::geom_text_repel(data = dat_avg |> 
              filter(temperature == 35),
            aes(label = gene),
            fontface = "bold",
            min.segment.length = 2,
            size = 5,
            hjust = 0,
            nudge_x = 0.5) +
  scale_color_manual(values = gene_colors) +
  scale_fill_manual(values = gene_colors) +
  scale_linetype_manual(values = c(2, 1)) +
  scale_y_continuous(limits = c(0, 1),
                     labels = scales::percent,
                     expand = c(0, 0),
                     name = "Hatching success") +
  scale_x_continuous(breaks = c(25, 35),
                     labels = c("25°C", "35°C"),
                     expand = c(0, 2),
                     name = "Temperature") +
  cowplot::theme_minimal_hgrid() +
  theme(legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_rect(fill = "grey95")) +
  facet_wrap(~gal4_sex, nrow = 1)


dat_avg_all |> 
  ggplot(aes(x = temperature,
             y = survival,
             fill = gene,
             color = gene)) +
  geom_line(aes(linetype = ctrl),
            size = 1) +
  geom_errorbar(aes(ymin = survival - se,
                    ymax = survival + se),
                width = 0.15) +
  geom_point(shape = 21,
             size = 4,
             alpha = 0.8,
             color = "grey80") +
  ggrepel::geom_text_repel(data = dat_avg_all |> 
                             filter(temperature == 35),
                           aes(label = gene),
                           fontface = "bold",
                           min.segment.length = 2,
                           direction = "y",
                           size = 5,
                           hjust = 0,
                           nudge_x = 0.5) +
  scale_color_manual(values = gene_colors) +
  scale_fill_manual(values = gene_colors) +
  scale_linetype_manual(values = c(2, 1)) +
  scale_y_continuous(limits = c(0, 1),
                     labels = scales::percent,
                     expand = c(0, 0),
                     name = "Hatching success") +
  scale_x_continuous(breaks = c(25, 35),
                     labels = c("25°C", "35°C"),
                     expand = c(0, 2),
                     name = "Temperature") +
  cowplot::theme_minimal_hgrid() +
  theme(legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_rect(fill = "grey95"))





