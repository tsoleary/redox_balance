# ------------------------------------------------------------------------------
# Natural Lines Survival Curves
# April 19, 2023
# TS O'Leary
# ------------------------------------------------------------------------------

# Description -----
# Count the number of observations of each species.

# Load libraries
library(tidyverse)
library(broom)
library(drc)

# Load data
dat <- read_csv("data/raw/fr_jp_hatch_2023-04-27.csv")

# Analyze data
df <- dat %>%
  group_by(Genotype) %>%
  filter(Genotype %in% c("FRMO-01", "FRMO-02", "JPSC")) %>%
  nest() %>%
  mutate(fit = map(data, ~ drm(N_hatched/N_eggs ~ Temperature, 
                               data = .x, 
                               weight = N_eggs,
                               fct = LL.3(names = c("slope", "upper limit", "LT50")),
                               type = "binomial")),
         lt_s = map(fit, ~ ED(.x, 
                              c(10, 50, 90), 
                              interval = "delta",
                              display = FALSE)))

df$lt_s

  
# Plot data
dat %>%
  #filter(Genotype %in% c("FRMO-01", "FRMO-02", "JPSC")) %>%
  ggplot(aes(x = Temperature, 
             y = N_hatched/N_eggs)) +
  geom_jitter(aes(fill = Genotype),
              color = "grey80",
              size = 2,
              alpha = 0.6,
              width = 0.1,
              height = 0,
              shape = 21) +
  geom_line(aes(color = Genotype),
            stat = "smooth",
            method = drm,
            alpha = 0.5,
            size = 1,
            method.args = list(fct = LL.3()),
            #linetype = 3,
            se = FALSE) +
  labs(x = "Temperature (°C)", 
       y = "Haching success") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                     labels = scales::percent_format(accuracy = 1)) +
  scale_x_continuous(breaks = c(25, 28, 30, 32, 34, 36, 38)) +
  cowplot::theme_minimal_grid(font_family = "Myriad Pro") +
  facet_wrap(~Genotype)

ggsave("~/Downloads/fr_jp.pdf", width = 10, height = 4)  
