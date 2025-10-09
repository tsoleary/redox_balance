# ------------------------------------------------------------------------------
# Colors for plotting and themes
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)

# Colors
colors_genotypes <- c("#9DC2D2", "#73A8BF", "#4B88A2",
                      "#C5A4CC", "#AC7CB6", "#92589D")

label_scientific_TSO <- function(x) {
  labs <- scales::scientific_format()(x)
  labs2 <- ifelse(
    x == 0,
    "0",
    gsub("([[:digit:]](?:[[:digit:]]|\\.)*)e\\+?(-?[[:digit:]]+)",
         "\\1*plain(x)*10^\\2",
         labs))
  parse(text = labs2)
}
