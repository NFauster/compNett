### Generating networks: Bianconi-Barabasi


# Packages ----
library(Rcpp)
library(magrittr)

sourceCpp(here::here("scripts/sample_BB.cpp"))

sample_BB(1e5, 2, 3) %>%
  matrix(ncol = 2, byrow = TRUE)
