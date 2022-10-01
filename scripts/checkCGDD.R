library(dplyr)
library(Rcpp)
library(here)


# Functions ----
source(here::here("scripts/functions.R"))
sourceCpp(here::here("scripts/countOrtmann_modSpeed.cpp"))


# matrix ----
(ni <- matrix(sample(1:10, 20, replace = T), ncol = 2))

for(spalte in 1: ncol(ni)){
  t_table <- table(ni[,spalte])/nrow(ni)
  x <- attr(t_table, "dimnames")[[1]] %>%
    as.numeric()
  x_px <- 0
  x2_px <- 0
  for(i in 1:length(t_table)){
    x_px <- x_px + x[i] * as.numeric(t_table[i])
    x2_px <- x2_px + x[i]*x[i] * as.numeric(t_table[i])
  }
  t_sd <- sqrt(x2_px-x_px)
  res <- matrix(nrow = 2, ncol = length(t_table))
  res[,1] <- c(x[1]/t_sd,as.numeric(t_table[1]))
  for(i in 2:length(t_table)){
    res[1,i] <- x[i]/t_sd
    res[2,i] <- res[2,i-1] + as.numeric(t_table[i])
  }
  print(res)
  print(t_sd)
}

CGDD_wo2(ni)
