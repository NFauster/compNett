library(Rcpp)
library(orca)
library(dplyr)

source(here::here("scripts/functions.R"))
sourceCpp(here::here("scripts/countOrtmann_modSpeed.cpp"))


count_and_GCM_complete <- function(edge_list){
  count4(matrix(as.integer(edge_list),ncol=2,byrow=F)) %>%
    GCM_wo2()
}


crt_name <- "BBexp-3e+05-5e-04_2.RData"

load(here::here(paste0("data/GDV_",
                       crt_name)))


crt_GCM <- GDV %>%
  GCM_wo2()

save(crt_GCM,
     file = here::here(paste0("data/GCM_complete_",
                              crt_name)))
rm(crt_GCM)
gc()



###########
crt_GCM <-  diag(nrow = 15, ncol = 15, names = FALSE)
GDV <- rbind(GDV, rep(1,15))

for(i in 1:14){
  for(j in (i+1):15){
    crt_GCM[i,j] <- spearman(GDV[,i], GDV[,j])
    crt_GCM[j,i] <- crt_GCM[i,j]
  }
}


