library(orca)
library(dplyr)

filenames <- c("ER-3e+05-0.005_1.RData",
               # "ER-3e+05-0.001_2.RData",
               
               "BA-3e+05-0.005_1.RData",
               
               "BBunif-3e+05-0.005_1.RData",
               # "BBunif-3e+05-0.001_2.RData",
               
               "BBexp-3e+05-0.005_1.RData")


for(crt_name in filenames){
  load(here::here(paste0("data/edgeList_",
                         crt_name)))
  
  GDV <- edge_list %>%
    as.integer() %>%
    matrix(ncol = 2, byrow = FALSE) %>%
    count4()

  print(paste0("done: ",
               crt_name))

  save(GDV,
       file = here::here(paste0("data/GDV_",
                                crt_name)))
  rm(GDV)
  gc()
}

