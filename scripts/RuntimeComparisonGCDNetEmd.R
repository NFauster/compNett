
# Packages ----
library(dplyr)
library(stringr)
library(ggplot2)
library(Rcpp)
library(here)
library(microbenchmark)


# Functions ----
source(here::here("scripts/functions.R"))
sourceCpp(here::here("scripts/countOrtmann_modSpeed.cpp"))


compare_GCD <- function(orbit_counts_1, orbit_counts_2){
  gcm1 <- GCM_wo2(orbit_counts_1)
  gcm2 <- GCM_wo2(orbit_counts_2)
  
  gcd12 <- GCD(gcm1, gcm2)
  
  return(gcd12)
}

compare_NetEmd <- function(orbit_counts_1, orbit_counts_2){
  cgdd1 <- CGDD_wo2(orbit_counts_1)
  cgdd2 <- CGDD_wo2(orbit_counts_2)
  
  netEmd12 <- NetEmd(cgdd1, cgdd2)
  
  return(netEmd12)
}



# Run Benchmark ----
benchmark <- data.frame("n_nodes" = numeric(),
                        "density" = numeric(),
                        "expr" = character(),
                        "time" = numeric())

network_names <- c("ER-1e+05-5e-04_1",
                   "ER-1e+05-5e-04_2"#,
                   
                   # "BA-1e+05-5e-04_1",
                   
                   # "BBunif-1e+05-5e-04_1",
                   # "BBunif-1e+05-5e-04_2",
                   
                   # "BBexp-1e+05-5e-04_1"
                   )
  
  
temp_network_names <- str_split(network_names,
                                "_",
                                simplify = T)[,1]

network_characteristics <- data.frame("Model" = str_split(temp_network_names, "-", simplify = T)[,1],
                                      "Nodes" = str_split(temp_network_names, "-", simplify = T)[,2],
                                      "Density" = str_split(temp_network_names, "-", simplify = T)[,-c(1,2)] %>%
                                        apply(MARGIN = 1, FUN = function(x){
                                          if(x[2] == "") return(x[1])
                                          else return(paste(x, collapse = "-"))
                                        }))


for(i in 1: (length(network_names) - 1)){
  for(j in (i + 1): length(network_names)){
    load(here(paste0("data/GDV_",
                     network_names[i],
                     ".RData")))
    orbit_counts_1 <- GDV
    
    load(here(paste0("data/GDV_",
                     network_names[j],
                     ".RData")))
    orbit_counts_2 <- GDV
    
    rm(GDV)
    gc()
    
    
    benchmark <- benchmark %>%
      rbind(microbenchmark("GCD" = compare_GCD(orbit_counts_1, orbit_counts_2),
                           "NetEmd" = compare_NetEmd(orbit_counts_1, orbit_counts_2),
                           times = 5) %>%
              as.data.frame() %>%
              mutate(n_nodes = network_characteristics$Nodes[i],
                     density = network_characteristics$Density[i]))
    
    print(paste(network_names[i], "--- and ---", network_names[j]))
  }
}

## Plot ====
benchmark %>%
  mutate(time = time /1e9,
         density = sapply(density, FUN = function(x) {
           paste("D =", format(x, scientific = FALSE))
         }) ) %>%
  
  ggplot(aes(x = n_nodes,
             y = time,
             colour = expr)) + 
  geom_point() + 
  
  scale_color_discrete(name = "Distance measure",
                       type = sjPlot::sjplot_pal(pal = "metro")) + 
  
  labs(x = "Number of nodes",
       y = "Execution time in s") + 
  
  facet_wrap(vars(density),
             scales = "free_y") + 
  
  sjPlot::theme_sjplot2() + 
  theme(legend.position = "bottom")

