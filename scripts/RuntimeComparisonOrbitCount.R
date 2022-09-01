# Runtime comparison orca, Ortmann (my implementation)

library(Rcpp)
library(orca)
library(igraph)
library(microbenchmark)
library(dplyr)
library(ggplot2)

source(here::here("scripts/countingGraphlets.R"))
sourceCpp(here::here("scripts/countOrtmann.cpp"))
sourceCpp(here::here("scripts/countOrtmann_modSpeed.cpp"))


# R vs Rcpp vs RcppParallel ----

n_nodes <- seq(1e2, 4e2, 1e2)
densities <- c(.0005, .001, .005, .01, .05, .1)


benchmark_R <- data.frame("n_nodes" = numeric(),
                        "density" = numeric(),
                        "expr" = character(),
                        "time" = numeric())

for(nodes in n_nodes){
  for(dens in densities){
    edges <- ceiling(dens * nodes * (nodes - 1)/2)
    
    edge_list <- as_edgelist(sample_gnm(nodes, edges))
    
    benchmark_R <- benchmark_R %>%
      rbind(microbenchmark("R" = countOrtmannR(edge_list),
                           "Rcpp" = countOrtmann(prep_network(edge_list)),
                           "RcppParallel" = parallelCountOrtmann(edge_list),
                           times = 5) %>%
              as.data.frame() %>%
              mutate(n_nodes = nodes,
                     density = dens))
    print(paste("ER", nodes, dens, sep = ", "))
  }
}

save(benchmark_R, file = here::here(paste0("data/benchmark_R-ER_",
                                         n,
                                         ".RData")))

benchmark_R %>%
  mutate(time = time /1e9,
         density = sapply(density, FUN = function(x) {
           paste("d =", format(x, scientific = FALSE))
         }) ) %>%
  
  ggplot(aes(x = n_nodes,
             y = time,
             colour = expr)) + 
  geom_point() + 
  
  scale_color_discrete(name = "Implementation") + 
  
  labs(x = "Number of nodes",
       y = "Execution time in s") + 
  
  facet_wrap(vars(density),
             scales = "free_y") + 
  
  sjPlot::theme_sjplot2() + 
  theme(legend.position = "bottom")



# RcppParallel vs Orca for ER different density ----

n_nodes <- seq(1e3, 1e4, 1e3)
densities <- c(.0001, .0005, .001, .005, .01, .05, .1, .2, .3)


benchmark <- data.frame("n_nodes" = numeric(),
                        "density" = numeric(),
                        "expr" = character(),
                        "time" = numeric())

for(nodes in n_nodes){
  for(dens in densities){
    edges <- round(dens * nodes * (nodes - 1)/2)
    
    edge_list <- as_edgelist(sample_gnm(nodes, edges))
    
    benchmark <- benchmark %>%
      rbind(microbenchmark("Orca" = count4(matrix(as.integer(edge_list),ncol=2,byrow=F)),
                           "Ortmann" = parallelCountOrtmann(edge_list),
                           times = 5) %>%
              as.data.frame() %>%
              mutate(n_nodes = nodes,
                     density = dens))
    print(paste("ER", nodes, dens, sep = ", "))
  }
}

save(benchmark, file = here::here(paste0("data/benchmark-ER_",
                                         "5000-02",
                                         ".RData")))

benchmark_R %>%
  mutate(time = time /1e9,
         density = sapply(density, FUN = function(x) {
           paste("d =", format(x, scientific = FALSE))
         }) ) %>%
  
  ggplot(aes(x = n_nodes,
             y = time,
             colour = expr)) + 
  geom_point() + 
  
  scale_color_discrete(name = "Algorithm") + 
  
  labs(x = "Number of nodes",
       y = "Execution time in s") + 
  
  facet_wrap(vars(density),
             scales = "free_y") + 
  
  sjPlot::theme_sjplot2() + 
  theme(legend.position = "bottom")




# Profiling Ortmann-code ----
n_nodes <- 3e2

## ER ====
edge_list <- igraph::as_edgelist(igraph::erdos.renyi.game(n_nodes, p.or.m = .3))

# microbenchmark("orca" = orca::count4(matrix(as.integer(edge_list),ncol=2,byrow=F)),
#                "ortmann" = parallelCountOrtmann(edge_list),
#                times = 1)
parallelCountOrtmann(edge_list)
profile_parallel

## BB ====
BBNetwork <- PAFit::generate_BB(N = n_nodes,
                                mode_f = "power_law")
edge_list <- BBNetwork$graph[,-3]

microbenchmark("orca" = orca::count4(matrix(as.integer(edge_list),ncol=2,byrow=F)),
               "ortmann" = parallelCountOrtmann(edge_list),
               times = 1)
profile_parallel


## BA ====
BANetwork <- PAFit::generate_BA(N = n_nodes)
edge_list <- BBNetwork$graph[,-3]

microbenchmark("orca" = orca::count4(matrix(as.integer(edge_list),ncol=2,byrow=F)),
               "ortmann" = parallelCountOrtmann(edge_list),
               times = 1)
profile_parallel


