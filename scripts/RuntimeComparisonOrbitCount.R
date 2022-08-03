# Runtime comparison orca, Ortmann (my implementation)

library(Rcpp)
library(RcppParallel)
library(orca)
library(microbenchmark)
library(dplyr)
library(ggplot2)


sourceCpp(here::here("scripts/countOrtmann_modSpeed.cpp"))


# Generate networks
# n_nodes <- c(seq(100, 1e3, 100),
#              seq(2e3, 1e4, 1e3),
#              seq(2e4, 1e5, 1e4))

n_nodes <- seq(100, 1e3, 100)


benchmark <- data.frame("n_nodes" = numeric(),
                        "network" = character(),
                        "expr" = character(),
                        "time" = numeric())

for(n in n_nodes){
  ## ER ====
  edge_list <- igraph::as_edgelist(igraph::erdos.renyi.game(n, p.or.m = .3))
  
  benchmark <- benchmark %>%
    rbind(microbenchmark("orca" = orca::count4(matrix(as.integer(edge_list),ncol=2,byrow=F)),
                         "ortmann" = parallelCountOrtmann(edge_list),
                         times = 10) %>%
            as.data.frame() %>%
            mutate(n_nodes = n,
                   network = "ER"))
  
  ## BB ====
  BBNetwork <- PAFit::generate_BB(N = n,
                                  mode_f = "power_law")
  edge_list <- BBNetwork$graph[,-3]
  
  benchmark <- benchmark %>%
    rbind(microbenchmark("orca" = orca::count4(matrix(as.integer(edge_list),ncol=2,byrow=F)),
                         "ortmann" = parallelCountOrtmann(edge_list),
                         times = 10) %>%
            as.data.frame() %>%
            mutate(n_nodes = n,
                   network = "BB"))
  
  ## BA ====
  BANetwork <- PAFit::generate_BA(N = n)
  edge_list <- BBNetwork$graph[,-3]
  
  benchmark <- benchmark %>%
    rbind(microbenchmark("orca" = orca::count4(matrix(as.integer(edge_list),ncol=2,byrow=F)),
                         "ortmann" = parallelCountOrtmann(edge_list),
                         times = 10) %>%
            as.data.frame() %>%
            mutate(n_nodes = n,
                   network = "BA"))
  
}


benchmark %>%
  ggplot(aes(x = n_nodes,
             y = time,
             colour = expr)) + 
  geom_point() + 
  facet_wrap(vars(network),
             scales = "free")
