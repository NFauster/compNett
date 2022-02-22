# Generate a graph with preferential attachment with choice of fitness

# Packages ----
library(magrittr)


# Functions ----
degree <- function(edge_list){
  nodes <- as.vector(edge_list) %>%
    unique() %>%
    na.omit()
  
  
  degrees <- sapply(nodes, function(node) {
    apply(edge_list, MARGIN = 1, FUN = function(edge_list) node %in% edge_list) %>%
      sum()
  })
  
  names(degrees) <- nodes
  
  return(degrees)
}


# Parameters ----
number_nodes <- 1e2
nodes_init_graph <- 1      ## Size of initial (seed) graph
init_graph <- matrix(ncol = 2)     ## Or specify inital graph (edge list)


# Generate graph ----
## RVs ====
fitness <- runif(number_nodes) %>%
  `names<-`(1:number_nodes)
R_nodes <- sample(1:1e1, number_nodes, replace = TRUE)


## Initial graph ====
## TODO: Implement inital graph with parameterized nodes or given seed graph
edge_list <- matrix(nrow = number_nodes, ncol = 3)
edge_list[1,] <- c(1,1,0)

n <- edge_list[,-3] %>%
  as.vector() %>%
  unique() %>% 
  na.omit() %>%
  length()    # number of nodes in the current (inital) graph


## Iteration ====
while(n < number_nodes){
  tictoc::tic()
  degrees <- degree(edge_list[,-3, drop = FALSE])
  P <- sample(1:n, 
              R_nodes[n], 
              replace = TRUE,
              prob = degrees[names(degrees) %in% 1:n]) %>%
    unique()
  
  fittest_node <- fitness[names(fitness) %in% P] %>%
    sort(decreasing = TRUE) %>%
    head(1) %>%
    names() %>%
    as.numeric()
  
  n <- n + 1                # one node is added each step
  
  edge_list[n,] <- c(fittest_node, n, n-1)
  
  tictoc::toc()
}






