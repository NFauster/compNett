# Generate a graph with preferential attachment with choice of fitness

# Packages ----
library(magrittr)


# Parameters ----
number_nodes <- 1e3


# Generate graph ----
## RVs ====
fitness <- runif(number_nodes)
R_nodes <- sample(1:1e2, number_nodes, replace = TRUE)


## Initial graph ====
edge_list <- matrix(c(1,1,0), ncol = 3)

n <- max(edge_list[,-3])    # number of nodes in the current (inital) graph


## Iteration ====
while(nrow(edge_list) < number_nodes){
  degrees <- degree(edge_list[,-3, drop = FALSE])
  P <- sample(1:n, 
              R_nodes[], 
              replace = TRUE,
              prob = degrees[names(degrees) == 1:n])
  
  
  n <- n + 1                # one node is added each step
}





degree <- function(edge_list){
  nodes <- as.vector(edge_list) %>%
    unique()
  
  degrees <- sapply(nodes, function(node) {
    apply(edge_list, MARGIN = 1, FUN = function(edge_list) node %in% edge_list) %>%
      sum()
    })
  
  names(degrees) <- nodes
  
  return(degrees)
  }
