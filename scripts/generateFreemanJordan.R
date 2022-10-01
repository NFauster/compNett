# Generate a graph with preferential attachment with choice of fitness

# Packages ----
library(magrittr)
library(here)

# Functions ----
source(here("scripts/functions.R"))


generate_cycle_graph <- function(n_node){
  # This function generates a network with the specified number of nodes.
  # The nodes of a the graph generated are linked to a circle, i.e. every node
  # has degree two (if undirected) or in- and out-degree both one (if directed).
  #
  # Arguments
  # n_nodes     numeric value specifying the number of nodes the graph should have
  #
  # Value
  # The generated network is given as an edge list, a two column matrix. 
  # The two columns indicate start- and end-node of the edge.
  
  if(!is.numeric(n_node)){
    stop("Argument n_node must be numeric, a whole number greater or equal than 1.")
  } else if (n_node < 1 | n_node %% 1){
    stop("Argument n_node must be numeric, a whole number greater or equal than 1.")
  } else if(n_node == 1){
    return(matrix(c(1,1),
                  ncol = 2))
  } 
  
  matrix(c(1, rep(2:n_node, rep(2, n_node-1)), 1),
         ncol = 2,
         byrow = TRUE)
}




# Parameters ----
number_nodes <- 1e2
nodes_init_graph <- 2      ## Size of initial (seed) graph
init_graph <- matrix(c(1,2,2,3,3,4,4,1,1,3), ncol = 2, byrow = T)     ## Or specify inital graph (edge list)
# init_graph <- NULL


## Check parameters ====
if(!is.numeric(number_nodes)){
  stop("Argument number_nodes must be numeric, a whole number greater or equal than 1.")
} else if (number_nodes < 1 | number_nodes %% 1){
  stop("Argument number_nodes must be numeric, a whole number greater or equal than 1.")
}


if(!is.numeric(nodes_init_graph)){
  stop("Argument nodes_init_graph must be numeric, a whole number greater or equal than 1.")
} else if (nodes_init_graph < 1 | nodes_init_graph %% 1){
  stop("Argument nodes_init_graph must be numeric, a whole number greater or equal than 1.")
} else if (nodes_init_graph == 1){
  message("nodes_init_grpah is chosen 1")
}


if(!is.null(init_graph)){
  if(!is.matrix(init_graph)){
    stop("Argument init_graph must be a two column matrix.")
  } else if (ncol(init_graph) !=2 & ncol(init_graph) >= 1){
    stop("Argument init_graph must be a two column matrix.")
  }
}
  

# Generate graph ----
## RVs ====
fitness <- runif(number_nodes) %>%
  `names<-`(1:number_nodes)
R_nodes <- sample(1:1e1, number_nodes, replace = TRUE)


## Initial graph ====
if(is.null(init_graph)){
  if(nodes_init_graph == 1){
    message("One node with self loop is chosen as seed graph.")
  }
  
  init_graph <- generate_cycle_graph(nodes_init_graph)
} else{
  init_graph <- recode_edge_list(init_graph)
}

crt_nodes <- max(init_graph) # number of nodes in the current (inital) graph
crt_edges <- nrow(init_graph)

edge_list <- matrix(nrow = number_nodes + crt_edges - crt_nodes,
                    ncol = 3)

edge_list[1:crt_edges,] <- cbind(init_graph,
                                        rep(0, crt_edges))

# # TODO: timer
# tictoc::tic.clearlog()


## Iteration ====

# tictoc::tic()
while(crt_nodes < number_nodes){
  # tictoc::tic()
  degrees <- degree(edge_list[,-3, drop = FALSE])
  # tictoc::toc(log = T, quiet = T)
  
  P <- sample(1:crt_nodes, 
              R_nodes[crt_nodes], 
              replace = TRUE,
              prob = degrees[names(degrees) %in% 1:crt_nodes]) %>%
    unique()
  
  fittest_node <- fitness[names(fitness) %in% P] %>%
    sort(decreasing = TRUE) %>%
    head(1) %>%
    names() %>%
    as.numeric()
  
  crt_nodes <- crt_nodes + 1                # one node is added each step
  crt_edges <- crt_edges + 1                # one edge is added each step
  
  edge_list[crt_edges,] <- c(fittest_node,
                             crt_nodes, 
                             edge_list[crt_edges-1,3]+1)
  
  
}
# tictoc::toc()


# t <- tictoc::tic.log(format = T) %>%
#   sapply(function(x) stringr::str_extract(x, "[0123456789.]*")) %>%
#   as.numeric() 
# 
# plot(t)




