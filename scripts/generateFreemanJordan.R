# Generate a graph with preferential attachment with choice of fitness

# Packages ----
library(magrittr)


# Functions ----
degree <- function(edge_list, directed = FALSE){
  # This function computes the degree for every node in the network given
  # by the edge list. If the network is directed, the in- and out-degree is
  # computed. 
  # Note: A self loop adds two to the degree of a node or one to each in-
  # and out-degree if the network is directed.
  #
  # Arguments
  # edge_list         a two-column matrix with one row for each edge in the
  #                   network. The two columns indicate start- and end-node
  #                   of the edge. If directed = TRUE, the first column is
  #                   interpreted as starting point.
  # directed          logical. If TRUE, the network is interpreted as directed,
  #                   otherwise undirected.
  #
  # Value
  # If directed is FALSE, the result is a named vector with the node degree
  # as value and the node-name as name.
  # If directed is TRUE, the result is a named list with two named vectors, one
  # for the in-degrees, one for the out-degrees.
  
  require(magrittr)
  
  # Check parameters ----
  if(!is.matrix(edge_list)){
    stop("Argument edge_list must be a two column matrix.")
  } else if (ncol(edge_list) !=2 & ncol(edge_list) >= 1){
    stop("Argument edge_list must be a two column matrix.")
  }
  
  if(!is.logical(directed)){
    stop("Argument directed must be a logical value.")
  }
  
  # Directed network
  if(directed){
    ## Out degrees
    from_freq <- as.vector(edge_list[,1]) %>%
      table(useNA = "no")
    
    out_degrees <- as.numeric(from_freq) %>%
      `names<-`(attr(from_freq, "dimnames")[[1]])
    
    ## In degrees
    to_freq <- as.vector(edge_list[,2]) %>%
      table(useNA = "no")
    
    in_degrees <- as.numeric(to_freq) %>%
      `names<-`(attr(to_freq, "dimnames")[[1]])
    
    ## Return
    return(list("in_deg" = in_degrees,
                "out_deg" = out_degrees))
  } else {
    # Undirected network
    node_freq <- as.vector(edge_list) %>%
      table(useNA = "no")
    
    degrees <- as.numeric(node_freq) %>%
      `names<-`(attr(node_freq, "dimnames")[[1]])
    
    ## Return
    return(degrees)
  }
}

recode_edge_list <- function(edge_list){
  # This function takes an edge_list of a network and renames the nodes
  # using all numbers from 1 to the number of nodes in the network. Thus,
  # node names have the same format (all numeric) and the highest node name
  # is equivalent to the node count.
  # One advantage is that when adding nodes to the network, there is no 
  # danger of using an already existing node name. Simply count upwards from
  # the highest node name.
  #
  # Arguments
  # edge_list         a two-column matrix with one row for each edge in the
  #                   network. The two columns indicate start- and end-node
  #                   of the edge.
  #
  # Value
  # A two-column matrix representing the recoded edge list.
  
  require(magrittr)
  require(dplyr)
  
  # Check parameters ----
  if(!is.matrix(edge_list)){
    stop("Argument edge_list must be a two column matrix.")
  } else if (ncol(edge_list) !=2 & ncol(edge_list) >= 1){
    stop("Argument edge_list must be a two column matrix.")
  }
  
  if(!is.logical(directed)){
    stop("Argument directed must be a logical value.")
  }
  
  # Existing node names
  old_nodes <- edge_list %>%
    as.vector() %>%
    unique()
  
  # Defining assignment old to new names with "codebook"
  matching_nodes <- data.frame("old" = old_nodes,
                               "new" = 1:length(old_nodes))
  
  # Renaming the nodes
  edge_list_rec <- edge_list %>%
    data.frame() %>%
    left_join(matching_nodes,
              by = c("X1" = "old")) %>%
    left_join(matching_nodes,
              by = c("X2" = "old")) %>%
    select(new.x, new.y) %>%
    as.matrix
  
  return(edge_list_rec)
}


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
number_nodes <- 1e3
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




