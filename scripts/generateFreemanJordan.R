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
  }else {
    # Undirected network
    node_freq <- as.vector(edge_list) %>%
      table(useNA = "no")
    
    degrees <- as.numeric(node_freq) %>%
      `names<-`(attr(node_freq, "dimnames")[[1]])
    
    ## Return
    return(degrees)
  }
}


# Parameters ----
number_nodes <- 1e3
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
  
  
}






