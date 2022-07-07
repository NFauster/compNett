########
# This script contains all functions needed. They are grouped for better overview.
######

# Editing networks ----
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
  
  # Main ----
  ## Existing node names
  old_nodes <- edge_list %>%
    as.vector() %>%
    unique()
  
  ## Defining assignment old to new names with "codebook"
  matching_nodes <- data.frame("old" = old_nodes,
                               "new" = 1:length(old_nodes))
  
  ## Renaming the nodes
  edge_list_rec <- edge_list %>%
    data.frame() %>%
    left_join(matching_nodes,
              by = c("X1" = "old")) %>%
    left_join(matching_nodes,
              by = c("X2" = "old")) %>%
    select(new.x, new.y) %>%
    as.matrix
  
  # Return ----
  return(edge_list_rec)
}


redirect_edge <- function(edge){
  # This function takes an edge of a network and redirects it such that it
  # starts at the node with lower degree and points to the node with higher
  # degree.
  # As indicator of which node has the higher degree, the name of the nodes
  # are used. The higher-named node is considered the node with higher degree.
  # Instead of using the degree, other criteria for redirecting may be used.
  # The information has just to be reported over the node names.
  #
  # Arguments
  # edge              a numeric vector with length 2, representing an edge
  #                   between two nodes.
  #
  # Value
  # A numeric vector with length 2, representing an edge starting at the node
  # with lower degree and pointing to the node with the higher degree.
  
  
  # Check parameters ----
  if(!is.vector(edge, mode = "numeric")){
    stop("Argument edge must be a numeric vector with two elements.")
  } else if (length(edge) !=2){
    stop("Argument edge must be a numeric vector with two elements.")
  }
  
  
  # Main ----
  
  names(edge) <- NULL
  if(edge[1] < edge[2]) {
    ## already correct direction
    result <- edge
  } else {
    result <- c(edge[2], edge[1])
  }
  
  return(result)
}


# Network statistics ----
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
  
  # Main ----
  ## Directed network ====
  if(directed){
    ### Out degrees
    from_freq <- as.vector(edge_list[,1]) %>%
      table(useNA = "no")
    
    out_degrees <- as.numeric(from_freq) %>%
      `names<-`(attr(from_freq, "dimnames")[[1]])
    
    ### In degrees
    to_freq <- as.vector(edge_list[,2]) %>%
      table(useNA = "no")
    
    in_degrees <- as.numeric(to_freq) %>%
      `names<-`(attr(to_freq, "dimnames")[[1]])
    
    ## Return
    return(list("in_deg" = in_degrees,
                "out_deg" = out_degrees))
  } else {
    
    ## Undirected network =====
    node_freq <- as.vector(edge_list) %>%
      table(useNA = "no")
    
    degrees <- as.numeric(node_freq) %>%
      `names<-`(attr(node_freq, "dimnames")[[1]])
    
    # Return ----
    return(degrees)
  }
}

list_neighbourhood <- function(edge_list, directed){
  # This function lists the neighbourhood for every node in the network given
  # by the edge list. If the network is directed, the incoming and and outgoing
  # neighbourhood is listed.
  # Note: The neighbourhood of a node "u" is the set of all adjacent nodes, i.e.
  # all nodes which are directly connected by an edge to the node "u". The node
  # "u" is not considered a neighbour to itself in the case of self-loops.
  # For directed networks, the incoming neighbourhood of "u" is the set of all 
  # adjacent nodes, where the connecting edges point towards "u". The outgoing 
  # neighbourhood is defined analogously for connecting edges point away from
  # the node "u" (Brandes, Ortmann 2017: 3f).
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
  # If directed is FALSE, the result is a list with one vector for each node in
  # the network. The names for the list entries are the nodes' names, each 
  # vector gives that node's neighbourhood.
  # If directed is TRUE, the result is a named list with thre list of the kind
  # just described for undirected case. On list is for the incoming neighbour-
  # hood (in_neighbour), one for the outgoing neighbourhood (out_neighbour), one
  # for the total neighbourhood.
  
  
  # Check parameters ----
  if(!is.matrix(edge_list)){
    stop("Argument edge_list must be a two column matrix.")
  } else if (ncol(edge_list) !=2 & ncol(edge_list) >= 1){
    stop("Argument edge_list must be a two column matrix.")
  }
  
  if(!is.logical(directed)){
    stop("Argument directed must be a logical value.")
  }
  
  
  # Main ----
  ## number of nodes in network
  n_node <- max(edge_list)
  
  ## initialise list for incoming and outgoing neighbourhood
  incoming <- vector(mode = "list", length = n_node)
  outgoing <- vector(mode = "list", length = n_node)
  total <- vector(mode = "list", length = n_node)
  
  
  for (i in 1:max(edge_list)) {
    ## look for edges ending at node i
    if (i %in% edge_list[,2]) incoming[[i]] <- edge_list[edge_list[,2] == i,1]
    
    ## look for edges starting at node i
    if (i %in% edge_list[,1]) outgoing[[i]] <- edge_list[edge_list[,1] == i,2]
    
    ## take the intersection of incoming and outgoing neigbourhood
    total[[i]] <- unique(c(incoming[[i]], outgoing[[i]]))
  }
  
  
  
  # Return ----
  if (directed) {
    result <- list("in_neighbourhood" = incoming,
                   "out_neighbourhood" = outgoing,
                   "total_neighbourhood" = total)
  }  else result <- total
  
  return(result)
}