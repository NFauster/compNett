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

l