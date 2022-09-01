#### Ortmann, Brandes algorithm for counting graphlets

# Packages ----
library(magrittr)
library(dplyr)

# Functions ----
source(here::here("scripts/functions.R"))


# Prepare the network ----
prep_network <- function(edge_list){
  ## Order nodes ====
  ## according to node degree
  order_nodes <- degree(edge_list, directed = FALSE) %>%
    sort() %>%
    names() %>%
    {
      data.frame(original = as.numeric((.)),
                 ordered = 1:length((.)))
    }
  
  ## Orient graph ====
  edge_list_rename <- edge_list %>%
    as.data.frame() %>%
    left_join(order_nodes,
              by = c("V1" = "original"))%>%
    left_join(order_nodes,
              by = c("V2" = "original"),
              suffix = c(".from", ".to")) %>%
    select(-c("V1", "V2")) %>%
    as.matrix()
  
  
  edge_list_rename_redirect <- apply(edge_list_rename, 
                                     MARGIN = 1, 
                                     FUN = redirect_edge,
                                     simplify = TRUE) %>%
    t()
  
  return(edge_list_rename_redirect)
}


# Counting algorithm ----
countOrtmannR <- function(edge_list_rename_redirect){
  ## Compute neighbourhoods ====
  neighbourhood <- list_neighbourhood(edge_list_rename_redirect,
                                      directed = TRUE)
  
  
  n_nodes <- max(edge_list_rename_redirect)     # number of nodes
  n_edges <- nrow(edge_list_rename_redirect)    # number of edges
  
  
  ## Initialise variables ====
  mark <- vector(mode = "numeric", length = n_nodes)
  visited <- vector(mode = "numeric", length = n_nodes)
  processed <- vector(mode = "numeric", length = n_nodes)
  k3 <- vector(mode = "numeric", length = n_nodes)
  k4 <- vector(mode = "numeric", length = n_nodes)
  c4 <- vector(mode = "numeric", length = n_nodes)
  
  
  ## Algorithm ====
  for(u in 2:n_nodes){
    u_N_in <- neighbourhood$in_neighbourhood[[u]]         # N-(u)
    
    mark[u_N_in] <- mark[u_N_in]+1
    
    for(v in u_N_in){
      mark[v] <- mark[v] - 1
      
      
      v_N <- neighbourhood$total_neighbourhood[[v]]     # N(v)
      v_N_sel <- v_N[v_N < u]                           # {w e N(v): w<u}
      
      visited[v_N_sel] <- visited[v_N_sel] + 1
      processed[v_N_sel] <- processed[v_N_sel] + 1
      
      ####
      v_N_out <- neighbourhood$out_neighbourhood[[v]]   # N+(v)
      v_N_out_sel <- v_N_out[v_N_out < u]               # {w e N+(v): w<u}
      
      mark[v_N_out_sel] <- mark[v_N_out_sel] + 2
      for(w in v_N_out_sel){
        mark[w] <- mark[w] - 2
        
        if(mark[w] != 0){
          k3[c(u,v,w)] <- k3[c(u,v,w)] + 1
          
          ###
          w_N_out <- neighbourhood$out_neighbourhood[[w]]   # N+(w)
          w_N_out_sel <- w_N_out[w_N_out < u]               # {x e N+(w): x<u}
          
          for(x in w_N_out_sel){
            if(mark[x] == 3) k4[c(x,w,v,u)] <- k4[c(x,w,v,u)] + 1
          }
        }
      }
    }
    
    
    for(v in u_N_in){
      v_N <- neighbourhood$total_neighbourhood[[v]]     # N(v)
      v_N_sel <- v_N[v_N < u]                           # {w e N(v): w<u}
      
      for(w in v_N_sel){
        processed[w] <- processed[w] - 1
        if(visited[w] > 0) c4[v] <- c4[v] + visited[w] - 1
        
        if(processed[w] == 0){
          c4[c(u,w)] <- c4[c(u,w)] + choose(visited[w], 2)
          visited[w] <- 0
        }
      }
    }
  }
  # tictoc::toc()
}


