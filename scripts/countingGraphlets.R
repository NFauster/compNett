#### Ortmann, Brandes algorithm for counting graphlets

# Packages ----
library(magrittr)
library(dplyr)

# Functions ----
source(here::here("scripts/functions.R"))




# Graph properties ----
n_nodes <- max(edge_list)     # number of nodes
n_edges <- nrow(edge_list)    # number of edges


# Initialise variables ----
mark <- vector(mode = "numeric", length = n_nodes)
visited <- vector(mode = "numeric", length = n_nodes)
processed <- vector(mode = "numeric", length = n_nodes)
k3 <- vector(mode = "numeric", length = n_nodes)
k4 <- vector(mode = "numeric", length = n_nodes)
c4 <- vector(mode = "numeric", length = n_nodes)


# Order nodes ----
## according to node degree
order_nodes <- degree(edge_list, directed = FALSE) %>%
  sort() %>%
  names() %>%
  {
    data.frame(original = as.numeric((.)),
               ordered = 1:length((.)))
  }


# Orient graph ----
edge_list_rename <- edge_list %>%
  as.data.frame() %>%
  left_join(order_nodes,
            by = c("V1" = "original"))%>%
  left_join(order_nodes,
            by = c("V2" = "original"),
            suffix = c(".from", ".to")) %>%
  select(-c("V1", "V2")) %>%
  as.matrix()


#
edge_list_rename_redirect <- apply(edge_list_rename, 
      MARGIN = 1, 
      FUN = redirect_edge,
      simplify = TRUE) %>%
  t()

# Compute neighbourhoods ----
neighbourhood <- list_neighbourhood(edge_list_rename_redirect,
                                    directed = TRUE)


# Algorithm ----
for(u in 2:n_nodes){
  u_N_in <- neighbourhood$in_neighbourhood[[u]] # N-(u)
  
  if (!is.null(u_N_in)){
    mark[u_N_in] <- mark[u_N_in]+1
    
    for(v in u_N_in){
      mark[v] <- mark[v] -1
      
      v_N <- neighbourhood$total_neighbourhood[[v]]   # N(v)
      if(!is.null(v_N)) {
        visited[v_N] <- visited[v_N] + 1
        processed[v_N] <- processed[v_N] + 1
      }
      
      v_N_out <- neighbourhood$out_neighbourhood[[v]]   # N+(v)
      if(!is.null(v_N_out)) {
        mark[v_N_out] <- mark[v_N_out] + 2
        
        for(w in v_N_out){
          mark[w] <- mark[w] - 2
          
          if(mark[w] != 0){
            k3[c(u,v,w)] <- k3[c(u,v,w)] + 1
            
            w_N_out <- neighbourhood$out_neighbourhood[[w]]   # N+(w)
            if(!is.null(w_N_out)){
              for(x in w_N_out){
                if(mark[x] == 3) k4[c(x,w,v,u)] <- k4[c(x,w,v,u)] + 1
              }
            }
          }
        }
      }
     }
  }
  
  
}