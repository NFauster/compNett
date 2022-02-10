### Generating networks: Bianconi-Barabasi


# Packages ----
library(PAFit)
library(ggraph)
library(igraph)
library(dplyr)


set.seed(7)

# Generating BB ----
BBNetwork <- generate_BB(N = 1e3,
                         mode_f = "power_law")


# Visualisation ----
BBNetwork$graph[,1:2] %>%
  graph_from_edgelist() %>%
  plot()
  
  ggraph()+
  geom_node_point() + geom_edge_link()  +  theme_graph()
