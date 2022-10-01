### Plot GCMs for networks

# Packages ----
library(Rcpp)
library(orca)
library(igraph)
library(dplyr)
library(ggplot2)
library(paletteer)

# Functions ----
source(here::here("scripts/functions.R"))
sourceCpp(here::here("scripts/sample_BB.cpp"))
sourceCpp(here::here("scripts/countOrtmann_modSpeed.cpp"))

count_and_GCM_complete <- function(edge_list){
  count4(matrix(as.integer(edge_list),ncol=2,byrow=F)) %>%
    GCM_wo2()
}
count_and_GCM_reduced <- function(edge_list){
  count4(matrix(as.integer(edge_list),ncol=2,byrow=F))[,-c(4,13:15)] %>%
    GCM_wo2()
}

orbit_names_complete <- paste0("O", 0:14)
orbit_names_reduced <- orbit_names_complete[-c(4,13:15)]



# Networks ----
## ER ====
n_nodes <- c(1e5, 3e5, 5e5, 8e5, 1e6)
densities <- c(.0005, .001, .005, .01)

for(i in 2:10){
  for(nodes in n_nodes){
    for(dens in densities){
      edges <- round(dens * nodes * (nodes - 1)/2)
      edge_list <- as_edgelist(sample_gnm(nodes, edges))
      
      ### complete ####
      # crt_GCM <-  edge_list %>%
      #   count_and_GCM_complete()
      # 
      # save(crt_GCM, 
      #      file = here::here(paste0("data/GCM_complete_ER-",
      #                               nodes, "-",
      #                               dens, "_",
      #                               i, ".RData")))
      save(edge_list,
           file = here::here(paste0("data/edgeList_ER-",
                                    nodes, "-",
                                    dens, "_",
                                    i,".RData")))
      # rm(crt_GCM)
      rm(edge_list)
      gc()
      
      
      ### reduced ####
      # crt_GCM <- crt_GCM[-c(4,13:15), -c(4,13:15)]
      
      # save(crt_GCM, 
      #      file = here::here(paste0("data/GCM_reduced_ER-",
      #                               nodes, "-",
      #                               dens, ".RData")))
      # rm(crt_GCM)
      gc()
      
      print(paste("ER", i, nodes, dens, sep = "-"))
    }
  }
}




## BA ====
n_nodes <- c(1e5, 3e5, 5e5, 8e5, 1e6)
densities <- c(.0005, .001, .005, .01, .05)


for(i in 1:10){
  for(nodes in n_nodes){
    for(dens in densities){
      edges <- round(dens * nodes * (nodes - 1)/2)
      edge_list <- as_edgelist(sample_pa(nodes,
                                         m = round(edges/nodes),
                                         directed = FALSE))
      
      ### complete ####
      # crt_GCM <- edge_list %>%
      #   count_and_GCM_complete()
      
      # save(crt_GCM, 
      #      file = here::here(paste0("data/GCM_complete_BA-",
      #                               nodes, "-",
      #                               dens, "_",
      #                               i, ".RData")))
      save(edge_list,
           file = here::here(paste0("data/edgeList_BA-",
                                    nodes, "-",
                                    dens, "_",
                                    i,".RData")))
      # rm(crt_GCM)
      rm(edge_list)
      gc()
      
      ### reduced ####
      # crt_GCM <- crt_GCM[-c(4,13:15), -c(4,13:15)]
      
      # save(crt_GCM, 
      #      file = here::here(paste0("data/GCM_reduced_BA-",
      #                               nodes, "-",
      #                               dens, ".RData")))
      # rm(crt_GCM)
      gc()
      
      print(paste("BA", i, nodes, dens, sep = "-"))
    }
  }
}




## BB unif ====
n_nodes <- c(1e5, 3e5, 5e5, 8e5, 1e6)
densities <- c(.0005, .001, .005, .01)

for(i in 2:10){
  for(nodes in n_nodes){
    for(dens in densities){
      edges <- round(dens * nodes * (nodes - 1)/2)
      edge_list <- sample_BB(nodes, 
                             round(edges/nodes), 
                             1) %>%
        matrix(ncol = 2, byrow = TRUE)
      
      ### complete ####
      # crt_GCM <- edge_list %>%
      #   count_and_GCM_complete()
      
      # save(crt_GCM, 
      #      file = here::here(paste0("data/GCM_complete_BBunif-",
      #                               nodes, "-",
      #                               dens, "_",
      #                               i,".RData")))
      save(edge_list,
           file = here::here(paste0("data/edgeList_BBunif-",
                                    nodes, "-",
                                    dens, "_",
                                    i,".RData")))
      # rm(crt_GCM)
      rm(edge_list)
      gc()
      
      ### reduced ####
      # crt_GCM <- crt_GCM[-c(4,13:15), -c(4,13:15)]
      
      # save(crt_GCM, 
      #      file = here::here(paste0("data/GCM_reduced_BB-",
      #                               nodes, "-",
      #                               pow, ".RData")))
      # rm(crt_GCM)
      gc()
      
      print(paste("BBunif", i, nodes, dens, sep = "-"))
    }
  }
}


## BB exp ====
n_nodes <- c(1e5, 3e5, 5e5, 8e5, 1e6)
densities <- c(.001)

for(i in 1:10){
  for(nodes in n_nodes){
    for(dens in densities){
      edges <- round(dens * nodes * (nodes - 1)/2)
      edge_list <- sample_BB(nodes, 
                             round(edges/nodes), 
                             2) %>%
        matrix(ncol = 2, byrow = TRUE)
      
      ### complete ####
      # crt_GCM <- edge_list %>%
      #   count_and_GCM_complete()
      
      # save(crt_GCM, 
      #      file = here::here(paste0("data/GCM_complete_BBexp-",
      #                               nodes, "-",
      #                               dens, "_",
      #                               i,".RData")))
      save(edge_list,
           file = here::here(paste0("data/edgeList_BBexp-",
                                    nodes, "-",
                                    dens, "_",
                                    i,".RData")))
      # rm(crt_GCM)
      rm(edge_list)
      gc()
      
      ### reduced ####
      # crt_GCM <- crt_GCM[-c(4,13:15), -c(4,13:15)]
      
      # save(crt_GCM, 
      #      file = here::here(paste0("data/GCM_reduced_BB-",
      #                               nodes, "-",
      #                               pow, ".RData")))
      # rm(crt_GCM)
      gc()
      
      print(paste("BBexp", i, nodes, dens, sep = "-"))
    }
  }
}


## BB beta 5-1 ====
n_nodes <- c(1e5, 3e5, 5e5, 8e5, 1e6)
densities <- c(.0005, .001, .005, .01, .05)

for(i in 1:10){
  for(nodes in n_nodes){
    for(dens in densities){
      edges <- round(dens * nodes * (nodes - 1)/2)
      edge_list <- sample_BB(nodes, 
                             round(edges/nodes), 
                             3,
                             shape1 = 5, shape2 = 1) %>%
        matrix(ncol = 2, byrow = TRUE)
      
      ### complete ####
      # crt_GCM <- edge_list %>%
      #   count_and_GCM_complete()
      
      # save(crt_GCM, 
      #      file = here::here(paste0("data/GCM_complete_BBbeta51-",
      #                               nodes, "-",
      #                               dens, "_",
      #                               i,".RData")))
      save(edge_list,
           file = here::here(paste0("data/edgeList_BBbeta51-",
                                    nodes, "-",
                                    dens, "_",
                                    i,".RData")))
      # rm(crt_GCM)
      rm(edge_list)
      gc()
      
      ### reduced ####
      # crt_GCM <- crt_GCM[-c(4,13:15), -c(4,13:15)]
      
      # save(crt_GCM, 
      #      file = here::here(paste0("data/GCM_reduced_BB-",
      #                               nodes, "-",
      #                               pow, ".RData")))
      # rm(crt_GCM)
      gc()
      
      print(paste("BBbeta51", i, nodes, dens, sep = "-"))
    }
  }
}


## BB beta 1-3 ====
n_nodes <- c(1e5, 3e5, 5e5, 8e5, 1e6)
densities <- c(.0005, .001, .005, .01, .05)

for(i in 1:10){
  for(nodes in n_nodes){
    for(dens in densities){
      edges <- round(dens * nodes * (nodes - 1)/2)
      edge_list <- sample_BB(nodes, 
                             round(edges/nodes), 
                             3,
                             shape1 = 1, shape2 = 3) %>%
        matrix(ncol = 2, byrow = TRUE)
      
      ### complete ####
      # crt_GCM <- edge_list %>%
      #   count_and_GCM_complete()
      
      # save(crt_GCM, 
      #      file = here::here(paste0("data/GCM_complete_BBbeta13-",
      #                               nodes, "-",
      #                               dens, "_",
      #                               i,".RData")))
      save(edge_list,
           file = here::here(paste0("data/edgeList_BBbeta13-",
                                    nodes, "-",
                                    dens, "_",
                                    i,".RData")))
      
      # rm(crt_GCM)
      rm(edge_list)
      gc()
      
      ### reduced ####
      # crt_GCM <- crt_GCM[-c(4,13:15), -c(4,13:15)]
      
      # save(crt_GCM, 
      #      file = here::here(paste0("data/GCM_reduced_BB-",
      #                               nodes, "-",
      #                               pow, ".RData")))
      # rm(crt_GCM)
      gc()
      
      print(paste("BBbeta13", i, nodes, dens, sep = "-"))
    }
  }
}




# Plot GCM ----
is_complete <- ncol(crt_GCM) == 15

crt_orbit_names <- if (is_complete) {
  orbit_names_complete
} else orbit_names_reduced

crt_GCM %>% 
  `colnames<-`(crt_orbit_names) %>%
  as.data.frame(row.names = crt_orbit_names) %>%
  
  tibble::rownames_to_column("f_id") %>%
  tidyr::pivot_longer(-c(f_id), 
                      names_to = "samples", 
                      values_to = "counts") %>%
  
  ggplot(aes(x = samples %>%
               forcats::fct_relevel(paste0("O", c(0,2,5,7,10,11,1,4,6,8,9))),
             y = f_id %>%
               forcats::fct_relevel(paste0("O", c(0,2,5,7,10,11,1,4,6,8,9))),
             fill = counts)) + 
  geom_raster() +
  
  scale_fill_gradientn(colours = paletteer_c("grDevices::RdYlBu", 30,
                                             direction = -1),
                       limits = c(-1, 1),
                       name = "") + 
  
  labs(x = "", y = "") + 
  
  theme_minimal() + 
  theme(legend.key.height= unit(2.5, 'cm'),
        legend.key.width= unit(.5, 'cm'))
