#### Comparing Networks

# Packages ----
library(Rcpp)
library(orca)
library(dplyr)
library(here)
library(stringr)
library(ggplot2)
library(cowplot)
library(dendextend)

# Functions ----
source(here::here("scripts/functions.R"))
sourceCpp(here::here("scripts/countOrtmann_modSpeed.cpp"))

GCD_comparison <- function(networks_list){
  count_GCM <- function(edge_list){
    count4(matrix(as.integer(edge_list),ncol=2,byrow=F))[,-c(4,13:15)] %>%
      GCM_wo2()
  }
  
  GCM_list <- lapply(networks_list, FUN = count_GCM)
  
  GCDs <- diag(length(networks_list)) %>%
    `diag<-`(0) %>%
    `colnames<-`(names(networks_list))
  
  for(i in 1: (ncol(GCDs) - 1)){
    for(j in (i + 1):ncol(GCDs)){
      GCDs[i,j] <- GCD(GCM_list[[i]], GCM_list[[j]])
      GCDs[j,i] <- GCDs[i,j]
    }
  }
  return(GCDs)
}

NetEmd_comparison <- function(GDV_list){
  count_CGDD <- function(edge_list){
    count4(matrix(as.integer(edge_list),ncol=2,byrow=F))[,-c(4,13:15)] %>%
      CGDD_wo()
  }
  
  CGDD_list <- lapply(GDV_list, FUN = CGDD_wo)
  gc()
  
  NetEmds <- diag(length(networks_list)) %>%
    `diag<-`(0) %>%
    `colnames<-`(names(networks_list))
  
  for(i in 1: (ncol(NetEmds) - 1)){
    for(j in (i + 1):ncol(NetEmds)){
      NetEmds[i,j] <- NetEmd(CGDD_list[[i]], CGDD_list[[j]])
      NetEmds[j,i] <- NetEmds[i,j]
    }
  }
  return(NetEmds)
}


# Networks ----
n_nodes <- rep(c(1e5), 3)
networks_list <- list()

## BB ====
for(n in n_nodes){
  BBNetwork <- PAFit::generate_BB(N = n,
                                  mode_f = "gamma")
  
  networks_list <- append(networks_list, list(BBNetwork$graph[,-3]) %>%
                            `names<-`(paste0("BB", n)))
}



## BA ====
for(n in n_nodes){
  BANetwork <- PAFit::generate_BA(N = n)
  
  networks_list <- append(networks_list, list(BANetwork$graph[,-3]) %>%
                            `names<-`(paste0("BA", n)))
}


## ER ====
for(n in n_nodes){
  networks_list <- append(networks_list, 
                          list(igraph::erdos.renyi.game(n, p.or.m = .3) %>%
                                 igraph::as_edgelist()) %>%
                            `names<-`(paste0("ER", n)))
}


# Network comparison ----
filenames_GCM <- c("ER-1e+05-5e-04_1",
                   "ER-1e+05-5e-04_2",
                   "ER-1e+05-0.001_1",
                   "ER-1e+05-0.001_2",
                   "ER-1e+05-0.005_1",
                   "ER-1e+05-0.005_2",
                   
                   "ER-3e+05-5e-04_1",
                   "ER-3e+05-5e-04_2",
                   "ER-3e+05-0.001_1",
                   
                   "BA-1e+05-5e-04_1",
                   "BA-1e+05-5e-04_2",
                   "BA-1e+05-0.001_1",
                   "BA-1e+05-0.001_2",
                   "BA-1e+05-0.005_1",
                   "BA-1e+05-0.005_2",

                   "BA-3e+05-5e-04_1",
                   "BA-3e+05-5e-04_2",
                   "BA-3e+05-0.001_1",

                   "BBunif-1e+05-5e-04_1",
                   "BBunif-1e+05-5e-04_2",
                   "BBunif-1e+05-0.001_1",
                   "BBunif-1e+05-0.001_2",
                   "BBunif-1e+05-0.005_1",
                   "BBunif-1e+05-0.005_2",
                   
                   "BBunif-3e+05-5e-04_1",
                   "BBunif-3e+05-5e-04_2",
                   "BBunif-3e+05-0.001_1",

                   "BBexp-1e+05-5e-04_1",
                   "BBexp-1e+05-5e-04_2",
                   "BBexp-1e+05-0.001_1",
                   "BBexp-1e+05-0.001_2",
                   "BBexp-1e+05-0.005_1",
                   "BBexp-1e+05-0.005_2",
                   
                   "BBexp-3e+05-5e-04_1",
                   "BBexp-3e+05-5e-04_2"#,
                   # "BBexp-3e+05-0.001_1"
                   )


## GCD ====
gcm_list <- list()

for(i in 1:length(filenames_GCM)){
  load(here(paste0("data/GCM_complete_", filenames_GCM[i], ".RData")))
  
  gcm_list[[i]] <- crt_GCM[-c(4,13:15),-c(4,13:15)]
}

GCDs <- diag(length(gcm_list)) %>%
  `diag<-`(0) %>%
  `colnames<-`(filenames_GCM) %>%
  `row.names<-`(filenames_GCM)

for(i in 1: (ncol(GCDs) - 1)){
  for(j in (i + 1):ncol(GCDs)){
    GCDs[i,j] <- GCD(gcm_list[[i]], gcm_list[[j]])
    GCDs[j,i] <- GCDs[i,j]
  }
}

save(GCDs, file = here::here("data/GCD.RData"))

### Dendrogram ####
temp_colnames <- str_split(colnames(GCDs),
                           "_",
                           simplify = T)[,1]

colour_coding_model <- data.frame("Model" = c("ER", "BA", "BBunif", "BBexp"),
                                  "Model_col" = sjPlot::sjplot_pal(pal = "metro")[1:4])
colour_coding_nodes <- data.frame("Nodes" = c(1e5, 3e5),
                                 "Nodes_col" = sjPlot::sjplot_pal(pal = "metro")[1:2])
colour_coding_density <- data.frame("Density" = c(5e-4, .001,.005,.01,.05),
                                    "Density_col" = sjPlot::sjplot_pal(pal = "metro")[1:5])

bars <- data.frame("Model" = str_split(temp_colnames, "-", simplify = T)[,1],
                   "Nodes" = str_split(temp_colnames, "-", simplify = T)[,2],
                   "Density" = str_split(temp_colnames, "-", simplify = T)[,-c(1,2)] %>%
                     apply(MARGIN = 1, FUN = function(x){
                       if(x[2] == "") return(x[1])
                       else return(paste(x, collapse = "-"))
                     })) %>%
  mutate(Nodes = as.numeric(Nodes),
         Density = as.numeric(Density)) %>%
  
  left_join(colour_coding_model) %>%
  left_join(colour_coding_nodes) %>%
  left_join(colour_coding_density)
  

crt_dend <- GCDs %>%
  `colnames<-`(NULL) %>%
  `row.names<-`(NULL) %>%
  as.dist() %>%
  hclust(method = "average") %>%
  as.dendrogram()


crt_dend %>%
  set("labels", "") %>%
  plot()
colored_bars(colors = bars[,6:4] %>%
               as.matrix(),
             dend = crt_dend,
             rowLabels = c("Density", "Nodes", "Model"),
             y_shift = -.5)
dendrogram.plot <- recordPlot()

plot_for_model <- ggplot(colour_coding_model) + 
  geom_point(aes(x = 1, y = 1, colour = "Model")) + 
  scale_color_manual(name = "Model",
                     values = colour_coding_model$Model_col %>%
                       `names<-`(colour_coding_model$Model),
                     labels = c("ER", "BA", "BB, unif", "BB, exp") %>%
                       `names<-`(colour_coding_model$Model)) + 
  sjPlot::theme_sjplot2() + 
  theme(legend.position = "bottom")

plot_for_nodes <- ggplot(colour_coding_nodes) + 
  geom_point(aes(x = 1, y = 1, colour = "Nodes")) + 
  scale_color_manual(name = "Nodes",
                     values = colour_coding_nodes$Nodes_col %>%
                       `names<-`(colour_coding_nodes$Nodes %>%
                                   as.character()),
                     labels = c("100 000", "300 000") %>%
                       `names<-`(colour_coding_nodes$Nodes %>%
                                   as.character())) + 
  sjPlot::theme_sjplot2() + 
  theme(legend.position = "bottom")

plot_for_density <- ggplot(colour_coding_density) + 
  geom_point(aes(x = 1, y = 1, colour = "Density")) + 
  scale_color_manual(name = "Density",
                     values = colour_coding_density$Density_col %>%
                       `names<-`(colour_coding_density$Density %>%
                                   as.character()),
                     labels = c("0.05%", "0.1%", "0.5%", "1%", "5%") %>%
                       `names<-`(colour_coding_density$Density %>%
                                   as.character())) + 
  sjPlot::theme_sjplot2() + 
  theme(legend.position = "bottom")

legend_model <- get_legend(plot_for_model)
legend_nodes <- get_legend(plot_for_nodes)
legend_density <- get_legend(plot_for_density)
plot_grid(dendrogram.plot, 
          legend_model,
          legend_nodes,
          legend_density,
          ncol = 1,
          rel_heights = c(20,1,1,1))



### Precision-Recall ####
pred_real <- data.frame("real" = logical(),
                        "dist" = numeric())


for(i in 1: (ncol(GCDs) - 1)){
  for(j in (i + 1):ncol(GCDs)){
    crt_colname <- colnames(GCDs)[j] %>%
      str_split("-")
    crt_rowname <- row.names(GCDs)[i] %>%
      str_split("-")
    
    crt_colname <- crt_colname[[1]][1]
    crt_rowname <- crt_rowname[[1]][1]
    
    pred_real <- pred_real %>%
      rbind(data.frame("real" = crt_colname == crt_rowname,
                       "dist" = GCDs[i,j]))
  }
}

PR <- data.frame("Precision" = 1,
                 "Recall" = 0)

for(decision in pred_real$dist){
  tp <- 0
  fp <- 0
  fn <- 0
  
  for(i in 1:nrow(pred_real)){
    if((pred_real$real[i] == TRUE) & (pred_real$dist[i] <= decision)) {
      tp <- tp +1
    } else if ((pred_real$real[i] == TRUE) & (pred_real$dist[i] > decision)) {
      fn <- fn + 1
    } else if ((pred_real$real[i] == FALSE) & (pred_real$dist[i] <= decision)) {
      fp <- fp + 1
    }
  }
  
  PR <- PR %>%
    rbind(data.frame("Precision" = tp/(tp+fp),
                     "Recall" = tp/(tp+fn)))
}

PR <- PR %>%
  group_by(Recall) %>%
  summarise(Precision = max(Precision)) %>%
  arrange(Recall)

save(PR, file = here::here("data/PR_GCD.RData"))


##### Plot #####
aupr <- 0

for(i in 2:nrow(PR)){
  aupr <- aupr + (PR$Recall[i] - PR$Recall[i - 1]) * PR$Precision[i - 1]
}
aupr <- aupr + (1 - PR$Recall[nrow(PR)]) * PR$Precision[nrow(PR)]


PR %>%
  ggplot(aes(x = Recall, y = Precision)) + 
  
  geom_line() + 
  annotate("text",
           x = .7, y = .9,
           label = paste("AUPR:", round(aupr, 3)))




## NetEmd ====
CGDD_list <- list()

for(i in 1:length(filenames_GCM)){
  load(here(paste0("data/GDV_", filenames_GCM[i], ".RData")))
  
  CGDD_list[[i]] <- GDV[,-c(4,13:15)] %>%
    CGDD_wo()
  print(i)
}

NetEmds <- diag(length(CGDD_list)) %>%
  `diag<-`(0) %>%
  `colnames<-`(filenames_GCM) %>%
  `row.names<-`(filenames_GCM)


for(i in 1: (ncol(NetEmds) - 1)){
  for(j in (i + 1):ncol(NetEmds)){
    NetEmds[i,j] <- NetEmd(CGDD_list[[i]], CGDD_list[[j]])
    NetEmds[j,i] <- NetEmds[i,j]
    print(paste(i,j))
  }
}


save(NetEmds, file = here::here("data/NetEmd.RData"))


### Dendrogram ####
temp_colnames <- str_split(colnames(NetEmds),
                           "_",
                           simplify = T)[,1]

colour_coding_model <- data.frame("Model" = c("ER", "BA", "BBunif", "BBexp"),
                                  "Model_col" = sjPlot::sjplot_pal(pal = "metro")[1:4])
colour_coding_nodes <- data.frame("Nodes" = c(1e5, 3e5),
                                  "Nodes_col" = sjPlot::sjplot_pal(pal = "metro")[1:2])
colour_coding_density <- data.frame("Density" = c(5e-4, .001,.005,.01,.05),
                                    "Density_col" = sjPlot::sjplot_pal(pal = "metro")[1:5])

bars <- data.frame("Model" = str_split(temp_colnames, "-", simplify = T)[,1],
                   "Nodes" = str_split(temp_colnames, "-", simplify = T)[,2],
                   "Density" = str_split(temp_colnames, "-", simplify = T)[,-c(1,2)] %>%
                     apply(MARGIN = 1, FUN = function(x){
                       if(x[2] == "") return(x[1])
                       else return(paste(x, collapse = "-"))
                     })) %>%
  mutate(Nodes = as.numeric(Nodes),
         Density = as.numeric(Density)) %>%
  
  left_join(colour_coding_model) %>%
  left_join(colour_coding_nodes) %>%
  left_join(colour_coding_density)


crt_dend <- NetEmds %>%
  `colnames<-`(NULL) %>%
  `row.names<-`(NULL) %>%
  as.dist() %>%
  hclust(method = "average") %>%
  as.dendrogram()


crt_dend %>%
  set("labels", "") %>%
  plot()
colored_bars(colors = bars[,6:4] %>%
               as.matrix(),
             dend = crt_dend,
             rowLabels = c("Density", "Nodes", "Model"))
dendrogram.plot <- recordPlot()

plot_for_model <- ggplot(colour_coding_model) + 
  geom_point(aes(x = 1, y = 1, colour = "Model")) + 
  scale_color_manual(name = "Model",
                     values = colour_coding_model$Model_col %>%
                       `names<-`(colour_coding_model$Model),
                     labels = c("ER", "BA", "BB, unif", "BB, exp") %>%
                       `names<-`(colour_coding_model$Model)) + 
  sjPlot::theme_sjplot2() + 
  theme(legend.position = "bottom")

plot_for_nodes <- ggplot(colour_coding_nodes) + 
  geom_point(aes(x = 1, y = 1, colour = "Nodes")) + 
  scale_color_manual(name = "Nodes",
                     values = colour_coding_nodes$Nodes_col %>%
                       `names<-`(colour_coding_nodes$Nodes %>%
                                   as.character()),
                     labels = c("100 000", "300 000") %>%
                       `names<-`(colour_coding_nodes$Nodes %>%
                                   as.character())) + 
  sjPlot::theme_sjplot2() + 
  theme(legend.position = "bottom")

plot_for_density <- ggplot(colour_coding_density) + 
  geom_point(aes(x = 1, y = 1, colour = "Density")) + 
  scale_color_manual(name = "Density",
                     values = colour_coding_density$Density_col %>%
                       `names<-`(colour_coding_density$Density %>%
                                   as.character()),
                     labels = c("0.05%", "0.1%", "0.5%", "1%", "5%") %>%
                       `names<-`(colour_coding_density$Density %>%
                                   as.character())) + 
  sjPlot::theme_sjplot2() + 
  theme(legend.position = "bottom")

legend_model <- get_legend(plot_for_model)
legend_nodes <- get_legend(plot_for_nodes)
legend_density <- get_legend(plot_for_density)
plot_grid(dendrogram.plot, 
          legend_model,
          legend_nodes,
          legend_density,
          ncol = 1,
          rel_heights = c(20,1,1,1))


### Precision-Recall ####
pred_real <- data.frame("real" = logical(),
                        "dist" = numeric())


for(i in 1: (ncol(NetEmds) - 1)){
  for(j in (i + 1):ncol(NetEmds)){
    crt_colname <- colnames(NetEmds)[j] %>%
      str_split("-")
    crt_rowname <- row.names(NetEmds)[i] %>%
      str_split("-")
    
    crt_colname <- crt_colname[[1]][1]
    crt_rowname <- crt_rowname[[1]][1]
    
    pred_real <- pred_real %>%
      rbind(data.frame("real" = crt_colname == crt_rowname,
                       "dist" = NetEmds[i,j]))
  }
}

PR <- data.frame("Precision" = 1,
                 "Recall" = 0)

for(decision in pred_real$dist){
  tp <- 0
  fp <- 0
  fn <- 0
  
  for(i in 1:nrow(pred_real)){
    if((pred_real$real[i] == TRUE) & (pred_real$dist[i] <= decision)) {
      tp <- tp +1
    } else if ((pred_real$real[i] == TRUE) & (pred_real$dist[i] > decision)) {
      fn <- fn + 1
    } else if ((pred_real$real[i] == FALSE) & (pred_real$dist[i] <= decision)) {
      fp <- fp + 1
    }
  }
  
  PR <- PR %>%
    rbind(data.frame("Precision" = tp/(tp+fp),
                     "Recall" = tp/(tp+fn)))
}

PR <- PR %>%
  group_by(Recall) %>%
  summarise(Precision = max(Precision)) %>%
  arrange(Recall)

save(PR, file = here::here("data/PR_NetEmd.RData"))


##### Plot #####
aupr <- 0

for(i in 2:nrow(PR)){
  aupr <- aupr + (PR$Recall[i] - PR$Recall[i - 1]) * PR$Precision[i - 1]
}
aupr <- aupr + (1 - PR$Recall[nrow(PR)]) * PR$Precision[nrow(PR)]


PR %>%
  ggplot(aes(x = Recall, y = Precision)) + 
  
  geom_line() + 
  annotate("text",
           x = .7, y = .9,
           label = paste("AUPR:", round(aupr, 3)))








###########
i <- 1
i <- i +1
GCM_list[[i]] <- GCM(networks_list[[i]])

j <- i + 1
j <- j +1
cor.test(gdvE[,1], gdvE[,j], method = "spearman")
