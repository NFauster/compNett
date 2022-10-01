## This script checks if Ortmann gives the same results as Orca

library(Rcpp)
library(orca)
library(dplyr)


sourceCpp(here::here("scripts/countOrtmann_modSpeed.cpp"))

n_nodes <- rep(c(1e3, 2e3, 3e3, 4e3, 5e3),1)


# Check methods ----
for(n in n_nodes){
  ## ER ====
  edge_list <- igraph::as_edgelist(igraph::erdos.renyi.game(n, p.or.m = .1))

  count_orca <- count4(matrix(as.integer(edge_list),ncol=2,byrow=F))
  count_ortmann <- parallelCountOrtmann(edge_list)
  print(paste0(n, " nodes, model ER"))
  for(i in 1:15){
    print(paste("Orbit", i-1, ":", identical(table(count_ortmann[,i+1]), 
                                             table(count_orca[,i]))))
  }
  print("")
  
  
  # ## BB ====
  # BBNetwork <- PAFit::generate_BB(N = n,
  #                                 mode_f = "gamma")
  # edge_list <- BBNetwork$graph[,-3]
  # 
  # count_orca <- count4(matrix(as.integer(edge_list),ncol=2,byrow=F))
  # count_ortmann <- parallelCountOrtmann(edge_list)
  # 
  # print(paste0(n, " nodes, model BB(gamma)"))
  # for(i in 1:15){
  #   print(paste("Orbit", i-1, ":", identical(table(count_ortmann[,i]), 
  #                                            table(count_orca[,i]))))
  # }
  # print("")
  # 
  # ## BA ====
  # BANetwork <- PAFit::generate_BA(N = n)
  # edge_list <- BANetwork$graph[,-3]
  # 
  # count_orca <- count4(matrix(as.integer(edge_list),ncol=2,byrow=F))
  # count_ortmann <- parallelCountOrtmann(edge_list)
  # print(paste0(n, " nodes, model BA"))
  # for(i in 1:15){
  #   print(paste("Orbit", i-1, ":", identical(table(count_ortmann[,i]), 
  #                                            table(count_orca[,i]))))
  # }
  # print("")
  
}
