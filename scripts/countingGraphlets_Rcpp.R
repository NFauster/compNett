library(Rcpp)
library(RcppParallel)
library(RcppClock)

# Main ----
sourceCpp(here::here("scripts/countOrtmann.cpp"))
sourceCpp(here::here("scripts/countOrtmann_modSpeed.cpp"))

tictoc::tic()
(count_rcpp_1 <- countOrtmann(edge_list_rename_redirect))
tictoc::toc()
profiling_countOrtmann
# profiling_non_induced

tictoc::tic()
(count_rcpp_modSpeed <- parallelCountOrtmann(edge_list_rename_redirect))
tictoc::toc()
tictoc::tic()
(count_rcpp_modSpeed1 <- parallelCountOrtmann(edge_list))
tictoc::toc()



microbenchmark::microbenchmark("Parallelold" = parallelCountOrtmann(edge_list_rename_redirect),
                               "Parallelnew" = parallelCountOrtmann_everything(edge_list),
                               "orca" = orca::count4(matrix(as.integer(edge_list),ncol=2,byrow=F)))

(gcm_modSpeed <- GCM(count_rcpp_modSpeed))
profile_parallel
# profiling_non_induced

identical(count_rcpp_1, count_rcpp_modSpeed)

(gcm_1 <- GCM(count_rcpp_1))

tictoc::tic()
GCD(gcm, gcm_1)
tictoc::toc()

## Heatmap ====
heatmap(gcm_1)

library(ggplot2)
i <- 1
i <- i+1
gcm %>% 
  as.data.frame() %>%
  tibble::rownames_to_column("f_id") %>%
  tidyr::pivot_longer(-c(f_id), names_to = "samples", values_to = "counts") %>%
  ggplot(aes(x=samples, y=f_id, fill=counts)) + 
  geom_raster() +
  scale_fill_viridis_c()


# Testing ----
sourceCpp(here::here("scripts/countOrtmann_1.cpp"))

# tictoc::tic()
countOrtmann_1(neighbourhood)
# tictoc::toc()



evalCpp("List::create(rep(0,5))")


sourceCpp(here::here("scripts/new_1.cpp"))
(edge_list <- igraph::as_edgelist(igraph::erdos.renyi.game(10, p.or.m = .3)))
sort(table(edge_list))
edge_list_rename_rcpp(edge_list)


microbenchmark::microbenchmark("R" = degree(edge_list, directed = FALSE) %>%
                                 sort() %>%
                                 names() %>%
                                 {
                                   data.frame(original = as.numeric((.)),
                                              ordered = 1:length((.)))
                                 },
                               "rcpp" = edge_list_rename_rcpp(edge_list))


sourceCpp(here::here("scripts/new_2.cpp"))
temp <- parallelCountOrtmann(edge_list_rename_redirect)
parallel


# Make package ----
RcppArmadillo::RcppArmadillo.package.skeleton("testPackage")
compileAttributes(here::here("testPackage/"),verbose=TRUE)
install.packages("testPackage", repos=NULL, type="source")
library(testPackage)

RcppArmadillo::RcppArmadillo.package.skeleton("testPackageParallel")
compileAttributes(here::here("testPackageParallel/"),verbose=TRUE)
install.packages("testPackageParallel", repos=NULL, type="source")
library(testPackageParallel)


