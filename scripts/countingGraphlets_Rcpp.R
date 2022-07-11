library(Rcpp)


sourceCpp(here::here("scripts/countOrtmann.cpp"))

tictoc::tic()
countOrtmann(edge_list_rename_redirect)
tictoc::toc()



sourceCpp(here::here("scripts/countOrtmann_1.cpp"))

# tictoc::tic()
countOrtmann_1(neighbourhood)
# tictoc::toc()



evalCpp("IntegerMatrix(2, 3)")
