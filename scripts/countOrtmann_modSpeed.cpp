#include <RcppArmadillo.h>
#include <RcppClock.h>

// [[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::depends(RcppClock)]]

using namespace Rcpp;
using namespace arma;
using namespace std;


double chooseC(double n, double k) {
	// as proposed by Dirk Eddelbuettel 
	// https://stackoverflow.com/questions/25005216/n-choose-k-function-crashes-rcpp
  return Rf_choose(n, k);
}

double spearman(IntegerVector x, IntegerVector y){
    
    Function spearman("spearman");  
	NumericVector result = spearman(x, y);
	
	LogicalVector check_NA = is_na(result);
	
    if(check_NA[0] == TRUE){
		return 0.0;
	}
	
	return result[0];
}

IntegerVector non_induced_orbits (unsigned int crt_node, unsigned int n_nodes, unsigned int n_edges, 
							  IntegerVector neighbourhood[][3], IntegerVector deg,
							  IntegerVector k3, IntegerVector c4, IntegerVector k4,
							  IntegerVector k3_edge[], IntegerVector completing_triangle[]){
								  
	IntegerVector nn(20,0);
	//IntegerVector ni(20,0);
	
	IntegerVector crt_N = neighbourhood[crt_node - 1][0];
	IntegerVector deg_N = deg[crt_N-1];
	IntegerVector crt_T = completing_triangle[crt_node - 1];
	
	unsigned int dv_2 = 0;
	for(unsigned int i = 0; i < n_nodes; i++){
		dv_2 += chooseC(deg[i], 2);
	}
	
	//add nn10
	unsigned int temp_nn10 = 0;
	for(unsigned int v:crt_N){
		IntegerVector v_N = neighbourhood[v - 1][0];
		temp_nn10 += sum(as<IntegerVector>(deg[v_N - 1])) - deg[v - 1];
	}
	
	nn[0] = chooseC(n_nodes - 1, 3);
	nn[1] = deg[crt_node - 1];
	nn[2] = (n_edges - deg[crt_node - 1])*(n_nodes - 3);
	nn[3] = crt_N.length()*n_edges - sum(deg_N) - deg[crt_node - 1]*(deg[crt_node - 1] - 1);
	nn[4] = chooseC(deg[crt_node - 1], 2);
	nn[5] = (sum(deg_N) - deg[crt_node - 1]);
	nn[6] = dv_2 - chooseC(deg[crt_node - 1], 2) - sum(deg_N) + deg[crt_node - 1];
	nn[7] = k3[crt_node - 1];
	nn[8] = sum(k3)/3 - k3[crt_node - 1];
	nn[9] = (deg[crt_node - 1] - 1)*(sum(deg_N) - crt_N.length()) - sum(k3_edge[crt_node - 1]);
	nn[10] = temp_nn10 - deg[crt_node - 1] * (deg[crt_node - 1] - 1) - 2*k3[crt_node - 1];
	nn[11] = chooseC(deg[crt_node - 1], 3);
	
	for(unsigned int v: crt_N){
		nn[12] += chooseC(deg[v - 1] - 1, 2);
	}
	
	nn[13] = k3[crt_node - 1]*(deg[crt_node - 1] - 2);
	
	for(unsigned int i = 0; i < crt_T.length(); i += 2){
		nn[14] += deg[crt_T[i] - 1] + deg[crt_T[i+1] - 1] - 4;
	}
	
	nn[15] = sum(as<IntegerVector>(k3[crt_N - 1])) - sum(k3_edge[crt_node - 1]);
	
	nn[16] = -chooseC(deg[crt_node - 1], 2);
	for(unsigned int i = 0; i < crt_N.length(); i++){
		IntegerVector v_N = neighbourhood[crt_N[i] - 1][0];
		for(unsigned int j = i+1; j < crt_N.length(); j++){
			IntegerVector w_N = neighbourhood[crt_N[j] - 1][0];
			nn[16] += intersect(v_N, w_N).length();
		}
	}
	
	nn[17] = -k3[crt_node - 1];
	for(unsigned int i = 0; i < crt_T.length(); i += 2){
		nn[17] += as<IntegerVector>(k3_edge[crt_T[i] - 1])[to_string(crt_T[i + 1])];
	}
	
	for(unsigned int t:k3_edge[crt_node - 1]){
		nn[18] += chooseC(t,2);
	}
	
	nn[19] = k4[crt_node - 1];
	
	return nn;
}

IntegerVector compute_induced_orbits(IntegerVector nn){
	arma::vec nn_arma = as<arma::vec>(wrap(nn));
	
	arma::mat LEM = {{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
					 {0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
					 {0, 0, 1, 1, 0, 1, 2, 1, 3, 1, 2, 0, 2, 1, 2, 3, 2, 3, 2, 3},
					 {0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 2, 2, 2, 3},
					 {0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
					 {0, 0, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
					 {0, 0, 0, 0, 0, 0, 1, 0, 3, 0, 1, 0, 1, 0, 1, 3, 1, 3, 1, 3},
					 {0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
					 {0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1},
					 {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 2, 1, 0, 2, 2, 4, 6},
					 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 2, 2, 4, 2, 6},
					 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1},
					 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 2, 1, 3},
					 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 2, 3},
					 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 2, 2, 6},
					 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 2, 0, 3},
					 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 3},
					 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 3},
					 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 3},
					 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}};
					 
	arma::vec ni_arma = solve(LEM, nn_arma);
	arma::uvec index_for_orbits = {1, 5, 4, 7, 10, 9, 12, 11, 16, 15, 14, 13, 17, 18, 19};
	ni_arma = ni_arma.elem(index_for_orbits);
	
	Rcpp::IntegerVector ni   = as<IntegerVector>(wrap(ni_arma));
	return ni;
}


// [[Rcpp::export]]
IntegerMatrix countOrtmann_modSpeed(IntegerMatrix edge_list){
	Rcpp::Clock clock;
	
	clock.tick("total_count");
	
	clock.tick("init");
	unsigned int n_nodes = max(edge_list);
	unsigned int n_edges = edge_list.nrow();
	
	// initialising variables
	IntegerVector k3(n_nodes,0);
	IntegerVector c4(n_nodes,0);
	IntegerVector k4(n_nodes,0);
	
	IntegerVector k3_edge[n_nodes];
	for(unsigned int i = 0; i < n_nodes; i++){
		k3_edge[i] = IntegerVector::create(Named("-1",0));
	}
	
	IntegerVector completing_triangle[n_nodes];
	/*for(unsigned int i = 0; i < n_nodes; i++){
		k3_edge[i] = IntegerVector::create(Named("-1",0));
	}*/
	
	/*IntegerVector k4_edge[n_nodes];
	for(unsigned int i = 0; i < n_nodes; i++){
		k3_edge[i] = IntegerVector::create(Named("-1",0));
	}
	
	IntegerVector c4_edge[n_nodes];
	for(unsigned int i = 0; i < n_nodes; i++){
		k3_edge[i] = IntegerVector::create(Named("-1",0));
	}*/
	
	
	IntegerVector mark(n_nodes,0);
	IntegerVector visited(n_nodes,0);
	IntegerVector processed(n_nodes,0);
	
	clock.tock("init");
	
	//call R function
	Function list_neighbourhood("list_neighbourhood");
	Function degree("degree");
	
	// Compute neighbourhoods
	clock.tick("neighbour");
	List neighbourhood_list = list_neighbourhood(edge_list, Named("directed") = true);
	IntegerVector neighbourhood[n_nodes][3];
	
	for(unsigned int i = 0; i < n_nodes; i++){
		neighbourhood[i][0] = as<List>(neighbourhood_list["total_neighbourhood"])[i];
		neighbourhood[i][1] = as<List>(neighbourhood_list["in_neighbourhood"])[i];
		neighbourhood[i][2] = as<List>(neighbourhood_list["out_neighbourhood"])[i];
	}
	clock.tock("neighbour");
									
	// Compute degree vector
	clock.tick("degree");
	IntegerVector deg = degree(edge_list, Named("directed") = false);
	clock.tock("degree");
	
	
	// Algorithm
	clock.tick("algorithm");
	for(unsigned int u = 2; u <= n_nodes; u++){
		//IntegerVector u_N_in = as<List>(neighbourhood["in_neighbourhood"])[u-1];
		
		clock.tick("mark_for");
		for(unsigned int v: neighbourhood[u-1][1]){
			mark[v-1] ++;
		}
		clock.tock("mark_for");
		
		clock.tick("count_k3");
		for(unsigned int v: neighbourhood[u-1][1]){
			mark[v-1] --;
			
			clock.tick("get_neighbour");
			//IntegerVector v_N = as<List>(neighbourhood["total_neighbourhood"])[v-1];
			IntegerVector v_N_sel(neighbourhood[v-1][0].begin(), lower_bound(neighbourhood[v-1][0].begin(), neighbourhood[v-1][0].end(), u));
			clock.tock("get_neighbour");
			
			for(unsigned int w: v_N_sel){
				visited[w-1] ++;
				processed[w-1] ++;
			}
			
			clock.tick("get_neighbour");
			//IntegerVector v_N_out = as<List>(neighbourhood["out_neighbourhood"])[v-1];
			IntegerVector v_N_out_sel(neighbourhood[v-1][2].begin(), lower_bound(neighbourhood[v-1][2].begin(), neighbourhood[v-1][2].end(), u));
			clock.tock("get_neighbour");
			
			for(unsigned int w: v_N_out_sel){
				mark[w-1] += 2;
			}
			
			clock.tick("w_for");
			for(unsigned int w: v_N_out_sel){
				mark[w-1] -= 2;
				
				if(mark[w-1] !=0){
					clock.tick("increment_k3");
					k3[u-1] ++;
					k3[v-1] ++;
					k3[w-1] ++;
					
					// increment (u,v)
					try{
						k3_edge[u-1][to_string(v)] = k3_edge[u-1][to_string(v)] + 1;
						k3_edge[v-1][to_string(u)] = k3_edge[v-1][to_string(u)] + 1;
						}
					catch(...){
						k3_edge[u-1].insert(k3_edge[u-1].end(), 1);
						as<CharacterVector>(k3_edge[u-1].names())[k3_edge[u-1].length()-1] = to_string(v);
						
						k3_edge[v-1].insert(k3_edge[v-1].end(), 1);
						as<CharacterVector>(k3_edge[v-1].names())[k3_edge[v-1].length()-1] = to_string(u);
					}
					
					// increment (w,v)
					try{
						k3_edge[w-1][to_string(v)] = k3_edge[w-1][to_string(v)] + 1;
						k3_edge[v-1][to_string(w)] = k3_edge[v-1][to_string(w)] + 1;
						}
					catch(...){
						k3_edge[w-1].insert(k3_edge[w-1].end(), 1);
						as<CharacterVector>(k3_edge[w-1].names())[k3_edge[w-1].length()-1] = to_string(v);
						
						k3_edge[v-1].insert(k3_edge[v-1].end(), 1);
						as<CharacterVector>(k3_edge[v-1].names())[k3_edge[v-1].length()-1] = to_string(w);
					}
					
					// increment (u,w)
					try{
						k3_edge[u-1][to_string(w)] = k3_edge[u-1][to_string(w)] + 1;
						k3_edge[w-1][to_string(u)] = k3_edge[w-1][to_string(u)] + 1;
						}
					catch(...){
						k3_edge[u-1].insert(k3_edge[u-1].end(), 1);
						as<CharacterVector>(k3_edge[u-1].names())[k3_edge[u-1].length()-1] = to_string(w);
						
						k3_edge[w-1].insert(k3_edge[w-1].end(), 1);
						as<CharacterVector>(k3_edge[w-1].names())[k3_edge[w-1].length()-1] = to_string(u);
					}
					
					// add {v,w} to T(u)
					completing_triangle[u-1].insert(completing_triangle[u-1].end(), v);
					completing_triangle[u-1].insert(completing_triangle[u-1].end(), w);
					
					// add {u,w} to T(v)
					completing_triangle[v-1].insert(completing_triangle[v-1].end(), u);
					completing_triangle[v-1].insert(completing_triangle[v-1].end(), w);
					
					// add {u,v} to T(w)
					completing_triangle[w-1].insert(completing_triangle[w-1].end(), u);
					completing_triangle[w-1].insert(completing_triangle[w-1].end(), v);
					clock.tock("increment_k3");
					
					clock.tick("get_neighbour");
					//IntegerVector w_N_out = as<List>(neighbourhood["out_neighbourhood"])[w-1];
					IntegerVector w_N_out_sel(neighbourhood[w-1][2].begin(), lower_bound(neighbourhood[w-1][2].begin(), neighbourhood[w-1][2].end(), u));
					clock.tock("get_neighbour");
					
					clock.tick("increment_k4");
					for(unsigned int x: w_N_out_sel){
						clock.tick("increment_k4_if");
						if(mark[x-1] == 3){
							clock.tick("increment_k4_inc");
							k4[u-1] ++;
							k4[v-1] ++;
							k4[w-1] ++;
							k4[x-1] ++;
							clock.tock("increment_k4_inc");
						}
						clock.tock("increment_k4_if");
					}
					clock.tock("increment_k4");
				}
				
			}
			clock.tock("w_for");
		}
		clock.tock("count_k3");
		
		clock.tick("count_c4");
		for(unsigned int v: neighbourhood[u-1][1]){
			//IntegerVector v_N = as<List>(neighbourhood["total_neighbourhood"])[v-1];
			IntegerVector v_N_sel(neighbourhood[v-1][0].begin(), lower_bound(neighbourhood[v-1][0].begin(), neighbourhood[v-1][0].end(), u));
			
			for(unsigned int w: v_N_sel){
				processed[w-1] --;
				
				if(visited[w-1] > 0){
					c4[v-1] += visited[w-1] - 1;
				}
				
				if(processed[w-1] == 0){
					c4[u-1] += chooseC(visited[w-1], 2);
					c4[w-1] += chooseC(visited(w-1), 2);
					
					visited[w-1] = 0;
				}
			}
		}
		clock.tock("count_c4");
		
		
	}
	clock.tock("algorithm");
	
	// solve system of equations
	IntegerMatrix all_non_induced_counts (n_nodes, 20);
	IntegerMatrix all_induced_counts (n_nodes, 15);
	
	for(unsigned int t = 1; t <= n_nodes; t++){
		all_non_induced_counts(t-1, _) = non_induced_orbits(t, n_nodes, n_edges, 
							  neighbourhood, deg,
							  k3, c4, k4, k3_edge, completing_triangle);
							  
		all_induced_counts(t-1, _) = compute_induced_orbits(all_non_induced_counts(t-1, _));
	}
	
	colnames(all_induced_counts) = CharacterVector::create("o0", "o1", "o2", "o3", "o4", "o5",
															   "o6", "o7", "o8", "o9", "o10", "o11",
															   "o12", "o13", "o14");
															   
	clock.tock("total_count");
	clock.stop("profiling_countOrtmann");

	return all_induced_counts;
}

// [[Rcpp::export]]
NumericMatrix GCM(IntegerMatrix induced_orbits){
	unsigned int n_orbits = induced_orbits.ncol();
	
	NumericMatrix GCM = NumericMatrix::diag(n_orbits, 1.0);
	
	for(unsigned int i = 0; i < n_orbits - 1; i++){
		for(unsigned int j = i + 1; j < n_orbits; j++){
			GCM(i, j) =  spearman(induced_orbits(_, i), induced_orbits(_, j));
			GCM(j, i) =  spearman(induced_orbits(_, i), induced_orbits(_, j));
		}
	}
	
	return GCM;
}

// [[Rcpp::export]]
double GCD(NumericMatrix X, NumericMatrix Y){
	unsigned int n_orbits = X.ncol();
	
	double GCD;
	
	for(unsigned int i = 0; i < n_orbits - 1; i++){
		for(unsigned int j = i + 1; j < n_orbits; j++){
			GCD += pow((X(i, j) - Y(i, j)), 2);
		}
	}
	
	return sqrt(GCD);
}