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

std::unordered_map<int, int> arrange_names(std::unordered_map<int, int>& M)
{
  
    // Declare a multimap
    std::multimap<int, int> MM;
  
    // Insert every (key-value) pairs from
    // map M to multimap MM as (value-key)
    // pairs
    for (auto& it : M) {
        MM.insert({ it.second, it.first });
    }
	int i = 1;
	for(auto& it : MM){
		M[it.second] = i;
		i++;
	}
  
    return M;
}

std::unordered_map<int, int> degree_rcpp(IntegerMatrix edge_list){
	//IntegerMatrix edge_list_orig = edge_list;
	//std::sort(edge_list.begin(), edge_list.end());
	
	/*int n_unique = std::distance(edge_list.begin(), 
					  std::unique(edge_list.begin(), edge_list.end()));*/
					  
	std::unordered_map<int, int> count;
	int n_edges = edge_list.nrow();
	
	for(int i = 0; i < n_edges; i++){
		if(*(edge_list.begin() + i) != *(edge_list.begin() + n_edges + i)){
			count[*(edge_list.begin() + i)] ++;
			count[*(edge_list.begin() + n_edges + i)] ++;
		}
		
	}
	
	return count;
}



void edge_list_rename_rcpp(IntegerMatrix edge_list, int& n_edges,
							std::vector<int> (*neighbourhood)[3],
							std::vector<int>& deg,
							std::vector<int>& os){
	n_edges = edge_list.nrow();
	
	// rename and redirect
	std::unordered_map<int, int> name_correspondance = degree_rcpp(edge_list);
	arrange_names(name_correspondance);
	
	int n_nodes = name_correspondance.size();
	
	deg.resize(n_nodes);
	os.resize(n_nodes);
	
	int temp_from;
	int temp_to;
	
	for(int i = 0; i < n_edges; i++){
		// change the original node names
		temp_from = min(name_correspondance[edge_list(i,0)], name_correspondance[edge_list(i,1)]);
		temp_to = max(name_correspondance[edge_list(i,0)], name_correspondance[edge_list(i,1)]);
		
		if(temp_from != temp_to){
			neighbourhood[temp_from - 1][0].push_back(temp_to);
			neighbourhood[temp_from - 1][2].push_back(temp_to);
				
			neighbourhood[temp_to - 1][0].push_back(temp_from);
			neighbourhood[temp_to - 1][1].push_back(temp_from);
			
			// the degree are not as defined but count without self-loops
			// (makes computation easier)
			deg[temp_from - 1] ++;
			deg[temp_to - 1] ++;
		} else{
			os[temp_from - 1] ++;
		}
		
	}

	
	// sort neighbourhood
	for(int i = 0; i < n_nodes; i++){
		for(int j = 0; j < 3; j++){
			sort(neighbourhood[i][j].begin(), neighbourhood[i][j].end());
		}
	}
}

IntegerVector non_induced_orbits (unsigned int crt_node, unsigned int n_nodes, unsigned int n_edges, 
							  IntegerVector neighbourhood[][3], IntegerVector deg,
							  IntegerVector k3, IntegerVector c4, IntegerVector k4,
							  IntegerVector k3_edge[], vector<unsigned int> completing_triangle[]){
	Rcpp::Clock clock;
	IntegerVector nn(20,0);
	//IntegerVector ni(20,0);
	
	IntegerVector crt_N = neighbourhood[crt_node - 1][0];
	IntegerVector deg_N = deg[crt_N-1];
	
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
	
	for(unsigned int i = 0; i < completing_triangle[crt_node - 1].size(); i += 2){
		nn[14] += deg[completing_triangle[crt_node - 1][i] - 1] + 
						deg[completing_triangle[crt_node - 1][i+1] - 1] - 4;
	}
	
	nn[15] = sum(as<IntegerVector>(k3[crt_N - 1])) - sum(k3_edge[crt_node - 1]);
	
	nn[16] = c4[crt_node - 1];
	
	nn[17] = -k3[crt_node - 1];
	for(unsigned int i = 0; i < completing_triangle[crt_node - 1].size(); i += 2){
		nn[17] += as<IntegerVector>(k3_edge[completing_triangle[crt_node - 1][i] - 1])[to_string(completing_triangle[crt_node - 1][i + 1])];
	}
	
	for(unsigned int t:k3_edge[crt_node - 1]){
		nn[18] += chooseC(t,2);
	}
	
	nn[19] = k4[crt_node - 1];
	
	//clock.stop("profiling_non_induced");
	
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
IntegerMatrix countOrtmann(IntegerMatrix edge_list){
	int n_nodes = max(edge_list);
	
	int n_edges;
	std::vector<int> temp_neighbourhood[n_nodes][3];
	std::vector<int> temp_deg(n_nodes, 0);
	std::vector<int> os(n_nodes, 0);
	
	edge_list_rename_rcpp(edge_list, n_edges, temp_neighbourhood, temp_deg, os);
	n_nodes = temp_deg.size();
	
	Rcpp::Clock clock;
	
	clock.tick("total_count");
	
	
	IntegerVector neighbourhood[n_nodes][3];
	
	for(unsigned int i = 0; i < n_nodes; i++){
		for(unsigned int j = 0; j < 3; j++){
			neighbourhood[i][j] = wrap(temp_neighbourhood[i][j]);
		}
	}
	
	IntegerVector deg = wrap(temp_deg);
	
	
	
	
	// Initialise variables
	clock.tick("init");
	IntegerVector k3(n_nodes,0);
	IntegerVector c4(n_nodes,0);
	IntegerVector k4(n_nodes,0);
	
	clock.tick("init_k3edge");
	IntegerVector k3_edge[n_nodes];
	for(unsigned int i = 1; i <= n_nodes; i++){
		k3_edge[i - 1] = IntegerVector(neighbourhood[i - 1][0].length());
		k3_edge[i - 1].names() = neighbourhood[i - 1][0];
	}
	clock.tock("init_k3edge");
	
	vector<unsigned int> completing_triangle[n_nodes];
	
	IntegerVector mark(n_nodes,0);
	IntegerVector visited(n_nodes,0);
	IntegerVector processed(n_nodes,0);
	clock.tock("init");
	
	
	// Algorithm
	clock.tick("algorithm");
	for(unsigned int u = 2; u <= n_nodes; u++){
		
		clock.tick("mark_for");
		for(unsigned int v: neighbourhood[u-1][1]){ //loop over N-(u)
			mark[v-1] ++;
		}
		clock.tock("mark_for");
		
		clock.tick("count_k3");
		for(unsigned int v: neighbourhood[u-1][1]){ //loop over N-(u)
			mark[v-1] --;
			
			clock.tick("get_neighbour");
			// {w e N(v):w < u}
			IntegerVector v_N_sel(neighbourhood[v-1][0].begin(), 
								  lower_bound(neighbourhood[v-1][0].begin(), 
								  neighbourhood[v-1][0].end(), u));
			clock.tock("get_neighbour");
			
			for(unsigned int w: v_N_sel){
				visited[w-1] ++;
				processed[w-1] ++;
			}
			
			clock.tick("get_neighbour");
			// {w e N+(v):w < u}
			IntegerVector v_N_out_sel(neighbourhood[v-1][2].begin(), 
									  lower_bound(neighbourhood[v-1][2].begin(),
									  neighbourhood[v-1][2].end(), u));
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
					
					clock.tick("increment_k3_edge");
					// increment (u,v)
					k3_edge[u-1][to_string(v)] = k3_edge[u-1][to_string(v)] + 1;
					k3_edge[v-1][to_string(u)] = k3_edge[v-1][to_string(u)] + 1;
					
					// increment (w,v)
					k3_edge[w-1][to_string(v)] = k3_edge[w-1][to_string(v)] + 1;
					k3_edge[v-1][to_string(w)] = k3_edge[v-1][to_string(w)] + 1;
					
					// increment (u,w)
					k3_edge[u-1][to_string(w)] = k3_edge[u-1][to_string(w)] + 1;
					k3_edge[w-1][to_string(u)] = k3_edge[w-1][to_string(u)] + 1;
					clock.tock("increment_k3_edge");
					
					clock.tick("increment_k3_compTriangle");
					// add {v,w} to T(u)
					vector<unsigned int> to_insert{v,w};
					completing_triangle[u-1].insert(completing_triangle[u-1].end(), 
													to_insert.begin(), 
													to_insert.end());
					
					// add {u,w} to T(v)
					to_insert = {u,w};
					completing_triangle[v-1].insert(completing_triangle[v-1].end(), 
													to_insert.begin(), 
													to_insert.end());
					
					// add {u,v} to T(w)
					to_insert = {u,v};
					completing_triangle[w-1].insert(completing_triangle[w-1].end(), 
													to_insert.begin(), 
													to_insert.end());
													
					clock.tock("increment_k3_compTriangle");
					clock.tock("increment_k3");
					
					clock.tick("get_neighbour");
					// {x e N+(w): x < u}
					IntegerVector w_N_out_sel(neighbourhood[w-1][2].begin(), 
											  lower_bound(neighbourhood[w-1][2].begin(), 
											  neighbourhood[w-1][2].end(), u));
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
		for(unsigned int v: neighbourhood[u-1][1]){ //loop over N-(u)
			// {w e N(v): w < u}
			IntegerVector v_N_sel(neighbourhood[v-1][0].begin(), 
								  lower_bound(neighbourhood[v-1][0].begin(), 
								  neighbourhood[v-1][0].end(), u));
			
			for(unsigned int w: v_N_sel){
				processed[w-1] --;
				
				clock.tick("count_c4_if");
				c4[v-1] += max(visited[w-1] - 1, 0);
				clock.tock("count_c4_if");
				
				if(processed[w-1] == 0){
					c4[u-1] += chooseC(visited[w-1], 2);
					c4[w-1] += chooseC(visited[w-1], 2);
					
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
	
	/*for(unsigned int t = 1; t <= n_nodes; t++){
		all_non_induced_counts(t-1, _) = non_induced_orbits(t, n_nodes, n_edges, 
							  neighbourhood, deg,
							  k3, c4, k4, k3_edge, completing_triangle);
							  
		all_induced_counts(t-1, _) = compute_induced_orbits(all_non_induced_counts(t-1, _));
	}
	
	colnames(all_induced_counts) = CharacterVector::create("o0", "o1", "o2", "o3", "o4", "o5",
															   "o6", "o7", "o8", "o9", "o10", "o11",
															   "o12", "o13", "o14");
															   
	clock.tock("total_count");
	clock.stop("profiling_countOrtmann");*/

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