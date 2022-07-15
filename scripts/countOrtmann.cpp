#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;


double chooseC(double n, double k) {
	// as proposed by Dirk Eddelbuettel 
	// https://stackoverflow.com/questions/25005216/n-choose-k-function-crashes-rcpp
  return Rf_choose(n, k);
}

IntegerVector induced_orbits (int crt_node, int n_nodes, int n_edges, 
							  List neighbourhood, IntegerVector deg, IntegerMatrix edge_list,
							  IntegerVector k3, IntegerVector c4, IntegerVector k4,
							  IntegerVector k3_edge[], IntegerVector completing_triangle[]){
								  
	IntegerVector nn(20,0);
	IntegerVector ni(20,0);
	
	IntegerVector crt_N = as<List>(neighbourhood["total_neighbourhood"])[crt_node - 1];
	IntegerVector deg_N = deg[crt_N-1];
	IntegerVector crt_T = completing_triangle[crt_node - 1];
	
	int dv_2 = 0;
	for(int i = 0; i < n_nodes; i++){
		dv_2 += chooseC(deg[i], 2);
	}
	
	//add nn10
	int temp_nn10 = 0;
	for(int v:crt_N){
		IntegerVector v_N = as<List>(neighbourhood["total_neighbourhood"])[v - 1];
		temp_nn10 += sum(as<IntegerVector>(deg[v_N - 1])) - deg[v - 1];
	}
	
	nn[0] = chooseC(n_nodes - 1, 3);
	nn[1] = chooseC(n_nodes - 2, 2) * deg[crt_node - 1];
	nn[2] = (n_edges - deg[crt_node - 1])*(n_nodes - 3);
	nn[3] = crt_N.length()*n_edges - sum(deg_N) - deg[crt_node - 1]*(deg[crt_node - 1] - 1);
	nn[4] = chooseC(deg[crt_node - 1], 2) * (n_nodes - 3);
	nn[5] = (n_nodes - 3)*(sum(deg_N) - deg[crt_node - 1]);
	nn[6] = dv_2 - chooseC(deg[crt_node - 1], 2) - sum(deg_N) + deg[crt_node - 1];
	nn[7] = k3[crt_node - 1]*(n_nodes - 3);
	nn[8] = sum(k3)/3 - k3[crt_node - 1];
	nn[9] = (deg[crt_node - 1] - 1)*(sum(deg_N) - crt_N.length()) - sum(k3_edge[crt_node - 1]);
	nn[10] = temp_nn10 - deg[crt_node - 1] * (deg[crt_node - 1] - 1) - 2*k3[crt_node - 1];
	nn[11] = chooseC(deg[crt_node - 1], 3);
	
	for(int v: crt_N){
		nn[12] += chooseC(deg[v - 1] - 1, 2);
	}
	
	nn[13] = k3[crt_node - 1]*(deg[crt_node - 1] - 2);
	
	for(int i = 0; i < crt_T.length(); i += 2){
		nn[14] += deg[crt_T[i] - 1] + deg[crt_T[i+1] - 1] - 4;
	}
	
	nn[15] = sum(as<IntegerVector>(k3[crt_N - 1])) - sum(k3_edge[crt_node - 1]);
	
	nn[16] = -chooseC(deg[crt_node - 1], 2);
	for(int i = 0; i < crt_N.length(); i++){
		IntegerVector v_N = as<List>(neighbourhood["total_neighbourhood"])[crt_N[i] - 1];
		for(int j = i+1; j < crt_N.length(); j++){
			IntegerVector w_N = as<List>(neighbourhood["total_neighbourhood"])[crt_N[j] - 1];
			nn[16] += intersect(v_N, w_N).length();
		}
	}
	
	nn[17] = -k3[crt_node - 1];
	for(int i = 0; i < crt_T.length(); i += 2){
		nn[17] += as<IntegerVector>(k3_edge[crt_T[i] - 1])[to_string(crt_T[i + 1])];
	}
	
	for(int t:k3_edge[crt_node - 1]){
		nn[18] += chooseC(t,2);
	}
	
	nn[19] = k4[crt_node - 1];
	
	return nn;
}

// [[Rcpp::export]]
List countOrtmann(IntegerMatrix edge_list){
	int n_nodes = max(edge_list);
	int n_edges = edge_list.nrow();
	
	// initialising variables
	IntegerVector k3(n_nodes,0);
	IntegerVector c4(n_nodes,0);
	IntegerVector k4(n_nodes,0);
	
	IntegerVector k3_edge[n_nodes];
	for(int i = 0; i < n_nodes; i++){
		k3_edge[i] = IntegerVector::create(Named("-1",0));
	}
	
	IntegerVector completing_triangle[n_nodes];
	/*for(int i = 0; i < n_nodes; i++){
		k3_edge[i] = IntegerVector::create(Named("-1",0));
	}*/
	
	/*IntegerVector k4_edge[n_nodes];
	for(int i = 0; i < n_nodes; i++){
		k3_edge[i] = IntegerVector::create(Named("-1",0));
	}
	
	IntegerVector c4_edge[n_nodes];
	for(int i = 0; i < n_nodes; i++){
		k3_edge[i] = IntegerVector::create(Named("-1",0));
	}*/
	
	
	IntegerVector mark(n_nodes,0);
	IntegerVector visited(n_nodes,0);
	IntegerVector processed(n_nodes,0);
	
	//call R function
	Function list_neighbourhood("list_neighbourhood");
	Function degree("degree");
	
	// Compute neighbourhoods
	List neighbourhood = list_neighbourhood(edge_list, Named("directed") = true);
									
	// Compute degree vector
	IntegerVector deg = degree(edge_list, Named("directed") = false);
	
	
	// Algorithm
	for(int u = 2; u <= n_nodes; u++){
		IntegerVector u_N_in = as<List>(neighbourhood["in_neighbourhood"])[u-1];
			
		for(int v: u_N_in){
			mark[v-1] ++;
		}
		
		for(int v: u_N_in){
			mark[v-1] --;
			
			IntegerVector v_N = as<List>(neighbourhood["total_neighbourhood"])[v-1];
			IntegerVector v_N_sel = v_N[v_N < u];
			
			for(int w: v_N_sel){
				visited[w-1] ++;
				processed[w-1] ++;
			}
			
			
			IntegerVector v_N_out = as<List>(neighbourhood["out_neighbourhood"])[v-1];
			IntegerVector v_N_out_sel = v_N_out[v_N_out < u];
			
			for(int w: v_N_out_sel){
				mark[w-1] += 2;
			}
			for(int w: v_N_out_sel){
				mark[w-1] -= 2;
				
				if(mark[w-1] !=0){
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
					
					
					IntegerVector w_N_out = as<List>(neighbourhood["out_neighbourhood"])[w-1];
					IntegerVector w_N_out_sel = w_N_out[w_N_out < u];
					
					for(int x: w_N_out_sel){
						if(mark[x-1] == 3){
							k4[u-1] ++;
							k4[v-1] ++;
							k4[w-1] ++;
							k4[x-1] ++;
						}
					}
				}
			}
		}
		
		for(int v: u_N_in){
			IntegerVector v_N = as<List>(neighbourhood["total_neighbourhood"])[v-1];
			IntegerVector v_N_sel = v_N[v_N < u];
			
			for(int w: v_N_sel){
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
		
		
	}
	
	
	// solve system of equations
	IntegerMatrix all_induced_counts (n_nodes, 20);
	
	for(int t = 1; t <= n_nodes; t++){
		all_induced_counts(t-1, _) = induced_orbits(t, n_nodes, n_edges, 
							  neighbourhood, deg, edge_list,
							  k3, c4, k4, k3_edge, completing_triangle);
	}
	
	

	return List::create(k3, c4, k4, all_induced_counts);
}


