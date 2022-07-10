#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;


double chooseC(double n, double k) {
	// as proposed by Dirk Eddelbuettel 
	// https://stackoverflow.com/questions/25005216/n-choose-k-function-crashes-rcpp
  return Rf_choose(n, k);
}

// [[Rcpp::export]]
List countOrtmann(IntegerMatrix edge_list){
	int n_nodes = max(edge_list);
	//int n_edges = nrow(edge_list);
	
	// initialising variables
	IntegerVector k3(n_nodes,0);
	IntegerVector c4(n_nodes,0);
	IntegerVector k4(n_nodes,0);
	
	IntegerVector mark(n_nodes,0);
	IntegerVector visited(n_nodes,0);
	IntegerVector processed(n_nodes,0);
	
	//call R function
	Function list_neighbourhood("list_neighbourhood");
	
	// Compute neighbourhoods
	List neighbourhood = list_neighbourhood(edge_list,
                                    Named("directed") = true);
	
	
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

	return List::create(k3, c4, k4);
}
