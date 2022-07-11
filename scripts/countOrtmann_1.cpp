#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List countOrtmann_1(List neighbourhood){
	int n_nodes = as<List>(neighbourhood["in_neighbourhood"]).length();
	
	IntegerVector t = as<List>(neighbourhood["in_neighbourhood"])[2];
	print(t);
	
	for(int x:t){
		print(t);
	}
	
	//int n_edges = nrow(edge_list);
	
	// initialising variables
	IntegerVector k3(n_nodes,0);
	IntegerVector c4(n_nodes,0);
	IntegerVector k4(n_nodes,0);
	
	IntegerVector mark(n_nodes,0);
	IntegerVector visited(n_nodes,0);
	IntegerVector processed(n_nodes,0);
	

	
	
	// Algorithm
	/*for(int u = 2; u <= n_nodes; u++){
		IntegerVector u_N_in;
		if(as<List>(neighbourhood["in_neighbourhood"])[u-1] != R_NilValue){
			u_N_in = as<List>(neighbourhood["in_neighbourhood"])[u-1];
		}
			
		for(int v: u_N_in){
			mark[v-1] ++;
		}
		
		for(int v: u_N_in){
			mark[v-1] --;
			
			IntegerVector v_N;
			if(as<List>(neighbourhood["total_neighbourhood"])[v-1] != R_NilValue){
				v_N = as<List>(neighbourhood["total_neighbourhood"])[v-1];
			}
			IntegerVector v_N_sel = v_N[v_N < u];
			
			for(int w: v_N_sel){
				visited[w-1] ++;
				processed[w-1] ++;
			}
			
			
			IntegerVector v_N_out;
			if(as<List>(neighbourhood["out_neighbourhood"])[v-1] != R_NilValue){
				v_N_out = as<List>(neighbourhood["out_neighbourhood"])[v-1];
			}
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
					
					IntegerVector w_N_out;
					if(as<List>(neighbourhood["out_neighbourhood"])[w-1] != R_NilValue){
						w_N_out = as<List>(neighbourhood["out_neighbourhood"])[w-1];
					}
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
		
	}*/

	return List::create(k3, c4, k4);
}
