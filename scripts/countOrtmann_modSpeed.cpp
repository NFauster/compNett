#include <RcppArmadillo.h>
#include <RcppParallel.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppParallel)]]

using namespace Rcpp;
using namespace arma;
using namespace RcppParallel;
using namespace std;


double chooseC(double n, double k) {
	// as proposed by Dirk Eddelbuettel 
	// https://stackoverflow.com/questions/25005216/n-choose-k-function-crashes-rcpp
  return Rf_choose(n, k);
}

struct CountOrtmann : public Worker
{
	// inputs
	const int n_nodes;
	const RVector<int> u_vec;
	const std::vector<int> (*neighbourhood)[3];
   
	// vector of counts
	std::vector<int> k3;
	std::vector<int> c4;
	std::vector<int> k4;
   
	// constructors
	CountOrtmann(int n_nodes_in,
				 const IntegerVector u_vec_in,
				 const std::vector<int> (*neighbourhood_in)[3]):
				 n_nodes(n_nodes_in),
				 u_vec(u_vec_in),
				 neighbourhood(neighbourhood_in),
				 k3(), c4(), k4()
				 {
					 k3.resize(n_nodes, 0);
					 c4.resize(n_nodes, 0);
					 k4.resize(n_nodes, 0);
					 }
	CountOrtmann(const CountOrtmann& countOrtmann, Split): 
		n_nodes(countOrtmann.n_nodes),
		u_vec(countOrtmann.u_vec),
		neighbourhood(countOrtmann.neighbourhood),
		k3(), c4(), k4(){
			k3.resize(n_nodes, 0);
			c4.resize(n_nodes, 0);
			k4.resize(n_nodes, 0);
			}

	// process just the elements of the range I've been asked to
	void operator()(std::size_t begin, std::size_t end) {
		std::vector<int> mark(n_nodes,0);
		std::vector<int> visited(n_nodes,0);
		std::vector<int> processed(n_nodes,0);
	   
		for(auto u = u_vec.begin() + begin; u != u_vec.begin() + end; ++u){
			
			for(unsigned int v: (*neighbourhood + 3*(*u-1))[1]){ //loop over N-(u)
				mark[v-1] ++;
			}
			
			for(unsigned int v: (*neighbourhood + 3*(*u-1))[1]){ //loop over N-(u)
				mark[v-1] --;
			
				// {w e N(v):w < u}
				std::vector<int> v_N_sel((*neighbourhood + 3*(v-1))[0].begin(), 
										 lower_bound((*neighbourhood + 3*(v-1))[0].begin(), 
										 (*neighbourhood + 3*(v-1))[0].end(), *u));
			
				for(unsigned int w: v_N_sel){
					visited[w-1] ++;
					processed[w-1] ++;
				}
			
				// {w e N+(v):w < u}
				std::vector<int> v_N_out_sel((*neighbourhood + 3*(v-1))[2].begin(), 
											 lower_bound((*neighbourhood + 3*(v-1))[2].begin(),
											 (*neighbourhood + 3*(v-1))[2].end(), *u));
			
				for(unsigned int w: v_N_out_sel){
					mark[w-1] += 2;
				}
			
				for(unsigned int w: v_N_out_sel){
					mark[w-1] -= 2;
				
					if(mark[w-1] !=0){;
						k3[*u-1] ++;
						k3[v-1] ++;
						k3[w-1] ++;
					
						// increment (u,v)
						/*k3_edge[*u-1][to_string(v)] = k3_edge[*u-1][to_string(v)] + 1;
						k3_edge[v-1][to_string(*u)] = k3_edge[v-1][to_string(*u)] + 1;

						// increment (w,v)
						k3_edge[w-1][to_string(v)] = k3_edge[w-1][to_string(v)] + 1;
						k3_edge[v-1][to_string(w)] = k3_edge[v-1][to_string(w)] + 1;
					
						// increment (u,w)
						k3_edge[*u-1][to_string(w)] = k3_edge[*u-1][to_string(w)] + 1;
						k3_edge[w-1][to_string(*u)] = k3_edge[w-1][to_string(*u)] + 1;*/
					
						// add {v,w} to T(u)
						/*vector<unsigned int> to_insert{v,w};
						completing_triangle[*u-1].insert(completing_triangle[*u-1].end(), 
													to_insert.begin(), 
													to_insert.end());
					
						// add {u,w} to T(v)
						to_insert = {*u,w};
						completing_triangle[v-1].insert(completing_triangle[v-1].end(), 
													to_insert.begin(), 
													to_insert.end());
					
						// add {u,v} to T(w)
						to_insert = {*u,v};
						completing_triangle[w-1].insert(completing_triangle[w-1].end(), 
													to_insert.begin(), 
													to_insert.end());*/
													
					
						// {x e N+(w): x < u}
						std::vector<int> w_N_out_sel((*neighbourhood + 3*(w-1))[2].begin(), 
													 lower_bound((*neighbourhood + 3*(w-1))[2].begin(), 
													 (*neighbourhood + 3*(w-1))[2].end(), *u));
					
						for(unsigned int x: w_N_out_sel){
							if(mark[x-1] == 3){
								k4[*u-1] ++;
								k4[v-1] ++;
								k4[w-1] ++;
								k4[x-1] ++;
							}
						}
					}
				}
			}
			
		   
			for(unsigned int v: (*neighbourhood + 3*(*u-1))[1]){ //loop over N-(u)
				// {w e N(v): w < u}
				std::vector<int> v_N_sel((*neighbourhood + 3*(v-1))[0].begin(), 
										 lower_bound((*neighbourhood + 3*(v-1))[0].begin(), 
										 (*neighbourhood + 3*(v-1))[0].end(), *u));
										 
				for(unsigned int w: v_N_sel){
					processed[w-1] --;
				
					c4[v-1] += max(visited[w-1] - 1, 0);
				
					if(processed[w-1] == 0){
						c4[*u-1] += chooseC(visited[w-1], 2);
						c4[w-1] += chooseC(visited[w-1], 2);
					
						visited[w-1] = 0;
					}
				}
			}
			
		}   
	}

	// join my value with that of another InnerProduct
	void join(const CountOrtmann& rhs) {
		for (size_t i = 0; i < k3.size(); i++) {
			k3[i] += rhs.k3[i];
			c4[i] += rhs.c4[i];
			k4[i] += rhs.k4[i];
        }
	}
};

// [[Rcpp::export]]
List parallelCountOrtmann(IntegerMatrix edge_list) {
		
	int n_nodes = max(edge_list);
	IntegerVector u_vec = seq(2, n_nodes);
	
	//call R function
	Function list_neighbourhood("list_neighbourhood");
	//Function degree("degree");
	
	// Compute neighbourhoods
	List neighbourhood_list = list_neighbourhood(edge_list, Named("directed") = true); //TODO in Rcpp
	std::vector<int> neighbourhood[n_nodes][3];
	
	for(int i = 0; i < n_nodes; i++){
		for(int j = 0; j < as<IntegerVector>(as<List>(neighbourhood_list["total_neighbourhood"])[i]).length(); j++)
		{
			neighbourhood[i][0].push_back(as<IntegerVector>(as<List>(neighbourhood_list["total_neighbourhood"])[i])[j]);
		}
		
		for(int j = 0; j < as<IntegerVector>(as<List>(neighbourhood_list["in_neighbourhood"])[i]).length(); j++)
		{
			neighbourhood[i][1].push_back(as<IntegerVector>(as<List>(neighbourhood_list["in_neighbourhood"])[i])[j]);
		}
		
		for(int j = 0; j < as<IntegerVector>(as<List>(neighbourhood_list["out_neighbourhood"])[i]).length(); j++)
		{
			neighbourhood[i][2].push_back(as<IntegerVector>(as<List>(neighbourhood_list["out_neighbourhood"])[i])[j]);
		}
	}
	
   // declare the InnerProduct instance that takes a pointer to the vector data
   CountOrtmann countOrtmann(n_nodes, u_vec, neighbourhood);

   // call paralleReduce to start the work
   parallelReduce(0, u_vec.length(), countOrtmann);

   // return the computed product
   return List::create(wrap(countOrtmann.k3), wrap(countOrtmann.c4), wrap(countOrtmann.k4));
}