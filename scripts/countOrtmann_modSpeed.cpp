#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <RcppClock.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppParallel)]]
//[[Rcpp::depends(RcppClock)]]

using namespace Rcpp;
using namespace arma;
using namespace RcppParallel;
using namespace std;


double chooseC(double n, double k) {
	// as proposed by Dirk Eddelbuettel 
	// https://stackoverflow.com/questions/25005216/n-choose-k-function-crashes-rcpp
  return Rf_choose(n, k);
}


struct non_induced_orbits_parallel : public Worker
{
	// source matrix
	const RVector<int> u_vec;
	const int n_nodes;
	const int n_edges;
	const std::vector<int> (*neighbourhood)[3];
	const RVector<int> deg;
	const RVector<int> k3;
	const RVector<int> c4;
	const RVector<int> k4;
	std::vector<std::vector<int> > k3_edge;
	std::vector<std::vector<int> > completing_triangle;
   
   
	// destination matrix
	RMatrix<int> nn;
   
	// initialize with source and destination
	non_induced_orbits_parallel(const IntegerVector u_vec,
								const int n_nodes, const int n_edges,
								const std::vector<int> (*neighbourhood)[3],
								const IntegerVector deg,
								const IntegerVector k3,
								const IntegerVector c4,
								const IntegerVector k4,
								std::vector<std::vector<int> > k3_edge,
								std::vector<std::vector<int> > completing_triangle,
								IntegerMatrix nn):
			u_vec(u_vec), n_nodes(n_nodes), n_edges(n_edges),
			neighbourhood(neighbourhood),
			deg(deg), k3(k3), c4(c4), k4(k4),
			k3_edge(k3_edge), completing_triangle(completing_triangle),
			nn(nn)
			{}
   
	// take the square root of the range of elements requested
	void operator()(std::size_t begin, std::size_t end) {
		for(auto u = u_vec.begin() + begin; u != u_vec.begin() + end; ++u){
			
			std::vector<int> crt_N = neighbourhood[*u - 1][0];
			/*std::vector<int> deg_N(crt_N.size());
			for(int n:crt_N){
				deg_N.push_back(deg[n - 1]);
			}*/
			
			//sum(deg_N)
			int deg_N_sum = 0;
			for(int n:crt_N){
				deg_N_sum += deg[n - 1];
			}
	
			unsigned int dv_2 = 0;
			for(int i = 0; i < n_nodes; i++){
				dv_2 += chooseC(deg[i], 2);
			}
	
			//add nn10
			unsigned int temp_nn10 = 0;
			for(int v:crt_N){
				std::vector<int> v_N = neighbourhood[v - 1][0];
				
				for(int n:v_N){
					temp_nn10 += deg[n - 1];
				}
				temp_nn10 -= deg[v - 1];
			}
	
			nn(*u - 1, 0) = chooseC(n_nodes - 1, 3);
			nn(*u - 1, 1) = deg[*u - 1];
			nn(*u - 1, 2) = (n_edges - deg[*u - 1])*(n_nodes - 3);
			nn(*u - 1, 3) = crt_N.size()*n_edges - deg_N_sum - deg[*u - 1]*(deg[*u - 1] - 1);
			nn(*u - 1, 4) = chooseC(deg[*u - 1], 2);
			nn(*u - 1, 5) = (deg_N_sum - deg[*u - 1]);
			nn(*u - 1, 6) = dv_2 - chooseC(deg[*u - 1], 2) - deg_N_sum + deg[*u - 1];
			nn(*u - 1, 7) = k3[*u - 1];
			
			nn(*u - 1, 8) = - k3[*u - 1];
			for(int k:k3){
				nn(*u - 1, 8) += k/3;
			}
			
			nn(*u - 1, 9) = (deg[*u - 1] - 1)*(deg_N_sum - crt_N.size());
			for(int k:k3_edge[*u - 1]){
				nn(*u - 1, 9) -= k;
			}

			nn(*u - 1, 10) = temp_nn10 - deg[*u - 1] * (deg[*u - 1] - 1) - 2*k3[*u - 1];
			nn(*u - 1, 11) = chooseC(deg[*u - 1], 3);
			
			nn(*u - 1, 12) = 0;
			for(int v: crt_N){
				nn(*u - 1, 12) += chooseC(deg[v - 1] - 1, 2);
			}
	
			nn(*u - 1, 13) = k3[*u - 1]*(deg[*u - 1] - 2);
	
			for(unsigned int i = 0; i < completing_triangle[*u - 1].size(); i += 2){
				nn(*u - 1, 14) += deg[completing_triangle[*u - 1][i] - 1] + 
						deg[completing_triangle[*u - 1][i+1] - 1] - 4;
			}
	
			nn(*u - 1, 15) = 0;
			for(int n:crt_N){
				nn(*u - 1, 15) += k3[n - 1];
			}
			for(int k:k3_edge[*u - 1]){
				nn(*u - 1, 15) -= k;
			}
	
			nn(*u - 1, 16) = -chooseC(deg[*u - 1], 2);
			for(unsigned int i = 0; i < crt_N.size(); i++){
				std::vector<int> v_N = neighbourhood[crt_N[i] - 1][0];
				
				for(unsigned int j = i+1; j < crt_N.size(); j++){
					std::vector<int> w_N = neighbourhood[crt_N[j] - 1][0];
					std::vector<int> intersection;
					std::set_intersection(v_N.begin(),v_N.end(),
															w_N.begin(), w_N.end(),
															std::back_inserter(intersection));
					nn(*u - 1, 16) += intersection.size();
				}
			}
	
			nn(*u - 1, 17) = -k3[*u - 1];
			for(unsigned int i = 0; i < completing_triangle[*u - 1].size(); i += 2){
				nn(*u - 1, 17) += k3_edge[completing_triangle[*u - 1][i] - 1]
									[distance(neighbourhood[completing_triangle[*u - 1][i]-1][0].begin(),
												lower_bound(neighbourhood[completing_triangle[*u - 1][i]-1][0].begin(), 
															neighbourhood[completing_triangle[*u - 1][i]-1][0].end(), 
															completing_triangle[*u - 1][i + 1])
											)
									];
			}
	
			for(unsigned int t:k3_edge[*u - 1]){
				nn(*u - 1, 18) += chooseC(t,2);
			}
	
			nn(*u - 1, 19) = k4[*u - 1];
		}
   }
};

struct compute_induced_orbits_parallel : public Worker
{
	// source matrix
	const RMatrix<int> nn;   
	const RVector<int> u_vec;
   
	// destination matrix
	RMatrix<int> ni;
   
	// initialize with source and destination
	compute_induced_orbits_parallel(const IntegerMatrix nn,
									IntegerMatrix ni,
									const IntegerVector u_vec):
			nn(nn), ni(ni), u_vec(u_vec)
			{}
   
	// take the square root of the range of elements requested
	void operator()(std::size_t begin, std::size_t end) {
		for(auto u = u_vec.begin() + begin; u != u_vec.begin() + end; ++u){
			
			RMatrix<int>::Row row1 = nn.row(*u - 1);
			arma::vec nn_arma(row1.length());
			for(int i=0; i < row1.length(); i++){
				nn_arma(i) = row1[i];
			}
			
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
			//ni_arma = ni_arma.elem(index_for_orbits);
			
			for(int i=0; i < ni.ncol(); i++){
				ni(*u - 1, i) = ni_arma[index_for_orbits[i]];
			}
		}
   }
};


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
	
	std::vector<std::vector<int> > k3_edge;
	std::vector<std::vector<int> > completing_triangle;
   
	// constructors
	CountOrtmann(int n_nodes_in,
				 const IntegerVector u_vec_in,
				 const std::vector<int> (*neighbourhood_in)[3]):
				 n_nodes(n_nodes_in),
				 u_vec(u_vec_in),
				 neighbourhood(neighbourhood_in),
				 k3(), c4(), k4(), k3_edge(),
				 completing_triangle()
				 {
					 k3.resize(n_nodes, 0);
					 c4.resize(n_nodes, 0);
					 k4.resize(n_nodes, 0);
					 
					 for(int i = 1; i <= n_nodes; i++){
						 k3_edge.push_back(std::vector<int>((*neighbourhood + 3*(i - 1))[0].size()));
						 //k3_edge[i - 1].names() = (*neighbourhood + 3*(i - 1))[0];
						 
						 completing_triangle.push_back(std::vector<int>(0));
					 }
				 }
	CountOrtmann(const CountOrtmann& countOrtmann, Split): 
		n_nodes(countOrtmann.n_nodes),
		u_vec(countOrtmann.u_vec),
		neighbourhood(countOrtmann.neighbourhood),
		k3(), c4(), k4(), k3_edge(),
		completing_triangle()
		{
			k3.resize(n_nodes, 0);
			c4.resize(n_nodes, 0);
			k4.resize(n_nodes, 0);
			
			for(int i = 1; i <= n_nodes; i++){
				k3_edge.push_back(std::vector<int>((*neighbourhood + 3*(i - 1))[0].size()));
				//k3_edge[i - 1].names() = (*neighbourhood + 3*(i - 1))[0];
				
				completing_triangle.push_back(std::vector<int>(0));
			}
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
						k3_edge[*u-1][distance((*neighbourhood + 3*(*u-1))[0].begin(),
												lower_bound((*neighbourhood + 3*(*u-1))[0].begin(), 
															(*neighbourhood + 3*(*u-1))[0].end(), 
															v))] ++;
						k3_edge[v-1][distance((*neighbourhood + 3*(v-1))[0].begin(),
												lower_bound((*neighbourhood + 3*(v-1))[0].begin(), 
															(*neighbourhood + 3*(v-1))[0].end(), 
															*u))] ++;

						// increment (w,v)
						//k3_edge[w-1][to_string(v)] ++;
						//k3_edge[v-1][to_string(w)] ++;
						k3_edge[w-1][distance((*neighbourhood + 3*(w-1))[0].begin(),
												lower_bound((*neighbourhood + 3*(w-1))[0].begin(), 
															(*neighbourhood + 3*(w-1))[0].end(), 
															v))] ++;
						k3_edge[v-1][distance((*neighbourhood + 3*(v-1))[0].begin(),
												lower_bound((*neighbourhood + 3*(v-1))[0].begin(), 
															(*neighbourhood + 3*(v-1))[0].end(), 
															w))] ++;
					
						// increment (u,w)
						//k3_edge[*u-1][to_string(w)] ++;
						//k3_edge[w-1][to_string(*u)] ++;
						k3_edge[*u-1][distance((*neighbourhood + 3*(*u-1))[0].begin(),
												lower_bound((*neighbourhood + 3*(*u-1))[0].begin(), 
															(*neighbourhood + 3*(*u-1))[0].end(), 
															w))] ++;
						k3_edge[w-1][distance((*neighbourhood + 3*(w-1))[0].begin(),
												lower_bound((*neighbourhood + 3*(w-1))[0].begin(), 
															(*neighbourhood + 3*(w-1))[0].end(), 
															*u))] ++;
					
						// add {v,w} to T(u)
						vector<unsigned int> to_insert{v,w};
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
													to_insert.end());
													
					
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
			
			for(size_t j = 0; j < k3_edge[i].size(); j++)
			{
				k3_edge[i][j] += rhs.k3_edge[i][j];
			}
			
			completing_triangle[i].insert(completing_triangle[i].end(),
											rhs.completing_triangle[i].begin(),
											rhs.completing_triangle[i].end());
        }
	}
};

// [[Rcpp::export]]
IntegerMatrix parallelCountOrtmann(IntegerMatrix edge_list) {
		Rcpp::Clock clock;
	int n_nodes = max(edge_list);
	IntegerVector u_vec = seq(1, n_nodes);
	
	//call R function
	Function list_neighbourhood("list_neighbourhood");
	
	
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
	
	clock.tick("total_count");
   // declare the InnerProduct instance that takes a pointer to the vector data
   CountOrtmann countOrtmann(n_nodes, u_vec, neighbourhood);

   // call paralleReduce to start the work
   parallelReduce(1, u_vec.length(), countOrtmann);
   
   clock.tock("total_count");
   
   
   clock.tick("total_compute");
   // ######
   	unsigned int n_edges = edge_list.nrow();
   Function degree("degree");
   IntegerVector deg = degree(edge_list, Named("directed") = false);
   
   // solve system of equations
	IntegerMatrix nn (n_nodes, 20);
	IntegerMatrix ni (n_nodes, 15);
	
	non_induced_orbits_parallel non_induced_orbits_parallel(u_vec,
								n_nodes, n_edges,
								neighbourhood,
								deg,
								wrap(countOrtmann.k3),
								wrap(countOrtmann.c4),
								wrap(countOrtmann.k4),
								countOrtmann.k3_edge,
								countOrtmann.completing_triangle,
								nn);
	parallelFor(0, u_vec.length(), non_induced_orbits_parallel);
	
	compute_induced_orbits_parallel compute_induced_orbits_parallel(nn,ni, u_vec);
	parallelFor(0, u_vec.length(), compute_induced_orbits_parallel);
	
	colnames(ni) = CharacterVector::create("o0", "o1", "o2", "o3", "o4", "o5",
															   "o6", "o7", "o8", "o9", "o10", "o11",
															   "o12", "o13", "o14");
   
   clock.tock("total_compute");
   clock.stop("profile_parallel");

   // return the computed product
   return ni;
}