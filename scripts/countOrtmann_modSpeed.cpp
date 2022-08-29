#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <RcppClock.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppClock)]]

using namespace Rcpp;
using namespace arma;
using namespace RcppParallel;
using namespace std;

// Recursive function to return
// gcd of a and b
// as proposed by https://www.geeksforgeeks.org/cpp-program-for-program-to-find-lcm-of-two-numbers/
long long gcd(long long int a,
              long long int b)
{
  if (b == 0)
    return a;
  return gcd(b, a % b);
}
 
// Function to return LCM of
// two numbers
// as proposed by https://www.geeksforgeeks.org/cpp-program-for-program-to-find-lcm-of-two-numbers/
long long lcm(int a, int b)
{
    return (a / gcd(a, b)) * b;
}

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

double sd(NumericVector x){
    
    Function sd("sd"); 
	NumericVector result = sd(x);
	
	return result[0];
}

double spearman2(IntegerVector x, IntegerVector y){
    
    Function spearman("spearman"); 
	x.push_back(1);
	y.push_back(1);
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

// [[Rcpp::export]]
std::unordered_map<int, int> degree_rcpp(IntegerMatrix edge_list){
	//IntegerMatrix edge_list_orig = edge_list;
	//std::sort(edge_list.begin(), edge_list.end());
	
	/*int n_unique = std::distance(edge_list.begin(), 
					  std::unique(edge_list.begin(), edge_list.end()));*/
					  
	std::unordered_map<int, int> count;
	std::unordered_map<int, int> count_self;
	int n_edges = edge_list.nrow();
	
	for(int i = 0; i < n_edges; i++){
		if(*(edge_list.begin() + i) != *(edge_list.begin() + n_edges + i)){
			count[*(edge_list.begin() + i)] ++;
			count[*(edge_list.begin() + n_edges + i)] ++;
		}
		else{
			count_self[*(edge_list.begin() + i)] ++;
		}
		
	}
	
	return count; //TODO:export count_self
}



void edge_list_rename_rcpp(IntegerMatrix edge_list, int& n_edges,
							std::vector<int> (*neighbourhood)[3],
							std::vector<int>& deg){
	n_edges = edge_list.nrow();
	
	// rename and redirect
	std::unordered_map<int, int> name_correspondance = degree_rcpp(edge_list);
	arrange_names(name_correspondance);
	
	std::vector<int> rename_redirect(n_edges * 2);
	
	for(int i = 0; i < n_edges; i++){
		if(name_correspondance[edge_list(i,0)] <= name_correspondance[edge_list(i,1)]){
			rename_redirect[2 * i + 0] = name_correspondance[edge_list(i,0)];
			rename_redirect[2 * i + 1] = name_correspondance[edge_list(i,1)];
		}
		else{
			rename_redirect[2 * i + 0] = name_correspondance[edge_list(i,1)];
			rename_redirect[2 * i + 1] = name_correspondance[edge_list(i,0)];
		}
		
		/*if(rename_redirect[2 * i + 0] == 0){
			rename_redirect[2 * i + 0] = R_NaN;
		}
		if(rename_redirect[2 * i + 1] == 0){
			rename_redirect[2 * i + 1] = R_NaN;
		}*/
	}
	
	
	// list neighbourhood
	int n_nodes = name_correspondance.size();

	//std::vector<int> neighbourhood[n_nodes][3];
		
	//This can be done with the previous loop.
	for(auto edge = rename_redirect.begin(); edge < rename_redirect.end(); edge += 2){
		if(*edge != *(edge + 1)){
				neighbourhood[(*edge - 1)][0].push_back(*(edge + 1));
				neighbourhood[(*edge - 1)][2].push_back(*(edge + 1));
				
				neighbourhood[*(edge + 1) - 1][0].push_back(*edge);
				neighbourhood[*(edge + 1) - 1][1].push_back(*edge);
			}
	}
	
	// sort neighbourhood
	for(int i = 0; i < n_nodes; i++){
		for(int j = 0; j < 3; j++){
			sort(neighbourhood[i][j].begin(), neighbourhood[i][j].end());
		}
	}
	
	
	// get degrees from neighbourhood
	//std::vector<int> deg(n_nodes);
	// this also can be done in the loop through rename_direct
	deg.resize(n_nodes);
	for(int i = 0; i < n_nodes; i++){
		deg[i] = neighbourhood[i][0].size();
	}
}

struct non_induced_orbits_parallel : public Worker
{
	// source matrix
	const RVector<int> u_vec;
	const int n_nodes;
	const int n_edges;
	const std::vector<int> (*neighbourhood)[3];
	const std::vector<int> deg;
	const std::vector<int> k3;
	const std::vector<int> c4;
	const std::vector<int> k4;
	std::vector<std::vector<int> > k3_edge;
	std::vector<std::vector<int> > completing_triangle;
   
   
	// destination matrix
	std::vector<int>& nn;
   
	// initialize with source and destination
	non_induced_orbits_parallel(const IntegerVector u_vec,
								const int n_nodes, 
								const int n_edges,
								const std::vector<int> (*neighbourhood)[3],
								const std::vector<int> deg,
								const std::vector<int> k3,
								const std::vector<int> c4,
								const std::vector<int> k4,
								std::vector<std::vector<int> > k3_edge,
								std::vector<std::vector<int> > completing_triangle,
								std::vector<int>& nn):
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
			
			//sum(deg_N)
			int deg_N_sum = 0;
			for(int n:crt_N){
				deg_N_sum += deg[n - 1];
			}
	
			int dv_2 = 0;
			for(int i = 0; i < n_nodes; i++){
				dv_2 += chooseC(deg[i], 2);
			}
	
			//add nn10
			int temp_nn10 = 0;
			for(int v:crt_N){
				std::vector<int> v_N = neighbourhood[v - 1][0];
				
				for(int n:v_N){
					temp_nn10 += deg[n - 1];
				}
				temp_nn10 -= deg[v - 1];
			}
	
			//nn[(*u - 1) * 20 + 0] = chooseC(n_nodes - 1, 3);
			nn[(*u - 1) * 20 + 1] = deg[*u - 1];
			//nn[(*u - 1) * 20 + 2] = (n_edges - deg[*u - 1])*(n_nodes - 3);
			//nn[(*u - 1) * 20 + 3] = crt_N.size()*n_edges - deg_N_sum - deg[*u - 1]*(deg[*u - 1] - 1);
			nn[(*u - 1) * 20 + 4] = chooseC(deg[*u - 1], 2);
			nn[(*u - 1) * 20 + 5] = (deg_N_sum - deg[*u - 1]);
			//nn[(*u - 1) * 20 + 6] = dv_2 - chooseC(deg[*u - 1], 2) - deg_N_sum + deg[*u - 1];
			nn[(*u - 1) * 20 + 7] = k3[*u - 1];
			
			//nn[(*u - 1) * 20 + 8] = accumulate(k3.begin(), k3.end(), 0)/3 - k3[*u - 1];
			
			nn[(*u - 1) * 20 + 9] = (deg[*u - 1] - 1)*(deg_N_sum - crt_N.size()) - 
									accumulate(k3_edge[*u - 1].begin(),
											   k3_edge[*u - 1].end(), 0);


			nn[(*u - 1) * 20 + 10] = temp_nn10 - deg[*u - 1] * (deg[*u - 1] - 1) - 2*k3[*u - 1];
			nn[(*u - 1) * 20 + 11] = chooseC(deg[*u - 1], 3);
			
			nn[(*u - 1) * 20 + 12] = 0;
			for(int v: crt_N){
				nn[(*u - 1) * 20 + 12] += chooseC(deg[v - 1] - 1, 2);
			}
	
			nn[(*u - 1) * 20 + 13] = k3[*u - 1]*(deg[*u - 1] - 2);
	
			for(int i = 0; i < completing_triangle[*u - 1].size(); i += 2){
				nn[(*u - 1) * 20 + 14] += deg[completing_triangle[*u - 1][i] - 1] + 
						deg[completing_triangle[*u - 1][i+1] - 1] - 4;
			}
	
			nn[(*u - 1) * 20 + 15] = - accumulate(k3_edge[*u - 1].begin(),
										  k3_edge[*u - 1].end(), 0);
			for(int n:crt_N){
				nn[(*u - 1) * 20 + 15] += k3[n - 1];
			}
			
			nn[(*u - 1) * 20 + 16] = c4[*u - 1];
	
			nn[(*u - 1) * 20 + 17] = -k3[*u - 1];
			for(int i = 0; i < completing_triangle[*u - 1].size(); i += 2){
				nn[(*u - 1) * 20 + 17] += k3_edge[completing_triangle[*u - 1][i] - 1]
									[distance(neighbourhood[completing_triangle[*u - 1][i]-1][0].begin(),
												lower_bound(neighbourhood[completing_triangle[*u - 1][i]-1][0].begin(), 
															neighbourhood[completing_triangle[*u - 1][i]-1][0].end(), 
															completing_triangle[*u - 1][i + 1])
											)
									];
			}
	
			for(int t:k3_edge[*u - 1]){
				nn[(*u - 1) * 20 + 18] += chooseC(t,2);
			}
	
			nn[(*u - 1) * 20 + 19] = k4[*u - 1];
		}
   }
};

struct compute_induced_orbits_parallel : public Worker
{
	// source matrix
	const std::vector<int> nn;   
	const RVector<int> u_vec;
   
	// destination matrix
	RMatrix<int> ni;
   
	// initialize with source and destination
	compute_induced_orbits_parallel(const std::vector<int> nn,
									IntegerMatrix ni,
									const IntegerVector u_vec):
			nn(nn), ni(ni), u_vec(u_vec)
			{}
   
	// take the square root of the range of elements requested
	void operator()(std::size_t begin, std::size_t end) {
		for(auto u = u_vec.begin() + begin; u != u_vec.begin() + end; ++u){
			
			arma::vec nn_arma(15);
			std::vector<int> index = {1,4,5,7,9,10,11,12,13,14,15,16,17,18,19};
			for(int i=0; i < nn_arma.size(); i++){
				nn_arma(i) = nn[(*u - 1) * 20 + index[i]];
			}
			
			arma::mat LEM = {
					 {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
					 {0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
					 {0, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
					 {0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
					 {0, 0, 0, 0, 1, 0, 0, 0, 2, 1, 0, 2, 2, 4, 6},
					 {0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 2, 2, 4, 2, 6},
					 {0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1},
					 {0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 2, 1, 3},
					 {0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 2, 3},
					 {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 2, 2, 6},
					 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 2, 0, 3},
					 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 3},
					 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 3},
					 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 3},
					 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}};
					 
			arma::vec ni_arma = solve(LEM, nn_arma);
			arma::uvec index_for_orbits = {0,2,1,3,5,4,7,6,11,10,9,8,12,13,14};
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
	const RVector< int> u_vec;
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
				
				completing_triangle.push_back(std::vector<int>(0));
			}
		}

	// process just the elements of the range I've been asked to
	void operator()(std::size_t begin, std::size_t end) {
		std::vector<int> mark(n_nodes,0);
		std::vector<int> visited(n_nodes,0);
		std::vector<int> processed(n_nodes,0);
	   
		for(auto u = u_vec.begin() + begin; u != u_vec.begin() + end; ++u){
			
			for(int v: (*neighbourhood + 3*(*u-1))[1]){ //loop over N-(u)
				mark[v-1] ++;
			}
			
			for(int v: (*neighbourhood + 3*(*u-1))[1]){ //loop over N-(u)
				mark[v-1] --;
			
				// {w e N(v):w < u}
				std::vector<int> v_N_sel((*neighbourhood + 3*(v-1))[0].begin(), 
										 lower_bound((*neighbourhood + 3*(v-1))[0].begin(), 
										 (*neighbourhood + 3*(v-1))[0].end(), *u));
			
				for(int w: v_N_sel){
					visited[w-1] ++;
					processed[w-1] ++;
				}
			
				// {w e N+(v):w < u}
				std::vector<int> v_N_out_sel((*neighbourhood + 3*(v-1))[2].begin(), 
											 lower_bound((*neighbourhood + 3*(v-1))[2].begin(),
											 (*neighbourhood + 3*(v-1))[2].end(), *u));
			
				for(int w: v_N_out_sel){
					mark[w-1] += 2;
				}
			
				for(int w: v_N_out_sel){
					mark[w-1] -= 2;
				
					if(mark[w-1] !=0){
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
						k3_edge[w-1][distance((*neighbourhood + 3*(w-1))[0].begin(),
												lower_bound((*neighbourhood + 3*(w-1))[0].begin(), 
															(*neighbourhood + 3*(w-1))[0].end(), 
															v))] ++;
						k3_edge[v-1][distance((*neighbourhood + 3*(v-1))[0].begin(),
												lower_bound((*neighbourhood + 3*(v-1))[0].begin(), 
															(*neighbourhood + 3*(v-1))[0].end(), 
															w))] ++;
					
						// increment (u,w)
						k3_edge[*u-1][distance((*neighbourhood + 3*(*u-1))[0].begin(),
												lower_bound((*neighbourhood + 3*(*u-1))[0].begin(), 
															(*neighbourhood + 3*(*u-1))[0].end(), 
															w))] ++;
						k3_edge[w-1][distance((*neighbourhood + 3*(w-1))[0].begin(),
												lower_bound((*neighbourhood + 3*(w-1))[0].begin(), 
															(*neighbourhood + 3*(w-1))[0].end(), 
															*u))] ++;
					
						// add {v,w} to T(u)
						vector<int> to_insert{v,w};
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
					
						for(int x: w_N_out_sel){
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
			
		   
			for(int v: (*neighbourhood + 3*(*u-1))[1]){ //loop over N-(u)
				// {w e N(v): w < u}
				std::vector<int> v_N_sel((*neighbourhood + 3*(v-1))[0].begin(), 
										 lower_bound((*neighbourhood + 3*(v-1))[0].begin(), 
										 (*neighbourhood + 3*(v-1))[0].end(), *u));
										 
				for(int w: v_N_sel){
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
	
	clock.tick("total_time");
	clock.tick("pre_proc");
	int n_nodes = max(edge_list);
	
	int n_edges;
	std::vector<int> neighbourhood[n_nodes][3];
	std::vector<int> deg(n_nodes);
	
	edge_list_rename_rcpp(edge_list, n_edges, neighbourhood, deg);
	n_nodes = deg.size();
	clock.tock("pre_proc");
	
	clock.tick("counting_graphlet");
	IntegerVector u_vec = seq(1, n_nodes);
	
	clock.tick("total_count");
   
   CountOrtmann countOrtmann(n_nodes, u_vec, neighbourhood);

   // call paralleReduce to start the work
   parallelReduce(1, u_vec.length(), countOrtmann);
   
   clock.tock("counting_graphlet");
   
   
   clock.tick("compute_nn");
   // ######
   
   // solve system of equations
	std::vector<int> nn(n_nodes * 20);
	IntegerMatrix ni (n_nodes, 15);
	
	non_induced_orbits_parallel non_induced_orbits_parallel(u_vec,
								n_nodes, n_edges,
								neighbourhood,
								deg,
								countOrtmann.k3,
								countOrtmann.c4,
								countOrtmann.k4,
								countOrtmann.k3_edge,
								countOrtmann.completing_triangle,
								nn);
	parallelFor(0, u_vec.length(), non_induced_orbits_parallel);
	clock.tock("compute_nn");
	
	for(int i = 0; i <100; i++){
		Rprintf("%i\n", nn[i]);
		if(i%5 == 4) Rprintf("\n");
	}
	
	clock.tick("compute_ni");
	compute_induced_orbits_parallel compute_induced_orbits_parallel(nn,ni, u_vec);
	parallelFor(0, u_vec.length(), compute_induced_orbits_parallel);
	
	colnames(ni) = CharacterVector::create("o0", "o1", "o2", "o3", "o4", "o5",
															   "o6", "o7", "o8", "o9", "o10", "o11",
															   "o12", "o13", "o14");
   clock.tock("compute_ni");
   clock.tock("total_time");
   clock.stop("profile_parallel");

   // return the computed product
   return ni;
}



// [[Rcpp::export]]
NumericMatrix GCM(IntegerMatrix edge_list){
	IntegerMatrix ni = parallelCountOrtmann(edge_list);
	
	unsigned int n_orbits = ni.ncol();
	
	NumericMatrix GCM = NumericMatrix::diag(n_orbits, 1.0);
	
	for(unsigned int i = 0; i < n_orbits - 1; i++){
		for(unsigned int j = i + 1; j < n_orbits; j++){
			GCM(i, j) =  spearman(ni(_, i), ni(_, j));
			GCM(j, i) =  spearman(ni(_, i), ni(_, j));
		}
	}
	
	return GCM;
}

// [[Rcpp::export]]
NumericMatrix GCM_wo(IntegerMatrix ni){
	unsigned int n_orbits = ni.ncol();
	
	NumericMatrix GCM = NumericMatrix::diag(n_orbits, 1.0);
	
	for(unsigned int i = 0; i < n_orbits - 1; i++){
		for(unsigned int j = i + 1; j < n_orbits; j++){
			GCM(i, j) =  spearman(ni(_, i), ni(_, j));
			GCM(j, i) =  spearman(ni(_, i), ni(_, j));
		}
	}
	
	return GCM;
}

// [[Rcpp::export]]
NumericMatrix GCM_wo2(IntegerMatrix ni){
	unsigned int n_orbits = ni.ncol();
	
	NumericMatrix GCM = NumericMatrix::diag(n_orbits, 1.0);
	
	for(unsigned int i = 0; i < n_orbits - 1; i++){
		for(unsigned int j = i + 1; j < n_orbits; j++){
			GCM(i, j) =  spearman2(ni(_, i), ni(_, j));
			GCM(j, i) =  spearman2(ni(_, i), ni(_, j));
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

// [[Rcpp::export]]
NumericMatrix GDD(IntegerMatrix edge_list){
	IntegerMatrix ni = parallelCountOrtmann(edge_list);
	
	unsigned int n_orbits = ni.ncol();
	unsigned int n_nodes = ni.nrow();
	int k_max = max(ni);
	
	NumericMatrix GDD = NumericMatrix(k_max + 1, n_orbits);
	
	for(unsigned int orbit = 0; orbit < n_orbits; orbit++){
		for(unsigned int n = 0; n < n_nodes; n++){
			GDD(ni(n, orbit), orbit) += 1.0/n_nodes;
		}
		
		GDD(_, orbit) = GDD(_, orbit) / sd(GDD(_, orbit));
	}
	
	return GDD;
}

// [[Rcpp::export]]
NumericMatrix GDD_wo(IntegerMatrix ni){
	
	unsigned int n_orbits = ni.ncol();
	unsigned int n_nodes = ni.nrow();
	int k_max = max(ni);
	
	NumericMatrix GDD = NumericMatrix(k_max + 1, n_orbits);
	
	for(unsigned int orbit = 0; orbit < n_orbits; orbit++){
		for(unsigned int n = 0; n < n_nodes; n++){
			GDD(ni(n, orbit), orbit) += 1.0/n_nodes;
		}
		
		GDD(_, orbit) = GDD(_, orbit) / sd(GDD(_, orbit));
	}
	
	
	return GDD;
}


// [[Rcpp::export]]
NumericMatrix CGDD(IntegerMatrix edge_list){
	NumericMatrix gdd = GDD(edge_list);
	
	unsigned int n_orbits = gdd.ncol();
	int ks = gdd.nrow();
	
	NumericMatrix cgdd_temp = gdd;
	
	for(unsigned int orbit = 0; orbit < n_orbits; orbit++){
		for(unsigned int k = 1; k < ks; k++){
			cgdd_temp(k, orbit) += cgdd_temp(k - 1, orbit);
		}
	}
	
	int sample_freq = floor(gdd.nrow()/4000) + 1;
	int new_nrow = floor(gdd.nrow()/sample_freq) + 1;
	
	NumericMatrix CGDD(new_nrow , gdd.ncol());
	
	for(unsigned int orbit = 0; orbit < n_orbits; orbit++){
		for(unsigned int k = 0; k < new_nrow - 1; k++){
			CGDD(k, orbit) = cgdd_temp(sample_freq * k, orbit);
		}
		CGDD(new_nrow - 1, orbit) = sample_freq;
	}
	
	return CGDD;
}

// [[Rcpp::export]]
NumericMatrix CGDD_wo(IntegerMatrix ni){
	NumericMatrix gdd = GDD_wo(ni);
	
	unsigned int n_orbits = gdd.ncol();
	int ks = gdd.nrow();
	
	NumericMatrix cgdd_temp = gdd;
	
	for(unsigned int orbit = 0; orbit < n_orbits; orbit++){
		for(unsigned int k = 1; k < ks; k++){
			cgdd_temp(k, orbit) += cgdd_temp(k - 1, orbit);
		}
	}
	
	int sample_freq = floor(gdd.nrow()/4000) + 1;
	int new_nrow = floor(gdd.nrow()/sample_freq) + 1;
	
	NumericMatrix CGDD(new_nrow , gdd.ncol());
	
	for(unsigned int orbit = 0; orbit < n_orbits; orbit++){
		for(unsigned int k = 0; k < new_nrow - 1; k++){
			CGDD(k, orbit) = cgdd_temp(sample_freq * k, orbit);
		}
		CGDD(new_nrow - 1, orbit) = sample_freq;
	}
	
	return CGDD;
}
// [[Rcpp::export]]
double emd_star(NumericVector X, NumericVector Y){
	NumericVector* shorter;
	NumericVector* longer;
	
	int l_short = std::min(X.length(), Y.length()) - 1;
	int l_long = std::max(X.length(), Y.length()) - 1;
	
	//int sample_freq_X = X[X.length() - 1];
	//int sample_freq_Y = Y[Y.length() - 1];
	
	if(X.length() < Y.length()){
		shorter = &X;
		longer = &Y;
	} else{
		shorter = &Y;
		longer = &X;
	}
	
	int sample_freq_shorter = (*shorter)[l_short];
	int sample_freq_longer = (*longer)[l_long];
	
	int sample_factor_shorter = lcm(sample_freq_shorter, sample_freq_longer)/sample_freq_shorter;
	int sample_factor_longer = lcm(sample_freq_shorter, sample_freq_longer)/sample_freq_longer;
	
	l_short = floor(l_short / sample_factor_shorter);
	l_long = floor(l_long / sample_factor_longer);
	
	
	double emd = l_long + l_short;
	double prop_emd = emd;
	int c = - l_short + 1;
	
	while(prop_emd <= emd){
		emd = prop_emd;
		prop_emd = 0;
		// sum shorter vector standing out left
		for(int i = 0; i < -c; i++){
			prop_emd += (*shorter)[i * sample_factor_shorter];
		}
		// sum both vectors when overlapping
		for(int i = std::max(0, -c); i < std::min(l_long-c, l_short); i++){
			prop_emd += abs((*shorter)[i * sample_factor_shorter]-
							(*longer)[i * sample_factor_longer+c]);
		}
		// sum longer vector not covered by shorter vector from left
		for(int i = 0; i < c; i++){
			prop_emd += (*longer)[i * sample_factor_longer];
		}
		// sum longer vector not covered by shorter vector from right
		for(int i = l_long-1; i > l_short+c - 1; i--){
			prop_emd += 1-(*longer)[i * sample_factor_longer];
		}
		
		// sum shorter vector standing out right
		for(int i = l_short-1; i > l_long-c - 1; i--){
			prop_emd += 1-(*shorter)[i * sample_factor_shorter];
		}
		
		c ++;
	}
	return emd;
}

// [[Rcpp::export]]
double NetEmd(NumericMatrix X, NumericMatrix Y){
	unsigned int n_orbits = X.ncol();
	double netemd = 0;
	
	for(unsigned int orbit = 0; orbit < n_orbits; orbit++){
		netemd =+ emd_star(X(_, orbit), Y(_, orbit));
	}
	
	return netemd/n_orbits;
}

