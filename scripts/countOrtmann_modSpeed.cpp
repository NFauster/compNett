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

double sd(std::vector<double> x, int* nonZeroElements){
    
    double x_px;    // holds sum(X*p(X))
	double x2_px;   // holds sum(X^2*p(X))
	*nonZeroElements = 0;
	
	for(int i = 0; i< x.size(); i++){
		x_px += (i*x[i]);
		x2_px += (pow(i,2)*x[i]);
		
		if(x[i] > 0){
			(*nonZeroElements) ++;
		}
	}
	
	return sqrt(x2_px-x_px);
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
			int temp_nn4 = 0;
			for(int v:crt_N){
				std::vector<int> v_N = neighbourhood[v - 1][0];
				
				for(int n:v_N){
					temp_nn4 += deg[n - 1];
				}
				temp_nn4 -= deg[v - 1];
			}
	
			nn[(*u - 1) * 15 + 0] = deg[*u - 1];
			nn[(*u - 1) * 15 + 1] = (deg_N_sum - deg[*u - 1]);
			nn[(*u - 1) * 15 + 2] = chooseC(deg[*u - 1], 2);
			nn[(*u - 1) * 15 + 3] = k3[*u - 1];
			nn[(*u - 1) * 15 + 4] = temp_nn4 - deg[*u - 1] * (deg[*u - 1] - 1) - 2*k3[*u - 1];
			
			nn[(*u - 1) * 15 + 5] = (deg[*u - 1] - 1)*(deg_N_sum - crt_N.size()) - 
									accumulate(k3_edge[*u - 1].begin(),
											   k3_edge[*u - 1].end(), 0);
			
			nn[(*u - 1) * 15 + 6] = 0;
			for(int v: crt_N){
				nn[(*u - 1) * 15 + 6] += chooseC(deg[v - 1] - 1, 2);
			}

			nn[(*u - 1) * 15 + 7] = chooseC(deg[*u - 1], 3);
			nn[(*u - 1) * 15 + 8] = c4[*u - 1];
			
			nn[(*u - 1) * 15 + 9] = - accumulate(k3_edge[*u - 1].begin(),
										  k3_edge[*u - 1].end(), 0);
			for(int n:crt_N){
				nn[(*u - 1) * 15 + 9] += k3[n - 1];
			}
	
			for(int i = 0; i < completing_triangle[*u - 1].size(); i += 2){
				nn[(*u - 1) * 15 + 10] += deg[completing_triangle[*u - 1][i] - 1] + 
						deg[completing_triangle[*u - 1][i+1] - 1] - 4;
			}
			
			nn[(*u - 1) * 15 + 11] = k3[*u - 1]*(deg[*u - 1] - 2);
	
			nn[(*u - 1) * 15 + 12] = -k3[*u - 1];
			for(int i = 0; i < completing_triangle[*u - 1].size(); i += 2){
				nn[(*u - 1) * 15 + 12] += k3_edge[completing_triangle[*u - 1][i] - 1]
									[distance(neighbourhood[completing_triangle[*u - 1][i]-1][0].begin(),
												lower_bound(neighbourhood[completing_triangle[*u - 1][i]-1][0].begin(), 
															neighbourhood[completing_triangle[*u - 1][i]-1][0].end(), 
															completing_triangle[*u - 1][i + 1])
											)
									];
			}
	
			for(int t:k3_edge[*u - 1]){
				nn[(*u - 1) * 15 + 13] += chooseC(t,2);
			}
	
			nn[(*u - 1) * 15 + 14] = k4[*u - 1];
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
			for(int i=0; i < nn_arma.size(); i++){
				nn_arma(i) = nn[(*u - 1) * 15 + i];
			}
			
			arma::mat LEM = {
					 {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
					 {0, 1, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
					 {0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
					 {0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
					 {0, 0, 0, 0, 1, 0, 0, 0, 2, 2, 1, 0, 4, 2, 6}, 
					 {0, 0, 0, 0, 0, 1, 0, 0, 2, 0, 1, 2, 2, 4, 6}, 
					 {0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 2, 1, 3}, 
					 {0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1}, 
					 {0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 3}, 
					 {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 2, 0, 3}, 
					 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 2, 2, 6}, 
					 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 2, 3}, 
					 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 3}, 
					 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 3}, 
					 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}};
					 
			arma::vec ni_arma = solve(LEM,  nn_arma);
			
			for(int i=0; i < ni.ncol(); i++){
				ni(*u - 1, i) = ni_arma[i];
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
	std::vector<int> nn(n_nodes * 15);
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
	
	/*for(int i = 0; i <100; i++){
		Rprintf("%i\n", nn[i]);
		if(i%5 == 4) Rprintf("\n");
	}*/
	
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
	//ni ... Matrix of graphlet degree vectors
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
List CGDD_wo(IntegerMatrix ni){
	//GDD
	unsigned int n_orbits = ni.ncol();
	unsigned int n_nodes = ni.nrow();
	
	std::vector<double> temp_gdd;
	List CGDD(n_orbits);
	int nonZeroElements;
	
	for(unsigned int orbit = 0; orbit < n_orbits; orbit++){
		temp_gdd.clear();
		temp_gdd.resize(max(ni(_, orbit)) + 1);
		
		for(unsigned int n = 0; n < n_nodes; n++){
			temp_gdd[ni(n, orbit)] += 1.0/n_nodes;
		}
		
		double gdd_sd = sd(temp_gdd, &nonZeroElements);
		if(gdd_sd == 0) gdd_sd = 1.0;
		
		
		CGDD[orbit] = NumericMatrix(2,nonZeroElements);
		int count_index = 0;
		for(int i = 0; i < temp_gdd.size(); i++)
		{
			if(temp_gdd[i] > 0.0 & count_index == 0) {
				as<NumericMatrix>(CGDD[orbit])(0,count_index) = i/gdd_sd;
				as<NumericMatrix>(CGDD[orbit])(1,count_index) = temp_gdd[i];
				count_index ++ ;
			} else if(temp_gdd[i] > 0.0) {
				as<NumericMatrix>(CGDD[orbit])(0,count_index) = i/gdd_sd;
				as<NumericMatrix>(CGDD[orbit])(1,count_index) = temp_gdd[i] + as<NumericMatrix>(CGDD[orbit])(1,count_index-1);
				count_index ++ ;
			}
		}
	}
	
	return CGDD;
}

// [[Rcpp::export]]
double EMD(NumericMatrix X, NumericMatrix Y, double c){
	double emd = 0;
	int index_x = 0;
	int index_y = 0;
	
	// shift distribution
	NumericMatrix X_shift = clone(X);
	X_shift(0,_) = X(0,_) + c;
	
	NumericMatrix* first_start;
	NumericMatrix* later_start;
	
	int* index_later;
	
	// determine which distribution has the first step
	if(X_shift(0,0) < Y(0,0)){
		first_start = &X_shift;
		later_start = &Y;
		
		index_later = &index_y;
	} else{
		first_start = &Y;
		later_start = &X_shift;
		
		index_later = &index_x;
	}
	
	double prev_xCoord_later = 0.0;
	double prev_value_later = 0.0;
	
	double x_upper;
	double x_lower;
	double y_diff;
	
	for(int index_first = 0; index_first < ((*first_start).ncol() - 1); 
				index_first ++){
		
		if((*first_start)(0, index_first + 1) < (*later_start)(0, 0)){
			// integrate the earlier starting distribution before the later one
			// starts
			
			x_upper = (*first_start)(0, index_first + 1);
			x_lower = (*first_start)(0, index_first);
			y_diff = (*first_start)(1, index_first);
			
			emd += (x_upper - x_lower) * y_diff;
			
		} else if(((*first_start)(0, index_first) >= 
					(*later_start)(0, (*later_start).ncol() - 1))){
			// integrate the difference of the distributions after the 
			// later starting distributions stops
			
			x_upper = (*first_start)(0, index_first + 1);
			x_lower = (*first_start)(0, index_first);
			y_diff = 1 - (*first_start)(1, index_first);
			
			emd += (x_upper - x_lower) * y_diff;
		} else {
			// integrate the difference of the distributions in the 
			// overlapping section
			
			while((*index_later) < (*later_start).ncol() & 
					((*later_start)(0, *index_later) <= 
					(*first_start)(0, index_first + 1))){
				if((*index_later) > 0){
					prev_xCoord_later = (*later_start)(0, (*index_later) - 1);
					prev_value_later = (*later_start)(1, (*index_later) - 1);
				}
			
				x_upper = min(max((*later_start)(0, *index_later), 
								(*first_start)(0, index_first)),
							(*first_start)(0, index_first + 1));
			
				x_lower = max(prev_xCoord_later, (*first_start)(0, index_first));
			
				y_diff = abs((*first_start)(1, index_first) - prev_value_later);
			
				emd += (x_upper - x_lower) * y_diff;
			
				(*index_later) ++;
			}
		
		
			x_upper = (*first_start)(0, index_first + 1);
			x_lower = max((*later_start)(0, (*index_later - 1)), 
						(*first_start)(0, index_first));	
			y_diff = abs((*first_start)(1, index_first) - 
							(*later_start)(1, (*index_later - 1)));
			emd += (x_upper - x_lower) * y_diff;
		}
	}
	
	// integrate the difference of the distributions after the 
	// earlier starting distributions stops
	while((*index_later) < (*later_start).ncol()){
		if((*index_later) > 0){
			prev_xCoord_later = (*later_start)(0, (*index_later) - 1);
			prev_value_later = (*later_start)(1, (*index_later) - 1);
		}
		x_upper = (*later_start)(0, *index_later);
		x_lower = max((*first_start)(0, (*first_start).ncol() - 1),
							prev_xCoord_later);
		
		y_diff = 1 - prev_value_later;
		
		emd += (x_upper - x_lower) * y_diff;
	
		(*index_later) ++;
	}

	return emd;
}


// [[Rcpp::export]]
double emd_star(NumericMatrix X, NumericMatrix Y){
	double emd = X(0, X.ncol()-1) + Y(0, Y.ncol()-1);
	double prop_emd = emd;
	
	int l_long = std::max(X(0, X.ncol()-1), Y(0, Y.ncol()-1));
	double c = - l_long + 1;
	int counter = 0;
	bool check_run = true;
	
	while(check_run){
		prop_emd = EMD(X, Y, c);
		c += .5;
		
		if(prop_emd <= emd){
			emd = prop_emd;
			counter = 0;
			check_run = true;
		} else if (counter <= 200){
			counter ++;
			check_run = true;
		} else {
			check_run = false;
		}
	}
	
	return emd;
}

// [[Rcpp::export]]
double NetEmd(List X, List Y){
	unsigned int n_orbits = X.length();
	double netemd = 0;
	
	for(unsigned int orbit = 0; orbit < n_orbits; orbit++){
		netemd += emd_star(X[orbit], Y[orbit]);
	}
	
	return netemd/n_orbits;
}

