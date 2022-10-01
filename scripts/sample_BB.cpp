#include <Rcpp.h>

// [[Rcpp::depends(Rcpp)]]

using namespace Rcpp;
using namespace std;

NumericVector unif(int n, int min, int max){
    
    Function runif("runif"); 
	NumericVector result = runif(n, min, max);
	
	return result;
}

NumericVector exp(int n, int rate){
    
    Function rexp("rexp"); 
	NumericVector result = rexp(n, rate);
	
	return result;
}

NumericVector beta(int n, int shape1, int shape2, int ncp){
    
    Function rbeta("rbeta"); 
	NumericVector result = rbeta(n, shape1, shape2, Named("ncp") = ncp);
	
	return result;
}

NumericVector gamma(int n, int shape, int scale){
    
    Function rgamma("rgamma"); 
	NumericVector result = rgamma(n, shape, Named("scale") = scale);
	
	return result;
}

IntegerVector sample(IntegerVector x, int size, NumericVector prob){
    
    Function sample("sample"); 
	IntegerVector result = sample(x, size, Named("prob") = prob);
	
	return result;
}

// [[Rcpp::export]]
IntegerVector sample_BB(int n_nodes, int m, int mode,
						int min = 0, int max = 1, //parameters for mode 1
						int rate = 1, // parameter for mode 2
						int shape1 = 2, int shape2 = 5, int ncp = 0, //mode 3
						int shape = 5, int scale = 1 // mode 4
						){
	std::vector<int> edge_list;
	std::vector<int> degree(n_nodes,0);		//BB call it connectivity k
	std::vector<double> temp_prob(n_nodes,0);		
	IntegerVector sampled_nodes;
	NumericVector fitness(n_nodes);	
	int crt_nodes = 0;
	
	//seed graph: cycle graph, 5 nodes
	vector<int> to_insert{1,2,2,3,3,4,4,5,5,1};
	edge_list.insert(edge_list.end(), to_insert.begin(), to_insert.end());
	
	//sample fitness
	if(mode == 0){
		// scale-free: all fitness equal
		for(int i = 0; i < n_nodes; i++){
			fitness[i] = 2;
		}
	} else if(mode == 1){
		// uniform distribution
		fitness = unif(n_nodes, min, max);
	} else if(mode == 2){
		// exponential distribution
		fitness = exp(n_nodes, rate);
	} else if(mode == 3){
		// beta distribution
		fitness = beta(n_nodes, shape1, shape2, ncp);
	} else if(mode == 4){
		// gamma distribution
		fitness = gamma(n_nodes, shape, scale);
	}
	
	
	for(int i = 1; i <= 5; i++){
		degree[i - 1] = 2;
		crt_nodes ++;
		
		temp_prob[i - 1] = fitness[i - 1] * degree[i - 1];
	}
	
	
	// generate graph
	while(crt_nodes < n_nodes){
		std:vector<double> to_pass(temp_prob.begin(),
									(temp_prob.begin() + crt_nodes));
									
		sampled_nodes = sample(seq(1, crt_nodes), std::min(m, crt_nodes), wrap(to_pass));
		
		for(int s: sampled_nodes){
			edge_list.push_back(crt_nodes + 1);
			edge_list.push_back(s);
			
			degree[s - 1] ++;
			
			temp_prob[s - 1] = fitness[s - 1] * degree[s - 1];
		}
		degree[crt_nodes] = m;
		temp_prob[crt_nodes] = fitness[crt_nodes] * degree[crt_nodes];
		crt_nodes ++;
	}
	
	return wrap(edge_list);
}