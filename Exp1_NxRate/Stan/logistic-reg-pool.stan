/*
 * Simple Hierarchical Hierarchical Logistic Regression Model 
 * -- Linear predictor w/ nSperm
 *
 */

data {
	int<lower=0> N;                    // number of individuals in data set
	int<lower=0> nEggs[N];             // total number of eggs counted for each individual (trials)
	int<lower=0> nFert[N];             // number of fertilized eggs counted for each individual (successes)
	vector[N] nSperm;                  // covariate of interest
}

parameters {
vector<lower=0, upper=1>[2] theta;    // Coefficients for linear predictor
}

model {
	nFert  ~  binomial_logit(nEggs, theta[1] + theta[2] * nSperm);  // Likelihood
}

generated quantities {
	int<lower=0> nFert_rep[N];      // replications for existing items 
	for (n in 1:N) 
    	nFert_rep[n] <- binomial_rng(nEggs[n], inv_logit(theta[1] + theta[2] * nSperm[n])); 
}