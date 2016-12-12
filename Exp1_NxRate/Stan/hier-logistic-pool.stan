/*
 * Simple Hierarchical Hierarchical Logistic Regression Model 
 * -- single mean
 * -- complete pooling 
 *
 */
data {
	int<lower=0> N;                    // number of individuals in data set
	int<lower=0> nEggs[N];             // total number of eggs counted for each individual (trials)
	int<lower=0> nFert[N];             // number of fertilized eggs counted for each individual (successes)
	vector[N] nSperm;                  // covariate of interest
}
parameters {
	real<lower=0,upper=1> theta;    // probability of fertilization for each individual
}
model {
	nFert  ~  binomial(nEggs, theta);  // Likelihood
}
