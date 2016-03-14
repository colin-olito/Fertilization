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
vector<lower=0, upper=1>[2] theta;    // probability of fertilization for each individual
}
model {
	nFert  ~  binomial(nEggs, theta[1] + theta[2] * nSperm);  // Likelihood
}