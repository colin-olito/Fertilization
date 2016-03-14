/*
 * Simple Hierarchical Hierarchical Logistic Regression Model 
 * -- Partial POOLING
 * -- Beta-Binomial parameterization
 *
 */

data {
	int<lower=0> N;                    // number of individuals in data set
	int<lower=0> nEggs[N];             // total number of eggs counted for each individual (trials)
	int<lower=0> nFert[N];             // number of fertilized eggs counted for each individual (successes)
	vector[N] nSperm;                  // covariate of interest
}

parameters {
  real<lower=0, upper=1> phi;         // population chance of successful fertilization
  real<lower=1> kappa;                // population concentration
  vector<lower=0, upper=1>[N] theta;    // probability of fertilization for each individual
}

model {
	kappa ~ pareto(1, 1.5);                        // hyperprior
	theta ~ beta(phi * kappa, (1 - phi) * kappa);  // prior
  	nFert  ~  binomial(nEggs, theta);  // Likelihood
}