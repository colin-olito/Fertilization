/*
 * Simple Hierarchical Hierarchical Logistic Regression Model 
 * -- Linear predictor w/ nSperm
 *
 */

data {
	int<lower=0> N;                    // number of individuals in data set
	int<lower=0> nEggs[N];             // total number of eggs counted for each individual (trials)
	int<lower=0> nFert[N];             // number of fertilized eggs counted for each individual (successes)
	int<lower=0> Run[N];               // Grouping variable for Run
	vector[N] nSperm;                  // covariate of interest
}

parameters {
  vector[8] a;
  real<lower=0,upper=1> theta;
  real<lower=0,upper=100> sigma_a;
  real<lower=0,upper=100> sigma_y;
  real mu_a;
} 

transformed parameters {
  vector[N] y_hat;

  for (i in 1:N)
    y_hat[i]  <-  a[Run[i]] + theta * nSperm[i];
}

model {
	mu_a ~ normal(0, 1);                     // Hyperprior
	a ~ normal (mu_a, sigma_a);              // Prior
	nFert  ~  binomial_logit(nEggs, y_hat);  // Likelihood
}

generated quantities {
	int<lower=0> nFert_rep[N];      // Posterior draws for prediction
	for (i in 1:N) 
    	nFert_rep[i] <- binomial_rng(nEggs[i], inv_logit(a[Run[i]] + theta * nSperm[i])); 
}