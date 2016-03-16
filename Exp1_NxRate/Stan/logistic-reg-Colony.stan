/*
 *  Hierarchical Logistic Regression Model 
 *  -- Random intercept ~ Colony
 *  -- Linear predictor ~ nSperm
 *
 */

data {
	int<lower=0> N;                    // number of individuals in data set
	int<lower=0> nEggs[N];             // total number of eggs counted for each individual (trials)
	int<lower=0> nFert[N];             // number of fertilized eggs counted for each individual (successes)
	int<lower=0> Colony[N];            // Grouping variable for Colony
	vector[N] nSperm;                  // covariate of interest
}

parameters {
  vector[3] a;
  real theta;
  real<lower=0> sigma_a;
  real mu_a;
} 

transformed parameters {
  vector[N] y_hat;

  for (i in 1:N)
    y_hat[i]  <-  a[Colony[i]] + theta * nSperm[i];
}

model {
	mu_a         ~  normal(0, 1);                   // Hyperprior for among-intercept mean
	sigma_a      ~  cauchy(0,5);                    // Hyperprior for among-intercept variance
	a            ~  normal (mu_a, sigma_a);         // Prior
	nFert        ~  binomial_logit(nEggs, y_hat);  // Likelihood
}

generated quantities {
	int<lower=0> nFert_rep[N];      // Posterior draws for prediction
	for (i in 1:N) 
    	nFert_rep[i] <- binomial_rng(nEggs[i], inv_logit(a[Colony[i]] + theta * nSperm[i])); 
}