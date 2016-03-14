/*
 * Poisson model with offset
 */

data {
	int<lower=0> N;                  // number of individuals
	int nFert[N];           // number of fertilized eggs (successes)
	vector[N] nEggs;           // number of eggs counted (trials)
	vector[N] nSperm;                  // covariate of interest
}

transformed data {
	vector[N] log_nEggs;

	log_nEggs <- log(nEggs);
}

parameters {
	vector[2] beta;
}

model {
	nFert ~ poisson_log(log_nEggs + beta[1] + beta[2] * nSperm);
}