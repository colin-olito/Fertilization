/*
 * Random intercept heirarchical model using
 * matrix notation. Modeled on Logan Murray's
 * advice from Monash stats course
 *
 */

data {
   int<lower=1> N;               // number of observations
   int<lower=1> P;               // number of main effect parameters (thetas)
   int<lower=1> K;               // number of runs (groups)
   matrix[N,P] X;               // Model matrix for main effects (by run)
   matrix[N,K] Z;               // Model matrix for run effects
   vector[N] y;                  // Response variable
}

parameters {
   vector [P] theta;                            // Vector of main effect estimates
   vector [K] gamma;                          // Vector of run effects
   real<lower=0,upper=10> sigma_gamma;   // Hyperprior for among-run variance
   real<lower=0,upper=10> sigma_y;         // error scale
}

model {

   // Priors
   gamma ~ normal(0, sigma_gamma);
   theta ~ normal(0,1);

   // Likelihood
   y ~ normal(X*theta + Z*gamma, sigma_y);
}