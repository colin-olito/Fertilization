/*
 * Random Effects Logistic Regression/ANOVA model 
 * -- all predictors 'random', 
 * -- using cellmean specification
 */

data {
   int<lower=1> N;      // Number of observations
   int<lower=1> K0;     // Number of group-level predictors
   int<lower=1> K1;     // Number of group-level predictors
   int<lower=0> nT[N];  // Number of trials
   int<lower=0> nS[N];  // Number of successes
   matrix[N,K0] Z0;      // Model matrix for run intercepts
   matrix[N,K1] Z1;      // Model matrix for run slopes
}

parameters {
   vector [K0] gamma0;           // Vector of random intercept estimates
   vector [K1] gamma1;           // Vector of random slopes estimates
   real mu_gamma0;              // Prior for among-Run mean intercept
   real mu_gamma1;              // Prior for among-Run mean slope
   real<lower=0> sigma_gamma0;  // Prior for among-Run variance in intercepts
   real<lower=0> sigma_gamma1;  // Prior for among-Run variance in slopes
}

model {

   // Hyperpriors
   mu_gamma0  ~  normal(0,2);
   mu_gamma1  ~  normal(0,2);
   sigma_gamma0  ~  cauchy(0,5);
   sigma_gamma1  ~  cauchy(0,5);

   // Priors
   gamma0 ~ normal(mu_gamma0, sigma_gamma0);
   gamma1 ~ normal(mu_gamma1, sigma_gamma1);

   // Likelihood
   nS  ~  binomial_logit(nT, Z0*gamma0 + Z1*gamma1);  // Likelihood
}

generated quantities {
   vector[N] log_lik;
 
   for (i in 1:N) {
      log_lik[i]  =  binomial_logit_lpmf(nS[i] | nT[i], Z0[i]*gamma0 + Z1[i]*gamma1);
   }
}
