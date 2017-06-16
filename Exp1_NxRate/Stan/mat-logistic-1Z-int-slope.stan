/*
 * Basic, flexible Mixed Effects Regression/ANOVA model 
 * using design matrix to define model structure. 
 */

data {
   int<lower=1> N;      // Number of observations
   int<lower=0> P;      // Number of predictors
   int<lower=0> K;      // Number of group-level predictors
   int<lower=0> nT[N];  // Number of trials
   int<lower=0> nS[N];  // Number of successes
   matrix [N,P] X;      //  Predictor matrix
   matrix[N,K] Z;   // Model matrix for run effects
}

parameters {
   vector [P] beta;            // Vector of fixed effect estimates
   vector [K] gamma;           // Vector of random effect estimates
   real<lower=0> sigma_gamma0;  // Prior for among-Run intercept variance
   real<lower=0> sigma_gamma1;  // Prior for among-Run slope variance
}

model {

   // Hyperpriors
   sigma_gamma0  ~  cauchy(0,5);
   sigma_gamma1  ~  cauchy(0,5);

   // Priors
   for (i in 1:K) {
      if(i <= 8)
         gamma[i] ~ normal(0, sigma_gamma0);
      if(i > 8)
         gamma[i] ~ normal(0, sigma_gamma1);
   }
   beta ~ normal(0,3);

   // Likelihood
   nS  ~  binomial_logit(nT, X*beta + Z*gamma);  // Likelihood
}

generated quantities {
   vector[N] log_lik;
 
   for (i in 1:N) {
      log_lik[i]  =  binomial_logit_log(nS[i], nT[i], X[i]*beta + Z[i]*gamma);
   }
}
