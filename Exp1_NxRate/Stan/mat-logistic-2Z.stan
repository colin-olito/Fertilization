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
   matrix[N,K] Z0;   // Model matrix for run-specific intercepts
   matrix[N,K] Z1;   // Model matrix for run-specific slopes
}

parameters {
   vector [P] beta;            // Vector of fixed effect estimates
   vector [K] gamma0;           // Vector of random effect estimates
   vector [K] gamma1;           // Vector of random effect estimates
   real<lower=0> sigma_gamma0;  // Prior for among-Run intercept variance
   real<lower=0> sigma_gamma1;  // Prior for among-Run slope variance
}

transformed parameters {
   vector[N] mu;
   mu  =  X*beta + Z0*gamma0 + Z1*gamma1;
}

model {
   // Hyperpriors
   sigma_gamma0  ~  cauchy(0,5);
   sigma_gamma1  ~  cauchy(0,5);

   // Priors
   gamma0 ~ normal(0, sigma_gamma0);
   gamma1 ~ normal(0, sigma_gamma1);
   beta ~ normal(0,3);

   // Likelihood
   nS  ~  binomial_logit(nT, mu);  // Likelihood
}

generated quantities {
   vector[N] log_lik;
 
   for (i in 1:N) {
      log_lik[i]  =  binomial_logit_log(nS[i], nT[i], X[i]*beta + Z0[i]*gamma0 + Z1[i]*gamma1);
   }
}
