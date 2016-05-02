/*
 * Basic, flexible logistic ME Regression model 
 * using design matrices to define model structure.
 * Covariance matrices estimated. 
 */

data {
   int<lower=1> N;               // Number of observations
   int<lower=0> P;               // Number of predictors
   int<lower=0> J;               // Number of groups (Runs)
   int<lower=0> K;               // Number of group-level predictors
   int<lower=1,upper=J> grp[N];  // Group index for each observation
   int<lower=0> nT[N];           // Number of trials
   int<lower=0> nS[N];           // Number of successes
   matrix [N,P] X;               // Predictor matrix
   row_vector[K] Z[N];           // Group predictor matrix (row_vector for efficiency)
}

parameters {
   cholesky_factor_corr[K] L_run;   // Cholesky factor of Run r.e. correlations
   vector<lower=0>[K] tau_run;      // standard deviations of unconditional Run r.e. dist
   vector[K] u[J];                  // spherical Run random effects
   vector [P] beta;                 // Vector of fixed effect estimates
   real<lower=0> sigma_y;           // error scale
}

transformed parameters {
  matrix[K,K] corrs;              
  vector[N] mu;                   
  matrix[K,K] Lambda_run; 
  vector[K] gamma[J];
//  vector[K] gamma;

  corrs <- tcrossprod(L_run);       // For monitoring correlations
  Lambda_run  <-  diag_pre_multiply(tau_run,L_run);
  for (i in 1:J)
    gamma[i] <- Lambda_run * u[i];  // Random-effect coefficients
  
  for (i in 1:N)                    // Loop to vectorize likelihood in model block
    mu[i]  <-  X[i] * beta + Z[i] * gamma[grp[i]];
}

model {
  tau_run  ~  cauchy(0,5);          // Priors for covariance matrix
  L_run    ~  lkj_corr_cholesky(2);

  for (i in 1:J)                    // Priors on spherical Random effects
    u[i] ~ normal(0,1);

  beta     ~  normal(0,5);          // Priors on fixed effects
  sigma_y  ~  cauchy(0,5);          // Prior on error scale

  nS  ~  binomial_logit(nT, mu);    // Likelihood
}

generated quantities {
  vector[N] log_lik;  // Log-likelihood

  for (i in 1:N) {
    log_lik[i]  <-  binomial_logit_log(nS[i], nT[i], mu[i]);
  }
}