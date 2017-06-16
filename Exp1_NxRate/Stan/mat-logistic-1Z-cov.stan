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

transformed data {
  real min_y;   // minimum successes
  real max_y;   // maximum successes
  real mean_y;  // sample mean successes
  real sd_y;    // sample std dev successes

  min_y   =  min(nS);
  max_y   =  max(nS);
  mean_y  =  mean(to_vector(nS));
  sd_y    =  sd(to_vector(nS));
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

  corrs = tcrossprod(L_run);       // For monitoring correlations
  Lambda_run  =  diag_pre_multiply(tau_run,L_run);
  for (i in 1:J)
    gamma[i] = Lambda_run * u[i];  // Random-effect coefficients
  
  for (i in 1:N)                    // Loop to vectorize likelihood in model block
    mu[i]  =  X[i] * beta + Z[i] * gamma[grp[i]];
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

   // Containers
   vector[N] log_lik;            // log-likelihood for LOO calculations
   vector[N] y_rep;              // vector for replicated data for posterior predictive checks
   real<lower=0> min_y_rep;      // posterior predictive min replicated successes
   real<lower=0> max_y_rep;      // posterior predictive max replicated successes
   real<lower=0> mean_y_rep;     // posterior predictive sample mean replicated successes
   real<lower=0> sd_y_rep;       // posterior predictive sample std dev replicated successes
   int<lower=0, upper=1> p_min;  // posterior predictive p-values
   int<lower=0, upper=1> p_max;  // ...
   int<lower=0, upper=1> p_mean; // ...
   int<lower=0, upper=1> p_sd;   // ...
   vector[N] y_hat;              // predicted values
   vector[N] X2_data;            // Chi-squared discrepancy between real data and prediction line
   vector[N] X2_rep;             // Chi-squared discrepancy between simulated data and prediction line
   real<lower=0> fit_data;       // ...
   real<lower=0> fit_rep;        // ...
 
   for (i in 1:N) {
      log_lik[i]  =  binomial_logit_lpmf(nS[i] | nT[i], mu[i]);
      y_rep[i]    =  binomial_rng(nT[i], inv_logit(mu[i]));      
      y_hat[i]    =  inv_logit(mu[i]);
      X2_data[i]  =  ((y_hat[i] - (to_vector(nS)[i]/to_vector(nT)[i]))^2)/(to_vector(nS)[i]/to_vector(nT)[i]);
      X2_rep[i]   =  ((y_hat[i] - (y_rep[i]/to_vector(nT)[i]))^2)/(y_rep[i]/to_vector(nT)[i]);
   }

   min_y_rep   =  min(y_rep);
   max_y_rep   =  max(y_rep);
   mean_y_rep  =  mean(to_vector(y_rep));
   sd_y_rep    =  sd(to_vector(y_rep));

   p_min   =  (min_y_rep >= min_y);
   p_max   =  (max_y_rep >= max_y);
   p_mean  =  (mean_y_rep >= mean_y);
   p_sd    =  (sd_y_rep >= sd_y);

   fit_data  =  sum(X2_data);
   fit_rep   =  sum(X2_rep);

}
