/*
 * Basic, flexible Mixed Effects Exponential Regression/ANOVA model 
 * using design matrix to define model structure. 
 */

data {
   int<lower=1> N;      // Number of observations
   int<lower=1> P;      // Number of predictors
   int<lower=1> K;      // Number of group-level predictors
   vector[N]    Y;      // Response variable of trials
   matrix [N,P] X;      //  Predictor matrix
   matrix [N,K] Z;      // Model matrix for random effects
}

transformed data {
  real min_y;   // minimum response
  real max_y;   // maximum response
  real mean_y;  // sample mean response
  real sd_y;    // sample std dev response

  min_y   =  min(Y);
  max_y   =  max(Y);
  mean_y  =  mean(to_vector(Y));
  sd_y    =  sd(to_vector(Y));
}

parameters {
   vector [P] beta;            // Vector of fixed effect estimates
   vector [K] gamma;           // Vector of random effect estimates
   real<lower=0> sigma_gamma;  // Prior for among-Run variance
}

transformed parameters {
  vector[N] mu;

  for(i in 1:N){
     mu[i]  =  (X[i]*beta + Z[i]*gamma)^-1;
  } 
}

model {

   // Hyperpriors
   sigma_gamma  ~  cauchy(0,100);
//   sigma_y      ~  cauchy(0,100);

   // Priors
   gamma ~ normal(0, sigma_gamma);
   beta  ~ normal(mean_y,100);

   // Likelihood
   Y  ~  exponential(mu);  
}

generated quantities {

   // Containers
   vector[N] log_lik;            // log-likelihood for LOO calculations
   vector[N] y_rep;              // vector for replicated data for posterior predictive checks
   real min_y_rep;      // posterior predictive min replicated successes
   real max_y_rep;      // posterior predictive max replicated successes
   real mean_y_rep;     // posterior predictive sample mean replicated successes
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
      log_lik[i]  =  exponential_lpdf(Y[i] | mu[i]);
      y_rep[i]    =  exponential_rng(mu[i]);      
      y_hat[i]    =  mu[i]^-1;
      X2_data[i]  =  ((y_hat[i] - to_vector(Y)[i])^2) / fabs(to_vector(Y)[i]);
      X2_rep[i]   =  ((y_hat[i] - y_rep[i])^2) / fabs(y_rep[i]);
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
