/*
 * Basic, flexible Mixed Effects Regression/ANOVA model 
 * using design matrix to define model structure.
 * Covariance matrices estimated. 
 */

data {
   int<lower=1> N;               // number of observations
   int<lower=1> P;               // number of main effect parameters (betas)
   int<lower=0> J;               // Number of groups (Plots)
   int<lower=1> K;               // number of group predictors (r.e. predictors; ncol(Z))
   int<lower=1,upper=J> grp[N];  // Group index for each individual
   matrix[N,P] X;                // Model matrix for main effects (by run)
   row_vector[K] Z[N];           // Group predictors
   vector[N] y;                  // Response variable
}

parameters {
   cholesky_factor_corr[K] Lp;   // Cholesky factor of Plot r.e. correlations
   vector<lower=0>[K] tauP;      // standard deviations of unconditional Plot r.e. dist
   vector[K] u[J];               // spherical Plot random effects
   vector [P] beta;              // Vector of fixed effect estimates
   real<lower=0> sigma_y;        // error scale
}

transformed parameters {
   matrix[K,K] corrs;
   corrs <- tcrossprod(Lp);      // for monitoring Plot correlations
}

 model {

   // Local variables
   vector[N] mu;
   matrix[K,K] LambdaP; 
   vector[K] gamma[J];

   // Hyperpriors
   tauP  ~  cauchy(0,2.5);
   Lp    ~  lkj_corr_cholesky(2);
   LambdaP  <-  diag_pre_multiply(tauP,Lp);
   for (i in 1:J) {
      u[i] ~ normal(0,1);
      gamma[i] <- LambdaP * u[i];
   }

   // Priors
   beta ~ normal(0,3);
   sigma_y  ~  cauchy(0,5);

   // Likelihood
  for (i in 1:N) {
      mu[i]  <-  X[i] * beta + Z[i] * gamma[grp[i]];
  }
   y ~ normal(mu, sigma_y);
}

generated quantities {
   vector[N] log_lik;
   vector[N] mu;
   vector[K] gamma[J];

   for (i in 1:N) {
      log_lik[i]  <-  normal_log(y[i], X[i] * beta + Z[i] * gamma[grp[i]], sigma_y);
   }
}