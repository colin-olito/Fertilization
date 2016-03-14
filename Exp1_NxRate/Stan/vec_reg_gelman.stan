/*
 * Vectorized regression model, as written by
 * the STAN developers in their User's Manual
 *
 *
 */

data {

   int<lower=0> N;             // Number of observations
   int<lower=0> P;             // Number of predictors
   matrix[N,P] X;              // predictor matrix
   vector[N] y;                // response variable
}

parameters {

   real alpha;                 // intercept
   vector[P] beta;             // coefficients for predictors
   real<lower=0> sigma;        // error scale
}

model {
   beta ~ normal(0,100);       // Prior for regression coefficients
   sigma ~ normal(0,100);      // Prior for residual error
   y ~ normal(X*beta, sigma);  // likelihood
}


