/*
 * Basic, flexible Regression/ANOVA model using design
 * matrix to define model structure. Basic design taken
 * from Murray Logan's JAGS models from his Stats Course
 * and vectorized as in STAN User Manual 2.5.0
 */

data {
   int<lower=1> N;	 // number of observations
   int<lower=1> P;	 // number of predictors
   matrix [N,P] X;	 // Predictor matrix
   vector[N] nEggs;	 // Number of eggs counted
   vector [N] y;	 // Response (number of fertilized eggs)
}

parameters {
   vector[P] beta;		 //  Coefficients
   real<lower=0> sigma;  //  var
}

model {
   y  ~  binomial_logit(nEggs, X * beta);   // Likelihood
}
