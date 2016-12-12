/*
 * Flexible Grouped Logistic Regression/ANOVA model using design
 * matrix to define model structure. 
 */

data {
   int<lower=1> N;      // Number of observations
   int<lower=1> P;      // Number of predictors
   int<lower=0> nT[N];  // Number of trials
   int<lower=0> nS[N];  // Number of successes
   matrix [N,P] X;      //  Predictor matrix
}

parameters {
   vector[P] beta;   //  Coefficients
}

model {
   beta  ~  normal(0,2);  // Prior for coefficients

   nS  ~  binomial_logit(nT, X * beta);  // Likelihood
}

generated quantities {
 	vector[N] log_lik;
 
 	for (i in 1:N) {
 		log_lik[i]  =  binomial_logit_log(nS[i], nT[i], X[i] * beta);
 	}
}
