/*
 * Basic, flexible Mixed Effects Beta-binomial 
 * Regression/ANOVA model using design matrices
 * to define model structure. 
 */

data {
   int<lower=1> N;      // Number of observations
   int<lower=0> P;      // Number of predictors
   int<lower=0> K;      // Number of group-level predictors
   int<lower=0> nT[N];  // Number of trials
   int<lower=0> nS[N];  // Number of successes
   matrix [N,P] X;      //  Predictor matrix
   matrix[N,K] Z;       // Model matrix for random effects
}

parameters {
   vector [P] beta;            // Vector of fixed effect estimates
   vector [K] gamma;           // Vector of random effect estimates
   real<lower=0> sigma_gamma;  // Prior for among-Run variance
   real<lower=0.1> kappa;      // prior count
}

transformed parameters {
   vector[N] a;             // Prior success count (alpha parameter in BBin dist.)
   vector[N] b;             // Prior failure count (beta parameter in BBin dist.)
   vector[N] lambda;        // mean chance of success

   for (n in 1:N)
      lambda[n] = inv_logit(X[n]*beta + Z[n]*gamma); //using logit link
   
   a  = lambda * kappa;
   b  = (1 - lambda) * kappa;
}


model {
   // Hyperpriors
   // Implicit Unif[0.1,Inf] hyperprior on kappa (see parameters{} block)
   sigma_gamma  ~  cauchy(0,5);

   // Priors
   gamma  ~  normal(0, sigma_gamma);
   beta   ~  normal(0,3);

   // Likelihood
   nS  ~  beta_binomial(nT, a, b);  
}

generated quantities {
   vector[N] log_lik;
 
   for (i in 1:N) {
      log_lik[i]  =  beta_binomial_log(nS[i], nT[i], a[i], b[i]);
   }
}
