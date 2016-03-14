/*
 *  Simple beta-binomial example
 */

data {
  int<lower=0> N;                  // number of individuals
  int<lower=0> nFert[N];            // number of fertilized eggs (successes)
  int<lower=0> nEggs[N];            // number of eggs counted (trials)
  real nSperm[N];                // covariate of interest
}

parameters {
    real<lower = 0> phi;
    real alpha;
    real beta;
}

model {
    vector[N] lambda;
    vector[N] a;
    vector[N] b;

    // Priors
    phi ~ pareto(0.1,1.5);
    alpha ~ uniform(0,1);
    beta ~ uniform(0,1);

    // Likelihood
    for(i in 1:N) {
        lambda[i] <- inv_logit(alpha + beta * nSperm[i]);
        a[i] <- lambda[i] * phi;
        b[i] <- (1 - lambda[i]) * phi;
    }

    nFert ~ beta_binomial(nEggs, a, b);
}

