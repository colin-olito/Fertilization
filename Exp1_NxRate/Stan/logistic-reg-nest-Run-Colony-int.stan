/*
 *  Hierarchical Logistic Regression Model 
 *	-- Nested random intercepts: Run w/in Colony
 *  -- Linear predictor ~ nSperm
 *
 */

data {
	int<lower=0> N;                    // number of individuals in data set
	int<lower=0> Nr;                   // number of runs in data set
	int<lower=0> Nc;                   // number of colonies in data set
	int<lower=0> nEggs[N];             // total number of eggs counted for each individual (trials)
	int<lower=0> nFert[N];             // number of fertilized eggs counted for each individual (successes)
	int<lower=0> Run[N];               // Grouping variable for Run
	int<lower=0> Colony[Nr];           // Grouping variable for Colony
	vector[N] nSperm;                  // covariate of interest
}

parameters {
	vector[Nr] a;			  //  Run-specific intercepts 
	vector[Nc] b;			  //  Colony-specific intercepts 
	real theta;				  //  overall slope ~ nSperm
	real mu_col;			  //  among-colony mean
	real<lower=0> sigma_col;  //  among-colony variance
	real mu_run;			  //  among-run mean
	real<lower=0> sigma_run;  //  among-run variance
} 

model {
	mu_col       ~  normal(0,5);				// Hyper-prior for Among-Colony mean
	sigma_col    ~  cauchy(0,5);				// Hyper-prior for among-colony sd
	mu_run       ~  normal(mu_col, sigma_col);  // Hyperprior for among-run mean
	sigma_run    ~  cauchy(0,5);                // Hyperprior for among-run variance
	b            ~  normal(mu_col,sigma_col);	// Prior for among-Colony variation
	a            ~  normal(mu_run, sigma_run);  // Prior: among-Run variation
	for(i in 1:N) {
		nFert[i]  ~  binomial_logit(nEggs[i], a[Run[i]] + theta * nSperm[i]);	// Likelihood
	}
	for (j in 1:Nr) {
		a[j]   ~  normal(b[Colony[j]], sigma_run);
	}
}
