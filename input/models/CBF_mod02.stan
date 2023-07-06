data {
  int<lower=1> n;           // number of observations
  int<lower=1> n_population;   // number of population
  int<lower=1, upper=n_population> population_id[n];  // population ID for each observation
  vector[n] y;              // size data
}
parameters {
  
  // normal parameters
  vector<lower = 0.0, upper = 5.0>[n_population] ln_mu;       // mean size for each population
  vector<lower= -4.0,upper=2.0>[n_population]  ln_cv;     // log-population coefficient of variation
  
	// log-normal parameters
  vector<lower= -1.0,upper=2.0>[n_population]  ln_meanlog;     // log of meanlog
  vector<lower= -10,upper=4.0>[n_population]  ln_sdlog;  // log of sdlog
}
model {
  
  // Values for normal distribution
  real sigma;
  real mu;
  
  // Values for Lognormal distribution
  real sdlog;
  real meanlog;
  
  // Normal priors
  ln_mu ~ normal(0, 10);
  ln_cv ~ cauchy(0, 1);
  
  // Lognormal priors
  ln_meanlog ~ normal(0, 10);
  ln_sdlog ~ cauchy(0, 0.1);

  
  // Likelihood for each observation
  for (i in 1:n) {
    mu     = exp(ln_mu[population_id[i]]);
    sigma  = mu * exp(ln_cv[population_id[i]]);
    y[i] ~ normal(mu, sigma);
    
    meanlog     = exp(ln_meanlog[population_id[i]]);
    sdlog  = mu * exp(ln_sdlog[population_id[i]]);
    y[i] ~ lognormal(meanlog, sdlog);

  }
}

generated quantities {
  vector<lower=0>[n_population] sigma; 
  vector<lower=0>[n_population] mu;  
  vector<lower=0>[n_population] cv;  
    vector<lower=0>[n_population] meanlog;  
  vector<lower=0>[n_population] sdlog;  

  for (s in 1:n_population) {
    mu[s]     = exp(ln_mu[s]);
    cv[s] =  exp(ln_cv[s]);
    sigma[s]  = mu[s]*cv[s];
    meanlog[s] = exp(ln_meanlog[s]);
    sdlog[s] = exp(ln_sdlog[s]);
  }
}
