data {
  int<lower=1>         B;       // number of length bins
  int<lower=1>         N;       // total number of fish
  real<lower = 0.0>   b_upr[B];
  real<lower = 0.0>   b_lwr[B];
  int<lower=1>         n[B];    // number of fish in the bin
  real<lower = 0.0> low_bound; // absolute lowest cutoff
}

parameters {
  real               mu;     // population meanlog
  real<lower=0.001>  sigma;  // population sdlog
  real<lower=0.001>  theta0; // variance parameter
}

model {
  real nconst; // normalising constant (for RLS data this is prob below 2.5cm)
  real p_upr; // cumulative probability of upper size bin boundary
  real p_lwr; // cumulative probability of lower size bin boundary
  real p;     // probability of being in bin
  real theta_i;
  real ln_pr;
  
  theta0 ~ exponential(0.01);
  sigma  ~ exponential(0.02);
  mu     ~ normal(10.0, 5.0);
  
  nconst = 1.0 - normal_cdf(low_bound, mu, sigma); 
  
  // for each size bin i

  target += lgamma(theta0) + lgamma(N+1) - lgamma(N+theta0); 

  for (i in 1:B) { 
    p_upr = normal_cdf(b_upr[i], mu, sigma);
    p_lwr = normal_cdf(b_lwr[i], mu, sigma);
    p = (p_upr - p_lwr)/ nconst;
    
    theta_i = theta0*p;
    
    target += lgamma(n[i] + theta_i) - lgamma(theta_i) - lgamma(n[i]+1);
  }
}

// generated quantities {
//     real rho;
//     
//    rho = 1/(sqrt(1.0 + theta0));
// }
