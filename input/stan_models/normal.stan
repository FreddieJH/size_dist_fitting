data {
  int<lower=1>         B;       // number of length bins
  real<lower = 0.0>   b_upr[B];
  real<lower = 0.0>   b_lwr[B];
  int<lower=1>         n[B];    // number of fish in the bin
  real<lower = 0.0> low_bound; // absolute lowest cutoff
}

parameters {
  real             mu;    // population meanlog
  real<lower=0.0>  sigma;     // population sdlog
}

model {
  real nconst; // normalising constant (for RLS data this is prob below 2.5cm)
  real p_upr; // cumulative probability of upper size bin boundary
  real p_lwr; // cumulative probability of lower size bin boundary
  real p;     // probability of being in bin
  
  nconst = 1.0 - normal_cdf(low_bound, mu, sigma); 
  
  // for each size bin i
  for (i in 1:B) { 
    p_upr = normal_cdf(b_upr[i], mu, sigma);
    p_lwr = normal_cdf(b_lwr[i], mu, sigma);
    p = (p_upr - p_lwr)/ nconst;
    target += n[i]*log(p); 
  }
}
