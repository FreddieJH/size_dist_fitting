data {
  int<lower=1>         B;       // number of length bins
  real<lower = 0.0>   b_upr[B];
  real<lower = 0.0>   b_lwr[B];
  int<lower=1>         n[B];    // number of fish in the bin
  real<lower = 0.0> low_bound; // absolute lowest cutoff (bottom of smallest bin in all data)
}

parameters {
  real             meanlog;    // population meanlog
  real<lower=0.0>  sdlog;     // population sdlog
}

model {
  
  real nconst; // normalising constant (probability > 1.25cm)
  
  real p_upr;      // cumulative probability of upper size bin boundary
  real p_lwr;      // cumulative probability of lower size bin boundary
  real p;          // probability of being in bin
  
  // RLS census data assumed nothing below 1.25cm (bottom of 2.5cm bin)
  nconst = 1.0 - lognormal_cdf(low_bound, meanlog, sdlog); 
  
  // for each size bin i
  for (i in 1:B) { 
    
    p_upr = lognormal_cdf(b_upr[i], meanlog, sdlog);
    p_lwr = lognormal_cdf(b_lwr[i], meanlog, sdlog);
    
    p = (p_upr - p_lwr)/ nconst;
    
    target += n[i]*log(p); 
  }
}
