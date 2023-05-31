data {
  int<lower=1> I;                // number of observations
  int<lower=1> S;                // number of species
  int<lower=1> B;                // number of length bins
  real<lower=0.0>      l[B+1];        // boundaries of each bin
  int<lower=1>         N_species[S];  // total number observed
  int<lower=1,upper=I> i_min[S]; // min obs index for each species
  int<lower=1,upper=I> i_max[S]; // max obs index for each species
  int<lower=1,upper=S> s[I];     // species ID
  int<lower=1,upper=B> b[I];     // bin
  int<lower=0>         n[I];     // number of fish in the bin
}
parameters {

  real<lower = 0, upper = 6> meanlog[S]; // mean size for each species
  
  real<lower = -3, upper = 3> beta_0;
  real<lower = -0.1, upper = 0.1> beta_1;
}

model {
  real bin_prob;
  real sdlog;
  real normalisation_const;
  // real out_prob;
  // Prior distributions for mu and sigma
  meanlog ~ normal(3, 1);
  
  beta_0 ~ normal(1, 0.5);
  beta_1 ~ normal(0, 0.5);

  // Likelihood for each observation
  
  // loop over each species
  
  for (sp in 1:S) { 
     // might be better to normalise meanlog to center around zero
    sdlog = beta_0 + (beta_1*meanlog[sp]);
    
    for (i in i_min[sp]:i_max[sp]) { // for each species
      
    // probability of NOT being in the first bin (i.e. less than 1.25cm)
    normalisation_const = 1 - lognormal_cdf(l[1], meanlog[sp], sdlog); 
        
      bin_prob = (lognormal_cdf(l[b[i]+1], meanlog[sp], sdlog) - 
        lognormal_cdf(l[b[i]], meanlog[sp], sdlog))/normalisation_const; 
        
      // this is questionable
      // shane will send through code to replace this part
      bin_prob = fmax(1e-8, bin_prob);

      target += n[i]*log(bin_prob);
      
    }
  }

}

generated quantities {
  
  real sdlog[S];
  
  for (sp in 1:S) { 
    sdlog[sp] = beta_0 + (beta_1*meanlog[sp]);
  }
  
}
