data {
  int<lower=1> I;                // number of observations
  int<lower=1> S;                // number of species
  int<lower=1> B;                // number of length bins
  real<lower=0.0> l[B+1];        // boundaries of each bin
  int<lower=1>    N_species[S];  // total number observed
  int<lower=1,upper=I> i_min[S]; // min obs index for each species
  int<lower=1,upper=I> i_max[S]; // max obs index for each species
  int<lower=1,upper=S> s[I];     // species ID
  int<lower=1,upper=B> b[I];     // bin
  int<lower=0>         n[I];     // fish in the bin
}
parameters {
  real<lower = -10, upper = 10> meanlog[S]; // mean size for each species
  real<lower = 0, upper = 5> sdlog[S]; // standard deviation for each species
}

model {
  real bin_prob;
  real normalisation_const;
  // Prior distributions for mu and sigma
  meanlog ~ normal(3, 5);
  sdlog ~ normal(0, 4);
  

  
  // Likelihood for each observation
  
  // loop over each species
  
  for (sp in 1:S) { 

    // within a species, all bins have to come from one or the other dist:
    for (i in i_min[sp]:i_max[sp]) { // for each species
    
    // probability of NOT being in the first bin (i.e. less than 1.25cm)
    normalisation_const = 1 - normal_cdf(l[1], meanlog[sp], sdlog[sp]); 
        
      bin_prob = (lognormal_cdf(l[b[i]+1], meanlog[sp], sdlog[sp]) - 
        lognormal_cdf(l[b[i]], meanlog[sp], sdlog[sp]))/normalisation_const; 
        
        bin_prob = fmax(1e-8, bin_prob);

      target += binomial_lpmf(n[i] | N_species[sp], bin_prob);
    }

    
  }

}
