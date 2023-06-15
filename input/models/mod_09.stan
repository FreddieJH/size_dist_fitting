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
   real<lower = 0, upper = 5> ln_mu[S]; // mean size for each species
  real<lower = 0, upper = 10> ln_sigma[S]; // standard deviation for each species
}

model {
  real bin_prob;
  real normalisation_const;
  real sigma;
  real mu;
  
  real total_binprob;
  // Prior distributions for mu and sigma
 ln_mu ~ normal(2.4, 0.8);
  ln_sigma ~ normal(0, 4);
  

  
  // Likelihood for each observation
  
  // loop over each species
  for (sp in 1:S) { 
    total_binprob = 0.0;
    mu = exp(ln_mu[sp]);
      sigma = exp(ln_sigma[sp]);
    // within a species, all bins have to come from one or the other dist:
    for (i in i_min[sp]:i_max[sp]) { // for each species
    
    // probability of NOT being in the first bin (i.e. less than 1.25cm)
    normalisation_const = 1 - normal_cdf(l[1], mu, sigma); 
        
      bin_prob = (normal_cdf(l[b[i]+1], mu, sigma) - 
        normal_cdf(l[b[i]], mu, sigma))/normalisation_const; 
        
        bin_prob = fmax(1e-200, bin_prob);
        
        total_binprob += bin_prob;
        target += n[i]*log(bin_prob);
    }

    print(total_binprob);
  }

}

generated quantities {
  real sigma[S];
  real mu[S];
  for (sp in 1:S) {
    mu[sp] = exp(ln_mu[sp]);
    sigma[sp] = exp(ln_sigma[sp]);
  }
}