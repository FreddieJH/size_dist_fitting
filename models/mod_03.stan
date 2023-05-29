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
  real<lower = 0, upper = 5> ln_sigma[S]; // standard deviation for each species
}

model {
  real bin_prob;
  real sigma;
  real mu;
  // Prior distributions for mu and sigma
  ln_mu ~ normal(2.4, 0.8);
  ln_sigma ~ normal(1.5, 0.8);
  
  // Likelihood for each observation
  
  // loop over each species
  for (sp in 1:S) { 
    
    sigma = exp(ln_sigma[sp]);
    mu = exp(ln_mu[sp]);
    
    // within a species, all bins have to come from one or the other dist:
    for (i in i_min[sp]:i_max[sp]) { // for each species
        
      bin_prob = (normal_cdf(l[b[i]+1], mu, sigma) - 
        normal_cdf(l[b[i]], mu, sigma)); 
        
        
        bin_prob = fmax(1e-8, bin_prob);

      // target += n[i]*log(bin_prob); // this is shanes likelihood
      target += binomial_lpmf(n[i] | N_species[sp], bin_prob); // my likelihood
    }

    
  }

}

generated quantities {
  real sigma[S];
  real mu[S];
  for (sp in 1:S) { 
    sigma[sp] = exp(ln_sigma[sp]);
    mu[sp] = exp(ln_mu[sp]);
  }
}
