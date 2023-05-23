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
  real<lower = 0, upper = 30> mu[S]; // mean size for each species
  real<lower = 0, upper = 10> ln_sigma[S]; // standard deviation for each species
  
 real<lower = -10, upper = 10> meanlog[S]; // mean size for each species
  real<lower = 0, upper = 5> sdlog[S]; // standard deviation for each species

}

model {
  real bin_prob_norm;
  real bin_prob_lnorm;
  real total_prob_norm;
  real total_prob_lnorm;
  real normalisation_const_norm;
   real normalisation_const_lnorm;
  real sigma;
  real sdlog;
  // Prior distributions for mu and sigma
  mu ~ normal(4, 5);


  // Likelihood for each observation
  
  // loop over each species
  for (sp in 1:S) { 
    total_prob_norm = 0; // reset for each species
    sigma = exp(ln_sigma[sp]);
    // within a species, all bins have to come from one or the other dist:
    for (i in i_min[sp]:i_max[sp]) { // for each species
    
    // probability of NOT being in the first bin (i.e. less than 1.25cm)
    normalisation_const_norm = 1 - normal_cdf(l[1], mu[sp], sigma); 
        
      bin_prob_norm = (normal_cdf(l[b[i]+1], mu[sp], sigma) - 
        normal_cdf(l[b[i]], mu[sp], sigma))/normalisation_const_norm; 
        
        bin_prob_norm = fmax(1e-8, bin_prob_norm);

      // target += n[i]*log(bin_prob); // this is shanes likelihood
      total_prob_norm = total_prob_norm + binomial_lpmf(n[i] | N_species[sp], bin_prob_norm); // my likelihood
    }

    target += total_prob_norm;
  }

}

generated quantities {

}
