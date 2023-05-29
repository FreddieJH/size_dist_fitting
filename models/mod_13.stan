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
  // normal parameters
  real<lower = 0, upper = 5> ln_mu[S]; // mean size for each species
  real<lower = 0, upper = 5> ln_sigma[S]; // standard deviation for each species
  // lognormal parameters
  real<lower = 0, upper = 6> meanlog[S]; // mean size for each species
  real<lower = 0, upper = 3> sdlog[S]; // standard deviation for each species

  vector<lower= -2.0, upper=2.0>[S] logit_q;   // logit prob of normal
}

model {
  real bin_prob_norm;
  real bin_prob_lnorm;
  
  real total_prob_norm;
  real total_prob_lnorm;
  
  real normalisation_const_norm;
  real normalisation_const_lnorm;
  
  real sigma;
  real mu;
  
  real max_prob;
  real normal_prob;
  real lognormal_prob;
  real q;
  
  // Prior distributions
  ln_mu ~ normal(2.4, 0.8);
  ln_sigma ~ normal(1.5, 0.8);

  meanlog ~ normal(2.4, 0.8);
  sdlog ~ normal(0.4, 0.2);
  
  
  logit_q ~ normal(0, 1);

  // Likelihood for each observation
  
  // loop over each species
  for (sp in 1:S) { 
    total_prob_norm = 0; // reset for each species
    total_prob_lnorm = 0; // reset for each species
    sigma = exp(ln_sigma[sp]);
    mu = exp(ln_mu[sp]);
    // unlogit
    q = exp(logit_q[sp]) / (1.0 + exp(logit_q[sp]));
    
    
    // within a species, all bins have to come from one or the other dist:
    for (i in i_min[sp]:i_max[sp]) { // for each species
    
    // probability of NOT being in the first bin (i.e. less than 1.25cm)
    normalisation_const_norm = 1 - normal_cdf(l[1], mu, sigma); 
    normalisation_const_lnorm = 1 - lognormal_cdf(l[1], meanlog, sdlog); 
        
      bin_prob_norm = (normal_cdf(l[b[i]+1], mu, sigma) - 
        normal_cdf(l[b[i]], mu, sigma))/normalisation_const_norm; 
      bin_prob_lnorm = (lognormal_cdf(l[b[i]+1], meanlog, sdlog) - 
        lognormal_cdf(l[b[i]], meanlog, sdlog))/normalisation_const_lnorm; 
        
        bin_prob_norm = fmax(1e-8, bin_prob_norm);
        bin_prob_lnorm = fmax(1e-8, bin_prob_lnorm);

      // target += n[i]*log(bin_prob); // this is shanes likelihood
      total_prob_norm += binomial_lpmf(n[i] | N_species[sp], bin_prob_norm); // my likelihood
      total_prob_lnorm += binomial_lpmf(n[i] | N_species[sp], bin_prob_lnorm); // my likelihood
    
    }
    
    max_prob = fmax(total_prob_norm, total_prob_lnorm);
    normal_prob = exp(total_prob_norm - max_prob);
    lognormal_prob = exp(total_prob_lnorm - max_prob);

    target += log(q*exp(normal_prob) + (1.0 - q)*exp(lognormal_prob)); // add species log-likelihood term

  }
}

generated quantities {
  
  real sigma[S];   // logit prob of normal
  real mu[S];   // logit prob of normal
  real q[S];
  
  for (sp in 1:S) {
    q[sp] = exp(logit_q[sp]) / (1.0 + exp(logit_q[sp]));
    sigma[sp] = exp(ln_sigma[sp]);
    mu[sp] = exp(ln_mu[sp]);
  }

}
