// Shane.Richards@utas.edu.au or freddie.heather@utas.edu.au
// June 2023
// Fit normal and log-normal distributions to fish-length distributions
// Independently estimate population mean and coefficient of variation
// Inlcuding a selectivity term for the detectibility of size classes

data {
  int<lower=1> I;                  // number of observations
  int<lower=1> S;                  // number of populations
  int<lower=1> B;                  // number of length bins
  real<lower=0.0> l[B+1];          // boundaries of each bin
  int<lower=1,upper=I> i_min[S];   // min obs index for each population
  int<lower=1,upper=I> i_max[S];   // max obs index for each population
  int<lower=1,upper=S> s[I];       // population ID (not used here)
  int<lower=1,upper=B> b[I];       // bin
  int<lower=1>         n[I];       // fish in the bin
  real<lower=0.0> meansize[S];      // mean body size of the population (this is not currently used)
}

parameters {
  // normal parameters
  vector<lower= -1.0,upper=6.0>[S]  ln_mu;       // log-population mean length
  vector<lower= -4.0,upper=2.0>[S]  ln_cv;       // log-population coefficient of variation
  vector<lower= -6.0,upper=-4.0>[S] logit_eps_N; // logistic-population misclassification
  
  // log-normal parameters
  vector<lower= -1.0,upper=2.0>[S]  ln_meanlog;   // log-population mean length
  vector<lower= -4,upper=2.0>[S]  ln_sdlog;       // log-population sigma
  vector<lower= -6.0,upper=-4.0>[S] logit_eps_LN; // logistic-population misclassification
  
  // detectability parameters
  real<lower=-4.0, upper=10.0> logit_detectability_1; // logistic-detectability of 2.5cm bin
  real<lower=-4.0, upper=10.0> logit_detectability_2; // logistic-detectability of 5cm bin
  real<lower=-4.0, upper=10.0> logit_detectability_3; // logistic-detectability of 7.5cm bin
}

model {
  
  vector[B] f;     // probability in bin when randomly allocated
  
  // normal variables
  real mu;       // population mean
  real sigma;    // population standard deviation
  real norm_c_N;   // normalising constant (less than 1)
  real p_N;        // probability of being in observed bin
  vector[S] eps_N; // probability randomly allocated for population
  
  // log-normal variables
  real meanlog;  // population meanlog 
  real sdlog;    // population sdlog
  real norm_c_LN;   // normalising constant (less than 1)
  real p_LN;        // probability of being in observed bin
  vector[S] eps_LN; // probability randomly allocated for population
  
  // detectibility variables
  real detectability_1; // detectability of 2.5cm bin (bin #1)
  real detectability_2; // detectability of 5cm bin (bin #2)
  real detectability_3; // detectability of 7.5cm bin (bin #3)
  real est_N;           // estimated N (accounting for detectibility)
  
  // priors on model parameters
  // Normal priors
  ln_mu     ~ normal( 2.0, 1.0); // Log prior for species mean lengths
  ln_cv     ~ normal(-0.6, 1.0); // Log prior on cv 
  logit_eps_N ~ normal(-4.0, 1.0); // Logistic prior on misclassification 
  
  // Lognormal priors
  ln_meanlog ~ normal(1.0, 0.5); // Log prior for species mean lengths
  ln_sdlog  ~ normal(-0.6, 1.0); // Log prior on cv 
  logit_eps_LN ~ normal(-4.0, 1.0); // Logistic prior on misclassification 
  
  // detectability priors
  logit_detectability_1 ~ normal(-1, 1);
  logit_detectability_2 ~ normal(-1, 1);
  logit_detectability_3 ~ normal(-1, 1);
  
  detectability_1 = exp(logit_detectability_1) / (1.0 + exp(logit_detectability_1));
  detectability_2 = exp(logit_detectability_2) / (1.0 + exp(logit_detectability_2));
  detectability_3 = exp(logit_detectability_3) / (1.0 + exp(logit_detectability_3));
  
  for (i in 1:B) {
    f[i] = (l[i+1] - l[i]) / (l[B+1] - l[1]);	// calculate relative bin widths
  }
  
  // calculate misclassification probabilities
  for (i in 1:S) {
    eps_N[i]  = exp(logit_eps_N[i]) / (1.0 + exp(logit_eps_N[i]));	
    eps_LN[i] = exp(logit_eps_LN[i]) / (1.0 + exp(logit_eps_LN[i]));	
  }
  
  for (i in 1:S) { // for each population i
  
  // normal
  mu     = exp(ln_mu[i]);
  sigma  = mu * exp(ln_cv[i]);
  norm_c_N = 1.0 - normal_cdf(l[1], mu, sigma); // normalising constant
  
  // log-normal
  meanlog     = exp(ln_meanlog[i]);
  sdlog  = exp(ln_sdlog[i]);
  norm_c_LN = 1.0 - lognormal_cdf(l[1], meanlog, sdlog); // normalising constant
  
  for (j in i_min[i]:i_max[i]) { // observation j
  
  est_N = n[j]; // est abundance is the abundance in the bin bin otherwise...
  
  if(b[j]==1){ // if in first bin (2.5cm)
  est_N = n[j]*(1/detectability_1);
  }
  if(b[j]==2){ // if in second bin (5cm)
  est_N = n[j]*(1/detectability_2);
  }
  if(b[j]==3){ // if in third bin (7.5cm)
  est_N = n[j]*(1/detectability_3);
  }
  
  
  // probability of being in bin (prior to misclassification)
  p_N = (normal_cdf(l[b[j]+1], mu, sigma) - 
  normal_cdf(l[b[j]], mu, sigma)); 
  p_N = p_N/norm_c_N;
  
  // add misclassification probability (ensures non-zero p)
  p_N = (1.0 - eps_N[i])*p_N + eps_N[i]*f[b[j]]; 
  
  // probability of being in bin (prior to misclassification)
  p_LN = (lognormal_cdf(l[b[j]+1], meanlog, sdlog) - 
  lognormal_cdf(l[b[j]], meanlog, sdlog));
  p_LN = p_LN/norm_c_LN;
  // add misclassification probabily (ensures non-zero p)
  p_LN = (1.0 - eps_LN[i])*p_LN + eps_LN[i]*f[b[j]]; 
  
  target += est_N*log(p_N) + est_N*log(p_LN); // add log-likelihood term
  }
  }
}

generated quantities {
  real mu[S];   // population mean
  real cv[S];   // population cv
  real sigma[S];   // population cv
  real eps_N[S];  // population misclassification
  
  real meanlog[S];    // population mean
  real sdlog[S]; // population sigma
  real eps_LN[S];   // population misclassification
  
  real detectability_1 = exp(logit_detectability_1)/(1.0 + exp(logit_detectability_1));
  real detectability_2 = exp(logit_detectability_2)/(1.0 + exp(logit_detectability_2));
  real detectability_3 = exp(logit_detectability_3)/(1.0 + exp(logit_detectability_3));
  
  for (i in 1:S) {
    mu[i]  = exp(ln_mu[i]);
    cv[i]  = exp(ln_cv[i]);
    sigma[i] = cv[i]*mu[i];
    eps_N[i] = exp(logit_eps_N[i]) / (1.0 + exp(logit_eps_N[i]));
    
    meanlog[i]   = exp(ln_meanlog[i]);
    sdlog[i] = exp(ln_sdlog[i]);
    eps_LN[i]   = exp(logit_eps_LN[i]) / (1.0 + exp(logit_eps_LN[i]));
  }
}
