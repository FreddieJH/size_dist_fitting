// Step_2_Norm_LNorm.stan
// Shane.Richards@utas.edu.au
// 31/05/2023
// Fit normal and log-normal distributions to fish-length distributions
// Independently estimate population mean and coefficient of variation
 
data {
  int<lower=1> I;                  // number of observations
  int<lower=1> S;                  // number of species
  int<lower=1> B;                  // number of length bins
  real<lower=0.0> l[B+1];          // boundaries of each bin
  int<lower=1,upper=I> i_min[S];   // min obs index for each species
  int<lower=1,upper=I> i_max[S];   // max obs index for each species
  int<lower=1,upper=S> s[I];       // species ID (not used here)
  int<lower=1,upper=B> b[I];       // bin
  int<lower=1>         n[I];       // fish in the bin
}
 
parameters {
	// normal parameters
  vector<lower= -1.0,upper=4.0>[S]  ln_mu;     // log-species mean length
  vector<lower= -4,upper=4.0>[S]  ln_cv;     // log-species coefficient of variation
  vector<lower= -8.0,upper=-4.0>[S] logit_eps_N; // logistic-species misclassification

	// log-normal parameters
  vector<lower= -1.0,upper=4.0>[S]  ln_meanlog;     // log-species mean length
  vector<lower= -4,upper=4.0>[S]  ln_sdlog;  // log-species sihgma
  vector<lower= -8.0,upper=-4.0>[S] logit_eps_LN; // logistic-species misclassification
  
  real<lower = 0, upper = 20> select_mu;
  real<lower = 0, upper = 10> select_s;
}

model {
	vector[B] f;     // probability in bin when randomly allocated

  // normal variables
  real mu;       // species standard deviation
  real sigma;    // species standard deviation
	real norm_c_N;   // normalising constant (less than 1)
	real p_N;        // probability of being in observed bin
	vector[S] eps_N; // probability randomly allocated for species

  // log-normal variables
  real meanlog;
  real sdlog;    // species standard deviation
	real norm_c_LN;   // normalising constant (less than 1)
	real p_LN;        // probability of being in observed bin
	vector[S] eps_LN; // probability randomly allocated for species
	
	real selectivity;

	// priors on model parameters
	ln_mu     ~ normal( 2.0, 1.0); // Log prior for species mean lengths
	ln_cv     ~ normal(-0.6, 1.0); // Log prior on cv 
	logit_eps_N ~ normal(-6.0, 1.0); // Logistic prior on misclassification 

	ln_meanlog ~ normal( 2.0, 1.0); // Log prior for species mean lengths
	ln_sdlog  ~ normal(-0.6, 1.0); // Log prior on cv 
	logit_eps_LN ~ normal(-6.0, 1.0); // Logistic prior on misclassification 
	
	select_mu ~ normal(0, 2);
	select_s ~ normal(0, 2);
	

  for (i in 1:B) {
    f[i] = (l[i+1] - l[i]) / (l[B+1] - l[1]);	// calculate relative bin widths
  }
  
  // calculate misclassification probabilities
  for (i in 1:S) {
    eps_N[i]  = exp(logit_eps_N[i]) / (1.0 + exp(logit_eps_N[i]));	
    eps_LN[i] = exp(logit_eps_LN[i]) / (1.0 + exp(logit_eps_LN[i]));	
  }

  for (i in 1:S) { // for each species i
    // normal
    mu     = exp(ln_mu[i]);
    sigma  = mu * exp(ln_cv[i]);
   	norm_c_N = 1.0 - normal_cdf(l[1], mu, sigma); // normalising constant
   	// log-normal
   	meanlog     = exp(ln_meanlog[i]);
    sdlog  = exp(ln_sdlog[i]);
   	norm_c_LN = 1.0 - lognormal_cdf(l[1], meanlog, sdlog); // normalising constant

    // selectivity = 1/(1+exp(-((l[1]-select_mu)/select_s)));	 

   	for (j in i_min[i]:i_max[i]) { // observation j
   		// probability of being in bin (prior to misclassification)
      p_N = (normal_cdf(l[b[j]+1], mu, sigma) - 
        normal_cdf(l[b[j]], mu, sigma))*(1/(1+exp(-((((l[b[j]]+l[b[j]+1])/2)-select_mu)/select_s)))); 
      // add misclassification probability (ensures non-zero p)
      p_N = (1.0 - eps_N[i])*p_N + eps_N[i]*f[b[j]]; 
      
   		// probability of being in bin (prior to misclassification)
      p_LN = (lognormal_cdf(l[b[j]+1], meanlog, sdlog) - 
        lognormal_cdf(l[b[j]], meanlog, sdlog))*(1/(1+exp(-((((l[b[j]]+l[b[j]+1])/2)-select_mu)/select_s)))); 
      // add misclassification probabily (ensures non-zero p)
      p_LN = (1.0 - eps_LN[i])*p_LN + eps_LN[i]*f[b[j]]; 
      
      target += n[j]*log(p_N) + n[j]*log(p_LN); // add log-likelihood term
   	}
  }
}

generated quantities {
  real mu[S];   // species mean
  real cv[S];   // species cv
  real sigma[S];   // species cv
  real eps_N[S];  // species misclassification
  
  real meanlog[S];    // species mean
  real sdlog[S]; // species sigma
  real eps_LN[S];   // species misclassification
  
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
