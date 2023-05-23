data {
  int<lower=1> n;           // number of observations
  int<lower=1> n_species;   // number of species
  int<lower=1, upper=n_species> species_id[n];  // species ID for each observation
  vector[n] y;              // size data
}
parameters {
  real<lower=0> mu[n_species];       // mean size for each species
  real<lower=0> sigma[n_species]; // standard deviation for each species
}
model {
  // Prior distributions for mu and sigma
  mu ~ normal(0, 10);
  sigma ~ cauchy(0, 5);

  // Likelihood for each observation
  for (i in 1:n) {
    y[i] ~ normal(mu[species_id[i]], sigma[species_id[i]]);
  }
}
