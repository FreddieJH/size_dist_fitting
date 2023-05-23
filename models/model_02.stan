data {
  int<lower=1> n;           // number of observations
  int<lower=1> n_species;   // number of species
  int<lower=1, upper=n_species> species_id[n];  // species ID for each observation
  vector[n] y;              // size data
}
parameters {
  real<lower = -10, upper = 10> meanlog[n_species]; // mean size for each species
  real<lower = 0, upper = 10> sdlog[n_species]; // standard deviation for each species
}

model {
  // Prior distributions for mu and sigma
  meanlog ~ normal(0, 10);
  sdlog ~ cauchy(0, 5);

  // Likelihood for each observation
  for (i in 1:n) {
    y[i] ~ lognormal(meanlog[species_id[i]], sdlog[species_id[i]]);
  }
}
