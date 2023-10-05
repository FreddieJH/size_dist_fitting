data {
  int<lower=1> n;            // number of observations
  vector[n] y;               // body size data
}
parameters {
  real            mu;    
  real<lower= 0>  sigma;  
}
model {
  // for (i in 1:n) {
    y ~ normal(mu, sigma);
  // }
}
