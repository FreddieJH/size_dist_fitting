data {
  int<lower=1> n;            // number of observations
  vector[n] y;               // body size data
}
parameters {
  real            meanlog;    
  real<lower= 0>  sdlog;  
}
model {
  // for (i in 1:n) {
    y ~ lognormal(meanlog, sdlog);
  // }
}
