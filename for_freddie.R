library(rstan)
stan_data <- readRDS("C:/Users/MS23/Desktop/Temp/freddie/stan_data.rds")
st_mod <- stan_model("simple.stan")

fit <- sampling(st_mod,data=stan_data)
