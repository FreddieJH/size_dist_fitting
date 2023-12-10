// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
#include <RcppNumerical.h>
using namespace Numer;

// P(0.3 < X < 0.8), X ~ Beta(a, b)
class normPDFmean: public Func
{
private:
  double mu;
  double sigma;
  double nc;
public:
  normPDFmean(double mu_, double sigma_, double nc_) : mu(mu_), sigma(sigma_), nc(nc_) {}

  double operator()(const double& x) const
  {
    return R::dnorm4(x, mu, sigma, 0)*x/nc;
  }
};

// [[Rcpp::export]]
Rcpp::NumericVector normConst(Rcpp::NumericVector mu, Rcpp::NumericVector sigma)
{
  int n = mu.size();
  Rcpp::NumericVector true_val(n);
  const double lower = 1.25, upper = 1e4;
  
  for(int i = 0; i < n; ++i) {
    true_val[i] = R::pnorm5(upper, mu[i], sigma[i], 1, 0) -
      R::pnorm5(lower, mu[i], sigma[i], 1, 0);
  }
  return true_val;
  
};


Rcpp::NumericVector num_int_mean(Rcpp::NumericVector mu, Rcpp::NumericVector sigma, Rcpp::NumericVector nc)
{
  int n = mu.size();
  Rcpp::NumericVector res(n);
  const double lower = 1.25, upper = 1e4;
  ggg func(a, b);
  ggg f(mu, sigma);
  res = integrate(f, lower, upper);
  
  for(int i = 0; i < n; ++i) {
    normPDFmean f(mu[i], sigma[i], nc[i]);
    res[i] = integrate(f, lower, upper, mu=mu, sigma=sigma, nc=nc);
  }
  return res;
  
}