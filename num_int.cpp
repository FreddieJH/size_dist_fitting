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
    return Rcpp::dnorm(x, mu, sigma)*x/nc;
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

}


// f <- function(x, mu, sigma, nc) (x*dnorm(x, mu, sigma))/nc
// f_int <- function(mu, sigma, nc) {
//   integrate(f, lower = 1.25, upper = Inf, abs.tol = 0, mu = mu, sigma = sigma, nc = nc)$value
// }

// [[Rcpp::export]]


// P(0.3 < X < 0.8), X ~ Beta(a, b)
class ggg: public Func
{
private:
  double a;
  double b;
public:
  ggg(double a_, double b_) : a(a_), b(b_) {}
  
  double operator()(const double& x) const
  {
    return R::dbeta(x, a, b, 0);
  }
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

// [[Rcpp::export]]
Rcpp::List integrate_test(Rcpp::NumericVector mu, Rcpp::NumericVector sigma, Rcpp::NumericVector nc)
{
  const double a = 3, b = 10;
  const double lower = 0.3, upper = 0.8;
  const double true_val = R::pbeta(upper, a, b, 1, 0) -
    R::pbeta(lower, a, b, 1, 0);
  
  int n = mu.size();
  Rcpp::NumericVector res(n);
  
  ggg func(a, b);
  double err_est;
  int err_code;
  res = integrate(func, lower, upper, err_est, err_code);
  return Rcpp::List::create(
    Rcpp::Named("true") = true_val,
    Rcpp::Named("approximate") = res,
    Rcpp::Named("error_estimate") = err_est,
    Rcpp::Named("error_code") = err_code
  );
}


// // [[Rcpp::export]]
// Rcpp::List integrate_test(Rcpp::NumericVector mu, Rcpp::NumericVector sigma)
// {
//   int n = ys.size();
//   const double lower = 1.25, upper = 1e4;
//   for(int i = 0; i < n; ++i) {
//     out[i] = sqrt(pow(ys[i] - x, 2.0));
//   }
//   const double true_val = R::pnorm5(upper, mu, sigma, 1, 0) -
//     R::pnorm5(lower, mu, sigma, 1, 0);
// 
//   normPDF f(mu, sigma);
//   double err_est;
//   int err_code;
//   const double res = integrate(f, lower, upper, err_est, err_code);
//   return Rcpp::List::create(
//     Rcpp::Named("true") = true_val,
//     Rcpp::Named("approximate") = res,
//     Rcpp::Named("error_estimate") = err_est,
//     Rcpp::Named("error_code") = err_code
//   );
// }