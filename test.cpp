#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double meanC(NumericVector x) {
  int n = x.size();
  double total = 0;
  
  for(int i = 0; i < n; ++i) {
    total += x[i];
  }
  return total / n;
};

// [[Rcpp::export]]
NumericVector pdistC(double x, NumericVector ys) {
  int n = ys.size();
  NumericVector out(n);
  
  for(int i = 0; i < n; ++i) {
    out[i] = sqrt(pow(ys[i] - x, 2.0));
  }
  return out;
};

// [[Rcpp::export]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
class Func
{
public:
  virtual double operator()(const double& x) const = 0;
  virtual void eval(double* x, const int n) const
  {
    for(int i = 0; i < n; i++)
      x[i] = this->operator()(x[i]);
  }
  
  virtual ~Func() {}
};

