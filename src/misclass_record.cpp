#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
NumericVector cpprbinom(int n, double size, arma::vec prob) {
  arma::vec v = arma::vec(n);
  std::transform(prob.begin(), prob.end(), v.begin(), [=](double p){return R::rbinom(size, p);});
  return(wrap(v));
}
