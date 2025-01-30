#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector cpprbinom(int n, double size, NumericVector prob) {
  NumericVector v = no_init(n);
  std::transform(prob.begin(), prob.end(), v.begin(), [=](double p){return R::rbinom(size, p);});
  return(v);
}
