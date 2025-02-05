#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::vec cpprbinom(int n, double size, arma::vec prob) {
  arma::vec v = arma::vec(n);
  std::transform(prob.begin(), prob.end(), v.begin(), [=](double p){return R::rbinom(size, p);});
  return(v);
}

// [[Rcpp::export]]
NumericMatrix calc_RRexpo(int iter, arma::mat& obs_mat, arma::mat& draws) {
  arma::vec p(obs_mat.n_rows);
  arma::mat test_mat = arma::mat(obs_mat.n_rows, iter);
  for (int i = 0; i < iter; i++) {
    p = obs_mat.col(2) * draws(i, 11) +
      obs_mat.col(3) * draws(i, 12) +
      obs_mat.col(4) * (1 - draws(i, 13)) +
      obs_mat.col(5) * (1 - draws(i, 14));
    test_mat.col(i) = p;
  }

  return wrap(test_mat);
}
