#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

//' Set the RNG Seed from within Rcpp
//'
//' Within Rcpp, one can set the R session seed without triggering
//' the CRAN rng modifier check.
//' @param seed A \code{unsigned int} that is the seed one wishes to use.
//' @return A set RNG scope.
//' @examples
//' set.seed(10)
//' x = rnorm(5,0,1)
//' set_seed(10)
//' y = rnorm(5,0,1)
//' all.equal(x,y, check.attributes = F)
// [[Rcpp::export]]
void set_seed(unsigned int seed) {
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(seed);
}


// [[Rcpp::export]]
arma::vec cpprbinom(int n, double size, arma::vec prob) {
  arma::vec v = arma::vec(n);
  std::transform(prob.begin(), prob.end(), v.begin(), [=](double p){return R::rbinom(size, p);});
  return(v);
}

// [[Rcpp::export]]
List calc_RRexpo(int iter, arma::mat& obs_mat, arma::mat& draws) {
  //Obtaining namespace of fastglm package
  Environment pkg = Environment::namespace_env("fastglm");
  Function f = pkg["fastglm"];

  arma::vec p(obs_mat.n_rows);
  arma::mat ptest_mat = arma::mat(obs_mat.n_rows, iter);  // TEMP: to check we've got the right p
  arma::vec e(obs_mat.n_rows);
  arma::vec d = obs_mat.col(1);
  arma::mat etest_mat = arma::mat(obs_mat.n_rows, iter);  // TEMP: to check we've got the right e
  arma::vec coef(iter);
  arma::vec mod_se(iter);
  Rcpp::List mod_pois;
  //  set_seed(seed);
  for (int i = 0; i < iter; i++) {
    p = obs_mat.col(2) * draws(i, 11) +
      obs_mat.col(3) * draws(i, 12) +
      obs_mat.col(4) * (1 - draws(i, 13)) +
      obs_mat.col(5) * (1 - draws(i, 14));
    ptest_mat.col(i) = p;  // TEMP: to check we've got the right p
    e = cpprbinom(obs_mat.n_rows, 1, p);
    etest_mat.col(i) = e;  // TEMP: to check we've got the right e
    //    if (all(is_na(wrap(e))).is_true()) {
    //      coef = NA_REAL;
    //      mod_se = NA_REAL;
      //    } else {
      mod_pois = f(Named("x", e), Named("y", d), Named("family", "poisson"));
      //    }
  }

  return Rcpp::List::create(Rcpp::Named("ptest") = wrap(ptest_mat),
			    Rcpp::Named("etest") = wrap(etest_mat));
}
