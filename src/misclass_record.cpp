#define ARMA_WARN_LEVEL 1
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>

// [[Rcpp::export]]
NumericVector cpprbinom(int n, double size, NumericVector prob) {
  NumericVector v = no_init(n);
  std::transform(prob.begin(), prob.end(), v.begin(), [=](double p){return R::rbinom(size, p);});
  return(v);
}

// [[Rcpp::export]]
List define_estar(int iter, arma::mat& obs_mat, arma::mat& draws) {
  // Obtain environment containing function
  //Rcpp::Environment base("package:stats");
  // Obtaining namespace of fastglm package
  Environment pkg = Environment::namespace_env("fastglm");
  Function f = pkg["fastglm"];
  Function s = pkg["summary.fastglm"];

  //  Rcpp::NumericVector d = obs_mat(_, 1);
  arma::vec p(obs_mat.n_rows);
  arma::mat e = arma::mat(obs_mat.n_rows, 2);
  e.ones();
  arma::rowvec outcome = obs_mat.col(1);
  arma::mat e_syst = arma::mat(obs_mat.n_rows, 2);
  e_syst.col(0) = obs_mat.col(1);
  arma::mat d1;
  arma::mat d0;
  double at;
  double ct;
  double dt;
  double bt;
  double coef_syst;
  Rcpp::List mod_pois;
  Rcpp::List mod_sum;
  arma::vec mod_df(obs_mat.n_rows);
  Rcpp::NumericVector mod_coef = no_init(d.length());
  Rcpp::NumericVector coef = no_init(iter);
  for (int i = 0; i < iter; i++) {
    p = obs_mat(_, 2) * draws(i, 11) +
      obs_mat(_, 3) * draws(i, 12) +
      obs_mat(_, 4) * (1 - draws(i, 13)) +
      obs_mat(_, 5) * (1 - draws(i, 14));
    e.col(1) = cpprbinom(obs_mat.n_rows, 1, p);
    e_syst.col(1) = obs_mat(_, 2) * draws(i, 17) +
      obs_mat(_, 3) * draws(i, 18) +
      obs_mat(_, 4) * (1 - draws(i, 19)) +
      obs_mat(_, 5) * (1 - draws(i, 20));
    d1 = e_syst.rows(find(outcome == 1));
    d0 = e_syst.rows(find(outcome == 0));
    at = as<double>(wrap(arma::accu(d1.col(1))));
    ct = as<double>(wrap(arma::accu(d0.col(1))));
    dt = ct - as<double>(wrap(d0.n_rows));
    bt = at - as<double>(wrap(d1.n_rows));
    // For systematic error only
    if (all(is_na(e_syst.col(1))).is_true()) {
      coef_syst = NA_REAL;
    } else {
      if ((at * bt * ct * dt) <= 0) {
	coef_syst = NA_REAL;
      } else {
	coef_syst = log((at * (bt + dt)) / (bt * (at + ct)));
      }
    }
    // Systematic and random error
    if (all(is_na(e.col(1))).is_true()) {
      coef = NA_REAL;
      mod_se = NA_REAL;
    } else {
      mod_pois = f(Named("x", e), Named("y", d), Named("family", "poisson"));
      mod_sum = s(Named("object", mod_pois));
      mod_df = mod_sum["df"];
      mod_df = {mod_df[0], mod_df[1]};
      mod_coef = mod_pois["coefficients"];
      coef = {mod_coef[1]};
      mod_res = mod_pois["residuals"];
      mod_wght = mod_pois["weights"];
      // Bread
      // Get unscaled covariance matrix
      mod_fit = mod_pois["fitted.values"];
      mod_fit2 = 1 - mod_fit;
      W = arma::diagmat(as<arma::vec>(mod_fit));
      unsccov = arma::solve(arma::trans(X) * W * X, I,
			    arma::solve_opts::allow_ugly);
      // bread as unscaled cov matrix * (ncol + nrow) * dispersion (but dispersion = 1)
      bread = unsccov * arma::accu(as<arma::vec>(mod_df));
      // Meat
      wres = mod_res * mod_wght;
      res = X.each_col() % as<arma::vec>(wres);
      res = res.col(0);
      // omega is res^2, but then we have to squared it for rval so no need to compute it
      rval = X.each_col() % res;
      meat = (rval.t() * rval) / d.length();
      sandwich = bread * meat * bread;
      sandwich = sandwich.each_col() % inv_nrow;
      se = sqrt(sandwich.diag());
      se = se(1);
      mod_se = as<double>(wrap(se));
    }
  }
  return Rcpp::List::create(Rcpp::Named("e") = e,
			    Rcpp::Named("coef") = coef);
}


// [[Rcpp::export]]
List expo_adjRR(int iter, NumericMatrix obs_mat, NumericMatrix draws,
		bool display_progress=true) {
  // Obtain environment containing function
  //Rcpp::Environment base("package:stats");
  // Obtaining namespace of fastglm package
  Environment pkg = Environment::namespace_env("fastglm");
  Function f = pkg["fastglm"];
  Function s = pkg["summary.fastglm"];

  Rcpp::NumericVector d = obs_mat(_, 1);
  //  Rcpp::NumericVector pr = no_init(d.length());  // obs_mat(_, 6)
  //  Rcpp::NumericMatrix e(d.length(), iter);  // obs_mat(_, 7)
  Rcpp::NumericVector e_syst = no_init(d.length());  // obs_mat(_, 8)
  arma::vec darma = obs_mat(_, 1);
  arma::mat de_obs = arma::mat(d.length(), 2);
  arma::rowvec outcome = as<arma::rowvec>(d);
  arma::mat d1;
  arma::mat d0;
  double at;
  double ct;
  double dt;
  double bt;
  //  Rcpp::NumericMatrix e_mat(d.length(), iter);
  double coef_syst;
  Rcpp::NumericMatrix ze_syst(iter, 2);

  Rcpp::NumericMatrix res_mat(iter, 2);
  Rcpp::NumericMatrix ematrix(d.length(), 2);
  Rcpp::NumericVector v(d.length(), 1);
  ematrix(_, 0) = v;
  Rcpp::List mod_pois;
  Rcpp::List mod_log;
  Rcpp::NumericVector mod_coef = no_init(d.length());
  double coef;
  Rcpp::List mod_sum;
  Rcpp::NumericVector mod_df = no_init(d.length());
  Rcpp::NumericVector mod_res = no_init(d.length());
  Rcpp::NumericVector mod_wght = no_init(d.length());
  Rcpp::NumericVector mod_fit = no_init(d.length());
  Rcpp::NumericVector mod_fit2 = no_init(d.length());
  arma::mat W(d.length(), d.length());
  arma::mat X(d.length(), 2);  // model matrix
  arma::mat I;
  I.eye(2, 2);
  arma::mat unsccov;
  arma::mat bread;
  Rcpp::NumericVector wres = no_init(d.length());
  arma::mat res(d.length(), 2);
  arma::mat rval(d.length(), 2);
  arma::mat meat;
  arma::mat sandwich;
  arma::vec inv_nrow(2);
  inv_nrow.fill(1.0/d.length());
  arma::vec se(2);
  double mod_se;
  Progress p(iter*iter, display_progress);
  for (int i = 0; i < iter; i++) {
    p.increment(); // update progress
    obs_mat(_, 6) = obs_mat(_, 2) * draws(i, 11) +
      obs_mat(_, 3) * draws(i, 12) +
      obs_mat(_, 4) * (1 - draws(i, 13)) +
      obs_mat(_, 5) * (1 - draws(i, 14));
    obs_mat(_, 7) = cpprbinom(d.length(), 1, obs_mat(_, 6));
    //    e_mat(_, i) = e;
    obs_mat(_, 8) = obs_mat(_, 2) * draws(i, 17) +
      obs_mat(_, 3) * draws(i, 18) +
      obs_mat(_, 4) * (1 - draws(i, 19)) +
      obs_mat(_, 5) * (1 - draws(i, 20));
    e_syst = obs_mat(_, 8);
    arma::vec e_systarma = as<arma::vec>(e_syst);
    de_obs.col(0) = darma;
    de_obs.col(1) = e_systarma;
    d1 = de_obs.rows(find(outcome == 1));
    d0 = de_obs.rows(find(outcome == 0));
    at = as<double>(wrap(arma::accu(d1.col(1))));
    ct = as<double>(wrap(arma::accu(d0.col(1))));
    dt = ct - as<double>(wrap(d0.n_rows));
    bt = at - as<double>(wrap(d1.n_rows));
    // For systematic error only
    if (all(is_na(e_syst)).is_true()) {
      coef_syst = NA_REAL;
    } else {
      if ((at * bt * ct * dt) <= 0) {
	coef_syst = NA_REAL;
      } else {
	coef_syst = log((at * (bt + dt)) / (bt * (at + ct)));
      }
    }

    ze_syst(i, 0) = R::rnorm(0, 1);
    ze_syst(i, 1) = coef_syst;

    ematrix(_, 1) = obs_mat(_, 7);
    X = as<arma::mat>(ematrix);
    // Systematic and random error
    if (all(is_na(obs_mat(_, 7))).is_true()) {
      coef = NA_REAL;
      mod_se = NA_REAL;
    } else {
      mod_pois = f(Named("x", ematrix), Named("y", d), Named("family", "poisson"));
      mod_sum = s(Named("object", mod_pois));
      mod_df = mod_sum["df"];
      mod_df = {mod_df[0], mod_df[1]};
      mod_coef = mod_pois["coefficients"];
      coef = {mod_coef[1]};
      mod_res = mod_pois["residuals"];
      mod_wght = mod_pois["weights"];
      // Bread
      // Get unscaled covariance matrix
      mod_fit = mod_pois["fitted.values"];
      mod_fit2 = 1 - mod_fit;
      W = arma::diagmat(as<arma::vec>(mod_fit));
      unsccov = arma::solve(arma::trans(X) * W * X, I,
			    arma::solve_opts::allow_ugly);
      // bread as unscaled cov matrix * (ncol + nrow) * dispersion (but dispersion = 1)
      bread = unsccov * arma::accu(as<arma::vec>(mod_df));
      // Meat
      wres = mod_res * mod_wght;
      res = X.each_col() % as<arma::vec>(wres);
      res = res.col(0);
      // omega is res^2, but then we have to squared it for rval so no need to compute it
      rval = X.each_col() % res;
      meat = (rval.t() * rval) / d.length();
      sandwich = bread * meat * bread;
      sandwich = sandwich.each_col() % inv_nrow;
      se = sqrt(sandwich.diag());
      se = se(1);
      mod_se = as<double>(wrap(se));
    }

    res_mat(i, 0) = coef;
    res_mat(i, 1) = mod_se;
  }

    return Rcpp::List::create(Rcpp::Named("res_mat") = res_mat,
			      Rcpp::Named("ze_syst") = ze_syst);
}


// [[Rcpp::export]]
List expo_adjOR(int iter, NumericMatrix obs_mat, NumericMatrix draws,
		bool display_progress=true) {
  // Obtain environment containing function
  //Rcpp::Environment base("package:stats");
  // Obtaining namespace of fastglm package
  Environment pkg = Environment::namespace_env("fastglm");
  Function f = pkg["fastglm"];
  Function s = pkg["summary.fastglm"];

  Rcpp::NumericVector d = obs_mat(_, 1);
  Rcpp::NumericVector pr = no_init(d.length());
  Rcpp::NumericMatrix e(d.length(), iter);
  Rcpp::NumericVector e_syst = no_init(d.length());
  arma::vec darma = obs_mat(_, 1);
  arma::mat de_obs = arma::mat(d.length(), 2);
  arma::rowvec outcome = as<arma::rowvec>(d);
  arma::mat d1;
  arma::mat d0;
  double at;
  double ct;
  double dt;
  double bt;
  //  Rcpp::NumericMatrix e_mat(d.length(), iter);
  double coef_syst;
  Rcpp::NumericMatrix ze_syst(iter, 2);

  Rcpp::NumericMatrix res_mat(iter, 2);
  Rcpp::NumericMatrix ematrix(d.length(), 2);
  Rcpp::NumericVector v(d.length(), 1);
  ematrix(_, 0) = v;
  Rcpp::List mod_pois;
  Rcpp::List mod_log;
  Rcpp::NumericVector mod_coef = no_init(d.length());
  double coef;
  Rcpp::List mod_sum;
  Rcpp::NumericVector mod_df = no_init(d.length());
  Rcpp::NumericVector mod_res = no_init(d.length());
  Rcpp::NumericVector mod_wght = no_init(d.length());
  Rcpp::NumericVector mod_fit = no_init(d.length());
  Rcpp::NumericVector mod_fit2 = no_init(d.length());
  arma::mat W(d.length(), d.length());
  arma::mat X(d.length(), 2);  // model matrix
  arma::mat I;
  I.eye(2, 2);
  arma::mat unsccov;
  arma::mat bread;
  Rcpp::NumericVector wres = no_init(d.length());
  arma::mat res(d.length(), 2);
  arma::mat rval(d.length(), 2);
  arma::mat meat;
  arma::mat sandwich;
  arma::vec inv_nrow(2);
  inv_nrow.fill(1.0/d.length());
  arma::vec se(2);
  double mod_se;
  Progress p(iter*iter, display_progress);
  for (int i = 0; i < iter; i++) {
    p.increment(); // update progress
    pr = obs_mat(_, 2) * draws(i, 11) +
      obs_mat(_, 3) * draws(i, 12) +
      obs_mat(_, 4) * (1 - draws(i, 13)) +
      obs_mat(_, 5) * (1 - draws(i, 14));
    e(_, i) = cpprbinom(d.length(), 1, pr);
    //    e_mat(_, i) = e;
    e_syst = obs_mat(_, 2) * draws(i, 17) +
      obs_mat(_, 3) * draws(i, 18) +
      obs_mat(_, 4) * (1 - draws(i, 19)) +
      obs_mat(_, 5) * (1 - draws(i, 20));
    arma::vec e_systarma = as<arma::vec>(e_syst);
    de_obs.col(0) = darma;
    de_obs.col(1) = e_systarma;
    d1 = de_obs.rows(find(outcome == 1));
    d0 = de_obs.rows(find(outcome == 0));
    at = as<double>(wrap(arma::accu(d1.col(1))));
    ct = as<double>(wrap(arma::accu(d0.col(1))));
    dt = ct - as<double>(wrap(d0.n_rows));
    bt = at - as<double>(wrap(d1.n_rows));
    // For systematic error only
    if (all(is_na(e_syst)).is_true()) {
      coef_syst = NA_REAL;
    } else {
      if ((at * bt * ct * dt) <= 0) {
	coef_syst = NA_REAL;
      } else {
	coef_syst = log(at * dt / bt / ct);
      }
    }

    ze_syst(i, 0) = R::rnorm(0, 1);
    ze_syst(i, 1) = coef_syst;
  }

  for (int j = 0; j < iter; j++) {
    ematrix(_, 1) = e(_, j);
    X = as<arma::mat>(ematrix);
    // Systematic and random error
    if (all(is_na(e(_, j))).is_true()) {
      coef = NA_REAL;
      mod_se = NA_REAL;
    } else {
      mod_log = f(Named("x", ematrix), Named("y", d), Named("family", "binomial"));
      mod_sum = s(Named("object", mod_log));
      mod_df = mod_sum["df"];
      mod_df = {mod_df[0], mod_df[1]};
      mod_coef = mod_log["coefficients"];
      coef = {mod_coef[1]};
      mod_res = mod_log["residuals"];
      mod_wght = mod_log["weights"];
      // Bread
      // Get unscaled covariance matrix
      mod_fit = mod_log["fitted.values"];
      mod_fit2 = 1 - mod_fit;
      W = arma::diagmat(as<arma::vec>(mod_fit) % as<arma::vec>(mod_fit2));
      unsccov = arma::solve(arma::trans(X) * W * X, I,
			    arma::solve_opts::allow_ugly);
      // bread as unscaled cov matrix * (ncol + nrow) * dispersion (but dispersion = 1)
      bread = unsccov * arma::accu(as<arma::vec>(mod_df));
      // Meat
      wres = mod_res * mod_wght;
      res = X.each_col() % as<arma::vec>(wres);
      res = res.col(0);
      // omega is res^2, but then we have to squared it for rval so no need to compute it
      rval = X.each_col() % res;
      meat = (rval.t() * rval) / d.length();
      sandwich = bread * meat * bread;
      sandwich = sandwich.each_col() % inv_nrow;
      se = sqrt(sandwich.diag());
      se = se(1);
      mod_se = as<double>(wrap(se));
    }

    res_mat(j, 0) = coef;
    res_mat(j, 1) = mod_se;
  }

  return Rcpp::List::create(Rcpp::Named("res_mat") = res_mat,
			    Rcpp::Named("ze_syst") = ze_syst);
}
