#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>
// [[Rcpp::export]]
NumericVector loop_misclass(int reps, int n, CharacterVector measure,
			    NumericMatrix obs_mat, NumericMatrix draws,
			    Formula formula,
			    bool display_progress=true) {
  // Obtain environment containing function
  Rcpp::Environment base("package:stats");

  // Picking up glm() function from base stats
  Rcpp::Function glm_r = base["glm"];
  Rcpp::Function sum_r = base["summary.glm"];

  Rcpp::NumericVector d = obs_mat(_, 1);
  Rcpp::NumericVector p = no_init(n);
  Rcpp::NumericVector e = no_init(n);
  Rcpp::NumericVector mod_coef = no_init(n);
  //Rcpp::Formula f("d ~ e");
  //  Rcpp::NumericVector mod_se = no_init(n);
  Progress b(reps*reps, display_progress);
  for (int i = 0; i < reps; i++) {
    b.increment(); // update progress
    obs_mat(_, 6) = obs_mat(_, 2) * draws(i, 11) +
      obs_mat(_, 3) * draws(i, 12) +
      obs_mat(_, 4) * (1 - draws(i, 13)) +
      obs_mat(_, 5) * (1 - draws(i, 14));
    p = obs_mat(_, 6);
    std::transform(p.begin(), p.end(), e.begin(), [=](double p){return R::rbinom(1, p);});
    e = as<IntegerVector>(e);
    obs_mat(_, 7) = e;
    obs_mat(_, 8) = obs_mat(_, 2) * draws(i, 17) +
      obs_mat(_, 3) * draws(i, 18) +
      obs_mat(_, 4) * (1 - draws(i, 19)) +
      obs_mat(_, 5) * (1 - draws(i, 20));
    //DataFrame df = DataFrame::create(Named("d") = d, _["e"] = e);
    //Rcpp::NumericMatrix M(d.length(), 2);
    //M(_, 0) = d;
    //M(_, 1) = e;
    if (all(is_na(e)).is_true()) {
      NumericVector mod_coef(n, NumericVector::get_na());
      NumericVector mod_se(n, NumericVector::get_na());
    } else {
      if (measure = "RR") {
	Rcpp::List mod_pois = glm_r(_["formula"] = formula,
				    _["family"] = "poisson");
	Rcpp::List mod_sum = sum_r(mod_pois);
	Rcpp::NumericMatrix M_coef = mod_sum[12];
	mod_coef = M_coef(2, 1);
	//	mod_cov = R::sandwich::vcovHC(mod_pois, type = "HC0");
	//mod_se = R::sqrt(diag(mod_cov))[2];
      }
      if (measure = "OR") {
	mod_coef = 100;
      }
    }
  }
  return mod_coef;
}
