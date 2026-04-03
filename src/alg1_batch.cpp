// =============================================================================
// alg1_batch.cpp
//
// Algorithm 1 (Props 1-3) batch processing — handles M new units at once
//
// Mathematical equivalence:
//   Produces algebraically exact results as applying alg1_new_unit() to
//   each of the M new units sequentially.
//
// Single-pass structure (Pass 1 only):
//   Compute S_i, s_i, SSy_i, M_ss_i, vecS_i per unit, then aggregate:
//     ZtZ_new   = ZtZ_old + Σ S_i
//     inv_new   = solve(ZtZ_new)
//     theta_new = inv_new * (s_old + Σ s_i)
//     sigma2_new
//     M_ss_add  = Σ s_i s_i'   ← used for Vcr via formula
//
// Vcr is computed on the R side in otwfe_finalize():
//   M = M_ss - A_N(θ⊗I) - [A_N(θ⊗I)]' + (θ'⊗I) B_N (θ⊗I)
//   Vcr = inv_new %*% M %*% inv_new
// =============================================================================

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//' Algorithm 1 batch: process M new units at once
//'
//' @param inv_dotZtZ  p × p, pre-update (\dot{Z}'\dot{Z})^{-1}
//' @param theta_hat   length p, pre-update coefficient vector
//' @param sigma2_hat  pre-update sigma^2
//' @param N_old       pre-update number of units
//' @param n_old       pre-update total observation count
//' @param x_mat       n_stream × k, covariates of new units (stacked in unit order)
//' @param time_vec    length n_stream, re-indexed calendar time (1..T_support)
//' @param y_vec       length n_stream, dependent variable
//' @param unit_lens   length M, number of observations per unit (sum = n_stream)
//' @param T_support   size of the current calendar time support
//' @param baseline_time dropped baseline time dummy
//' @return list: inv_dotZtZ, theta_hat, sigma2_hat, M_ss_add,
//'               A_N_add, B_N_add, N_add, n_add
// [[Rcpp::export]]
List alg1_batch_cpp(
    const arma::mat&  inv_dotZtZ,
    const arma::vec&  theta_hat,
    double            sigma2_hat,
    int               N_old,
    int               n_old,
    const arma::mat&  x_mat,
    const arma::ivec& time_vec,
    const arma::vec&  y_vec,
    const arma::ivec& unit_lens,
    int               T_support,
    int               baseline_time
) {
  const int p  = (int)inv_dotZtZ.n_rows;
  const int k  = (int)x_mat.n_cols;
  const int M  = (int)unit_lens.n_elem;
  const int p2 = p * p;

  // Time dummy indices excluding the baseline time (ascending, 1-indexed)
  std::vector<int> dummy_times;
  dummy_times.reserve(T_support - 1);
  for (int t = 1; t <= T_support; t++) {
    if (t != baseline_time) dummy_times.push_back(t);
  }
  const int n_dummy = (int)dummy_times.size();  // = p - k

  // --------------------------------------------------------------------------
  // Pass 1: compute S_i, s_i, SSy_i, M_ss_i, vecS_i per unit → aggregate
  // --------------------------------------------------------------------------
  arma::mat total_S(p, p, fill::zeros);
  arma::vec total_s(p, fill::zeros);
  double    total_SSy = 0.0;
  arma::mat M_ss_add(p, p, fill::zeros);   // Σ s_i s_i'  (for Vcr formula)
  arma::mat A_N_add(p, p2, fill::zeros);   // Σ s_i vec(S_i)'
  arma::mat B_N_add(p2, p2, fill::zeros);  // Σ vec(S_i) vec(S_i)'

  int row_start = 0;
  for (int i = 0; i < M; i++) {
    const int Ti = unit_lens(i);

    // Build Z_raw_i (Ti × p) and y_i (Ti)
    arma::mat Z_raw_i(Ti, p);
    arma::vec y_i(Ti);

    for (int t = 0; t < Ti; t++) {
      const int obs   = row_start + t;
      const int t_val = time_vec(obs);
      y_i(t) = y_vec(obs);

      // Covariate columns
      for (int j = 0; j < k; j++)
        Z_raw_i(t, j) = x_mat(obs, j);

      // Time dummy columns
      for (int j = 0; j < n_dummy; j++)
        Z_raw_i(t, k + j) = (t_val == dummy_times[j]) ? 1.0 : 0.0;
    }

    // Within transformation: dotZ_i = Z_raw_i - barZ_i,  dotY_i = y_i - barY_i
    const arma::rowvec barZ_i = arma::mean(Z_raw_i, 0);
    const double       barY_i = arma::mean(y_i);
    const arma::mat    dotZ_i = Z_raw_i.each_row() - barZ_i;
    const arma::vec    dotY_i = y_i - barY_i;

    const arma::mat S_i   = dotZ_i.t() * dotZ_i;
    const arma::vec s_i   = dotZ_i.t() * dotY_i;
    const double    SSy_i = arma::dot(dotY_i, dotY_i);

    total_S   += S_i;
    total_s   += s_i;
    total_SSy += SSy_i;

    // M_ss: s_i s_i'  (for Vcr formula)
    M_ss_add += s_i * s_i.t();

    // A_N, B_N contribution: vec(S_i)
    const arma::vec vecS_i = arma::vectorise(S_i);
    A_N_add += s_i * vecS_i.t();
    B_N_add += vecS_i * vecS_i.t();

    row_start += Ti;
  }
  const int n_new = row_start;  // total streaming observation count

  // --------------------------------------------------------------------------
  // Compute theta_new, inv_new, sigma2_new
  // --------------------------------------------------------------------------
  // ZtZ_old = inv_dotZtZ^{-1},  s_old = ZtZ_old * theta_hat
  const arma::mat ZtZ_old = arma::inv_sympd(inv_dotZtZ);
  const arma::vec s_old   = ZtZ_old * theta_hat;

  // ZtZ_new = ZtZ_old + Σ S_i,  inv_new = ZtZ_new^{-1}
  const arma::mat ZtZ_new   = ZtZ_old + total_S;
  const arma::mat inv_new   = arma::inv_sympd(ZtZ_new);
  const arma::vec s_all     = s_old + total_s;
  const arma::vec theta_new = inv_new * s_all;

  // sigma2: recover SSy_warmup = sigma2_old * df_old + theta_old' ZtZ_old theta_old
  const int    df_old     = n_old - N_old - p;
  const int    df_new     = (n_old + n_new) - (N_old + M) - p;
  const double SSy_warmup = sigma2_hat * df_old + arma::dot(theta_hat, ZtZ_old * theta_hat);
  const double SSy_all    = SSy_warmup + total_SSy;
  const double RSS_new    = SSy_all - arma::dot(theta_new, s_all);
  const double sigma2_new = RSS_new / (double)df_new;

  // --------------------------------------------------------------------------
  // Return
  // --------------------------------------------------------------------------
  return List::create(
    Named("inv_dotZtZ")  = inv_new,
    Named("theta_hat")   = theta_new,
    Named("sigma2_hat")  = sigma2_new,
    Named("M_ss_add")    = M_ss_add,
    Named("A_N_add")     = A_N_add,
    Named("B_N_add")     = B_N_add,
    Named("N_add")       = M,
    Named("n_add")       = n_new
  );
}
