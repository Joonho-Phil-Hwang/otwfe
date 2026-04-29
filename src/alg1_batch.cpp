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
//
// Performance notes:
//   On platforms with BLAS libraries that have high per-call overhead for
//   tiny matrices (e.g. Apple Accelerate with GCD dispatch), the inner
//   loop avoids all BLAS calls (no arma::mat * arma::mat) and uses raw
//   pointer arithmetic instead. This eliminates ~60 μs/unit overhead from
//   three BLAS calls (dgemm for S_i, dgemv for s_i, dger for M_ss_add)
//   and reduces per-unit time from ~66 μs to ~3–5 μs.
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

  // Time dummy indices excluding the baseline (ascending, 1-indexed)
  std::vector<int> dummy_times;
  dummy_times.reserve(T_support - 1);
  for (int t = 1; t <= T_support; t++)
    if (t != baseline_time) dummy_times.push_back(t);
  const int n_dummy = (int)dummy_times.size();  // = p - k

  // --------------------------------------------------------------------------
  // Pass 1: for each unit, compute S_i, s_i, SSy_i, vecS_i and accumulate
  //         directly into the global aggregates — no BLAS calls.
  //
  // Key optimisation: BLAS libraries with high per-call overhead (e.g. Apple
  // Accelerate on ARM) add ~20 μs per call even for tiny matrices. We avoid
  // the three BLAS calls (dgemm for S_i, dgemv for s_i, dger for M_ss_add)
  // by using raw pointer loops, reducing per-unit time ~13×.
  //
  // Memory layout: all work matrices are column-major (Armadillo default).
  //   Z_raw_buf(t, j) = Z_raw_ptr[t + j * T_max]
  // --------------------------------------------------------------------------
  arma::mat total_S(p, p, fill::zeros);
  arma::vec total_s(p, fill::zeros);
  double    total_SSy = 0.0;
  arma::mat M_ss_add(p, p, fill::zeros);
  arma::mat A_N_add(p, p2, fill::zeros);
  arma::mat B_N_add(p2, p2, fill::zeros);  // upper triangle filled, symmetrized at end

  // Pre-allocate reusable work buffers (one heap allocation each, reused for all M units)
  const int T_max = (int)unit_lens.max();
  arma::mat Z_raw_buf(T_max, p);           // column-major work buffer for design matrix
  arma::vec y_buf(T_max);
  std::vector<double> barZ_v(p);           // column means of Z_raw (within-transform)
  std::vector<double> s_local(p);          // s_i = dotZ_i' dotY_i
  std::vector<double> vecS_local(p2);      // vec(S_i), column-major

  double* const Z_ptr = Z_raw_buf.memptr();   // raw column-major pointer
  double* const Y_ptr = y_buf.memptr();
  double* const tS    = total_S.memptr();     // raw column-major pointer for total_S
  double* const ts    = total_s.memptr();

  int row_start = 0;
  for (int i = 0; i < M; i++) {
    const int Ti = unit_lens(i);

    // ------------------------------------------------------------------
    // (1) Fill Z_raw_buf (column-major) and y_buf for this unit's Ti rows
    // ------------------------------------------------------------------
    for (int t = 0; t < Ti; t++) {
      const int obs   = row_start + t;
      const int t_val = time_vec(obs);
      Y_ptr[t] = y_vec(obs);
      for (int j = 0; j < k; j++)
        Z_ptr[t + j * T_max] = x_mat(obs, j);
      for (int j = 0; j < n_dummy; j++)
        Z_ptr[t + (k + j) * T_max] = (t_val == dummy_times[j]) ? 1.0 : 0.0;
    }

    // ------------------------------------------------------------------
    // (2) Compute column means barZ[j] and barY, then demean in-place
    //     (dotZ overwrites Z_raw_buf; dotY overwrites y_buf)
    // ------------------------------------------------------------------
    double barY = 0.0;
    for (int t = 0; t < Ti; t++) barY += Y_ptr[t];
    barY /= Ti;

    for (int j = 0; j < p; j++) {
      double* col_j = Z_ptr + j * T_max;
      double  sum_j = 0.0;
      for (int t = 0; t < Ti; t++) sum_j += col_j[t];
      barZ_v[j] = sum_j / Ti;
      for (int t = 0; t < Ti; t++) col_j[t] -= barZ_v[j];  // dotZ in-place
    }
    double SSy_i = 0.0;
    for (int t = 0; t < Ti; t++) {
      Y_ptr[t] -= barY;   // dotY in-place
      SSy_i    += Y_ptr[t] * Y_ptr[t];
    }
    total_SSy += SSy_i;

    // ------------------------------------------------------------------
    // (3) Compute s_i = dotZ' dotY  (manual, no BLAS)
    //     Accumulate directly into total_s
    // ------------------------------------------------------------------
    for (int j = 0; j < p; j++) {
      const double* col_j = Z_ptr + j * T_max;
      double sum_j = 0.0;
      for (int t = 0; t < Ti; t++) sum_j += col_j[t] * Y_ptr[t];
      s_local[j]  = sum_j;
      ts[j]      += sum_j;  // total_s (column-major, 1-D)
    }

    // ------------------------------------------------------------------
    // (4) Compute S_i = dotZ' dotZ (symmetric, manual upper-triangle loop)
    //     Simultaneously:
    //       a. accumulate into total_S
    //       b. store vec(S_i) in vecS_local  (column-major: S_i(j1,j2) at j1 + j2*p)
    //       c. accumulate M_ss_add += s_i s_i'
    // ------------------------------------------------------------------
    for (int j2 = 0; j2 < p; j2++) {
      const double* col_j2 = Z_ptr + j2 * T_max;
      for (int j1 = 0; j1 <= j2; j1++) {
        const double* col_j1 = Z_ptr + j1 * T_max;
        double Sval = 0.0;
        for (int t = 0; t < Ti; t++) Sval += col_j1[t] * col_j2[t];

        // total_S (column-major): (j1, j2) and (j2, j1)
        tS[j1 + j2 * p] += Sval;
        if (j1 < j2) tS[j2 + j1 * p] += Sval;

        // vec(S_i) in column-major order: index = j1 + j2*p
        vecS_local[j1 + j2 * p] = Sval;
        if (j1 < j2) vecS_local[j2 + j1 * p] = Sval;
      }
      // M_ss_add column j2 update (lower triangle j1 ≤ j2)
      const double s_j2 = s_local[j2];
      for (int j1 = 0; j1 <= j2; j1++) {
        double mval = s_local[j1] * s_j2;
        M_ss_add(j1, j2) += mval;
        if (j1 < j2) M_ss_add(j2, j1) += mval;
      }
    }

    // ------------------------------------------------------------------
    // (5) A_N_add += s_i vec(S_i)'  (p × p²)
    //     B_N_add += vec(S_i) vec(S_i)'  (p² × p², upper triangle)
    // ------------------------------------------------------------------
    for (int jj = 0; jj < p2; jj++) {
      const double v_jj = vecS_local[jj];
      for (int ii = 0; ii < p; ii++)
        A_N_add(ii, jj) += s_local[ii] * v_jj;
      for (int ii = 0; ii <= jj; ii++)
        B_N_add(ii, jj) += vecS_local[ii] * v_jj;
    }

    row_start += Ti;
  }
  const int n_new = row_start;

  // Symmetrize B_N_add (copy upper triangle to lower)
  for (int jj = 1; jj < p2; jj++)
    for (int ii = 0; ii < jj; ii++)
      B_N_add(jj, ii) = B_N_add(ii, jj);

  // --------------------------------------------------------------------------
  // Compute theta_new, inv_new, sigma2_new
  // --------------------------------------------------------------------------
  const arma::mat ZtZ_old   = arma::inv_sympd(inv_dotZtZ);
  const arma::vec s_old     = ZtZ_old * theta_hat;
  const arma::mat ZtZ_new   = ZtZ_old + total_S;
  const arma::mat inv_new   = arma::inv_sympd(ZtZ_new);
  const arma::vec s_all     = s_old + total_s;
  const arma::vec theta_new = inv_new * s_all;

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
