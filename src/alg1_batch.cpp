// =============================================================================
// alg1_batch.cpp
//
// Algorithm 1 (Prop 1-3) 배치 처리 — 새 unit M개를 한 번에 처리
//
// 수학적 등치:
//   M개의 새 unit을 순차적으로 alg1_new_unit() 적용한 결과와
//   algebraically exact하게 동일.
//
// 단일 패스 구조 (Pass 1만):
//   각 unit의 S_i, s_i, SSy_i, M_ss_i, vecS_i 계산 후 합산
//     ZtZ_new  = ZtZ_old + Σ S_i
//     inv_new  = solve(ZtZ_new)
//     theta_new = inv_new * (s_old + Σ s_i)
//     sigma2_new
//     M_ss_add = Σ s_i s_i'   ← Vcr 계산에 사용 (formula 방식)
//
// Vcr 계산은 R-side otwfe_finalize()에서 수행:
//   M = M_ss - A_N(θ⊗I) - [A_N(θ⊗I)]' + (θ'⊗I) B_N (θ⊗I)
//   Vcr = inv_new %*% M %*% inv_new
// =============================================================================

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//' Algorithm 1 배치: 새 unit M개를 한 번에 처리
//'
//' @param inv_dotZtZ  p × p, 업데이트 전 (\dot{Z}'\dot{Z})^{-1}
//' @param theta_hat   p, 업데이트 전 계수 벡터
//' @param sigma2_hat  업데이트 전 sigma^2
//' @param N_old       업데이트 전 unit 수
//' @param n_old       업데이트 전 총 관측치 수
//' @param x_mat       n_stream × k, 새 unit들의 공변량 (unit 순서대로 연결)
//' @param time_vec    n_stream, 재인덱싱된 calendar time (1..T_support)
//' @param y_vec       n_stream, 종속변수
//' @param unit_lens   M, 각 unit의 관측치 수 (합 = n_stream)
//' @param T_support   현재 calendar time support 크기
//' @param baseline_time 제거한 기준 time dummy
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

  // 기준 time을 제외한 time dummy 인덱스 (오름차순, 1-indexed)
  std::vector<int> dummy_times;
  dummy_times.reserve(T_support - 1);
  for (int t = 1; t <= T_support; t++) {
    if (t != baseline_time) dummy_times.push_back(t);
  }
  const int n_dummy = (int)dummy_times.size();  // = p - k

  // --------------------------------------------------------------------------
  // Pass 1: S_i, s_i, SSy_i, M_ss_i, vecS_i 계산 → 집계
  // --------------------------------------------------------------------------
  arma::mat total_S(p, p, fill::zeros);
  arma::vec total_s(p, fill::zeros);
  double    total_SSy = 0.0;
  arma::mat M_ss_add(p, p, fill::zeros);   // Σ s_i s_i'  (Vcr formula용)
  arma::mat A_N_add(p, p2, fill::zeros);   // Σ s_i vec(S_i)'
  arma::mat B_N_add(p2, p2, fill::zeros);  // Σ vec(S_i) vec(S_i)'

  int row_start = 0;
  for (int i = 0; i < M; i++) {
    const int Ti = unit_lens(i);

    // Z_raw_i (Ti × p), y_i (Ti) 구성
    arma::mat Z_raw_i(Ti, p);
    arma::vec y_i(Ti);

    for (int t = 0; t < Ti; t++) {
      const int obs   = row_start + t;
      const int t_val = time_vec(obs);
      y_i(t) = y_vec(obs);

      // 공변량 열
      for (int j = 0; j < k; j++)
        Z_raw_i(t, j) = x_mat(obs, j);

      // time dummy 열
      for (int j = 0; j < n_dummy; j++)
        Z_raw_i(t, k + j) = (t_val == dummy_times[j]) ? 1.0 : 0.0;
    }

    // Within 변환: dotZ_i = Z_raw_i - barZ_i,  dotY_i = y_i - barY_i
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

    // M_ss: s_i s_i'  (Vcr formula용)
    M_ss_add += s_i * s_i.t();

    // A_N, B_N 기여: vec(S_i)
    const arma::vec vecS_i = arma::vectorise(S_i);
    A_N_add += s_i * vecS_i.t();
    B_N_add += vecS_i * vecS_i.t();

    row_start += Ti;
  }
  const int n_new = row_start;  // 총 스트리밍 관측치 수

  // --------------------------------------------------------------------------
  // theta_new, inv_new, sigma2_new 계산
  // --------------------------------------------------------------------------
  // ZtZ_old = inv_dotZtZ^{-1},  s_old = ZtZ_old * theta_hat
  const arma::mat ZtZ_old = arma::inv_sympd(inv_dotZtZ);
  const arma::vec s_old   = ZtZ_old * theta_hat;

  // ZtZ_new = ZtZ_old + Σ S_i,  inv_new = ZtZ_new^{-1}
  const arma::mat ZtZ_new   = ZtZ_old + total_S;
  const arma::mat inv_new   = arma::inv_sympd(ZtZ_new);
  const arma::vec s_all     = s_old + total_s;
  const arma::vec theta_new = inv_new * s_all;

  // sigma2: SSy_warmup = sigma2_old * df_old + theta_old' ZtZ_old theta_old
  const int    df_old     = n_old - N_old - p;
  const int    df_new     = (n_old + n_new) - (N_old + M) - p;
  const double SSy_warmup = sigma2_hat * df_old + arma::dot(theta_hat, ZtZ_old * theta_hat);
  const double SSy_all    = SSy_warmup + total_SSy;
  const double RSS_new    = SSy_all - arma::dot(theta_new, s_all);
  const double sigma2_new = RSS_new / (double)df_new;

  // --------------------------------------------------------------------------
  // 반환
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
