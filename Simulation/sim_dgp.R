# =============================================================================
# sim_dgp.R
#
# 논문: "Online Updating for Linear Panel Regressions" (Hwang & Lee, 2026)
# DGP: Subsection 5.2 Simulation
#
# 모형: Y_it = X_it' β + α_i + λ_t + u_it
#   β = (1.3, -0.7)'
#   α_i  ~ N(0, σ_α²),  σ_α = 1
#   λ_t  ~ N(0, σ_λ²),  σ_λ = 1
#   u_it ~ N(0, σ_u²),  σ_u = 1
#   X_it ~ N(0, I_2),   k = 2
#   T_i  ~ Uniform{T_min, ..., T},  T = 5, T_min = 2
#   S_i  ⊆ {1,...,T}: 크기 T_i의 비복원 랜덤 부분집합
# 기술적 조건: 모든 t ∈ {1,...,T}가 전체 표본에 최소 1회 등장
# =============================================================================

#' 패널 데이터 생성 — Section 5.2 DGP
#'
#' @param N       unit 수
#' @param T       calendar time 지지(support) 크기  (default 5)
#' @param T_min   unit당 최소 관측 time 수          (default 2)
#' @param sigma_a 개인 고정효과 표준편차             (default 1)
#' @param sigma_l 시간 고정효과 표준편차             (default 1)
#' @param sigma_u 오차항 표준편차                   (default 1)
#' @param seed    랜덤 시드
#' @param verbose 진행 상황 출력 여부
#' @return list:
#'   \item{df}        data.frame (id, time, y, x1, x2)
#'   \item{t_gen}     데이터 생성 소요 시간 (초)
#'   \item{N_obs}     총 관측치 수
#'   \item{T_i_vec}   unit별 관측 time 수 벡터 (길이 N)
#'   \item{beta_true} 참 계수 벡터 c(1.3, -0.7)
#'   \item{lambda_t}  시간 고정효과 실현값 (길이 T)
generate_panel_52 <- function(N,
                               T       = 5L,
                               T_min   = 2L,
                               sigma_a = 1,
                               sigma_l = 1,
                               sigma_u = 1,
                               seed    = 42L,
                               verbose = TRUE) {

  if (verbose)
    cat(sprintf("  DGP 생성: N = %s,  T = %d ...",
                format(N, big.mark = ","), T))
  t0 <- proc.time()

  set.seed(seed)

  beta     <- c(1.3, -0.7)                              # 참 계수
  alpha_i  <- rnorm(N,   sd = sigma_a)                  # α_i ~ N(0, σ_α²)
  lambda_t <- rnorm(T,   sd = sigma_l)                  # λ_t ~ N(0, σ_λ²)
  T_i_vec  <- sample(T_min:T, N, replace = TRUE)        # T_i ~ Uniform{T_min,...,T}
  total_obs <- sum(T_i_vec)

  # ------------------------------------------------------------------
  # id_vec 벡터화 (O(n))
  # ------------------------------------------------------------------
  id_vec <- rep(seq_len(N), T_i_vec)

  # ------------------------------------------------------------------
  # time_vec: 각 unit i에서 {1,...,T}의 크기 T_i 비복원 부분집합 S_i
  # sample.int(T, Ti) 호출이 N번 필요 → loop 불가피, 사전 할당으로 최적화
  # ------------------------------------------------------------------
  time_vec <- integer(total_obs)
  pos <- 1L
  for (i in seq_len(N)) {
    Ti  <- T_i_vec[i]
    end <- pos + Ti - 1L
    time_vec[pos:end] <- sample.int(T, Ti)   # 비복원 추출
    pos <- end + 1L
  }

  # ------------------------------------------------------------------
  # 기술적 조건: 모든 t ∈ {1,...,T}가 최소 1회 등장하도록 보장
  # (N이 크면 사실상 자동 만족되지만 명시적으로 확인·수정)
  # ------------------------------------------------------------------
  missing_t <- setdiff(seq_len(T), unique(time_vec))
  if (length(missing_t) > 0L) {
    warning(sprintf(
      "누락된 calendar time (%s) 강제 커버 — 매우 소규모 N에서만 발생",
      paste(missing_t, collapse = ",")))
    cum_end <- cumsum(T_i_vec)
    cum_beg <- c(1L, cum_end[-N] + 1L)
    for (t_miss in missing_t) {
      # T_i >= 2인 unit 중 처음으로 나오는 unit의 첫 관측치를 t_miss로 교체
      for (ii in seq_len(N)) {
        if (T_i_vec[ii] >= 2L) {
          time_vec[cum_beg[ii]] <- t_miss
          break
        }
      }
    }
  }

  # ------------------------------------------------------------------
  # X, u 벡터화 (fast)
  # ------------------------------------------------------------------
  x1_vec <- rnorm(total_obs)
  x2_vec <- rnorm(total_obs)
  u_vec  <- rnorm(total_obs, sd = sigma_u)

  # ------------------------------------------------------------------
  # y 벡터화: Y_it = α_i + λ_t + 1.3 x1 - 0.7 x2 + u
  # ------------------------------------------------------------------
  y_vec <- rep(alpha_i, T_i_vec) +
           lambda_t[time_vec]    +
           beta[1L] * x1_vec    +
           beta[2L] * x2_vec    +
           u_vec

  t_gen <- (proc.time() - t0)[["elapsed"]]

  if (verbose)
    cat(sprintf("  완료 %.1f초  (관측치 %s,  %.0f MB)\n",
                t_gen,
                format(total_obs, big.mark = ","),
                object.size(data.frame(
                  id=id_vec, time=time_vec, y=y_vec, x1=x1_vec, x2=x2_vec
                )) / 1e6))

  list(
    df        = data.frame(id   = id_vec,
                           time = time_vec,
                           y    = y_vec,
                           x1   = x1_vec,
                           x2   = x2_vec,
                           stringsAsFactors = FALSE),
    t_gen     = t_gen,
    N_obs     = total_obs,
    T_i_vec   = T_i_vec,
    beta_true = beta,
    lambda_t  = lambda_t
  )
}
