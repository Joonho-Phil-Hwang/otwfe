# =============================================================================
# online_twfe_core.R
#
# 논문: "Online Updating for Linear Panel Regressions" (Hwang & Lee, 2026)
# 핵심 알고리즘 구현 — 패키지 공개용 일반 함수 모음
#
# 구성:
#   Section 1 : 유틸리티
#   Section 2 : 설계 행렬 & within 변환
#   Section 3 : State 초기화 및 관리
#   Section 4 : Algorithm 1 / 2 / 3  (Props 1-9)
#   Section 5 : 메인 API  (otwfe)
#   Section 6 : Output  (tidy / summary / print)
# =============================================================================

# Rcpp 배치 알고리즘 (Algorithm 1 벡터화) 로드
# — src/alg1_batch.cpp 컴파일 실패 시 pure-R fallback 사용
# — 현재 디렉토리 및 상위 디렉토리까지 탐색 (Simulation/, uni/ 등 하위 폴더에서 실행 시 대응)
.alg1_batch_rcpp_available <- FALSE
tryCatch({
  .cpp_candidates <- c(
    file.path(getwd(), "src", "alg1_batch.cpp"),              # 현재 디렉토리
    file.path(dirname(getwd()), "src", "alg1_batch.cpp")      # 한 단계 위
  )
  .cpp_found <- Filter(file.exists, .cpp_candidates)
  if (length(.cpp_found) == 0L) stop("src/alg1_batch.cpp를 찾을 수 없습니다.")
  Rcpp::sourceCpp(.cpp_found[[1L]])
  .alg1_batch_rcpp_available <- TRUE
}, error = function(e) {
  message("[otwfe] Rcpp 배치 모듈 로드 실패 — pure-R fallback 사용\n  ", conditionMessage(e))
})


# =============================================================================
# Section 1: 유틸리티
# =============================================================================

#' 행렬을 열 우선(column-major)으로 벡터화
vec <- function(M) as.vector(M)

#' vec() 인덱스 임베딩: p_old × p_old 행렬이 (p_old+1) × (p_old+1)로 확장될 때
#' vec(M_old)의 각 원소가 vec(M_new)의 몇 번째 위치로 가는지 반환
#'
#' vec()은 열 우선 순서이므로 열 c, 행 r의 원소는
#' 확장된 행렬에서  r + (c-1) * p_new  번째 위치로 이동
vec_embed_index <- function(p_old) {
  p_new <- p_old + 1L
  idx   <- integer(p_old * p_old)
  k     <- 0L
  for (c in seq_len(p_old)) {
    for (r in seq_len(p_old)) {
      k      <- k + 1L
      idx[k] <- r + (c - 1L) * p_new
    }
  }
  idx
}

#' 리스트 내 객체들의 메모리 합계 (bytes)
sum_object_size <- function(x) {
  sum(vapply(x, function(z) as.numeric(object.size(z)), numeric(1)))
}

#' 표현식을 실행하고 소요 시간과 결과를 함께 반환
safe_system_time <- function(expr, label = NULL, verbose_error = FALSE) {
  res <- NULL
  tt  <- system.time({
    res <- tryCatch(
      list(ok = TRUE,  value = eval.parent(substitute(expr)), error = NULL),
      error = function(e) list(ok = FALSE, value = NULL, error = e)
    )
  })
  out <- list(
    ok    = isTRUE(res$ok),
    time  = unname(tt[["elapsed"]]),
    value = res$value,
    error = res$error,
    label = label
  )
  if (!isTRUE(out$ok) && isTRUE(verbose_error)) {
    msg <- conditionMessage(out$error)
    cat(if (!is.null(label)) sprintf("  [ERROR] %s: %s\n", label, msg)
        else                 sprintf("  [ERROR] %s\n", msg))
  }
  out
}

# 진행 표시줄 (대형 N에서 integer overflow 방지를 위해 numeric 사용)
make_progress_bar <- function(min_val, max_val, style = 3) {
  txtProgressBar(min = as.numeric(min_val), max = as.numeric(max_val), style = style)
}
safe_set_progress <- function(pb, value) {
  if (is.null(pb)) return(invisible(NULL))
  tryCatch(suppressWarnings(setTxtProgressBar(pb, as.numeric(value))),
           error = function(e) invisible(NULL))
  invisible(NULL)
}


# =============================================================================
# Section 2: 설계 행렬 & within 변환
# =============================================================================

#' TWFE 설계 행렬 Z_it = [X_it, time dummies] 구성
#'
#' @param x_mat   n × k 공변량 행렬 (이미 이름이 붙어 있어야 함)
#' @param time_vec 길이 n의 정수형 calendar time 벡터
#' @param T_support 현재 sample의 최대 calendar time
#' @param baseline_time 제거할 기준 time dummy (default = 1)
#' @return n × (k + T_support - 1) 행렬
build_Z_raw <- function(x_mat, time_vec, T_support, baseline_time = 1L) {
  # baseline_time을 제외한 모든 calendar time에 대해 dummy 자동 생성
  dummy_times  <- setdiff(seq_len(T_support), baseline_time)
  dummy_names  <- paste0("factor(time)", dummy_times)

  # outer()로 한 번에 dummy 행렬 생성
  D_mat        <- outer(time_vec, dummy_times, `==`) + 0L
  colnames(D_mat) <- dummy_names

  cbind(x_mat, D_mat)
}

#' individual i의 within 변환
#'
#' \dot{Z}_{it} = Z_{it} - \bar{Z}_i,  \dot{Y}_{it} = Y_{it} - \bar{Y}_i
#'
#' @param Z_raw_i T_i × p 행렬
#' @param y_i     길이 T_i 벡터
#' @return list(dotZ_i, dotY_i, barZ_i, barY_i)
within_transform_unit <- function(Z_raw_i, y_i) {
  barZ_i <- colMeans(Z_raw_i)
  barY_i <- mean(y_i)
  dotZ_i <- sweep(Z_raw_i, 2, barZ_i, "-")
  dotY_i <- y_i - barY_i
  list(dotZ_i = dotZ_i, dotY_i = dotY_i, barZ_i = barZ_i, barY_i = barY_i)
}


# =============================================================================
# Section 3: State 초기화 및 관리
# =============================================================================

#' unit별 요약 객체 생성 (내부용)
#'
#' @param id      unit 식별자
#' @param T_i     관측 시점 수
#' @param barZ_i  within 평균 벡터 (길이 p)
#' @param barY_i  within 평균 스칼라
#' @param S_i     \dot{Z}_i' \dot{Z}_i  (p × p)
#' @param s_i     \dot{Z}_i' \dot{Y}_i  (p × 1 행렬)
make_unit_state <- function(id, T_i, barZ_i, barY_i, S_i, s_i) {
  list(id = id, T_i = T_i, barZ_i = barZ_i, barY_i = barY_i, S_i = S_i, s_i = s_i)
}

#' warm-up 데이터프레임으로부터 초기 state 구성
#'
#' @param df_warm     warm-up용 data.frame
#' @param id_col      unit ID 컬럼명
#' @param time_col    calendar time 컬럼명
#' @param y_col       종속변수 컬럼명
#' @param x_cols      설명변수 컬럼명 벡터
#' @param baseline_time 제거할 기준 time dummy (NULL이면 데이터의 최솟값)
#' @return list(S_pre = state 객체,  units = named list of unit 요약)
state_init <- function(df_warm, id_col, time_col, y_col, x_cols,
                       baseline_time = NULL) {

  id_vec   <- df_warm[[id_col]]
  time_vec <- as.integer(df_warm[[time_col]])
  y_vec    <- df_warm[[y_col]]
  x_mat    <- as.matrix(df_warm[, x_cols, drop = FALSE])

  # baseline_time 자동 결정
  if (is.null(baseline_time)) baseline_time <- min(time_vec)
  T_support <- max(time_vec)

  n <- nrow(df_warm)
  Z_raw <- build_Z_raw(x_mat, time_vec, T_support, baseline_time)
  p     <- ncol(Z_raw)

  dotZ  <- matrix(0, nrow = n, ncol = p, dimnames = list(NULL, colnames(Z_raw)))
  dotY  <- numeric(n)

  uid_vec  <- unique(id_vec)
  uid_ch_vec <- as.character(uid_vec)
  units    <- vector("list", length(uid_vec))
  names(units) <- uid_ch_vec

  # 인덱스를 unit별로 한 번에 분할 (O(n)) — which() 반복 루프 O(n²) 방지
  idx_list <- split(seq_along(id_vec), as.character(id_vec))

  for (uid_ch in uid_ch_vec) {
    idx_i <- idx_list[[uid_ch]]
    wt_i  <- within_transform_unit(Z_raw[idx_i, , drop = FALSE], y_vec[idx_i])

    dotZ[idx_i, ] <- wt_i$dotZ_i
    dotY[idx_i]   <- wt_i$dotY_i

    S_i <- crossprod(wt_i$dotZ_i)
    s_i <- crossprod(wt_i$dotZ_i, matrix(wt_i$dotY_i, ncol = 1))

    units[[uid_ch]] <- make_unit_state(
      id     = uid_ch,
      T_i    = length(idx_i),
      barZ_i = wt_i$barZ_i,
      barY_i = wt_i$barY_i,
      S_i    = S_i,
      s_i    = s_i
    )
  }

  # Alg 1/2에서 cluster-robust 업데이트에 필요한 집계 행렬
  # A_N = \sum_i \dot{Z}_i'\dot{Y}_i  vec(\dot{Z}_i'\dot{Z}_i)'   (p × p^2)
  # B_N = \sum_i vec(\dot{Z}_i'\dot{Z}_i) vec(\dot{Z}_i'\dot{Z}_i)'  (p^2 × p^2)
  A_N <- matrix(0, nrow = p, ncol = p^2)
  B_N <- matrix(0, nrow = p^2, ncol = p^2)
  for (u in units) {
    vecS_i <- vec(u$S_i)
    A_N    <- A_N + u$s_i %*% vecS_i
    B_N    <- B_N + tcrossprod(vecS_i)
  }

  dotZtZ     <- crossprod(dotZ)
  inv_dotZtZ <- solve(dotZtZ)
  theta_hat  <- drop(inv_dotZtZ %*% crossprod(dotZ, dotY))
  names(theta_hat) <- colnames(Z_raw)

  eps_hat    <- drop(dotY - dotZ %*% theta_hat)
  N_pre      <- length(uid_vec)
  df_pre     <- n - N_pre - p
  sigma2_hat <- sum(eps_hat^2) / df_pre

  # M = \sum_i r_i r_i',  r_i = \dot{Z}_i' \hat{e}_i  (cluster score 외적 합)
  # idx_list 재사용 — O(n²) 루프 방지
  M_hat <- matrix(0, nrow = p, ncol = p,
                  dimnames = list(colnames(Z_raw), colnames(Z_raw)))
  for (uid_ch in uid_ch_vec) {
    idx_i <- idx_list[[uid_ch]]
    r_i   <- crossprod(dotZ[idx_i, , drop = FALSE],
                       matrix(eps_hat[idx_i], ncol = 1))
    M_hat <- M_hat + tcrossprod(r_i)
  }
  Vcr_hat <- inv_dotZtZ %*% M_hat %*% inv_dotZtZ

  S_pre <- list(
    # 메타데이터
    x_cols       = x_cols,
    baseline_time = baseline_time,
    # 차원 정보
    T_support    = T_support,
    p            = p,
    N            = N_pre,
    n            = n,
    # 핵심 추정량
    inv_dotZtZ   = inv_dotZtZ,
    theta_hat    = theta_hat,
    sigma2_hat   = sigma2_hat,
    # 분산 추정량
    Vcr_hat      = Vcr_hat,
    M_hat        = M_hat,          # Alg 3에서 사용; Alg 1/2 후엔 NULL로 초기화
    # cluster-robust 업데이트용 집계 행렬 (Prop 3/6)
    A_N          = A_N,
    B_N          = B_N,
    # 청크 방식 Vcr formula용: Σ s_i s_i' (스트리밍 unit 누적)
    M_ss         = matrix(0, nrow = p, ncol = p)
  )

  list(S_pre = S_pre, units = units)
}

#' 새로운 time dummy 도입 후 unit 요약 객체를 zero-extend
#'
#' Algorithm 3 실행 후 나머지 unit들의 S_i, s_i, barZ_i를
#' (p+1) 차원으로 확장 (새 dummy는 pre-update에서 항상 0)
#'
#' @param unit         unit 요약 객체
#' @param new_dummy_name 새로 추가되는 time dummy 이름
#' @return 차원이 확장된 unit 요약 객체
expand_unit_state_zero <- function(unit, new_dummy_name) {
  p_old <- nrow(unit$S_i)
  nms   <- c(rownames(unit$S_i), new_dummy_name)

  S_ext <- rbind(cbind(unit$S_i, 0), 0)
  s_ext <- rbind(unit$s_i, 0)
  dimnames(S_ext) <- list(nms, nms)
  rownames(s_ext) <- nms

  unit$barZ_i <- c(unit$barZ_i, 0)
  names(unit$barZ_i)[p_old + 1L] <- new_dummy_name
  unit$S_i <- S_ext
  unit$s_i <- s_ext
  unit
}


# =============================================================================
# Section 4: 알고리즘
# =============================================================================

# -----------------------------------------------------------------------------
# Algorithm 1 (Props 1–3): 새로운 unit 도착
#
# 입력  S_pre  : 업데이트 전 state
#       x_mat_i: T_i × k 공변량 행렬 (새 unit)
#       time_vec_i: 길이 T_i의 calendar time 벡터
#       y_i    : 길이 T_i의 종속변수 벡터
# 반환  S_post : 업데이트 후 state
# -----------------------------------------------------------------------------
alg1_new_unit <- function(S_pre, x_mat_i, time_vec_i, y_i) {

  stopifnot(all(time_vec_i %in% seq_len(S_pre$T_support)))

  Z_raw_i <- build_Z_raw(x_mat_i, as.integer(time_vec_i),
                         S_pre$T_support, S_pre$baseline_time)
  wt_i    <- within_transform_unit(Z_raw_i, y_i)
  dotZ_i  <- wt_i$dotZ_i
  dotY_i  <- wt_i$dotY_i

  inv_pre   <- S_pre$inv_dotZtZ
  theta_pre <- S_pre$theta_hat
  sig2_pre  <- S_pre$sigma2_hat
  Vcr_pre   <- S_pre$Vcr_hat
  A_pre     <- S_pre$A_N
  B_pre     <- S_pre$B_N
  p         <- S_pre$p
  T_i       <- length(y_i)

  # H_i = I_{T_i} + \dot{Z}_i (\dot{Z}'\dot{Z})^{-1} \dot{Z}_i'
  # solve(H_i, [dotZ_i, e_tilde]) 방식으로 H_i^{-1} 명시적 구성 회피
  H_i       <- diag(T_i) + dotZ_i %*% inv_pre %*% t(dotZ_i)
  e_tilde_i <- drop(dotY_i - dotZ_i %*% theta_pre)
  sol       <- solve(H_i, cbind(dotZ_i, e_tilde_i))
  W_i       <- sol[, seq_len(p), drop = FALSE]   # H_i^{-1} \dot{Z}_i
  v_i       <- sol[, p + 1L]                     # H_i^{-1} \tilde{e}_i

  # --- Prop 1: 계수 업데이트
  # \hat{\beta}^* = \hat{\beta} + (\dot{X}'\dot{X})^{-1} \dot{X}_{N+1}' H^{-1} \tilde{e}_{N+1}
  q_i            <- drop(crossprod(dotZ_i, v_i))   # \dot{Z}_i' H_i^{-1} \tilde{e}_i
  d_i            <- drop(inv_pre %*% q_i)
  theta_post     <- drop(theta_pre + d_i)
  names(theta_post) <- names(theta_pre)

  # --- Prop 2: (\dot{Z}'\dot{Z})^{-1} 업데이트 (Woodbury)
  G_i         <- crossprod(dotZ_i, W_i)            # \dot{Z}_i' H_i^{-1} \dot{Z}_i
  inv_post    <- inv_pre - inv_pre %*% G_i %*% inv_pre

  # sigma^2 업데이트
  # \hat{\sigma}^{2*} = [(n-N-k)/(n^*-(N+1)-k)] \hat{\sigma}^2
  #                   + \tilde{e}_{N+1}' H^{-1}_{N+1} \tilde{e}_{N+1} / [n^*-(N+1)-k]
  n_post   <- S_pre$n + T_i
  N_post   <- S_pre$N + 1L
  df_pre   <- S_pre$n - S_pre$N - p
  df_post  <- n_post - N_post - p
  sig2_post <- (df_pre / df_post) * sig2_pre +
    as.numeric(crossprod(e_tilde_i, v_i)) / df_post

  # --- Prop 3: cluster-robust 업데이트
  # d_i = (\dot{Z}'\dot{Z})^{-1} \dot{Z}_{N+1}' H^{-1}_{N+1} \tilde{e}_{N+1}
  # U_{N+1} = I_k - (\dot{Z}'\dot{Z})^{-1} \dot{Z}_{N+1}' H^{-1}_{N+1} \dot{Z}_{N+1}
  I_p   <- diag(p)
  U_i   <- I_p - inv_pre %*% G_i           # = I_p - inv_pre %*% dotZ_i' %*% W_i
  K_th  <- kronecker(t(theta_pre), I_p)    # (\hat{\beta}' \otimes I_k)
  K_d   <- kronecker(I_p, matrix(d_i, ncol = 1))  # (I_k \otimes d_i)

  # \hat{\Omega}^{(1)}_{N+1} = (A_N - K_\beta B_N)(I_k \otimes d_i)
  Omega1 <- -(A_pre - K_th %*% B_pre) %*% K_d
  Omega2 <- t(Omega1)
  # \hat{\Omega}^{(2)}_{N+1} = (d_i' \otimes I_k) B_N (I_k \otimes d_i)
  Omega3 <- kronecker(t(d_i), I_p) %*% B_pre %*% K_d
  # \hat{\Omega}^{(3)}_{N+1} = \dot{Z}_{N+1}' H^{-1}_{N+1} \tilde{e}_{N+1}
  #                             \tilde{e}_{N+1}' H^{-1}_{N+1} \dot{Z}_{N+1}
  Omega4 <- tcrossprod(q_i)

  Vcr_post <- U_i %*%
    (Vcr_pre + inv_pre %*% (Omega1 + Omega2 + Omega3 + Omega4) %*% inv_pre) %*%
    t(U_i)

  # A_N, B_N 업데이트 (새 unit의 기여 추가)
  S_i    <- crossprod(dotZ_i)
  s_i    <- crossprod(dotZ_i, matrix(dotY_i, ncol = 1))
  vecS_i <- matrix(vec(S_i), ncol = 1)

  S_post <- S_pre
  S_post$N          <- N_post
  S_post$n          <- n_post
  S_post$inv_dotZtZ <- inv_post
  S_post$theta_hat  <- theta_post
  S_post$sigma2_hat <- sig2_post
  S_post$Vcr_hat    <- Vcr_post
  S_post$M_hat      <- NULL  # Alg 1/2 후엔 M_hat 추적하지 않음 (Alg 3 진입 시 재구성)
  S_post$A_N        <- A_pre + s_i %*% t(vecS_i)
  S_post$B_N        <- B_pre + vecS_i %*% t(vecS_i)

  S_post
}


# -----------------------------------------------------------------------------
# Algorithm 2 (Props 4–6): 기존 unit의 기존 calendar time에 새 관측치 도착
#
# 입력  S_pre    : 업데이트 전 state
#       unit_pre : 해당 unit의 업데이트 전 요약 객체
#       x_new    : 길이 k 공변량 벡터 (새 관측치)
#       time_new : 새 관측치의 calendar time (기존 support 내)
#       y_new    : 새 관측치의 종속변수 값
# 반환  list(S_post, unit_post)
# -----------------------------------------------------------------------------
alg2_existing_unit <- function(S_pre, unit_pre, x_new, time_new, y_new) {

  stopifnot(as.integer(time_new) %in% seq_len(S_pre$T_support))

  # 새 관측치를 기존 T_support 기반으로 p차원 Z로 변환
  x_row   <- matrix(x_new, nrow = 1, dimnames = list(NULL, S_pre$x_cols))
  Z_new   <- build_Z_raw(x_row, as.integer(time_new),
                         S_pre$T_support, S_pre$baseline_time)

  # within 편차: \dot{z}_{1,T_1^*} = Z_{new} - \bar{Z}_1,  \dot{y}_{1,T_1^*} = y - \bar{Y}_1
  kappa    <- unit_pre$T_i / (unit_pre$T_i + 1L)
  doty_new <- as.numeric(y_new - unit_pre$barY_i)
  dotz_new <- drop(Z_new - unit_pre$barZ_i)
  dotz_col <- matrix(dotz_new, ncol = 1)

  inv_pre   <- S_pre$inv_dotZtZ
  theta_pre <- S_pre$theta_hat
  sig2_pre  <- S_pre$sigma2_hat
  Vcr_pre   <- S_pre$Vcr_hat
  A_pre     <- S_pre$A_N
  B_pre     <- S_pre$B_N
  p         <- S_pre$p

  # s_{1,T_1^*} = (kappa^{-1} + h_{1,T_1^*})^{-1},  h = \dot{z}' (\dot{Z}'\dot{Z})^{-1} \dot{z}
  h        <- as.numeric(t(dotz_new) %*% inv_pre %*% dotz_new)
  s_scalar <- as.numeric(kappa / (1 + kappa * h))
  e_tilde  <- as.numeric(doty_new - t(dotz_new) %*% theta_pre)

  # --- Prop 4: 계수 업데이트
  # \hat{\beta}^* = \hat{\beta} + s_{1,T_1^*} (\dot{Z}'\dot{Z})^{-1} \dot{z}_{1,T_1^*} \tilde{e}
  theta_post <- drop(theta_pre + s_scalar * inv_pre %*% dotz_col * e_tilde)
  names(theta_post) <- names(theta_pre)

  # --- Prop 5: (\dot{Z}'\dot{Z})^{-1} 및 sigma^2 업데이트 (Sherman-Morrison)
  # (\dot{X}^{*'}\dot{X}^*)^{-1} = U_{1,T_1^*} (\dot{X}'\dot{X})^{-1}
  # U_{1,T_1^*} = I_k - s (\dot{X}'\dot{X})^{-1} \dot{x} \dot{x}'
  inv_post <- inv_pre -
    s_scalar * (inv_pre %*% dotz_col %*% t(dotz_col) %*% inv_pre)

  n_post  <- S_pre$n + 1L
  N_post  <- S_pre$N
  df_pre  <- S_pre$n - S_pre$N - p
  df_post <- n_post - N_post - p
  # \hat{\sigma}^{2*} = [(n-N-k)/(n^*-N-k)] \hat{\sigma}^2
  #                   + [\tilde{e}^2 / (n^*-N-k)] * [kappa / (1 + kappa h)]
  sig2_post <- (df_pre / df_post) * sig2_pre +
    (e_tilde^2 / df_post) * (kappa / (1 + kappa * h))

  # --- Prop 6: cluster-robust 업데이트
  I_p       <- diag(p)
  U_rank1   <- I_p - s_scalar * (inv_pre %*% dotz_col %*% t(dotz_col))
  # 개인 cluster score: \dot{X}_i' \hat{e}_i = s_i - S_i \hat{\beta}
  Xi_e_pre  <- unit_pre$s_i - unit_pre$S_i %*% matrix(theta_pre, ncol = 1)

  Omega1 <- e_tilde *
    (A_pre - kronecker(t(theta_pre), I_p) %*% B_pre) %*%
    kronecker(I_p, s_scalar * (inv_pre %*% dotz_col))
  Omega2 <- e_tilde^2 *
    kronecker(s_scalar * t(dotz_col) %*% inv_pre, I_p) %*%
    B_pre %*%
    kronecker(I_p, s_scalar * (inv_pre %*% dotz_col))
  # \Omega^{(3)} = s \tilde{e} (\dot{X}_i' \hat{e}_i \dot{x}' + \dot{x} \hat{e}_i' \dot{X}_i)
  Omega3 <- s_scalar * e_tilde *
    (Xi_e_pre %*% t(dotz_col) + dotz_col %*% t(Xi_e_pre))
  # \Omega^{(4)} 항
  Omega4 <- s_scalar^2 * e_tilde^2 * (
    dotz_col %*% t(dotz_col) -
    dotz_col %*% t(dotz_col) %*% inv_pre %*% unit_pre$S_i -
    unit_pre$S_i %*% inv_pre %*% dotz_col %*% t(dotz_col)
  )
  Vcr_post <- U_rank1 %*%
    (Vcr_pre + inv_pre %*% (-Omega1 - t(Omega1) + Omega2 + Omega3 + Omega4) %*% inv_pre) %*%
    t(U_rank1)

  # 해당 unit 요약 업데이트
  G_old    <- unit_pre$S_i
  s_old    <- unit_pre$s_i

  unit_post        <- unit_pre
  unit_post$S_i    <- unit_pre$S_i + kappa * (dotz_col %*% t(dotz_col))
  unit_post$s_i    <- unit_pre$s_i + kappa * (dotz_col * doty_new)
  unit_post$barZ_i <- unit_pre$barZ_i + dotz_new / (unit_pre$T_i + 1L)
  unit_post$barY_i <- unit_pre$barY_i + doty_new  / (unit_pre$T_i + 1L)
  unit_post$T_i    <- unit_pre$T_i + 1L

  vecG_old <- matrix(vec(G_old), ncol = 1)
  vecG_new <- matrix(vec(unit_post$S_i), ncol = 1)

  S_post <- S_pre
  S_post$N          <- N_post
  S_post$n          <- n_post
  S_post$inv_dotZtZ <- inv_post
  S_post$theta_hat  <- theta_post
  S_post$sigma2_hat <- sig2_post
  S_post$Vcr_hat    <- Vcr_post
  S_post$M_hat      <- NULL
  S_post$A_N        <- A_pre - (s_old %*% t(vecG_old)) + (unit_post$s_i %*% t(vecG_new))
  S_post$B_N        <- B_pre - (vecG_old %*% t(vecG_old)) + (vecG_new %*% t(vecG_new))

  list(S_post = S_post, unit_post = unit_post)
}


# -----------------------------------------------------------------------------
# Algorithm 3 (Props 7–9): 기존 unit의 새로운 calendar time T+1 도착
#
# parameter dimension p → p+1 (새 time dummy 추가)
#
# 입력  S_pre    : 업데이트 전 state  (T_support = T)
#       unit_pre : 해당 unit의 업데이트 전 요약 객체
#       x_new    : 길이 k 공변량 벡터
#       time_new : 새 calendar time (반드시 T_support + 1)
#       y_new    : 새 관측치의 종속변수 값
# 반환  list(S_post, unit_post)
# -----------------------------------------------------------------------------
alg3_new_caltime <- function(S_pre, unit_pre, x_new, time_new, y_new) {

  time_new <- as.integer(time_new)
  stopifnot(time_new == S_pre$T_support + 1L)

  p_pre   <- S_pre$p
  p_post  <- p_pre + 1L
  inv_pre <- S_pre$inv_dotZtZ
  theta_pre <- S_pre$theta_hat
  A_pre   <- S_pre$A_N
  B_pre   <- S_pre$B_N

  new_dummy_name <- paste0("factor(time)", time_new)

  # Z_new를 기존 p차원 공간에서 계산
  # (새 time dummy = 1 이지만, 기존 support에 속하지 않으므로 old Z에선 0으로 처리)
  x_row    <- matrix(x_new, nrow = 1, dimnames = list(NULL, S_pre$x_cols))
  Z_new_p  <- build_Z_raw(x_row, as.integer(time_new),
                          S_pre$T_support, S_pre$baseline_time)
  # Z_new_p는 기존 p차원 (새 dummy column 없음)

  # within 편차 (기존 p차원 barZ_i 기준)
  kappa    <- unit_pre$T_i / (unit_pre$T_i + 1L)
  doty_new <- as.numeric(y_new - unit_pre$barY_i)
  g_old    <- drop(Z_new_p - unit_pre$barZ_i)   # 길이 p_pre

  # --- Prop 7: 확장된 역행렬 및 계수 업데이트
  # g = [g_old, 1] ∈ R^{p+1}  (새 time dummy의 within 편차 = 1 - 0 = 1)
  # \dot{Z}^{*'} \dot{Z}^* = diag(\dot{Z}'\dot{Z}, 0) + kappa * g g'  (4.3)
  # 역행렬 (4.4):
  #   [ (Z'Z)^{-1}         -q/g2 ]
  #   [ -q'/g2    (1+kappa*h)/(kappa*g2^2) ]
  # 여기서 g2 = 1 (새 time dummy 값), q = (Z'Z)^{-1} g_old, h = g_old' q
  h   <- as.numeric(t(g_old) %*% inv_pre %*% g_old)
  q   <- drop(inv_pre %*% g_old)   # 길이 p_pre

  nms_post <- c(names(theta_pre), new_dummy_name)

  inv_post <- rbind(
    cbind(inv_pre,           -matrix(q,    ncol = 1)),
    cbind(-matrix(q, nrow = 1), (1 + kappa * h) / kappa)
  )
  dimnames(inv_post) <- list(nms_post, nms_post)

  # \hat{v}_{T+1} = \dot{y}_{1,t_{new}} - g_old' \hat{\theta}  (4.5)
  v_new      <- as.numeric(doty_new - t(g_old) %*% theta_pre)
  theta_post <- c(theta_pre, v_new)
  names(theta_post) <- nms_post

  # --- Prop 8: sigma^2 불변
  # n^* - N - p^* = (n+1) - N - (p+1) = n - N - p  → df 불변, RSS 불변
  # \hat{sigma}^{2*} = \hat{sigma}^2
  sig2_post <- S_pre$sigma2_hat

  # A_N, B_N를 (p+1)차원으로 zero-extend 후 해당 unit 기여 교체
  idx_embed <- vec_embed_index(p_pre)
  A_ext     <- matrix(0, nrow = p_post, ncol = p_post^2)
  B_ext     <- matrix(0, nrow = p_post^2, ncol = p_post^2)
  A_ext[seq_len(p_pre), idx_embed] <- A_pre
  B_ext[idx_embed, idx_embed]      <- B_pre

  # M = \sum_i r_i r_i'를 재구성 (Alg 1/2에서 M_hat을 NULL로 관리하므로 재계산)
  if (is.null(S_pre$M_hat)) {
    dotZtZ_pre <- solve(inv_pre)
    M_pre      <- dotZtZ_pre %*% S_pre$Vcr_hat %*% dotZtZ_pre
  } else {
    M_pre <- S_pre$M_hat
  }
  dimnames(M_pre) <- list(names(theta_pre), names(theta_pre))

  # M을 (p+1)차원으로 zero-extend: \tilde{M}
  M_ext <- rbind(cbind(M_pre, 0), 0)
  dimnames(M_ext) <- list(nms_post, nms_post)

  # --- Prop 9: cluster score 외적 합 M^* 업데이트
  # r_1 = s_1 - S_1 \hat{\theta},  \tilde{r}_1 = [r_1, 0]
  r_pre   <- unit_pre$s_i - unit_pre$S_i %*% matrix(theta_pre, ncol = 1)
  r_tilde <- rbind(r_pre, 0)

  # 업데이트 후 unit 1의 \dot{Z}^{*'}\dot{Z}^*,  \dot{Z}^{*'}\dot{Y}^*
  g_col      <- matrix(g_old, ncol = 1)
  S1_kappa   <- unit_pre$S_i + kappa * (g_col %*% t(g_col))   # 기존 p 차원
  s1_kappa   <- unit_pre$s_i + kappa * (g_col * doty_new)

  # (p+1)차원으로 확장: 새 time dummy와의 외적 추가
  c_star  <- kappa * g_col                # \kappa g_{old} ∈ R^{p_pre}
  S_i_post <- rbind(
    cbind(S1_kappa,          c_star),
    cbind(t(c_star),         kappa)
  )
  s_i_post <- rbind(s1_kappa, kappa * doty_new)
  dimnames(S_i_post) <- list(nms_post, nms_post)
  rownames(s_i_post) <- nms_post

  # r_1^* = \tilde{r}_1 + \kappa g \tilde{e}_{1,t_{new}} - (\tilde{S}_1 + \kappa g g')(θ^* - \tilde{θ})
  # \tilde{θ} = [θ, 0] ∈ R^{p+1}
  theta_tilde <- c(theta_pre, 0)
  S_i_tilde   <- rbind(cbind(unit_pre$S_i, 0), 0)   # \tilde{S}_1

  e_tilde_new <- doty_new - t(g_old) %*% theta_pre   # \dot{y} - g' θ (= v_new 전 값)
  g_full      <- c(g_old, 1)                          # g ∈ R^{p+1}
  r_post      <- r_tilde +
    kappa * g_full * as.numeric(e_tilde_new) -
    (S_i_tilde + kappa * outer(g_full, g_full)) %*%
    matrix(theta_post - theta_tilde, ncol = 1)

  # M^* = \tilde{M} - \tilde{r}_1 \tilde{r}_1' + r_1^* r_1^{*'}
  M_post  <- M_ext - tcrossprod(r_tilde) + tcrossprod(r_post)
  Vcr_post <- inv_post %*% M_post %*% inv_post

  # A_N, B_N에서 unit 1의 기존 기여 제거 후 새 기여 추가
  S_i_old_ext <- rbind(cbind(unit_pre$S_i, 0), 0)
  s_i_old_ext <- rbind(unit_pre$s_i, 0)
  vec_old     <- matrix(vec(S_i_old_ext), ncol = 1)
  vec_new     <- matrix(vec(S_i_post),    ncol = 1)

  A_post <- A_ext - (s_i_old_ext %*% t(vec_old)) + (s_i_post %*% t(vec_new))
  B_post <- B_ext - (vec_old %*% t(vec_old))      + (vec_new %*% t(vec_new))

  # unit 요약 업데이트 (p+1 차원)
  unit_post        <- unit_pre
  barZ_upd         <- unit_pre$barZ_i + g_old / (unit_pre$T_i + 1L)
  unit_post$barZ_i <- c(barZ_upd, 1 / (unit_pre$T_i + 1L))
  names(unit_post$barZ_i) <- nms_post
  unit_post$barY_i <- unit_pre$barY_i + doty_new / (unit_pre$T_i + 1L)
  unit_post$T_i    <- unit_pre$T_i + 1L
  unit_post$S_i    <- S_i_post
  unit_post$s_i    <- s_i_post

  S_post <- S_pre
  S_post$T_support  <- time_new
  S_post$p          <- p_post
  S_post$n          <- S_pre$n + 1L
  S_post$inv_dotZtZ <- inv_post
  S_post$theta_hat  <- theta_post
  S_post$sigma2_hat <- sig2_post
  S_post$Vcr_hat    <- Vcr_post
  S_post$M_hat      <- NULL   # 다음 Alg 3 진입 시 재구성
  S_post$A_N        <- A_post
  S_post$B_N        <- B_post

  list(S_post = S_post, unit_post = unit_post)
}


# =============================================================================
# Section 5: 메인 API
# =============================================================================

#' 온라인 TWFE 추정
#'
#' data.frame을 unit × time 순으로 순차 처리하여 TWFE 계수 및 분산 추정.
#' 전체 데이터를 메모리에 올릴 필요 없이 unit 단위로 스트리밍 처리 가능.
#'
#' 알고리즘 선택 규칙 (자동):
#'   - 새로운 unit 도착     → Algorithm 1  (Props 1–3)
#'   - 기존 unit, 기존 time → Algorithm 2  (Props 4–6)  [warm-up unit만 해당]
#'   - 기존 unit, 새 time   → Algorithm 3  (Props 7–9)  [warm-up unit만 해당]
#'
#' @param data            data.frame (모든 관측치 포함)
#' @param id_col          unit ID 컬럼명  (character)
#' @param time_col        calendar time 컬럼명  (character, 정수형으로 변환됨)
#' @param y_col           종속변수 컬럼명  (character)
#' @param x_cols          설명변수 컬럼명 벡터  (character vector)
#' @param baseline_time   제거할 기준 time dummy.  NULL이면 min(time) 자동 사용
#' @param warmup_ids      warm-up으로 사용할 unit ID 벡터.  NULL이면 자동 선택
#' @param warmup_n        warm-up unit 최소 수 (warmup_ids가 NULL일 때 사용).
#'                        NULL이면 max(3*p, 30) 사용
#' @param track_all_units TRUE이면 스트리밍 unit도 per-unit state 보존.
#'                        FALSE(기본값)이면 warm-up unit 상태만 보존 → 메모리 O(n_warm * p^2)
#' @param verbose         진행 상황 출력 여부
#' @return class \code{"otwfe"} 객체
#' @export
# =============================================================================
# Vcr formula 헬퍼
# M = Σ r_i r_i',  r_i = s_i - S_i θ
# 수식: M = M_ss - A_N(θ⊗I) - [A_N(θ⊗I)]' + (θ'⊗I) B_N (θ⊗I)
# 여기서 M_ss = Σ s_i s_i',  A_N = Σ s_i vec(S_i)',  B_N = Σ vec(S_i) vec(S_i)'
# =============================================================================
.compute_vcr <- function(S, units) {
  p  <- S$p
  th <- S$theta_hat

  # M_ss_warm: warm-up unit들의 s_i s_i' 합 (units 리스트에서 계산)
  M_ss_warm <- matrix(0, p, p)
  for (u in units) M_ss_warm <- M_ss_warm + tcrossprod(u$s_i)
  M_ss_total <- M_ss_warm + S$M_ss

  # Kronecker 행렬: (θ'⊗I_p) = p × p², (θ⊗I_p) = p² × p
  KtI   <- kronecker(t(th), diag(p))   # p × p²
  KI    <- t(KtI)                       # p² × p
  cross <- S$A_N %*% KI                 # p × p

  M <- M_ss_total - cross - t(cross) + KtI %*% S$B_N %*% KI
  S$inv_dotZtZ %*% M %*% S$inv_dotZtZ
}

# =============================================================================
# df를 unit 경계에서 chunk_size unit 단위로 분할
# id 컬럼 기준 정렬된 데이터 가정
# =============================================================================
.split_df_by_unit <- function(df, id_col, chunk_size) {
  ids      <- df[[id_col]]
  rle_res  <- rle(ids)
  u_lens   <- rle_res$lengths          # 각 unit의 관측치 수
  n_units  <- length(u_lens)
  cum_rows <- cumsum(u_lens)
  row_beg  <- c(1L, cum_rows[-n_units] + 1L)

  n_chunks <- ceiling(n_units / chunk_size)
  chunks   <- vector("list", n_chunks)
  for (ci in seq_len(n_chunks)) {
    u_s <- (ci - 1L) * chunk_size + 1L
    u_e <- min(n_units, ci * chunk_size)
    chunks[[ci]] <- df[row_beg[u_s]:cum_rows[u_e], , drop = FALSE]
  }
  chunks
}

# =============================================================================
# otwfe_init(): 빈 handle 생성
# =============================================================================
#' Online TWFE 빈 handle 초기화
#'
#' @param x_cols        공변량 컬럼명 벡터
#' @param time_col      calendar time 컬럼명
#' @param id_col        unit id 컬럼명
#' @param y_col         종속변수 컬럼명
#' @param T_support     calendar time 지지 크기 (NULL이면 첫 청크에서 자동 감지)
#' @param baseline_time 기준 time dummy (기본값 1)
#' @param verbose       진행 상황 출력 여부
#' @return handle 리스트
#' @export
otwfe_init <- function(x_cols,
                       time_col      = "time",
                       id_col        = "id",
                       y_col         = "y",
                       T_support     = NULL,
                       baseline_time = 1L,
                       verbose       = TRUE) {
  list(
    x_cols        = x_cols,
    time_col      = time_col,
    id_col        = id_col,
    y_col         = y_col,
    T_support     = if (!is.null(T_support)) as.integer(T_support) else NULL,
    baseline_time = as.integer(baseline_time),
    verbose       = verbose,
    S             = NULL,
    units         = list(),
    warmup_done   = FALSE,
    time_map      = NULL,   # 원래 time 값 → 재인덱싱 매핑
    time_levels   = NULL    # 원래 time 값 보존 (출력용)
  )
}

# =============================================================================
# otwfe_update(): chunk_df를 handle에 반영
# =============================================================================
#' Online TWFE 청크 업데이트
#'
#' @param handle    otwfe_init() 또는 이전 otwfe_update()가 반환한 handle
#' @param chunk_df  처리할 데이터프레임 청크 (unit 경계에서 분할된 것)
#' @return 업데이트된 handle
#' @export
otwfe_update <- function(handle, chunk_df) {
  if (!handle$warmup_done) {
    handle <- .otwfe_first_chunk(handle, chunk_df)
  } else {
    handle <- .otwfe_stream_chunk(handle, chunk_df)
  }
  handle
}

# 첫 번째 청크: warm-up 선택 + state 초기화 + Rcpp 배치
.otwfe_first_chunk <- function(handle, chunk_df) {
  x_cols        <- handle$x_cols
  time_col      <- handle$time_col
  id_col        <- handle$id_col
  y_col         <- handle$y_col
  verbose       <- handle$verbose

  # 입력 검증
  for (col in c(id_col, time_col, y_col, x_cols)) {
    if (!col %in% names(chunk_df))
      stop(sprintf("컬럼을 찾을 수 없습니다: '%s'", col))
  }

  chunk_df[[time_col]] <- as.integer(chunk_df[[time_col]])
  all_times_chunk      <- sort(unique(chunk_df[[time_col]]))

  # T_support 결정 (미지정이면 첫 청크의 max time)
  if (is.null(handle$T_support)) {
    T_support <- max(all_times_chunk)
    if (verbose) message(sprintf("T_support 자동 감지 (첫 청크): %d", T_support))
  } else {
    T_support <- handle$T_support
  }

  # time 재인덱싱 (time_map 미설정이면 첫 청크에서 결정)
  if (is.null(handle$time_map)) {
    time_levels          <- sort(unique(c(seq_len(T_support), all_times_chunk)))
    time_map             <- setNames(match(time_levels, time_levels),
                                     as.character(time_levels))
    handle$time_map      <- time_map
    handle$time_levels   <- time_levels
  }
  chunk_df[[time_col]] <- handle$time_map[as.character(chunk_df[[time_col]])]

  baseline_time <- handle$baseline_time
  if (is.na(baseline_time) || baseline_time < 1L || baseline_time > T_support)
    baseline_time <- 1L

  # eligible units (관측치 >= 2)
  unit_counts <- tapply(chunk_df[[id_col]], chunk_df[[id_col]], length)
  eligible    <- names(unit_counts)[unit_counts >= 2L]
  if (length(eligible) == 0L)
    stop("첫 청크에 eligible unit(관측치 >= 2)이 없습니다.")

  # 그리디 warm-up 선택 (모든 calendar time 커버 보장)
  all_times_idx <- seq_len(T_support)
  p_guess  <- length(x_cols) + max(0L, T_support - 1L)
  min_warm <- min(max(3L * p_guess, 30L), length(eligible))

  id_ch_vec  <- as.character(chunk_df[[id_col]])
  unit_times <- tapply(as.integer(chunk_df[[time_col]]), id_ch_vec,
                       function(x) unique(x), simplify = FALSE)
  elig_times <- unit_times[eligible]
  cover_cnt  <- sapply(elig_times, function(ts) sum(all_times_idx %in% ts))
  elig_ord   <- names(sort(cover_cnt, decreasing = TRUE))

  uncovered  <- all_times_idx
  warmup_ids <- character(0)
  for (uid_ch in elig_ord) {
    if (length(uncovered) == 0L) break
    ts <- elig_times[[uid_ch]]
    if (any(uncovered %in% ts)) {
      warmup_ids <- c(warmup_ids, uid_ch)
      uncovered  <- setdiff(uncovered, ts)
    }
  }
  if (length(warmup_ids) < min_warm) {
    extra      <- setdiff(elig_ord, warmup_ids)
    warmup_ids <- c(warmup_ids, head(extra, min_warm - length(warmup_ids)))
  }

  if (verbose)
    cat(sprintf("초기 State 구성 중 (warm-up: %d units)...\n", length(warmup_ids)))

  df_warm         <- chunk_df[id_ch_vec %in% warmup_ids, ]
  df_stream_chunk <- chunk_df[!(id_ch_vec %in% warmup_ids), ]

  # state 초기화 (warm-up 배치 OLS)
  init      <- state_init(df_warm, id_col, time_col, y_col, x_cols, baseline_time)
  S_pre     <- init$S_pre
  units     <- init$units

  handle$S             <- S_pre
  handle$units         <- units
  handle$T_support     <- T_support
  handle$baseline_time <- baseline_time
  handle$warmup_done   <- TRUE

  # 첫 청크의 나머지 unit들을 Rcpp 배치로 처리
  if (nrow(df_stream_chunk) > 0L)
    handle <- .otwfe_stream_chunk(handle, df_stream_chunk)

  handle
}

# 이후 청크: 모두 새 unit → Rcpp 배치
.otwfe_stream_chunk <- function(handle, chunk_df) {
  S_pre    <- handle$S
  x_cols   <- handle$x_cols
  time_col <- handle$time_col
  id_col   <- handle$id_col
  y_col    <- handle$y_col

  # time 재인덱싱 (handle에 저장된 map 사용)
  if (!is.null(handle$time_map))
    chunk_df[[time_col]] <- handle$time_map[as.character(chunk_df[[time_col]])]

  if (!.alg1_batch_rcpp_available)
    stop("Rcpp 배치 모듈이 필요합니다. sourceCpp('src/alg1_batch.cpp')를 확인하세요.")

  stream_idx_list <- split(seq_len(nrow(chunk_df)),
                            as.character(chunk_df[[id_col]]))
  new_units <- names(stream_idx_list)
  if (length(new_units) == 0L) return(handle)

  new_all_idx  <- unlist(stream_idx_list, use.names = FALSE)
  t_all_new    <- as.integer(chunk_df[[time_col]])[new_all_idx]
  over_support <- t_all_new > S_pre$T_support

  if (any(over_support)) {
    warning(sprintf("%d개 관측치가 T_support(%d)를 초과합니다. 해당 관측치 제외.",
                    sum(over_support), S_pre$T_support))
  }

  keep_idx <- new_all_idx[!over_support]
  t_batch  <- t_all_new[!over_support]
  x_batch  <- as.matrix(chunk_df[keep_idx, x_cols, drop = FALSE])
  y_batch  <- as.numeric(chunk_df[[y_col]])[keep_idx]

  if (any(over_support)) {
    keep_flags    <- !over_support
    flags_by_unit <- split(keep_flags,
                            as.character(chunk_df[[id_col]])[new_all_idx])
    lens_batch <- sapply(flags_by_unit[new_units], sum)
  } else {
    lens_batch <- lengths(stream_idx_list)
  }
  lens_batch <- lens_batch[lens_batch > 0L]
  if (length(lens_batch) == 0L) return(handle)

  rcpp_res <- alg1_batch_cpp(
    inv_dotZtZ    = S_pre$inv_dotZtZ,
    theta_hat     = S_pre$theta_hat,
    sigma2_hat    = S_pre$sigma2_hat,
    N_old         = as.integer(S_pre$N),
    n_old         = as.integer(S_pre$n),
    x_mat         = x_batch,
    time_vec      = t_batch,
    y_vec         = y_batch,
    unit_lens     = as.integer(lens_batch),
    T_support     = as.integer(S_pre$T_support),
    baseline_time = as.integer(S_pre$baseline_time)
  )

  theta_new <- rcpp_res$theta_hat
  names(theta_new) <- names(S_pre$theta_hat)
  inv_new <- rcpp_res$inv_dotZtZ
  dimnames(inv_new) <- dimnames(S_pre$inv_dotZtZ)

  S_pre$inv_dotZtZ <- inv_new
  S_pre$theta_hat  <- theta_new
  S_pre$sigma2_hat <- rcpp_res$sigma2_hat
  S_pre$M_hat      <- NULL
  S_pre$N          <- S_pre$N + rcpp_res$N_add
  S_pre$n          <- S_pre$n + rcpp_res$n_add
  S_pre$M_ss       <- S_pre$M_ss + rcpp_res$M_ss_add
  S_pre$A_N        <- S_pre$A_N + rcpp_res$A_N_add
  S_pre$B_N        <- S_pre$B_N + rcpp_res$B_N_add

  handle$S <- S_pre
  handle
}

# =============================================================================
# otwfe_finalize(): 최종 결과 추출
# =============================================================================
#' Online TWFE 최종 결과 추출
#'
#' @param handle    otwfe_update() 완료 후 handle
#' @return class \code{"otwfe"} 객체
#' @export
otwfe_finalize <- function(handle) {
  if (is.null(handle$S))
    stop("otwfe_update()를 먼저 호출하세요.")

  S       <- handle$S
  nm_list <- dimnames(S$inv_dotZtZ)

  # Vcr 계산 (formula: M_ss, A_N, B_N, theta_final 사용)
  Vcr <- .compute_vcr(S, handle$units)
  dimnames(Vcr) <- list(nm_list[[1]], nm_list[[2]])
  S$Vcr_hat <- Vcr

  structure(
    list(
      state       = S,
      units       = handle$units,
      id_col      = handle$id_col,
      time_col    = handle$time_col,
      y_col       = handle$y_col,
      x_cols      = handle$x_cols,
      time_levels = handle$time_levels
    ),
    class = "otwfe"
  )
}

# =============================================================================
otwfe <- function(data,
                  id_col,
                  time_col,
                  y_col,
                  x_cols,
                  baseline_time   = NULL,
                  warmup_ids      = NULL,
                  warmup_n        = NULL,
                  track_all_units = FALSE,
                  verbose         = TRUE,
                  chunk_size      = NULL) {

  t_start <- proc.time()   # 전체 실행 시간 측정 시작

  # --------------------------------------------------
  # chunk_size 경로: otwfe_init/update/finalize 위임
  # --------------------------------------------------
  if (!is.null(chunk_size)) {
    chunk_size <- as.integer(chunk_size)
    if (chunk_size < 1L) stop("chunk_size는 양의 정수여야 합니다.")

    # 전체 데이터에서 time 재인덱싱 (청크 간 일관성 유지)
    stopifnot(is.data.frame(data))
    data[[time_col]] <- as.integer(data[[time_col]])
    all_times_cs     <- sort(unique(data[[time_col]]))
    T_support_cs     <- length(all_times_cs)
    time_levels_cs   <- all_times_cs
    time_map_cs      <- setNames(seq_along(time_levels_cs),
                                 as.character(time_levels_cs))
    data[[time_col]] <- time_map_cs[as.character(data[[time_col]])]

    bl_cs <- if (is.null(baseline_time)) {
      1L
    } else {
      as.integer(time_map_cs[as.character(as.integer(baseline_time))])
    }

    handle <- otwfe_init(
      x_cols        = x_cols,
      time_col      = time_col,
      id_col        = id_col,
      y_col         = y_col,
      T_support     = T_support_cs,
      baseline_time = bl_cs,
      verbose       = verbose
    )
    # time이 이미 재인덱싱됐으므로 identity map 설정
    handle$time_map    <- setNames(seq_len(T_support_cs),
                                   as.character(seq_len(T_support_cs)))
    handle$time_levels <- time_levels_cs

    # 청크 경계를 row 범위로 계산 (청크를 미리 복사하지 않음 → 메모리 절약)
    rle_cs   <- rle(data[[id_col]])
    u_lens   <- rle_cs$lengths
    n_units  <- length(u_lens)
    cum_rows <- cumsum(u_lens)
    n_chunks <- ceiling(n_units / chunk_size)

    if (verbose)
      cat(sprintf("청크 방식 처리: %d 청크 (chunk_size=%s units)\n",
                  n_chunks, format(chunk_size, big.mark = ",")))

    pb <- progress::progress_bar$new(
      format = "  |:bar| :percent",
      total  = n_chunks, clear = FALSE, width = 79,
      force  = verbose
    )
    prev_end <- 0L
    for (ci in seq_len(n_chunks)) {
      u_e      <- min(n_units, ci * chunk_size)
      row_e    <- cum_rows[u_e]
      chunk_df <- data[(prev_end + 1L):row_e, , drop = FALSE]
      handle   <- otwfe_update(handle, chunk_df)
      rm(chunk_df)
      invisible(gc())
      prev_end <- row_e
      if (verbose) pb$tick()
    }

    result <- otwfe_finalize(handle)
    t_elapsed <- (proc.time() - t_start)[["elapsed"]]
    S <- result$state
    if (verbose)
      cat(sprintf("\n완료: N = %s, n = %s, p = %d, T_support = %d | 소요 시간: %.2f 초\n",
                  format(S$N, big.mark = ","), format(S$n, big.mark = ","),
                  S$p, S$T_support, t_elapsed))
    return(result)
  }

  # --------------------------------------------------
  # 1. 입력 검증 및 전처리 (기존 단일 호출 경로)
  # --------------------------------------------------
  stopifnot(is.data.frame(data))
  for (col in c(id_col, time_col, y_col, x_cols)) {
    if (!col %in% names(data))
      stop(sprintf("컬럼을 찾을 수 없습니다: '%s'", col))
  }

  data[[time_col]] <- as.integer(data[[time_col]])

  all_times <- sort(unique(data[[time_col]]))

  if (length(all_times) < 2L)
    stop("최소 2개 이상의 calendar time이 필요합니다.")

  # calendar time을 1, 2, ..., T로 재인덱싱
  # seq_len(T_support) = 1..T가 실제 관측 time 수와 일치해야
  # time dummy 개수가 올바르게 결정됨 (예: 1935~1954 → 1~20)
  time_levels <- all_times                         # 원래 time 값 보존 (출력용)
  time_map    <- setNames(seq_along(time_levels),
                          as.character(time_levels))
  data[[time_col]] <- time_map[as.character(data[[time_col]])]
  all_times <- seq_along(time_levels)              # 1, 2, ..., T

  # baseline_time도 재인덱싱
  if (is.null(baseline_time)) {
    baseline_time <- 1L
  } else {
    baseline_time <- as.integer(time_map[as.character(as.integer(baseline_time))])
    if (is.na(baseline_time))
      stop("baseline_time이 데이터에 존재하지 않는 time 값입니다.")
  }

  # within 변환 최소 조건: 각 unit 관측치 >= 2
  unit_counts <- tapply(data[[id_col]], data[[id_col]], length)
  eligible    <- names(unit_counts)[unit_counts >= 2L]

  # --------------------------------------------------
  # 2. Warm-up unit 선택
  #    핵심 조건: warm-up 데이터가 모든 calendar time을 커버해야
  #    state_init() 후 T_support = T_max 가 보장됨.
  #    (Algorithm 3에 해당하는 time dummy 확장이 배치 OLS 안에서 처리됨)
  # --------------------------------------------------
  if (!is.null(warmup_ids)) {
    # 사용자가 직접 지정한 경우
    warmup_ids <- as.character(warmup_ids)
    warmup_ids <- intersect(warmup_ids, eligible)
    if (length(warmup_ids) < 1L)
      stop("지정한 warmup_ids 중 eligible한 unit이 없습니다.")

  } else {
    # 자동 선택: 모든 calendar time을 커버하는 unit들 우선 선택 후 최솟값 보장
    p_guess  <- length(x_cols) + max(0L, length(all_times) - 1L)
    min_warm <- min(max(3L * p_guess, 30L), length(eligible))
    if (!is.null(warmup_n)) min_warm <- min(warmup_n, length(eligible))

    # unit별 관측 time 집합 (eligible만)
    id_ch_vec  <- as.character(data[[id_col]])
    unit_times <- tapply(as.integer(data[[time_col]]), id_ch_vec,
                         function(x) unique(x), simplify = FALSE)
    elig_times <- unit_times[eligible]

    # 그리디: 커버 기여도(내림차순) 순으로 unit 선택, 모든 time이 커버될 때까지
    cover_cnt  <- sapply(elig_times, function(ts) sum(all_times %in% ts))
    elig_ord   <- names(sort(cover_cnt, decreasing = TRUE))

    uncovered  <- all_times
    warmup_ids <- character(0)
    for (uid_ch in elig_ord) {
      if (length(uncovered) == 0L) break
      ts <- elig_times[[uid_ch]]
      if (any(uncovered %in% ts)) {
        warmup_ids <- c(warmup_ids, uid_ch)
        uncovered  <- setdiff(uncovered, ts)
      }
    }

    # 최솟값(min_warm)까지 추가 확보
    if (length(warmup_ids) < min_warm) {
      extra      <- setdiff(elig_ord, warmup_ids)
      warmup_ids <- c(warmup_ids,
                      head(extra, min_warm - length(warmup_ids)))
    }
  }

  if (length(warmup_ids) < 1L)
    stop("warm-up 가능한 unit이 없습니다. 각 unit은 최소 2개의 관측치가 필요합니다.")

  df_warm   <- data[id_ch_vec %in% warmup_ids, ]
  df_stream <- data[!(id_ch_vec %in% warmup_ids), ]

  # --------------------------------------------------
  # 3. 초기 State 구성 (배치 OLS — Alg 1/2/3와 algebraically equivalent)
  # --------------------------------------------------
  if (verbose) {
    cat(sprintf("초기 State 구성 중 (warm-up: %d units)...\n", length(warmup_ids)))
  }

  init  <- state_init(df_warm, id_col, time_col, y_col, x_cols, baseline_time)
  S_pre <- init$S_pre
  units <- init$units   # warm-up unit 상태만 보존

  # --------------------------------------------------
  # 4. 스트리밍 처리: 나머지 unit들 Algorithm 1으로 처리
  #    track_all_units = FALSE(기본값)이면 per-unit state를 저장하지 않음
  #    → units 리스트는 warm-up unit 상태만 유지 (메모리 절약)
  # --------------------------------------------------
  if (nrow(df_stream) == 0L) {
    t_elapsed <- (proc.time() - t_start)[["elapsed"]]
    if (verbose)
      cat(sprintf("추가 스트리밍 데이터 없음. 초기 상태 반환. | 소요 시간: %.2f 초\n",
                  t_elapsed))
    return(structure(
      list(state = S_pre, units = units,
           id_col = id_col, time_col = time_col,
           y_col = y_col, x_cols = x_cols,
           warmup_ids = warmup_ids,
           time_levels = time_levels),
      class = "otwfe"
    ))
  }

  # unit × time 기준 정렬 후 unit별 인덱스 사전 분할 (O(n))
  ord       <- order(df_stream[[id_col]], df_stream[[time_col]])
  df_stream <- df_stream[ord, ]

  stream_idx_list <- split(seq_len(nrow(df_stream)),
                           as.character(df_stream[[id_col]]))
  stream_unit_chs <- names(stream_idx_list)
  n_stream        <- length(stream_unit_chs)

  # 스트리밍 unit을 두 그룹으로 분리:
  #   warm_unit_chs : warm-up unit (Algorithm 2 or 3)
  #   new_unit_chs  : 새 unit     (Algorithm 1 — Rcpp 배치 또는 R fallback)
  warm_in_stream <- stream_unit_chs[stream_unit_chs %in% names(units)]
  new_in_stream  <- stream_unit_chs[!stream_unit_chs %in% names(units)]

  if (verbose) {
    cat(sprintf("스트리밍 처리: %s개 unit (%d warm-up / %s new), %s개 관측치\n",
                format(n_stream, big.mark = ","),
                length(warm_in_stream),
                format(length(new_in_stream), big.mark = ","),
                format(nrow(df_stream), big.mark = ",")))
    pb <- make_progress_bar(0, n_stream)
  }

  # ------------------------------------------------------------------
  # STEP A: warm-up unit 스트리밍 처리 (Algorithm 2 / 3, R 루프)
  # ------------------------------------------------------------------
  for (u_idx in seq_along(warm_in_stream)) {
    uid_ch  <- warm_in_stream[u_idx]
    df_u    <- df_stream[stream_idx_list[[uid_ch]], ]
    times_u <- df_u[[time_col]]
    y_u     <- df_u[[y_col]]
    x_mat_u <- as.matrix(df_u[, x_cols, drop = FALSE])

    for (obs_idx in seq_len(nrow(df_u))) {
      t_obs <- times_u[obs_idx]
      y_obs <- y_u[obs_idx]
      x_obs <- x_mat_u[obs_idx, ]

      if (t_obs <= S_pre$T_support) {
        # 케이스 B: 기존 calendar time → Algorithm 2
        res             <- alg2_existing_unit(S_pre, units[[uid_ch]], x_obs, t_obs, y_obs)
        S_pre           <- res$S_post
        units[[uid_ch]] <- res$unit_post

      } else {
        # 케이스 C: 새 calendar time T+1 → Algorithm 3
        stopifnot(t_obs == S_pre$T_support + 1L)

        res             <- alg3_new_caltime(S_pre, units[[uid_ch]], x_obs, t_obs, y_obs)
        S_pre           <- res$S_post
        units[[uid_ch]] <- res$unit_post

        # 나머지 모든 기존 unit 요약을 zero-extend
        new_dummy <- paste0("factor(time)", t_obs)
        for (other_id in names(units)) {
          if (other_id != uid_ch)
            units[[other_id]] <- expand_unit_state_zero(units[[other_id]], new_dummy)
        }
      }
    }
    if (verbose) safe_set_progress(pb, u_idx)
  }

  # ------------------------------------------------------------------
  # STEP B: 새 unit 처리 (Algorithm 1)
  #   track_all_units = FALSE & Rcpp 사용 가능 → Rcpp 배치 (고속)
  #   그 외                                    → pure-R 루프 (fallback)
  # ------------------------------------------------------------------
  use_rcpp_batch <- (!track_all_units) && .alg1_batch_rcpp_available &&
                    (length(new_in_stream) > 0L)

  if (use_rcpp_batch) {
    # ---- Rcpp 배치 경로 ----
    # 새 unit 전체 인덱스를 한 번에 취합 — per-unit 루프 없이 벡터화 처리
    new_all_idx   <- unlist(stream_idx_list[new_in_stream], use.names = FALSE)
    t_all_new     <- as.integer(df_stream[[time_col]])[new_all_idx]
    over_support  <- t_all_new > S_pre$T_support

    if (any(over_support)) {
      warning(sprintf(
        "%d개 관측치가 T_support(%d)를 초과합니다. 해당 관측치 제외.",
        sum(over_support), S_pre$T_support))
    }

    # T_support 이하 관측치만 유지 — 필요한 행만 추출해 메모리 절약
    keep_idx   <- new_all_idx[!over_support]
    t_batch    <- t_all_new[!over_support]
    x_batch    <- as.matrix(df_stream[keep_idx, x_cols, drop = FALSE])
    y_batch    <- as.numeric(df_stream[[y_col]])[keep_idx]

    # unit_lens: T_support 필터링 후 각 unit의 실제 관측치 수
    # (필터링 없는 일반 경우엔 lengths()로 직접 계산해 빠름)
    if (any(over_support)) {
      # 필터링 있는 경우: unit별 keep 여부를 계산
      lens_raw  <- lengths(stream_idx_list[new_in_stream])
      # 각 unit의 관측치가 연속으로 쌓여있으므로 cumsum으로 단위 경계 파악
      cum_lens  <- c(0L, cumsum(lens_raw))
      lens_vec  <- integer(length(new_in_stream))
      for (j in seq_along(new_in_stream)) {
        obs_j     <- (cum_lens[j] + 1L) : cum_lens[j + 1L]
        lens_vec[j] <- sum(!over_support[obs_j])
      }
      # 관측치 2개 미만 unit 제거
      valid_j   <- lens_vec >= 2L
      if (!all(valid_j)) {
        # 해당 unit의 기여 구간을 keep_idx에서도 제거
        cum_keep  <- c(0L, cumsum(lens_vec))
        keep_rows <- unlist(lapply(which(valid_j),
                                   function(j) seq(cum_keep[j]+1L, cum_keep[j+1L])))
        x_batch   <- x_batch[keep_rows, , drop = FALSE]
        t_batch   <- t_batch[keep_rows]
        y_batch   <- y_batch[keep_rows]
        lens_vec  <- lens_vec[valid_j]
      }
    } else {
      lens_vec <- lengths(stream_idx_list[new_in_stream])
      valid_j  <- lens_vec >= 2L
      if (!all(valid_j)) {
        cum_lens  <- c(0L, cumsum(lens_vec))
        keep_rows <- unlist(lapply(which(valid_j),
                                   function(j) seq(cum_lens[j]+1L, cum_lens[j+1L])))
        x_batch   <- x_batch[keep_rows, , drop = FALSE]
        t_batch   <- t_batch[keep_rows]
        y_batch   <- y_batch[keep_rows]
        lens_vec  <- lens_vec[valid_j]
      }
    }
    lens_batch <- lens_vec

    vj <- seq_along(lens_batch)   # 유효 unit 인덱스
    if (length(vj) > 0L) {

      rcpp_res <- alg1_batch_cpp(
        inv_dotZtZ    = S_pre$inv_dotZtZ,
        theta_hat     = S_pre$theta_hat,
        sigma2_hat    = S_pre$sigma2_hat,
        N_old         = as.integer(S_pre$N),
        n_old         = as.integer(S_pre$n),
        x_mat         = x_batch,
        time_vec      = t_batch,
        y_vec         = y_batch,
        unit_lens     = as.integer(lens_batch),
        T_support     = as.integer(S_pre$T_support),
        baseline_time = as.integer(S_pre$baseline_time)
      )

      # 이름 복원 (Rcpp는 dimnames 없이 반환)
      theta_new <- rcpp_res$theta_hat
      names(theta_new) <- names(S_pre$theta_hat)
      inv_new   <- rcpp_res$inv_dotZtZ
      dimnames(inv_new) <- dimnames(S_pre$inv_dotZtZ)

      # 상태 업데이트: M_ss/A_N/B_N 누적 (Vcr는 otwfe_finalize에서 formula로 계산)
      S_pre$inv_dotZtZ <- inv_new
      S_pre$theta_hat  <- theta_new
      S_pre$sigma2_hat <- rcpp_res$sigma2_hat
      S_pre$M_hat      <- NULL
      S_pre$N          <- S_pre$N + rcpp_res$N_add
      S_pre$n          <- S_pre$n + rcpp_res$n_add
      S_pre$M_ss       <- S_pre$M_ss + rcpp_res$M_ss_add
      S_pre$A_N        <- S_pre$A_N + rcpp_res$A_N_add
      S_pre$B_N        <- S_pre$B_N + rcpp_res$B_N_add

      # Vcr 계산 (단일 호출 방식에서는 theta_new = theta_final이므로 즉시 계산)
      nm_list  <- dimnames(S_pre$inv_dotZtZ)
      Vcr_new  <- .compute_vcr(S_pre, units)
      dimnames(Vcr_new) <- list(nm_list[[1]], nm_list[[2]])
      S_pre$Vcr_hat <- Vcr_new
    }

    if (verbose) {
      safe_set_progress(pb, length(warm_in_stream) + length(new_in_stream))
    }

  } else {
    # ---- pure-R 루프 fallback (track_all_units = TRUE 또는 Rcpp 미사용) ----
    for (u_idx in seq_along(new_in_stream)) {
      uid_ch  <- new_in_stream[u_idx]
      df_u    <- df_stream[stream_idx_list[[uid_ch]], ]
      times_u <- df_u[[time_col]]
      y_u     <- df_u[[y_col]]
      x_mat_u <- as.matrix(df_u[, x_cols, drop = FALSE])

      if (any(times_u > S_pre$T_support)) {
        warning(sprintf(
          "Unit '%s'가 T_support(%d)를 초과하는 calendar time을 포함합니다. 해당 관측치 제외.",
          uid_ch, S_pre$T_support))
        keep    <- times_u <= S_pre$T_support
        times_u <- times_u[keep]
        y_u     <- y_u[keep]
        x_mat_u <- x_mat_u[keep, , drop = FALSE]
        if (length(times_u) == 0L) {
          if (verbose) safe_set_progress(pb, length(warm_in_stream) + u_idx)
          next
        }
      }

      S_pre <- alg1_new_unit(S_pre, x_mat_u, times_u, y_u)

      if (track_all_units) {
        Z_raw_u <- build_Z_raw(x_mat_u, as.integer(times_u),
                               S_pre$T_support, S_pre$baseline_time)
        wt_u    <- within_transform_unit(Z_raw_u, y_u)
        S_iu    <- crossprod(wt_u$dotZ_i)
        s_iu    <- crossprod(wt_u$dotZ_i, matrix(wt_u$dotY_i, ncol = 1))
        units[[uid_ch]] <- make_unit_state(
          id     = uid_ch,
          T_i    = length(times_u),
          barZ_i = wt_u$barZ_i,
          barY_i = wt_u$barY_i,
          S_i    = S_iu,
          s_i    = s_iu
        )
      }

      if (verbose) safe_set_progress(pb, length(warm_in_stream) + u_idx)
    }
  }

  t_elapsed <- (proc.time() - t_start)[["elapsed"]]

  if (verbose) {
    close(pb)
    cat(sprintf("\n완료: N = %s, n = %s, p = %d, T_support = %d | 소요 시간: %.2f 초\n",
                format(S_pre$N, big.mark = ","),
                format(S_pre$n, big.mark = ","),
                S_pre$p, S_pre$T_support,
                t_elapsed))
  }

  structure(
    list(
      state       = S_pre,
      units       = units,
      id_col      = id_col,
      time_col    = time_col,
      y_col       = y_col,
      x_cols      = x_cols,
      warmup_ids  = warmup_ids,
      time_levels = time_levels   # 원래 calendar time 값 (재인덱싱 역변환용)
    ),
    class = "otwfe"
  )
}


# =============================================================================
# Section 6: Output — tidy / summary / print
# =============================================================================

#' otwfe 결과에서 계수 테이블 추출
#'
#' @param x         \code{otwfe()} 반환 객체
#' @param vcov      분산 종류: \code{"cluster"} (default) 또는 \code{"classical"}
#' @param conf_level 신뢰 수준 (default 0.95)
#' @param include_time_fe time fixed effect 계수를 테이블에 포함할지 여부
#' @return \code{data.frame}: term, estimate, std.error, statistic, p.value,
#'         conf.low, conf.high
#' @export
tidy.otwfe <- function(x,
                             vcov            = c("cluster", "classical"),
                             conf_level      = 0.95,
                             include_time_fe = FALSE) {

  vcov   <- match.arg(vcov)
  S      <- x$state
  theta  <- S$theta_hat
  p      <- length(theta)

  # 분산-공분산 행렬 선택
  V <- if (vcov == "classical") {
    S$sigma2_hat * S$inv_dotZtZ
  } else {
    S$Vcr_hat
  }
  colnames(V) <- rownames(V) <- names(theta)

  se   <- sqrt(pmax(diag(V), 0))
  tval <- theta / se
  pval <- 2 * pt(abs(tval), df = S$n - S$N - p, lower.tail = FALSE)

  alpha <- 1 - conf_level
  z     <- qt(1 - alpha / 2, df = S$n - S$N - p)

  tab <- data.frame(
    term      = names(theta),
    estimate  = theta,
    std.error = se,
    statistic = tval,
    p.value   = pval,
    conf.low  = theta - z * se,
    conf.high = theta + z * se,
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  # time FE 제외 옵션
  if (!include_time_fe) {
    is_time_fe <- grepl("^factor\\(time\\)", tab$term)
    tab        <- tab[!is_time_fe, ]
  }

  tab
}

#' @export
print.otwfe <- function(x, ...) {
  S <- x$state
  cat("Online TWFE Estimator\n")
  cat(sprintf("  Units (N): %s\n",       format(S$N,         big.mark = ",")))
  cat(sprintf("  Observations (n): %s\n", format(S$n,         big.mark = ",")))
  cat(sprintf("  Covariates (k): %d\n",   length(x$x_cols)))
  cat(sprintf("  T_support: %d\n",        S$T_support))
  cat(sprintf("  Parameters (p): %d\n",   S$p))
  cat("\nCoefficients (cluster-robust SE):\n")
  tab <- tidy.otwfe(x, vcov = "cluster", include_time_fe = FALSE)
  coef_tab <- data.frame(
    Estimate  = round(tab$estimate,  6),
    `Std.Err` = round(tab$std.error, 6),
    `t value` = round(tab$statistic, 4),
    `Pr(>|t|)` = format.pval(tab$p.value, digits = 3),
    check.names = FALSE
  )
  rownames(coef_tab) <- tab$term
  print(coef_tab)
  invisible(x)
}

#' @export
summary.otwfe <- function(object, vcov = c("cluster", "classical"), ...) {
  vcov <- match.arg(vcov)
  S    <- object$state

  cat("==============================================\n")
  cat("Online TWFE (Hwang & Lee, 2026)\n")
  cat("==============================================\n")
  cat(sprintf("Dependent variable : %s\n", object$y_col))
  cat(sprintf("Covariates         : %s\n", paste(object$x_cols, collapse = ", ")))
  cat(sprintf("ID column          : %s\n", object$id_col))
  cat(sprintf("Time column        : %s\n", object$time_col))
  cat(sprintf("Baseline time      : %s\n", S$baseline_time))
  cat(sprintf("T_support          : %d\n", S$T_support))
  cat(sprintf("N (units)          : %s\n", format(S$N, big.mark = ",")))
  cat(sprintf("n (observations)   : %s\n", format(S$n, big.mark = ",")))
  cat(sprintf("Parameters (p)     : %d\n", S$p))
  cat(sprintf("Residual df        : %d\n", S$n - S$N - S$p))
  cat(sprintf("sigma^2            : %.6f\n", S$sigma2_hat))
  cat(sprintf("Variance type      : %s\n",
              if (vcov == "cluster") "Cluster-robust (HC0, Arellano)"
              else                   "Classical (homoskedastic)"))
  cat("----------------------------------------------\n")

  tab <- tidy.otwfe(object, vcov = vcov, include_time_fe = FALSE)
  coef_tab <- data.frame(
    Estimate  = round(tab$estimate,  6),
    `Std.Err` = round(tab$std.error, 6),
    `t value` = round(tab$statistic, 4),
    `Pr(>|t|)` = format.pval(tab$p.value, digits = 3),
    `2.5%`    = round(tab$conf.low,  6),
    `97.5%`   = round(tab$conf.high, 6),
    check.names = FALSE
  )
  rownames(coef_tab) <- tab$term
  print(coef_tab)
  cat("==============================================\n")
  invisible(object)
}
