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


