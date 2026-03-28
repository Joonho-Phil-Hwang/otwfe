# =============================================================================
# online_twfe_sim.R
#
# 논문 Section 5.2 Simulation 재현 코드
# online_twfe_core.R의 일반 함수를 사용하며,
# simulation 전용 로직(데이터 생성, 벤치마크)은 이 파일에 유지
# =============================================================================

# online_twfe_core.R를 같은 폴더에서 자동으로 찾아 source
# Rscript 실행, RStudio, source() 환경 모두 대응
.sim_core_path <- local({
  # 1순위: 이 파일 자체의 경로 (source()로 불릴 때)
  self <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NULL)
  # 2순위: Rscript --file= 인자
  args <- commandArgs(trailingOnly = FALSE)
  rscript <- tryCatch({
    f <- sub("--file=", "", args[grep("--file=", args)])
    if (length(f) == 1L && nchar(f) > 0L) normalizePath(f) else NULL
  }, error = function(e) NULL)
  # 3순위: 현재 작업 디렉토리
  candidate <- Filter(Negate(is.null), list(self, rscript))
  dir_path  <- if (length(candidate) > 0L) dirname(candidate[[1L]]) else "."
  file.path(dir_path, "..", "online_twfe_core.R")
})

if (!exists("vec", mode = "function"))   # 이미 로드돼 있으면 재로드 생략
  source(.sim_core_path, local = FALSE)

suppressPackageStartupMessages(library(plm))


# =============================================================================
# Simulation 전용: 고정 T=3, balanced panel 빠른 경로
# =============================================================================

#' Algorithm 1 fast path: T=3 balanced panel, p=4 고정
#'
#' kronecker() 구성을 피하고 solve(H_i, ...) 대신 행렬 직접 연산으로 최적화.
#' 알고리즘 결과는 alg1_new_unit()과 동일.
alg1_new_unit_fast_T3 <- function(S_pre, x1_i, x2_i, y_i) {
  # 일반 경로로 fallback하는 조건
  if (!isTRUE(S_pre$T_support == 3L) || !isTRUE(S_pre$p == 4L) ||
      length(x1_i) != 3L) {
    x_mat_i <- cbind(x1 = x1_i, x2 = x2_i)
    return(alg1_new_unit(S_pre, x_mat_i, 1:3, y_i))
  }

  # T=3 within 변환 (time = 1,2,3 balanced → dotD2, dotD3 고정)
  x1_dot <- x1_i - mean(x1_i);  x2_dot <- x2_i - mean(x2_i)
  y_dot  <- y_i  - mean(y_i)
  dotD2  <- c(-1/3,  2/3, -1/3);  dotD3 <- c(-1/3, -1/3,  2/3)
  dotZ_i <- cbind(x1 = x1_dot, x2 = x2_dot,
                  `factor(time)2` = dotD2, `factor(time)3` = dotD3)
  dotY_i <- y_dot

  inv_pre   <- S_pre$inv_dotZtZ
  theta_pre <- S_pre$theta_hat
  sig2_pre  <- S_pre$sigma2_hat
  Vcr_pre   <- S_pre$Vcr_hat
  A_pre     <- S_pre$A_N
  B_pre     <- S_pre$B_N
  p <- 4L;  T_i <- 3L

  H_i       <- diag(T_i) + dotZ_i %*% inv_pre %*% t(dotZ_i)
  e_tilde_i <- drop(dotY_i - dotZ_i %*% theta_pre)
  sol       <- solve(H_i, cbind(dotZ_i, e_tilde_i))
  W_i       <- sol[, 1:4, drop = FALSE]
  v_i       <- sol[, 5]
  G_i       <- crossprod(dotZ_i, W_i)
  q_i       <- drop(crossprod(dotZ_i, v_i))

  # Prop 1
  d_i        <- drop(inv_pre %*% q_i)
  theta_post <- drop(theta_pre + d_i)
  names(theta_post) <- names(theta_pre)

  # Prop 2
  inv_post <- inv_pre - inv_pre %*% G_i %*% inv_pre
  n_post   <- S_pre$n + T_i;  N_post <- S_pre$N + 1L
  df_pre   <- S_pre$n - S_pre$N - p
  df_post  <- n_post - N_post - p
  sig2_post <- (df_pre / df_post) * sig2_pre +
    as.numeric(crossprod(e_tilde_i, v_i)) / df_post

  # Prop 3 (p=4 고정 → kronecker() 없이 수동 인덱싱)
  I_p   <- diag(p)
  U_i   <- I_p - inv_pre %*% G_i
  th    <- as.numeric(theta_pre);  d <- as.numeric(d_i)

  KbetaB <- th[1]*B_pre[ 1:4, ] + th[2]*B_pre[ 5:8, ] +
             th[3]*B_pre[ 9:12,] + th[4]*B_pre[13:16,]
  X_mat  <- A_pre - KbetaB

  Omega1        <- matrix(0, p, p)
  Omega1[, 1]   <- -X_mat[,  1:4 ] %*% d
  Omega1[, 2]   <- -X_mat[,  5:8 ] %*% d
  Omega1[, 3]   <- -X_mat[,  9:12] %*% d
  Omega1[, 4]   <- -X_mat[, 13:16] %*% d
  Omega2        <- t(Omega1)

  Bd <- matrix(0, p^2, p)
  Bd[, 1] <- B_pre[,  1:4 ] %*% d;  Bd[, 2] <- B_pre[,  5:8 ] %*% d
  Bd[, 3] <- B_pre[,  9:12] %*% d;  Bd[, 4] <- B_pre[, 13:16] %*% d
  Omega3 <- d[1]*Bd[ 1:4,] + d[2]*Bd[ 5:8,] +
             d[3]*Bd[ 9:12,] + d[4]*Bd[13:16,]
  Omega4 <- tcrossprod(q_i)

  Vcr_post <- U_i %*%
    (Vcr_pre + inv_pre %*% (Omega1 + Omega2 + Omega3 + Omega4) %*% inv_pre) %*%
    t(U_i)

  S_i    <- crossprod(dotZ_i)
  s_i    <- crossprod(dotZ_i, matrix(dotY_i, ncol = 1))
  vecS_i <- matrix(as.vector(S_i), ncol = 1)

  S_post <- S_pre
  S_post$N          <- N_post
  S_post$n          <- n_post
  S_post$inv_dotZtZ <- inv_post
  S_post$theta_hat  <- theta_post
  S_post$sigma2_hat <- sig2_post
  S_post$Vcr_hat    <- Vcr_post
  S_post$M_hat      <- NULL
  S_post$A_N        <- A_pre + s_i %*% t(vecS_i)
  S_post$B_N        <- B_pre + vecS_i %*% t(vecS_i)
  S_post
}


# =============================================================================
# Simulation 전용: warm-up 초기화 (DGP 인자 포함)
# =============================================================================

#' Simulation용 warm-up 초기화 (state_init 래퍼)
#'
#' @param df_warm   warm-up data.frame (컬럼: id, time, y, x1, x2)
#' @param N_warm    warm-up unit 수
#' @param alpha_warm unit FE 벡터 (unit 요약 객체에 alpha 필드 추가용)
#' @param lambda_t   time FE 벡터 (state에 저장, 이후 데이터 생성에 사용)
init_state_sim <- function(df_warm, N_warm, alpha_warm, lambda_t) {
  stopifnot(all(df_warm$time %in% 1:2))
  stopifnot(length(alpha_warm) == N_warm)

  init  <- state_init(df_warm, "id", "time", "y", c("x1", "x2"),
                      baseline_time = 1L)
  S_pre <- init$S_pre
  units <- init$units

  # simulation 전용 필드 추가
  S_pre$lambda_t <- lambda_t

  # unit 객체에 alpha (true FE) 추가 (검증용)
  uid_vec <- unique(df_warm$id)
  for (k in seq_along(uid_vec)) {
    uid_ch <- as.character(uid_vec[k])
    units[[uid_ch]]$alpha <- alpha_warm[k]
  }

  list(S_pre = S_pre, units = units)
}


# =============================================================================
# run_one_N: N개 unit의 balanced TWFE simulation (T=3)
# =============================================================================

#' N개 unit, T=3 balanced panel에 대한 온라인 TWFE 벤치마크
#'
#' @param N                 총 unit 수
#' @param N_warm            warm-up unit 수 (default 5)
#' @param beta_true         true 계수 벡터 (길이 2)
#' @param sd_alpha          unit FE 표준편차
#' @param sd_lambda         time FE 표준편차
#' @param sd_u              오차항 표준편차
#' @param seed              random seed
#' @param do_offline        offline (plm) 비교 실행 여부
#' @param offline_maxN_attempt offline 시도 최대 N (이를 초과하면 skip)
#' @param verify_smallN     소규모에서 수치 검증 수행 여부
#' @param verify_maxN       검증 수행 최대 N
#' @param verbose           진행 출력 여부
#' @param alg1_batch_size   Algorithm 1 타이밍 배치 크기
#' @return 결과 리스트
run_one_N <- function(N,
                      N_warm = 5L,
                      beta_true = c(1.3, -0.7),
                      sd_alpha = 1.0,
                      sd_lambda = 1.0,
                      sd_u = 1.0,
                      seed = 123,
                      do_offline = TRUE,
                      offline_maxN_attempt = Inf,
                      verify_smallN = TRUE,
                      verify_maxN = 100000L,
                      verbose = TRUE,
                      alg1_batch_size = 50000L) {

  stopifnot(N >= N_warm)
  set.seed(seed)

  offline_requested <- isTRUE(do_offline) && (N <= offline_maxN_attempt)
  store_for_offline <- offline_requested
  store_for_verif   <- isTRUE(verify_smallN) && (N <= verify_maxN)
  store_data        <- store_for_offline || store_for_verif

  storage_failed <- FALSE;  storage_fail_msg <- NA_character_
  x1_store <- x2_store <- y_store <- filled <- NULL

  if (store_data) {
    ok_alloc <- tryCatch({
      x1_store <- numeric(3L * N);  x2_store <- numeric(3L * N)
      y_store  <- numeric(3L * N);  filled   <- logical(3L * N)
      TRUE
    }, error = function(e) {
      storage_failed   <<- TRUE;  storage_fail_msg <<- conditionMessage(e);  FALSE
    })
    if (!isTRUE(ok_alloc)) {
      store_for_offline <- FALSE;  store_for_verif <- FALSE;  store_data <- FALSE
    }
  }

  idx_fun <- function(i, t) (i - 1L) * 3L + t

  # ------------------------------------------------------------------
  # Warm-up 데이터 생성 (타이밍 제외)
  # ------------------------------------------------------------------
  lambda_t   <- rnorm(3L, sd = sd_lambda)
  alpha_warm <- rnorm(N_warm, sd = sd_alpha)

  id_w   <- rep(seq_len(N_warm), each = 2L)
  time_w <- rep.int(1:2, times = N_warm)
  n_w    <- length(id_w)
  x1_w   <- rnorm(n_w);  x2_w <- rnorm(n_w);  u_w <- rnorm(n_w, sd = sd_u)
  y_w    <- beta_true[1]*x1_w + beta_true[2]*x2_w +
            alpha_warm[id_w] + lambda_t[time_w] + u_w

  df_warm <- data.frame(id = id_w, time = time_w, y = y_w, x1 = x1_w, x2 = x2_w)

  if (store_data) {
    for (k in seq_len(n_w)) {
      idx            <- idx_fun(id_w[k], time_w[k])
      x1_store[idx]  <- x1_w[k];  x2_store[idx] <- x2_w[k]
      y_store[idx]   <- y_w[k];   filled[idx]   <- TRUE
    }
  }

  # ------------------------------------------------------------------
  # 초기 state 구성
  # ------------------------------------------------------------------
  t_init_build <- system.time({
    init  <- init_state_sim(df_warm, N_warm, alpha_warm, lambda_t)
    S_pre <- init$S_pre
    units <- init$units
  })[["elapsed"]]

  # ------------------------------------------------------------------
  # 스트리밍 알고리즘 업데이트 (타이밍 대상)
  # ------------------------------------------------------------------
  t_alg1 <- t_alg2 <- t_alg3 <- t_expand <- 0

  # Event 1: unit 1이 t=3을 받음 → Algorithm 3 (새 calendar time)
  x1_new <- rnorm(1);  x2_new <- rnorm(1);  u_new <- rnorm(1, sd = sd_u)
  y_new  <- beta_true[1]*x1_new + beta_true[2]*x2_new +
            units[["1"]]$alpha + lambda_t[3] + u_new

  if (store_data) {
    idx <- idx_fun(1L, 3L)
    x1_store[idx] <- x1_new;  x2_store[idx] <- x2_new
    y_store[idx]  <- y_new;   filled[idx]   <- TRUE
  }

  t_alg3 <- system.time({
    out3       <- alg3_new_caltime(S_pre, units[["1"]],
                                   c(x1 = x1_new, x2 = x2_new), 3L, y_new)
    S_pre      <- out3$S_post
    units[["1"]] <- out3$unit_post
  })[["elapsed"]]

  # warm-up unit 2..N_warm의 요약을 zero-extend
  new_dummy_name <- "factor(time)3"
  t_expand <- system.time({
    for (j in 2:N_warm) {
      uid_ch        <- as.character(j)
      units[[uid_ch]] <- expand_unit_state_zero(units[[uid_ch]], new_dummy_name)
    }
  })[["elapsed"]]

  # Event 2: unit 2..N_warm이 t=3을 받음 → Algorithm 2 (기존 calendar time)
  for (j in 2:N_warm) {
    uid_ch <- as.character(j)
    x1_new <- rnorm(1);  x2_new <- rnorm(1);  u_new <- rnorm(1, sd = sd_u)
    y_new  <- beta_true[1]*x1_new + beta_true[2]*x2_new +
              units[[uid_ch]]$alpha + lambda_t[3] + u_new

    if (store_data) {
      idx <- idx_fun(j, 3L)
      x1_store[idx] <- x1_new;  x2_store[idx] <- x2_new
      y_store[idx]  <- y_new;   filled[idx]   <- TRUE
    }

    t_alg2 <- t_alg2 + system.time({
      out2          <- alg2_existing_unit(S_pre, units[[uid_ch]],
                                          c(x1 = x1_new, x2 = x2_new), 3L, y_new)
      S_pre         <- out2$S_post
      units[[uid_ch]] <- out2$unit_post
    })[["elapsed"]]
  }

  # Event 3: 신규 unit (N_warm+1)..N 도착 → Algorithm 1 (배치 처리)
  if (N > N_warm) {
    n_new      <- N - N_warm
    batch_size <- as.integer(max(1L, alg1_batch_size))
    i_start    <- N_warm + 1L

    last_pct <- -1L
    pb_step  <- as.integer(max(1L, ceiling(as.numeric(n_new) / 200)))
    progress_tick <- function(i_done) {
      if (!isTRUE(verbose)) return(invisible(NULL))
      if ((i_done %% pb_step) == 0L || i_done == n_new) {
        pct <- as.integer(floor(100 * as.numeric(i_done) / as.numeric(n_new)))
        if (!is.na(pct) && pct != last_pct) {
          cat(sprintf("\r  |%3d%%|", pct));  flush.console();  last_pct <<- pct
        }
      }
      invisible(NULL)
    }

    while (i_start <= N) {
      i_end <- min(N, i_start + batch_size - 1L)
      m     <- i_end - i_start + 1L

      # 배치 데이터 생성 (타이밍 제외)
      alpha_b <- rnorm(m, sd = sd_alpha)
      X1 <- matrix(rnorm(3L * m), nrow = m, ncol = 3L, byrow = TRUE)
      X2 <- matrix(rnorm(3L * m), nrow = m, ncol = 3L, byrow = TRUE)
      U  <- matrix(rnorm(3L * m, sd = sd_u), nrow = m, ncol = 3L, byrow = TRUE)
      LAM <- matrix(lambda_t, nrow = m, ncol = 3L, byrow = TRUE)
      Y   <- beta_true[1]*X1 + beta_true[2]*X2 + alpha_b + LAM + U

      if (store_data) {
        base <- (i_start - 1L) * 3L + 1L
        idxs <- base:(base + 3L * m - 1L)
        x1_store[idxs] <- as.vector(t(X1));  x2_store[idxs] <- as.vector(t(X2))
        y_store[idxs]  <- as.vector(t(Y));   filled[idxs]   <- TRUE
      }

      # 배치 내 unit별 Algorithm 1 적용 (타이밍 대상)
      t_alg1 <- t_alg1 + system.time({
        for (k in seq_len(m))
          S_pre <- alg1_new_unit_fast_T3(S_pre, X1[k, ], X2[k, ], Y[k, ])
      })[["elapsed"]]

      progress_tick(i_end - N_warm)
      i_start <- i_end + 1L
    }
    if (verbose) cat("\n")
  }

  t_upd       <- t_alg1 + t_alg2 + t_alg3 + t_expand
  bytes_state <- sum_object_size(S_pre) + sum_object_size(units)

  # ------------------------------------------------------------------
  # Offline 참조값 계산 (동일 데이터에 plm 적용)
  # ------------------------------------------------------------------
  offline_ok <- offline_ok_fit <- offline_ok_V0 <- offline_ok_Vcr <- NA
  t_off_build_df <- t_off_pdata <- t_off_fit <- t_off_V0 <- t_off_Vcr <- t_off_total <- NA_real_
  bytes_off_df   <- NA_real_
  offline_fail_stage <- offline_fail_msg <- NA_character_
  coef_plm <- V0_plm <- Vcr_plm <- NULL

  if (isTRUE(offline_requested) && isTRUE(storage_failed)) {
    offline_ok         <- FALSE
    offline_fail_stage <- "allocate storage"
    offline_fail_msg   <- storage_fail_msg
  }

  if (isTRUE(store_for_offline) && !isTRUE(storage_failed)) {
    if (!all(filled)) stop("Storage error: not all (id,time) cells were filled.")
    offline_ok <- FALSE

    res_build_df <- safe_system_time({
      data.frame(
        id   = factor(rep(seq_len(N), each = 3L)),
        time = factor(rep.int(1:3, times = N), levels = 1:3),
        y = y_store, x1 = x1_store, x2 = x2_store
      )
    }, label = "build df_stream", verbose_error = TRUE)

    if (!isTRUE(res_build_df$ok)) {
      offline_fail_stage <- "build df_stream"
      offline_fail_msg   <- conditionMessage(res_build_df$error)
    } else {
      df_stream      <- res_build_df$value
      t_off_build_df <- res_build_df$time
      bytes_off_df   <- as.numeric(object.size(df_stream))

      res_pdata <- safe_system_time(
        pdata.frame(df_stream, index = c("id", "time")),
        label = "pdata.frame", verbose_error = TRUE
      )

      if (!isTRUE(res_pdata$ok)) {
        offline_fail_stage <- "pdata.frame"
        offline_fail_msg   <- conditionMessage(res_pdata$error)
      } else {
        df_p        <- res_pdata$value
        t_off_pdata <- res_pdata$time

        res_fit <- safe_system_time(
          plm(y ~ x1 + x2 + factor(time), data = df_p,
              model = "within", effect = "individual"),
          label = "plm fit", verbose_error = TRUE
        )

        if (!isTRUE(res_fit$ok)) {
          offline_ok_fit     <- FALSE
          offline_fail_stage <- "plm fit"
          offline_fail_msg   <- conditionMessage(res_fit$error)
        } else {
          m_fit          <- res_fit$value
          coef_plm       <- coef(m_fit)
          t_off_fit      <- res_fit$time
          offline_ok_fit <- TRUE

          res_V0 <- safe_system_time(vcov(m_fit), label = "vcov (V0)",
                                     verbose_error = TRUE)
          if (!isTRUE(res_V0$ok)) {
            offline_ok_V0      <- FALSE
            offline_fail_stage <- "vcov (V0)"
            offline_fail_msg   <- conditionMessage(res_V0$error)
          } else {
            V0_plm        <- res_V0$value
            t_off_V0      <- res_V0$time
            offline_ok_V0 <- TRUE

            res_Vcr <- safe_system_time(
              vcovHC(m_fit, type = "HC0", method = "arellano", cluster = "group"),
              label = "vcovHC (Vcr)", verbose_error = TRUE
            )
            if (!isTRUE(res_Vcr$ok)) {
              offline_ok_Vcr     <- FALSE
              offline_fail_stage <- "vcovHC (Vcr)"
              offline_fail_msg   <- conditionMessage(res_Vcr$error)
            } else {
              Vcr_plm        <- res_Vcr$value
              t_off_Vcr      <- res_Vcr$time
              offline_ok_Vcr <- TRUE
              offline_ok     <- TRUE
              t_off_total    <- t_off_build_df + t_off_pdata + t_off_fit +
                                t_off_V0 + t_off_Vcr
            }
          }
          rm(m_fit)
        }
        rm(df_p)
      }
      rm(df_stream)
    }
    gc()
  }

  # ------------------------------------------------------------------
  # 수치 검증: online vs plm
  # ------------------------------------------------------------------
  verif <- NULL
  if (store_for_verif) {
    if (!all(filled)) stop("Storage error: not all (id,time) cells were filled.")

    if (is.null(coef_plm) || is.null(V0_plm) || is.null(Vcr_plm)) {
      df_stream <- data.frame(
        id   = factor(rep(seq_len(N), each = 3L)),
        time = factor(rep.int(1:3, times = N), levels = 1:3),
        y = y_store, x1 = x1_store, x2 = x2_store
      )
      df_p     <- pdata.frame(df_stream, index = c("id", "time"))
      m_ref    <- plm(y ~ x1 + x2 + factor(time), data = df_p,
                      model = "within", effect = "individual")
      coef_plm <- coef(m_ref)
      V0_plm   <- vcov(m_ref)
      Vcr_plm  <- vcovHC(m_ref, type = "HC0", method = "arellano", cluster = "group")
      if (is.null(colnames(V0_plm)))
        colnames(V0_plm) <- rownames(V0_plm) <- names(coef_plm)
      if (is.null(colnames(Vcr_plm)))
        colnames(Vcr_plm) <- rownames(Vcr_plm) <- names(coef_plm)
    }

    theta_online <- S_pre$theta_hat
    nm           <- intersect(names(theta_online), names(coef_plm))
    theta_max_abs <- max(abs(theta_online[nm] - coef_plm[nm]))

    V0_online <- S_pre$sigma2_hat * S_pre$inv_dotZtZ
    colnames(V0_online) <- rownames(V0_online) <- names(theta_online)
    rn0 <- intersect(rownames(V0_online), rownames(V0_plm))
    cn0 <- intersect(colnames(V0_online), colnames(V0_plm))
    V0_max_abs  <- max(abs(V0_online[rn0, cn0] - V0_plm[rn0, cn0]))

    Vcr_online <- S_pre$Vcr_hat
    colnames(Vcr_online) <- rownames(Vcr_online) <- names(theta_online)
    rnc <- intersect(rownames(Vcr_online), rownames(Vcr_plm))
    cnc <- intersect(colnames(Vcr_online), colnames(Vcr_plm))
    Vcr_max_abs <- max(abs(Vcr_online[rnc, cnc] - Vcr_plm[rnc, cnc]))

    verif <- list(
      theta_max_abs = theta_max_abs,
      V0_max_abs    = V0_max_abs,
      Vcr_max_abs   = Vcr_max_abs,
      n_common      = length(nm)
    )
  }

  list(
    N      = N,
    online = list(
      ok           = TRUE,
      t_init_build = t_init_build,
      t_update     = t_upd,
      t_alg1       = t_alg1,
      t_alg2       = t_alg2,
      t_alg3       = t_alg3,
      t_expand     = t_expand,
      t_total      = t_upd,
      bytes_state  = bytes_state,
      theta_hat    = S_pre$theta_hat
    ),
    offline = list(
      attempted  = offline_requested,
      ok         = offline_ok,
      ok_fit     = offline_ok_fit,
      ok_V0      = offline_ok_V0,
      ok_Vcr     = offline_ok_Vcr,
      t_build_df = t_off_build_df,
      t_pdata    = t_off_pdata,
      t_fit      = t_off_fit,
      t_V0       = t_off_V0,
      t_Vcr      = t_off_Vcr,
      t_total    = if (isTRUE(offline_ok)) t_off_total else NA_real_,
      bytes_df   = bytes_off_df,
      fail_stage = offline_fail_stage,
      fail_msg   = offline_fail_msg
    ),
    verification = verif
  )
}


# =============================================================================
# run_benchmark_grid: N 값 그리드에 대한 벤치마크 실행
# =============================================================================

#' N 값 그리드 벤치마크
#'
#' @param N_grid              N 값 벡터
#' @param N_warm              warm-up unit 수
#' @param seed                random seed
#' @param do_offline          offline 비교 여부
#' @param offline_maxN_attempt offline 시도 최대 N
#' @param verify_smallN       소규모 검증 여부
#' @param verify_maxN         검증 최대 N
#' @param alg1_batch_size     Algorithm 1 배치 크기
#' @param verbose             출력 여부
#' @param out_csv             결과 저장 CSV 경로 (NULL이면 저장 안 함)
#' @return list(raw, table)
run_benchmark_grid <- function(N_grid,
                               N_warm = 5L,
                               seed = 123,
                               do_offline = TRUE,
                               offline_maxN_attempt = Inf,
                               verify_smallN = TRUE,
                               verify_maxN = 100000L,
                               alg1_batch_size = 50000L,
                               verbose = TRUE,
                               out_csv = NULL) {

  out <- vector("list", length(N_grid))

  for (g in seq_along(N_grid)) {
    N <- as.integer(N_grid[g])
    cat("\n====================================\n")
    cat("Running N =", format(N, big.mark = ","), "\n")
    cat("====================================\n")

    res <- run_one_N(
      N = N, N_warm = N_warm, seed = seed,
      do_offline = do_offline,
      offline_maxN_attempt = offline_maxN_attempt,
      verify_smallN = verify_smallN,
      verify_maxN = verify_maxN,
      verbose = verbose,
      alg1_batch_size = alg1_batch_size
    )
    out[[g]] <- res

    cat(sprintf(
      "Online (update only): %.3f sec (Alg1 %.3f / Alg2 %.3f / Alg3 %.3f + expand %.3f; init-build %.3f)\n",
      res$online$t_total, res$online$t_alg1, res$online$t_alg2,
      res$online$t_alg3, res$online$t_expand, res$online$t_init_build))
    cat(sprintf("Online state size: %.2f MB\n", res$online$bytes_state / 1024^2))

    if (res$offline$attempted) {
      if (isTRUE(res$offline$ok)) {
        cat(sprintf(
          "Offline (same data): builddf %.3f / pdata %.3f / fit %.3f / V0 %.3f / Vcr %.3f (total %.3f)\n",
          res$offline$t_build_df, res$offline$t_pdata, res$offline$t_fit,
          res$offline$t_V0, res$offline$t_Vcr, res$offline$t_total))
        cat(sprintf("Offline df size: %.2f MB\n", res$offline$bytes_df / 1024^2))
      } else {
        cat(sprintf("Offline (same data): FAILED (%s: %s)\n",
                    res$offline$fail_stage, res$offline$fail_msg))
      }
    } else {
      cat("Offline: SKIPPED (N too large)\n")
    }

    if (!is.null(res$verification)) {
      cat(sprintf(
        "Verification (max abs): theta %.2e / V0 %.2e / Vcr %.2e (n=%d)\n",
        res$verification$theta_max_abs, res$verification$V0_max_abs,
        res$verification$Vcr_max_abs, res$verification$n_common))
    }
  }

  tab <- data.frame(
    N                   = vapply(out, function(z) z$N, integer(1)),
    online_t_upd        = vapply(out, function(z) z$online$t_update,   numeric(1)),
    online_t_alg1       = vapply(out, function(z) z$online$t_alg1,     numeric(1)),
    online_t_alg2       = vapply(out, function(z) z$online$t_alg2,     numeric(1)),
    online_t_alg3       = vapply(out, function(z) z$online$t_alg3,     numeric(1)),
    online_t_expand     = vapply(out, function(z) z$online$t_expand,   numeric(1)),
    online_bytes_state  = vapply(out, function(z) z$online$bytes_state, numeric(1)),
    offline_attempted   = vapply(out, function(z) z$offline$attempted, logical(1)),
    offline_ok          = vapply(out, function(z) isTRUE(z$offline$ok), logical(1)),
    offline_t_build_df  = vapply(out, function(z) z$offline$t_build_df, numeric(1)),
    offline_t_pdata     = vapply(out, function(z) z$offline$t_pdata,   numeric(1)),
    offline_t_fit       = vapply(out, function(z) z$offline$t_fit,     numeric(1)),
    offline_t_V0        = vapply(out, function(z) z$offline$t_V0,      numeric(1)),
    offline_t_Vcr       = vapply(out, function(z) z$offline$t_Vcr,     numeric(1)),
    offline_t_total     = vapply(out, function(z) z$offline$t_total,   numeric(1)),
    offline_bytes_df    = vapply(out, function(z) z$offline$bytes_df,  numeric(1)),
    offline_fail_stage  = vapply(out, function(z)
                            if (is.null(z$offline$fail_stage)) NA_character_
                            else z$offline$fail_stage, character(1)),
    offline_fail_msg    = vapply(out, function(z)
                            if (is.null(z$offline$fail_msg)) NA_character_
                            else z$offline$fail_msg, character(1)),
    verif_theta_max_abs = vapply(out, function(z)
                            if (is.null(z$verification)) NA_real_
                            else z$verification$theta_max_abs, numeric(1)),
    verif_V0_max_abs    = vapply(out, function(z)
                            if (is.null(z$verification)) NA_real_
                            else z$verification$V0_max_abs, numeric(1)),
    verif_Vcr_max_abs   = vapply(out, function(z)
                            if (is.null(z$verification)) NA_real_
                            else z$verification$Vcr_max_abs, numeric(1))
  )

  if (!is.null(out_csv)) write.csv(tab, out_csv, row.names = FALSE)
  list(raw = out, table = tab)
}
