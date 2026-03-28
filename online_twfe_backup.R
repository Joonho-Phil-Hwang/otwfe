options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(plm)
})

# ============================================================
# Helpers
# ============================================================
vec <- function(M) as.vector(M)

# Embed indices when expanding vec() from p_old to p_new = p_old + 1
# vec() stacks columns, so old entries map to positions r + (c-1)*p_new in the new vector.
vec_embed_index <- function(p_old) {
  p_new <- p_old + 1L
  idx <- integer(p_old * p_old)
  k <- 0L
  for (c in seq_len(p_old)) {
    for (r in seq_len(p_old)) {
      k <- k + 1L
      idx[k] <- r + (c - 1L) * p_new
    }
  }
  idx
}

sum_object_size <- function(x) {
  sum(vapply(x, function(z) as.numeric(object.size(z)), numeric(1)))
}

safe_system_time <- function(expr, label = NULL, verbose_error = FALSE) {
  res <- NULL
  tt <- system.time({
    res <- tryCatch(
      list(ok = TRUE, value = eval.parent(substitute(expr)), error = NULL),
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
    if (!is.null(label)) {
      cat(sprintf("  [OFFLINE ERROR] %s: %s\n", label, msg))
    } else {
      cat(sprintf("  [OFFLINE ERROR] %s\n", msg))
    }
  }
  out
}

# Progress bar helpers: avoid integer overflow in setTxtProgressBar for huge N
make_progress_bar <- function(min_val, max_val, style = 3) {
  txtProgressBar(min = as.numeric(min_val), max = as.numeric(max_val), style = style)
}

safe_set_progress <- function(pb, value) {
  if (is.null(pb)) return(invisible(NULL))
  tryCatch({
    suppressWarnings(setTxtProgressBar(pb, as.numeric(value)))
  }, error = function(e) invisible(NULL))
  invisible(NULL)
}


# ============================================================
# Build raw design matrix Z_it = [x1_it, x2_it, time dummies]
# Baseline time (t=1) is dropped. Dummy names match plm:
#   factor(time)2, factor(time)3
# ============================================================
build_Z_raw_TWFE <- function(x1, x2, time, T_support) {
  stopifnot(T_support %in% c(2L, 3L))
  D2 <- as.numeric(time == 2L)
  if (T_support == 2L) {
    Z_raw <- cbind(x1 = x1, x2 = x2, `factor(time)2` = D2)
    return(Z_raw)
  }
  D3 <- as.numeric(time == 3L)
  Z_raw <- cbind(x1 = x1, x2 = x2, `factor(time)2` = D2, `factor(time)3` = D3)
  Z_raw
}

# ============================================================
# Within transform by unit i:
#   \bar{Z}_i = mean_t Z_it,  \bar{Y}_i = mean_t Y_it
#   \dot{Z}_it = Z_it - \bar{Z}_i,  \dot{Y}_it = Y_it - \bar{Y}_i
# ============================================================
within_transform_unit <- function(Z_raw_i, y_i) {
  barZ_i <- colMeans(Z_raw_i)
  barY_i <- mean(y_i)
  dotZ_i <- sweep(Z_raw_i, 2, barZ_i, "-")
  dotY_i <- y_i - barY_i
  list(dotZ_i = dotZ_i, dotY_i = dotY_i, barZ_i = barZ_i, barY_i = barY_i)
}

# ============================================================
# Initialize state from PROVIDED warm-up data:
#   warm-up: units i=1..N_warm with t=1,2 only (T_support=2)
# ============================================================
init_state_from_warmup_df <- function(df_warm, N_warm, alpha_warm, lambda_t) {
  T_support <- 2L

  stopifnot(all(df_warm$time %in% 1:2))
  stopifnot(length(alpha_warm) == N_warm)
  stopifnot(length(lambda_t) == 3L)

  id   <- df_warm$id
  time <- df_warm$time
  y    <- df_warm$y
  x1   <- df_warm$x1
  x2   <- df_warm$x2

  Z_raw <- build_Z_raw_TWFE(x1, x2, time, T_support = T_support)

  n <- nrow(df_warm)
  p <- ncol(Z_raw)

  dotZ <- matrix(0, nrow = n, ncol = p)
  colnames(dotZ) <- colnames(Z_raw)
  dotY <- numeric(n)

  units <- vector("list", N_warm)

  for (i in seq_len(N_warm)) {
    idx_i <- which(id == i)
    wt_i <- within_transform_unit(Z_raw[idx_i, , drop = FALSE], y[idx_i])

    dotZ[idx_i, ] <- wt_i$dotZ_i
    dotY[idx_i]   <- wt_i$dotY_i

    S_i <- crossprod(wt_i$dotZ_i)
    s_i <- crossprod(wt_i$dotZ_i, matrix(wt_i$dotY_i, ncol = 1))

    units[[i]] <- list(
      id     = i,
      alpha  = alpha_warm[i],
      T_i    = length(idx_i),
      barZ_i = wt_i$barZ_i,
      barY_i = wt_i$barY_i,
      S_i    = S_i,
      s_i    = s_i
    )
  }

  dotZtZ     <- crossprod(dotZ)
  inv_dotZtZ <- solve(dotZtZ)
  dotZtY     <- crossprod(dotZ, dotY)

  theta_hat <- drop(inv_dotZtZ %*% dotZtY)
  names(theta_hat) <- colnames(Z_raw)

  eps_hat <- drop(dotY - dotZ %*% theta_hat)

  N_pre  <- N_warm
  n_pre  <- n
  df_pre <- n_pre - N_pre - p
  sigma2_hat <- sum(eps_hat^2) / df_pre

  # Cluster-robust "cluster score" aggregation:
  #   r_i = \dot{Z}_i' \hat{e}_i,  M = sum_i r_i r_i'
  M_hat <- matrix(0, nrow = p, ncol = p, dimnames = list(colnames(Z_raw), colnames(Z_raw)))
  for (i in seq_len(N_warm)) {
    idx_i <- which(id == i)
    r_i   <- crossprod(dotZ[idx_i, , drop = FALSE], matrix(eps_hat[idx_i], ncol = 1))
    M_hat <- M_hat + tcrossprod(r_i)
  }
  Vcr_hat <- inv_dotZtZ %*% M_hat %*% inv_dotZtZ

  # Aggregates for Prop 3/6:
  #   A_N = sum_i (\dot{Z}_i'\dot{Y}_i) vec(\dot{Z}_i'\dot{Z}_i)'
  #   B_N = sum_i vec(\dot{Z}_i'\dot{Z}_i) vec(\dot{Z}_i'\dot{Z}_i)'
  A_N <- matrix(0, nrow = p, ncol = p^2)
  B_N <- matrix(0, nrow = p^2, ncol = p^2)
  for (i in seq_len(N_warm)) {
    vecG_i <- matrix(vec(units[[i]]$S_i), ncol = 1)
    A_N    <- A_N + units[[i]]$s_i %*% t(vecG_i)
    B_N    <- B_N + vecG_i %*% t(vecG_i)
  }

  S_pre <- list(
    T_support  = T_support,
    p          = p,
    N          = N_pre,
    n          = n_pre,
    inv_dotZtZ = inv_dotZtZ,
    theta_hat  = theta_hat,
    sigma2_hat = sigma2_hat,
    Vcr_hat    = Vcr_hat,
    M_hat      = NULL,  # not carried through Alg 1/2; reconstructed only in Alg 3
    A_N        = A_N,
    B_N        = B_N,
    lambda_t   = lambda_t
  )

  list(S_pre = S_pre, units = units)
}

# ============================================================
# (Utility) Zero-extend a stored unit summary after a NEW time dummy is introduced
# (Dimension expansion p -> p+1; the new dummy is identically 0 pre-update)
# ============================================================
expand_unit_summary_zero <- function(unit, new_dummy_name) {
  unit$barZ_i <- c(unit$barZ_i, 0)
  names(unit$barZ_i)[length(unit$barZ_i)] <- new_dummy_name

  p_old <- nrow(unit$S_i)
  S_ext <- rbind(
    cbind(unit$S_i, matrix(0, nrow = p_old, ncol = 1)),
    cbind(matrix(0, nrow = 1, ncol = p_old), 0)
  )
  s_ext <- rbind(unit$s_i, 0)

  nms <- c(names(unit$barZ_i)[-length(unit$barZ_i)], new_dummy_name)
  colnames(S_ext) <- nms
  rownames(S_ext) <- nms
  rownames(s_ext) <- nms

  unit$S_i <- S_ext
  unit$s_i <- s_ext
  unit
}

# ============================================================
# Algorithm 1 (Props 1--3): NEW unit arrives with a block of observations
# Inputs: S_pre, (time_i, x1_i, x2_i, y_i) with time_i subset of {1,...,T_support}
# ============================================================
alg1_new_unit_block <- function(S_pre, time_i, x1_i, x2_i, y_i) {
  stopifnot(all(time_i %in% 1:S_pre$T_support))

  Z_raw_i <- build_Z_raw_TWFE(x1_i, x2_i, time_i, T_support = S_pre$T_support)
  wt_i    <- within_transform_unit(Z_raw_i, y_i)
  dotZ_i  <- wt_i$dotZ_i
  dotY_i  <- wt_i$dotY_i

  inv_dotZtZ_pre <- S_pre$inv_dotZtZ
  theta_hat_pre  <- S_pre$theta_hat
  sigma2_pre     <- S_pre$sigma2_hat
  Vcr_pre        <- S_pre$Vcr_hat
  A_N_pre        <- S_pre$A_N
  B_N_pre        <- S_pre$B_N
  p  <- S_pre$p
  T_i <- length(time_i)

  # H_i = I_{T_i} + \dot{Z}_i (\dot{Z}'\dot{Z})^{-1} \dot{Z}_i'
  H_i     <- diag(T_i) + dotZ_i %*% inv_dotZtZ_pre %*% t(dotZ_i)
  H_i_inv <- solve(H_i)

  e_tilde_i <- drop(dotY_i - dotZ_i %*% theta_hat_pre)

  # --- Prop 1: coefficient update
  d_i <- drop(inv_dotZtZ_pre %*% t(dotZ_i) %*% H_i_inv %*% e_tilde_i)
  theta_hat_post <- drop(theta_hat_pre + d_i)
  names(theta_hat_post) <- names(theta_hat_pre)

  # --- Prop 2: inverse update
  inv_dotZtZ_post <- inv_dotZtZ_pre -
    inv_dotZtZ_pre %*% t(dotZ_i) %*% H_i_inv %*% dotZ_i %*% inv_dotZtZ_pre

  # --- Prop 2: sigma^2 update (block case)
  n_post  <- S_pre$n + T_i
  N_post  <- S_pre$N + 1L
  df_pre  <- S_pre$n - S_pre$N - p
  df_post <- n_post - N_post - p

  sigma2_post <- (df_pre / df_post) * sigma2_pre +
    as.numeric(t(e_tilde_i) %*% H_i_inv %*% e_tilde_i) / df_post

  # --- Prop 3: Cluster-robust vcov update (uses A_N, B_N)
  I_p <- diag(p)
  U_new  <- I_p - inv_dotZtZ_pre %*% t(dotZ_i) %*% H_i_inv %*% dotZ_i
  K_beta <- kronecker(t(theta_hat_pre), I_p)
  K_Id   <- kronecker(I_p, matrix(d_i, ncol = 1))

  Omega1 <- -(A_N_pre - K_beta %*% B_N_pre) %*% K_Id
  Omega2 <- t(Omega1)
  Omega3 <- kronecker(t(d_i), I_p) %*% B_N_pre %*% K_Id
  Omega4 <- t(dotZ_i) %*% H_i_inv %*% tcrossprod(e_tilde_i) %*% H_i_inv %*% dotZ_i

  Vcr_post <- U_new %*%
    (Vcr_pre + inv_dotZtZ_pre %*% (Omega1 + Omega2 + Omega3 + Omega4) %*% inv_dotZtZ_pre) %*%
    t(U_new)

  # --- Update aggregates A_N, B_N with the NEW unit contribution
  S_i    <- crossprod(dotZ_i)
  s_i    <- crossprod(dotZ_i, matrix(dotY_i, ncol = 1))
  vecS_i <- matrix(vec(S_i), ncol = 1)

  S_post <- S_pre
  S_post$inv_dotZtZ <- inv_dotZtZ_post
  S_post$theta_hat  <- theta_hat_post
  S_post$sigma2_hat <- sigma2_post
  S_post$Vcr_hat    <- Vcr_post
  S_post$M_hat      <- NULL
  S_post$A_N        <- A_N_pre + s_i %*% t(vecS_i)
  S_post$B_N        <- B_N_pre + vecS_i %*% t(vecS_i)
  S_post$n          <- n_post
  S_post$N          <- N_post
  S_post
}

# Fast path for the big benchmark (T_support = 3, balanced, p = 4).
# Avoids repeated kronecker() construction and explicit H_i^{-1}.
# Algorithmic result is identical to alg1_new_unit_block().
alg1_new_unit_block_fast_T3 <- function(S_pre, x1_i, x2_i, y_i) {
  if (!isTRUE(S_pre$T_support == 3L) || !isTRUE(S_pre$p == 4L) ||
      length(x1_i) != 3L || length(x2_i) != 3L || length(y_i) != 3L) {
    return(alg1_new_unit_block(S_pre, 1:3, x1_i, x2_i, y_i))
  }

  # Fast within-transform (T=3, time=1:3)
  x1_dot <- x1_i - mean(x1_i)
  x2_dot <- x2_i - mean(x2_i)
  y_dot  <- y_i  - mean(y_i)

  dotD2 <- c(-1/3,  2/3, -1/3)
  dotD3 <- c(-1/3, -1/3,  2/3)

  dotZ_i <- cbind(x1_dot, x2_dot, dotD2, dotD3)
  dotY_i <- y_dot

  inv_dotZtZ_pre <- S_pre$inv_dotZtZ
  theta_hat_pre  <- S_pre$theta_hat
  sigma2_pre     <- S_pre$sigma2_hat
  Vcr_pre        <- S_pre$Vcr_hat
  A_N_pre        <- S_pre$A_N
  B_N_pre        <- S_pre$B_N
  p  <- 4L
  T_i <- 3L

  H_i       <- diag(T_i) + dotZ_i %*% inv_dotZtZ_pre %*% t(dotZ_i)
  e_tilde_i <- drop(dotY_i - dotZ_i %*% theta_hat_pre)

  # Solve H_i * [dotZ_i, e_tilde_i] to avoid forming H_i^{-1} explicitly
  sol <- solve(H_i, cbind(dotZ_i, e_tilde_i))
  W   <- sol[, 1:4, drop = FALSE]   # H^{-1} dotZ_i
  v   <- sol[, 5]                   # H^{-1} e_tilde_i

  G <- crossprod(dotZ_i, W)         # dotZ_i' H^{-1} dotZ_i
  q <- drop(crossprod(dotZ_i, v))   # dotZ_i' H^{-1} e_tilde_i

  # --- Prop 1
  d_i            <- drop(inv_dotZtZ_pre %*% q)
  theta_hat_post <- drop(theta_hat_pre + d_i)
  names(theta_hat_post) <- names(theta_hat_pre)

  # --- Prop 2
  inv_dotZtZ_post <- inv_dotZtZ_pre - inv_dotZtZ_pre %*% G %*% inv_dotZtZ_pre

  n_post  <- S_pre$n + T_i
  N_post  <- S_pre$N + 1L
  df_pre  <- S_pre$n - S_pre$N - p
  df_post <- n_post - N_post - p

  sigma2_post <- (df_pre / df_post) * sigma2_pre +
    as.numeric(crossprod(e_tilde_i, v)) / df_post

  # --- Prop 3: avoid kronecker() for fixed p = 4
  I_p   <- diag(p)
  U_new <- I_p - inv_dotZtZ_pre %*% G
  th    <- as.numeric(theta_hat_pre)
  d     <- as.numeric(d_i)

  # K_beta %*% B_N_pre  (sum over k of th[k] * corresponding block rows)
  KbetaB <- th[1]*B_N_pre[ 1:4,  ] + th[2]*B_N_pre[ 5:8,  ] +
            th[3]*B_N_pre[ 9:12, ] + th[4]*B_N_pre[13:16, ]

  X <- A_N_pre - KbetaB

  Omega1 <- matrix(0, p, p)
  Omega1[, 1] <- -X[,  1:4 ] %*% d
  Omega1[, 2] <- -X[,  5:8 ] %*% d
  Omega1[, 3] <- -X[,  9:12] %*% d
  Omega1[, 4] <- -X[, 13:16] %*% d
  Omega2 <- t(Omega1)

  Bd <- matrix(0, p^2, p)
  Bd[, 1] <- B_N_pre[,  1:4 ] %*% d
  Bd[, 2] <- B_N_pre[,  5:8 ] %*% d
  Bd[, 3] <- B_N_pre[,  9:12] %*% d
  Bd[, 4] <- B_N_pre[, 13:16] %*% d

  Omega3 <- d[1]*Bd[ 1:4, ] + d[2]*Bd[ 5:8, ] +
            d[3]*Bd[ 9:12,] + d[4]*Bd[13:16, ]
  Omega4 <- tcrossprod(q)

  Vcr_post <- U_new %*%
    (Vcr_pre + inv_dotZtZ_pre %*% (Omega1 + Omega2 + Omega3 + Omega4) %*% inv_dotZtZ_pre) %*%
    t(U_new)

  S_i    <- crossprod(dotZ_i)
  s_i    <- crossprod(dotZ_i, matrix(dotY_i, ncol = 1))
  vecS_i <- matrix(as.vector(S_i), ncol = 1)

  S_post <- S_pre
  S_post$inv_dotZtZ <- inv_dotZtZ_post
  S_post$theta_hat  <- theta_hat_post
  S_post$sigma2_hat <- sigma2_post
  S_post$Vcr_hat    <- Vcr_post
  S_post$M_hat      <- NULL
  S_post$A_N        <- A_N_pre + s_i %*% t(vecS_i)
  S_post$B_N        <- B_N_pre + vecS_i %*% t(vecS_i)
  S_post$n          <- n_post
  S_post$N          <- N_post
  S_post
}


# ============================================================
# Algorithm 2 (Props 4--6): EXISTING unit gets one NEW observation at an EXISTING time
# ============================================================
alg2_existing_unit_existing_time <- function(S_pre, unit_pre, time_new, x1_new, x2_new, y_new) {
  stopifnot(time_new %in% 1:S_pre$T_support)

  Z_new <- build_Z_raw_TWFE(x1_new, x2_new, time_new, T_support = S_pre$T_support)

  kappa    <- unit_pre$T_i / (unit_pre$T_i + 1)
  doty_new <- as.numeric(y_new - unit_pre$barY_i)
  dotx_new <- drop(Z_new - unit_pre$barZ_i)

  inv_dotZtZ_pre <- S_pre$inv_dotZtZ
  theta_hat_pre  <- S_pre$theta_hat
  sigma2_pre     <- S_pre$sigma2_hat
  Vcr_pre        <- S_pre$Vcr_hat
  A_N_pre        <- S_pre$A_N
  B_N_pre        <- S_pre$B_N
  p <- S_pre$p

  h        <- as.numeric(t(dotx_new) %*% inv_dotZtZ_pre %*% dotx_new)
  s_scalar <- as.numeric(kappa / (1 + kappa * h))
  e_tilde  <- as.numeric(doty_new - t(dotx_new) %*% theta_hat_pre)
  dotx_col <- matrix(dotx_new, ncol = 1)

  # --- Prop 4: coefficient update
  theta_hat_post <- drop(theta_hat_pre + s_scalar * (inv_dotZtZ_pre %*% dotx_col) * e_tilde)
  names(theta_hat_post) <- names(theta_hat_pre)

  # --- Prop 5: rank-one inverse and sigma^2 update
  inv_dotZtZ_post <- inv_dotZtZ_pre -
    s_scalar * (inv_dotZtZ_pre %*% dotx_col %*% t(dotx_col) %*% inv_dotZtZ_pre)

  n_post  <- S_pre$n + 1L
  N_post  <- S_pre$N
  df_pre  <- S_pre$n - S_pre$N - p
  df_post <- n_post - N_post - p

  sigma2_post <- (df_pre / df_post) * sigma2_pre +
    (e_tilde^2 / df_post) * (kappa / (1 + kappa * h))

  # --- Prop 6: cluster-robust vcov update
  I_p       <- diag(p)
  X_i_e_pre <- unit_pre$s_i - unit_pre$S_i %*% matrix(theta_hat_pre, ncol = 1)

  Omega1 <- e_tilde *
    (A_N_pre - kronecker(t(theta_hat_pre), I_p) %*% B_N_pre) %*%
    kronecker(I_p, s_scalar * (inv_dotZtZ_pre %*% dotx_col))

  Omega2 <- e_tilde^2 *
    kronecker(s_scalar * t(dotx_col) %*% inv_dotZtZ_pre, I_p) %*%
    B_N_pre %*%
    kronecker(I_p, s_scalar * (inv_dotZtZ_pre %*% dotx_col))

  Omega3 <- s_scalar * e_tilde *
    (X_i_e_pre %*% t(dotx_col) + dotx_col %*% t(X_i_e_pre))

  Omega4 <- s_scalar^2 * e_tilde^2 * (
    dotx_col %*% t(dotx_col) -
    dotx_col %*% t(dotx_col) %*% inv_dotZtZ_pre %*% unit_pre$S_i -
    unit_pre$S_i %*% inv_dotZtZ_pre %*% dotx_col %*% t(dotx_col)
  )

  U_rank1 <- I_p - s_scalar * (inv_dotZtZ_pre %*% dotx_col %*% t(dotx_col))

  Vcr_post <- U_rank1 %*%
    (Vcr_pre + inv_dotZtZ_pre %*% (-Omega1 - t(Omega1) + Omega2 + Omega3 + Omega4) %*% inv_dotZtZ_pre) %*%
    t(U_rank1)

  # --- Update this unit's summary
  G_old <- unit_pre$S_i
  s_old <- unit_pre$s_i

  unit_post        <- unit_pre
  unit_post$S_i    <- unit_pre$S_i + kappa * (dotx_col %*% t(dotx_col))
  unit_post$s_i    <- unit_pre$s_i + kappa * (dotx_col * doty_new)
  unit_post$barZ_i <- unit_pre$barZ_i + dotx_new / (unit_pre$T_i + 1)
  unit_post$barY_i <- unit_pre$barY_i + doty_new  / (unit_pre$T_i + 1)
  unit_post$T_i    <- unit_pre$T_i + 1L

  vecG_old <- matrix(vec(G_old), ncol = 1)
  vecG_new <- matrix(vec(unit_post$S_i), ncol = 1)

  S_post <- S_pre
  S_post$inv_dotZtZ <- inv_dotZtZ_post
  S_post$theta_hat  <- theta_hat_post
  S_post$sigma2_hat <- sigma2_post
  S_post$Vcr_hat    <- Vcr_post
  S_post$M_hat      <- NULL
  S_post$A_N        <- A_N_pre - (s_old %*% t(vecG_old)) + (unit_post$s_i %*% t(vecG_new))
  S_post$B_N        <- B_N_pre - (vecG_old %*% t(vecG_old)) + (vecG_new %*% t(vecG_new))
  S_post$n          <- n_post
  S_post$N          <- N_post

  list(S_post = S_post, unit_post = unit_post)
}

# ============================================================
# Algorithm 3 (Props 7--9): EXISTING unit gets a NEW observation at a NEW calendar time
# Expands p by 1. Benchmark-specific: support 2 -> 3 only.
# ============================================================
alg3_existing_unit_new_time <- function(S_pre, unit_pre, time_new, x1_new, x2_new, y_new) {
  stopifnot(time_new == S_pre$T_support + 1L)
  stopifnot(S_pre$T_support == 2L)

  Z_new_old <- build_Z_raw_TWFE(x1_new, x2_new, time_new, T_support = S_pre$T_support)

  kappa    <- unit_pre$T_i / (unit_pre$T_i + 1)
  doty_new <- as.numeric(y_new - unit_pre$barY_i)
  g_old    <- drop(Z_new_old - unit_pre$barZ_i)

  inv_dotZtZ_pre <- S_pre$inv_dotZtZ
  theta_hat_pre  <- S_pre$theta_hat
  Vcr_pre        <- S_pre$Vcr_hat
  A_N_pre        <- S_pre$A_N
  B_N_pre        <- S_pre$B_N
  p_pre <- S_pre$p

  h <- as.numeric(t(g_old) %*% inv_dotZtZ_pre %*% g_old)
  q <- drop(inv_dotZtZ_pre %*% g_old)

  new_dummy_name <- paste0("factor(time)", time_new)

  # --- Prop 7: expanded inverse and coefficients
  inv_dotZtZ_post <- rbind(
    cbind(inv_dotZtZ_pre, -matrix(q, ncol = 1)),
    cbind(-matrix(q, nrow = 1), (1 + kappa * h) / kappa)
  )
  colnames(inv_dotZtZ_post) <- c(names(theta_hat_pre), new_dummy_name)
  rownames(inv_dotZtZ_post) <- c(names(theta_hat_pre), new_dummy_name)

  v_new          <- as.numeric(doty_new - t(g_old) %*% theta_hat_pre)
  theta_hat_post <- c(theta_hat_pre, v_new)
  names(theta_hat_post) <- c(names(theta_hat_pre), new_dummy_name)

  # --- Prop 8: sigma^2 unchanged in this expansion step
  sigma2_post <- S_pre$sigma2_hat

  # Zero-extend A_N, B_N to new dimension
  p_post    <- p_pre + 1L
  idx_embed <- vec_embed_index(p_pre)
  A_N_ext   <- matrix(0, nrow = p_post, ncol = p_post^2)
  B_N_ext   <- matrix(0, nrow = p_post^2, ncol = p_post^2)
  A_N_ext[1:p_pre, idx_embed] <- A_N_pre
  B_N_ext[idx_embed, idx_embed] <- B_N_pre

  # Reconstruct M_pre for Prop 9 (not tracked by Alg 1/2)
  dotZtZ_pre <- solve(inv_dotZtZ_pre)
  M_pre      <- dotZtZ_pre %*% Vcr_pre %*% dotZtZ_pre
  dimnames(M_pre) <- list(names(theta_hat_pre), names(theta_hat_pre))

  M_ext <- rbind(
    cbind(M_pre, matrix(0, nrow = p_pre, ncol = 1)),
    cbind(matrix(0, nrow = 1, ncol = p_pre), 0)
  )
  dimnames(M_ext) <- list(c(names(theta_hat_pre), new_dummy_name),
                           c(names(theta_hat_pre), new_dummy_name))

  # Cluster score vectors for this unit
  r_pre   <- unit_pre$s_i - unit_pre$S_i %*% matrix(theta_hat_pre, ncol = 1)
  r_tilde <- rbind(r_pre, 0)

  # Updated cross-products for this unit in expanded dim
  g_col      <- matrix(g_old, ncol = 1)
  S_old_star <- unit_pre$S_i + kappa * (g_col %*% t(g_col))
  s_old_star <- unit_pre$s_i + kappa * (g_col * doty_new)
  c_star     <- kappa * g_col

  S_i_post <- rbind(cbind(S_old_star, c_star), cbind(t(c_star), kappa))
  s_i_post <- rbind(s_old_star, kappa * doty_new)
  dimnames(S_i_post) <- list(c(names(theta_hat_pre), new_dummy_name),
                              c(names(theta_hat_pre), new_dummy_name))
  rownames(s_i_post) <- c(names(theta_hat_pre), new_dummy_name)

  # --- Prop 9: outer-product replacement for this cluster
  r_post   <- s_i_post - S_i_post %*% matrix(theta_hat_post, ncol = 1)
  M_post   <- M_ext - tcrossprod(r_tilde) + tcrossprod(r_post)
  Vcr_post <- inv_dotZtZ_post %*% M_post %*% inv_dotZtZ_post

  # Replace this unit's contribution in A_N, B_N
  S_i_old_ext <- rbind(
    cbind(unit_pre$S_i, matrix(0, nrow = p_pre, ncol = 1)),
    cbind(matrix(0, nrow = 1, ncol = p_pre), 0)
  )
  s_i_old_ext <- rbind(unit_pre$s_i, 0)

  vec_old <- matrix(vec(S_i_old_ext), ncol = 1)
  vec_new <- matrix(vec(S_i_post), ncol = 1)

  A_N_post <- A_N_ext - (s_i_old_ext %*% t(vec_old)) + (s_i_post %*% t(vec_new))
  B_N_post <- B_N_ext - (vec_old %*% t(vec_old)) + (vec_new %*% t(vec_new))

  # Update unit means (expanded)
  unit_post <- unit_pre
  barZ_old_updated <- unit_pre$barZ_i + g_old / (unit_pre$T_i + 1)
  unit_post$barZ_i <- c(barZ_old_updated, 1 / (unit_pre$T_i + 1))
  names(unit_post$barZ_i) <- c(names(theta_hat_pre), new_dummy_name)
  unit_post$barY_i <- unit_pre$barY_i + doty_new / (unit_pre$T_i + 1)
  unit_post$T_i    <- unit_pre$T_i + 1L
  unit_post$S_i    <- S_i_post
  unit_post$s_i    <- s_i_post

  S_post <- S_pre
  S_post$T_support  <- S_pre$T_support + 1L
  S_post$p          <- p_post
  S_post$n          <- S_pre$n + 1L
  S_post$inv_dotZtZ <- inv_dotZtZ_post
  S_post$theta_hat  <- theta_hat_post
  S_post$sigma2_hat <- sigma2_post
  S_post$Vcr_hat    <- Vcr_post
  S_post$M_hat      <- NULL
  S_post$A_N        <- A_N_post
  S_post$B_N        <- B_N_post

  list(S_post = S_post, unit_post = unit_post)
}

# ============================================================
# One benchmark run for given N (balanced final panel, T=3)
#   - Online time: algorithm update time ONLY (no data-generation time)
#   - Offline uses the SAME data as online (built from stored arrays)
#   - Offline timing split: pdata.frame vs plm fit vs V0 vs Vcr
#   - Verification: max_abs only
# ============================================================
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

  storage_failed   <- FALSE
  storage_fail_msg <- NA_character_
  x1_store <- x2_store <- y_store <- filled <- NULL

  if (store_data) {
    ok_alloc <- tryCatch({
      x1_store <- numeric(3L * N)
      x2_store <- numeric(3L * N)
      y_store  <- numeric(3L * N)
      filled   <- logical(3L * N)
      TRUE
    }, error = function(e) {
      storage_failed   <<- TRUE
      storage_fail_msg <<- conditionMessage(e)
      FALSE
    })

    if (!isTRUE(ok_alloc)) {
      store_for_offline <- FALSE
      store_for_verif   <- FALSE
      store_data        <- FALSE
    }
  }

  idx_fun <- function(i, t) (i - 1L) * 3L + t

  # ----------------------------------------------------------
  # Warm-up data generation (not counted in algorithm update time)
  # ----------------------------------------------------------
  lambda_t   <- rnorm(3L, sd = sd_lambda)
  alpha_warm <- rnorm(N_warm, sd = sd_alpha)

  id_w   <- rep(seq_len(N_warm), each = 2L)
  time_w <- rep.int(1:2, times = N_warm)
  n_w    <- length(id_w)

  x1_w <- rnorm(n_w)
  x2_w <- rnorm(n_w)
  u_w  <- rnorm(n_w, sd = sd_u)
  y_w  <- beta_true[1]*x1_w + beta_true[2]*x2_w + alpha_warm[id_w] + lambda_t[time_w] + u_w

  df_warm <- data.frame(id = id_w, time = time_w, y = y_w, x1 = x1_w, x2 = x2_w)

  if (store_data) {
    for (k in seq_len(n_w)) {
      idx <- idx_fun(id_w[k], time_w[k])
      x1_store[idx] <- x1_w[k]
      x2_store[idx] <- x2_w[k]
      y_store[idx]  <- y_w[k]
      filled[idx]   <- TRUE
    }
  }

  # ----------------------------------------------------------
  # Initial state construction
  # ----------------------------------------------------------
  t_init_build <- system.time({
    init  <- init_state_from_warmup_df(df_warm, N_warm, alpha_warm, lambda_t)
    S_pre <- init$S_pre
    units <- init$units
  })[["elapsed"]]

  # ----------------------------------------------------------
  # Streaming algorithm updates (THIS is the target time)
  # ----------------------------------------------------------
  t_alg1 <- t_alg2 <- t_alg3 <- t_expand <- 0

  # Event 1: unit 1 gets new time t=3 (Algorithm 3)
  x1_new <- rnorm(1); x2_new <- rnorm(1); u_new <- rnorm(1, sd = sd_u)
  y_new  <- beta_true[1]*x1_new + beta_true[2]*x2_new + units[[1]]$alpha + lambda_t[3] + u_new

  if (store_data) {
    idx <- idx_fun(1L, 3L)
    x1_store[idx] <- x1_new; x2_store[idx] <- x2_new
    y_store[idx]  <- y_new;  filled[idx]   <- TRUE
  }

  t_alg3 <- system.time({
    out3      <- alg3_existing_unit_new_time(S_pre, units[[1]], time_new = 3L, x1_new, x2_new, y_new)
    S_pre     <- out3$S_post
    units[[1]] <- out3$unit_post
  })[["elapsed"]]

  # Zero-extend warm-up units 2..N_warm after new time dummy introduced
  new_dummy_name <- "factor(time)3"
  t_expand <- system.time({
    for (j in 2:N_warm) units[[j]] <- expand_unit_summary_zero(units[[j]], new_dummy_name)
  })[["elapsed"]]

  # Event 2: units 2..N_warm get t=3 (Algorithm 2)
  for (j in 2:N_warm) {
    x1_new <- rnorm(1); x2_new <- rnorm(1); u_new <- rnorm(1, sd = sd_u)
    y_new  <- beta_true[1]*x1_new + beta_true[2]*x2_new + units[[j]]$alpha + lambda_t[3] + u_new

    if (store_data) {
      idx <- idx_fun(j, 3L)
      x1_store[idx] <- x1_new; x2_store[idx] <- x2_new
      y_store[idx]  <- y_new;  filled[idx]   <- TRUE
    }

    t_alg2 <- t_alg2 + system.time({
      out2       <- alg2_existing_unit_existing_time(S_pre, units[[j]], time_new = 3L, x1_new, x2_new, y_new)
      S_pre      <- out2$S_post
      units[[j]] <- out2$unit_post
    })[["elapsed"]]
  }

  # Event 3: new units (N_warm+1)..N arrive with t=1,2,3 (Algorithm 1, batched)
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
          cat(sprintf("\r  |%3d%%|", pct))
          flush.console()
          last_pct <<- pct
        }
      }
      invisible(NULL)
    }

    while (i_start <= N) {
      i_end <- min(N, i_start + batch_size - 1L)
      m     <- i_end - i_start + 1L

      # Generate batch data (NOT timed)
      alpha_b <- rnorm(m, sd = sd_alpha)
      X1 <- matrix(rnorm(3L * m), nrow = m, ncol = 3L, byrow = TRUE)
      X2 <- matrix(rnorm(3L * m), nrow = m, ncol = 3L, byrow = TRUE)
      U  <- matrix(rnorm(3L * m, sd = sd_u), nrow = m, ncol = 3L, byrow = TRUE)
      LAM <- matrix(lambda_t, nrow = m, ncol = 3L, byrow = TRUE)
      Y  <- beta_true[1]*X1 + beta_true[2]*X2 + alpha_b + LAM + U

      if (store_data) {
        base <- (i_start - 1L) * 3L + 1L
        idxs <- base:(base + 3L * m - 1L)
        x1_store[idxs] <- as.vector(t(X1))
        x2_store[idxs] <- as.vector(t(X2))
        y_store[idxs]  <- as.vector(t(Y))
        filled[idxs]   <- TRUE
      }

      t_alg1 <- t_alg1 + system.time({
        for (k in seq_len(m))
          S_pre <- alg1_new_unit_block_fast_T3(S_pre, X1[k, ], X2[k, ], Y[k, ])
      })[["elapsed"]]

      progress_tick(i_end - N_warm)
      i_start <- i_end + 1L
    }

    if (verbose) cat("\n")
  }

  t_upd       <- t_alg1 + t_alg2 + t_alg3 + t_expand
  bytes_state <- sum_object_size(S_pre) + sum_object_size(units)

  # ----------------------------------------------------------
  # Offline reference on the same data
  # ----------------------------------------------------------
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
        y    = y_store, x1 = x1_store, x2 = x2_store
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
          plm(y ~ x1 + x2 + factor(time), data = df_p, model = "within", effect = "individual"),
          label = "plm fit", verbose_error = TRUE
        )

        if (!isTRUE(res_fit$ok)) {
          offline_ok_fit     <- FALSE
          offline_fail_stage <- "plm fit"
          offline_fail_msg   <- conditionMessage(res_fit$error)
        } else {
          m_fit       <- res_fit$value
          coef_plm    <- coef(m_fit)
          t_off_fit   <- res_fit$time
          offline_ok_fit <- TRUE

          res_V0 <- safe_system_time(vcov(m_fit), label = "vcov (V0)", verbose_error = TRUE)

          if (!isTRUE(res_V0$ok)) {
            offline_ok_V0      <- FALSE
            offline_fail_stage <- "vcov (V0)"
            offline_fail_msg   <- conditionMessage(res_V0$error)
          } else {
            V0_plm    <- res_V0$value
            t_off_V0  <- res_V0$time
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
              t_off_total    <- t_off_build_df + t_off_pdata + t_off_fit + t_off_V0 + t_off_Vcr
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

  # ----------------------------------------------------------
  # Verification: max absolute deviation vs plm on same data
  # ----------------------------------------------------------
  verif <- NULL
  if (store_for_verif) {
    if (!all(filled)) stop("Storage error: not all (id,time) cells were filled.")

    if (is.null(coef_plm) || is.null(V0_plm) || is.null(Vcr_plm)) {
      df_stream <- data.frame(
        id   = factor(rep(seq_len(N), each = 3L)),
        time = factor(rep.int(1:3, times = N), levels = 1:3),
        y    = y_store, x1 = x1_store, x2 = x2_store
      )
      df_p     <- pdata.frame(df_stream, index = c("id", "time"))
      m_ref    <- plm(y ~ x1 + x2 + factor(time), data = df_p, model = "within", effect = "individual")
      coef_plm <- coef(m_ref)
      V0_plm   <- vcov(m_ref)
      Vcr_plm  <- vcovHC(m_ref, type = "HC0", method = "arellano", cluster = "group")
      if (is.null(colnames(V0_plm)))  { colnames(V0_plm)  <- rownames(V0_plm)  <- names(coef_plm) }
      if (is.null(colnames(Vcr_plm))) { colnames(Vcr_plm) <- rownames(Vcr_plm) <- names(coef_plm) }
    }

    theta_online <- S_pre$theta_hat
    nm <- intersect(names(theta_online), names(coef_plm))
    theta_max_abs <- max(abs(theta_online[nm] - coef_plm[nm]))

    V0_online <- S_pre$sigma2_hat * S_pre$inv_dotZtZ
    colnames(V0_online) <- rownames(V0_online) <- names(theta_online)
    rn0 <- intersect(rownames(V0_online), rownames(V0_plm))
    cn0 <- intersect(colnames(V0_online), colnames(V0_plm))
    V0_max_abs <- max(abs(V0_online[rn0, cn0] - V0_plm[rn0, cn0]))

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

# ============================================================
# Grid runner (N values) -> returns a results data.frame
# ============================================================
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

    cat(sprintf("Online (update only): %.3f sec (Alg1 %.3f / Alg2 %.3f / Alg3 %.3f + expand %.3f; init-build %.3f)\n",
                res$online$t_total, res$online$t_alg1, res$online$t_alg2,
                res$online$t_alg3, res$online$t_expand, res$online$t_init_build))
    cat(sprintf("Online state size: %.2f MB\n", res$online$bytes_state / 1024^2))

    if (res$offline$attempted) {
      if (isTRUE(res$offline$ok)) {
        cat(sprintf("Offline (same data): builddf %.3f / pdata %.3f / fit %.3f / V0 %.3f / Vcr %.3f (total %.3f)\n",
                    res$offline$t_build_df, res$offline$t_pdata, res$offline$t_fit,
                    res$offline$t_V0, res$offline$t_Vcr, res$offline$t_total))
        cat(sprintf("Offline df size: %.2f MB\n", res$offline$bytes_df / 1024^2))
      } else {
        if (!is.na(res$offline$fail_stage) || !is.na(res$offline$fail_msg)) {
          cat(sprintf("Offline (same data): FAILED (%s: %s)\n",
                      res$offline$fail_stage, res$offline$fail_msg))
        } else {
          cat("Offline (same data): FAILED\n")
        }
      }
    } else {
      cat("Offline: SKIPPED (N too large for offline attempt)\n")
    }

    if (!is.null(res$verification)) {
      cat(sprintf("Verification (max abs): theta %.2e / V0 %.2e / Vcr %.2e (n=%d)\n",
                  res$verification$theta_max_abs, res$verification$V0_max_abs,
                  res$verification$Vcr_max_abs, res$verification$n_common))
    }
  }

  tab <- data.frame(
    N                   = vapply(out, function(z) z$N, integer(1)),
    online_t_upd        = vapply(out, function(z) z$online$t_update, numeric(1)),
    online_t_alg1       = vapply(out, function(z) z$online$t_alg1, numeric(1)),
    online_t_alg2       = vapply(out, function(z) z$online$t_alg2, numeric(1)),
    online_t_alg3       = vapply(out, function(z) z$online$t_alg3, numeric(1)),
    online_t_expand     = vapply(out, function(z) z$online$t_expand, numeric(1)),
    online_bytes_state  = vapply(out, function(z) z$online$bytes_state, numeric(1)),
    offline_attempted   = vapply(out, function(z) z$offline$attempted, logical(1)),
    offline_ok          = vapply(out, function(z) isTRUE(z$offline$ok), logical(1)),
    offline_t_build_df  = vapply(out, function(z) z$offline$t_build_df, numeric(1)),
    offline_t_pdata     = vapply(out, function(z) z$offline$t_pdata, numeric(1)),
    offline_t_fit       = vapply(out, function(z) z$offline$t_fit, numeric(1)),
    offline_t_V0        = vapply(out, function(z) z$offline$t_V0, numeric(1)),
    offline_t_Vcr       = vapply(out, function(z) z$offline$t_Vcr, numeric(1)),
    offline_t_total     = vapply(out, function(z) z$offline$t_total, numeric(1)),
    offline_bytes_df    = vapply(out, function(z) z$offline$bytes_df, numeric(1)),
    offline_fail_stage  = vapply(out, function(z) if (is.null(z$offline$fail_stage)) NA_character_ else z$offline$fail_stage, character(1)),
    offline_fail_msg    = vapply(out, function(z) if (is.null(z$offline$fail_msg))   NA_character_ else z$offline$fail_msg,   character(1)),
    verif_theta_max_abs = vapply(out, function(z) if (is.null(z$verification)) NA_real_ else z$verification$theta_max_abs, numeric(1)),
    verif_V0_max_abs    = vapply(out, function(z) if (is.null(z$verification)) NA_real_ else z$verification$V0_max_abs,    numeric(1)),
    verif_Vcr_max_abs   = vapply(out, function(z) if (is.null(z$verification)) NA_real_ else z$verification$Vcr_max_abs,   numeric(1))
  )

  if (!is.null(out_csv)) write.csv(tab, out_csv, row.names = FALSE)
  list(raw = out, table = tab)
}
