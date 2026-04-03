# =============================================================================
# Section 1: Utilities
# =============================================================================

#' Vectorize a matrix in column-major order
vec <- function(M) as.vector(M)

#' Index embedding for vec(): when a p_old × p_old matrix is expanded to
#' (p_old+1) × (p_old+1), returns the position in vec(M_new) where each
#' element of vec(M_old) lands.
#'
#' vec() uses column-major order, so element at row r, column c maps to
#' position r + (c-1) * p_new in the expanded matrix.
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

#' Total memory size of objects in a list (bytes)
sum_object_size <- function(x) {
  sum(vapply(x, function(z) as.numeric(object.size(z)), numeric(1)))
}

#' Evaluate an expression and return the result together with elapsed time
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

# Progress bar (uses numeric to avoid integer overflow for large N)
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
# Section 2: Design matrix & within transformation
# =============================================================================

#' Build the TWFE design matrix Z_it = [X_it, time dummies]
#'
#' @param x_mat       n × k covariate matrix (must already have column names)
#' @param time_vec    integer vector of calendar times, length n
#' @param T_support   size of the current calendar time support
#' @param baseline_time time dummy to drop (default = 1)
#' @return n × (k + T_support - 1) matrix
build_Z_raw <- function(x_mat, time_vec, T_support, baseline_time = 1L) {
  # Generate time dummies for all calendar times except baseline_time
  dummy_times  <- setdiff(seq_len(T_support), baseline_time)
  dummy_names  <- paste0("factor(time)", dummy_times)

  # Build the entire dummy matrix at once using outer()
  D_mat        <- outer(time_vec, dummy_times, `==`) + 0L
  colnames(D_mat) <- dummy_names

  cbind(x_mat, D_mat)
}

#' Within transformation for individual i
#'
#' \dot{Z}_{it} = Z_{it} - \bar{Z}_i,  \dot{Y}_{it} = Y_{it} - \bar{Y}_i
#'
#' @param Z_raw_i T_i × p matrix
#' @param y_i     vector of length T_i
#' @return list(dotZ_i, dotY_i, barZ_i, barY_i)
within_transform_unit <- function(Z_raw_i, y_i) {
  barZ_i <- colMeans(Z_raw_i)
  barY_i <- mean(y_i)
  dotZ_i <- sweep(Z_raw_i, 2, barZ_i, "-")
  dotY_i <- y_i - barY_i
  list(dotZ_i = dotZ_i, dotY_i = dotY_i, barZ_i = barZ_i, barY_i = barY_i)
}


# =============================================================================
# Section 3: State initialization and management
# =============================================================================

#' Create a per-unit summary object (internal)
#'
#' @param id      unit identifier
#' @param T_i     number of observed time periods
#' @param barZ_i  within mean vector (length p)
#' @param barY_i  within mean scalar
#' @param S_i     \dot{Z}_i' \dot{Z}_i  (p × p)
#' @param s_i     \dot{Z}_i' \dot{Y}_i  (p × 1 matrix)
make_unit_state <- function(id, T_i, barZ_i, barY_i, S_i, s_i) {
  list(id = id, T_i = T_i, barZ_i = barZ_i, barY_i = barY_i, S_i = S_i, s_i = s_i)
}

#' Initialize the state object from a warm-up data frame
#'
#' @param df_warm     data.frame used for warm-up
#' @param id_col      column name for unit ID
#' @param time_col    column name for calendar time
#' @param y_col       column name for the dependent variable
#' @param x_cols      character vector of covariate column names
#' @param baseline_time time dummy to drop (NULL = minimum time in data)
#' @return list(S_pre = state object,  units = named list of unit summaries)
state_init <- function(df_warm, id_col, time_col, y_col, x_cols,
                       baseline_time = NULL) {

  id_vec   <- df_warm[[id_col]]
  time_vec <- as.integer(df_warm[[time_col]])
  y_vec    <- df_warm[[y_col]]
  x_mat    <- as.matrix(df_warm[, x_cols, drop = FALSE])

  # Determine baseline_time automatically if not specified
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

  # Split indices by unit at once O(n) — avoids O(n²) repeated which() calls
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

  # Aggregate matrices required for cluster-robust updating in Alg 1/2
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

  # M = \sum_i r_i r_i',  r_i = \dot{Z}_i' \hat{e}_i  (sum of cluster score outer products)
  # Reuse idx_list — avoids O(n²) loop
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
    # Metadata
    x_cols       = x_cols,
    baseline_time = baseline_time,
    # Dimensions
    T_support    = T_support,
    p            = p,
    N            = N_pre,
    n            = n,
    # Core estimates
    inv_dotZtZ   = inv_dotZtZ,
    theta_hat    = theta_hat,
    sigma2_hat   = sigma2_hat,
    # Variance estimates
    Vcr_hat      = Vcr_hat,
    M_hat        = M_hat,          # used in Alg 3; set to NULL after Alg 1/2
    # Aggregate matrices for cluster-robust updating (Prop 3/6)
    A_N          = A_N,
    B_N          = B_N,
    # Vcr formula accumulator for chunk mode: Σ s_i s_i' (streaming units)
    M_ss         = matrix(0, nrow = p, ncol = p)
  )

  list(S_pre = S_pre, units = units)
}

#' Zero-extend a unit summary object after a new time dummy is introduced
#'
#' After Algorithm 3, expand S_i, s_i, barZ_i of remaining units to
#' (p+1) dimensions (the new dummy is always 0 in the pre-update sample)
#'
#' @param unit           unit summary object
#' @param new_dummy_name name of the newly added time dummy
#' @return unit summary object with expanded dimensions
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
# Section 4: Algorithms
# =============================================================================

# -----------------------------------------------------------------------------
# Algorithm 1 (Props 1–3): arrival of a new unit
#
# Input  S_pre     : pre-update state
#        x_mat_i   : T_i × k covariate matrix (new unit)
#        time_vec_i: calendar time vector of length T_i
#        y_i       : dependent variable vector of length T_i
# Output S_post    : post-update state
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
  # Solve H_i jointly for dotZ_i and e_tilde — avoids explicit inversion of H_i
  H_i       <- diag(T_i) + dotZ_i %*% inv_pre %*% t(dotZ_i)
  e_tilde_i <- drop(dotY_i - dotZ_i %*% theta_pre)
  sol       <- solve(H_i, cbind(dotZ_i, e_tilde_i))
  W_i       <- sol[, seq_len(p), drop = FALSE]   # H_i^{-1} \dot{Z}_i
  v_i       <- sol[, p + 1L]                     # H_i^{-1} \tilde{e}_i

  # --- Prop 1: coefficient update
  # \hat{\beta}^* = \hat{\beta} + (\dot{X}'\dot{X})^{-1} \dot{X}_{N+1}' H^{-1} \tilde{e}_{N+1}
  q_i            <- drop(crossprod(dotZ_i, v_i))   # \dot{Z}_i' H_i^{-1} \tilde{e}_i
  d_i            <- drop(inv_pre %*% q_i)
  theta_post     <- drop(theta_pre + d_i)
  names(theta_post) <- names(theta_pre)

  # --- Prop 2: (\dot{Z}'\dot{Z})^{-1} update (Woodbury)
  G_i         <- crossprod(dotZ_i, W_i)            # \dot{Z}_i' H_i^{-1} \dot{Z}_i
  inv_post    <- inv_pre - inv_pre %*% G_i %*% inv_pre

  # sigma^2 update
  # \hat{\sigma}^{2*} = [(n-N-k)/(n^*-(N+1)-k)] \hat{\sigma}^2
  #                   + \tilde{e}_{N+1}' H^{-1}_{N+1} \tilde{e}_{N+1} / [n^*-(N+1)-k]
  n_post   <- S_pre$n + T_i
  N_post   <- S_pre$N + 1L
  df_pre   <- S_pre$n - S_pre$N - p
  df_post  <- n_post - N_post - p
  sig2_post <- (df_pre / df_post) * sig2_pre +
    as.numeric(crossprod(e_tilde_i, v_i)) / df_post

  # --- Prop 3: cluster-robust update
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

  # A_N, B_N update (add new unit's contribution)
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
  S_post$M_hat      <- NULL  # M_hat not tracked after Alg 1/2 (reconstructed when entering Alg 3)
  S_post$A_N        <- A_pre + s_i %*% t(vecS_i)
  S_post$B_N        <- B_pre + vecS_i %*% t(vecS_i)

  S_post
}


# -----------------------------------------------------------------------------
# Algorithm 2 (Props 4–6): new observation for an existing unit at an existing calendar time
#
# Input  S_pre    : pre-update state
#        unit_pre : pre-update summary object for the unit
#        x_new    : covariate vector of length k (new observation)
#        time_new : calendar time of the new observation (within current support)
#        y_new    : dependent variable value of the new observation
# Output list(S_post, unit_post)
# -----------------------------------------------------------------------------
alg2_existing_unit <- function(S_pre, unit_pre, x_new, time_new, y_new) {

  stopifnot(as.integer(time_new) %in% seq_len(S_pre$T_support))

  # Map new observation to p-dimensional Z using existing T_support
  x_row   <- matrix(x_new, nrow = 1, dimnames = list(NULL, S_pre$x_cols))
  Z_new   <- build_Z_raw(x_row, as.integer(time_new),
                         S_pre$T_support, S_pre$baseline_time)

  # Within deviations: \dot{z}_{1,T_1^*} = Z_{new} - \bar{Z}_1,  \dot{y}_{1,T_1^*} = y - \bar{Y}_1
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

  # --- Prop 4: coefficient update
  # \hat{\beta}^* = \hat{\beta} + s_{1,T_1^*} (\dot{Z}'\dot{Z})^{-1} \dot{z}_{1,T_1^*} \tilde{e}
  theta_post <- drop(theta_pre + s_scalar * inv_pre %*% dotz_col * e_tilde)
  names(theta_post) <- names(theta_pre)

  # --- Prop 5: (\dot{Z}'\dot{Z})^{-1} and sigma^2 update (Sherman-Morrison)
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

  # --- Prop 6: cluster-robust update
  I_p       <- diag(p)
  U_rank1   <- I_p - s_scalar * (inv_pre %*% dotz_col %*% t(dotz_col))
  # Individual cluster score: \dot{X}_i' \hat{e}_i = s_i - S_i \hat{\beta}
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
  # \Omega^{(4)} term
  Omega4 <- s_scalar^2 * e_tilde^2 * (
    dotz_col %*% t(dotz_col) -
    dotz_col %*% t(dotz_col) %*% inv_pre %*% unit_pre$S_i -
    unit_pre$S_i %*% inv_pre %*% dotz_col %*% t(dotz_col)
  )
  Vcr_post <- U_rank1 %*%
    (Vcr_pre + inv_pre %*% (-Omega1 - t(Omega1) + Omega2 + Omega3 + Omega4) %*% inv_pre) %*%
    t(U_rank1)

  # Update the unit summary
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
# Algorithm 3 (Props 7–9): new calendar time T+1 for an existing unit
#
# Parameter dimension increases p → p+1 (new time dummy added)
#
# Input  S_pre    : pre-update state  (T_support = T)
#        unit_pre : pre-update unit summary object
#        x_new    : covariate vector of length k
#        time_new : new calendar time (must equal T_support + 1)
#        y_new    : dependent variable value of the new observation
# Output list(S_post, unit_post)
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

  # Compute Z_new in the existing p-dimensional space
  # (new time dummy = 1, but not in current support, so treated as 0 in old Z)
  x_row    <- matrix(x_new, nrow = 1, dimnames = list(NULL, S_pre$x_cols))
  Z_new_p  <- build_Z_raw(x_row, as.integer(time_new),
                          S_pre$T_support, S_pre$baseline_time)
  # Z_new_p is in the existing p dimensions (no new dummy column yet)

  # Within deviations (relative to existing p-dimensional barZ_i)
  kappa    <- unit_pre$T_i / (unit_pre$T_i + 1L)
  doty_new <- as.numeric(y_new - unit_pre$barY_i)
  g_old    <- drop(Z_new_p - unit_pre$barZ_i)   # length p_pre

  # --- Prop 7: expanded inverse and coefficient update
  # g = [g_old, 1] ∈ R^{p+1}  (within deviation of new time dummy = 1 - 0 = 1)
  # \dot{Z}^{*'} \dot{Z}^* = diag(\dot{Z}'\dot{Z}, 0) + kappa * g g'  (4.3)
  # Inverse (4.4):
  #   [ (Z'Z)^{-1}         -q/g2 ]
  #   [ -q'/g2    (1+kappa*h)/(kappa*g2^2) ]
  # where g2 = 1 (new time dummy value), q = (Z'Z)^{-1} g_old, h = g_old' q
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

  # --- Prop 8: sigma^2 unchanged
  # n^* - N - p^* = (n+1) - N - (p+1) = n - N - p  → df and RSS unchanged
  # \hat{sigma}^{2*} = \hat{sigma}^2
  sig2_post <- S_pre$sigma2_hat

  # Zero-extend A_N, B_N to (p+1) dimensions, then replace the unit's contribution
  idx_embed <- vec_embed_index(p_pre)
  A_ext     <- matrix(0, nrow = p_post, ncol = p_post^2)
  B_ext     <- matrix(0, nrow = p_post^2, ncol = p_post^2)
  A_ext[seq_len(p_pre), idx_embed] <- A_pre
  B_ext[idx_embed, idx_embed]      <- B_pre

  # Reconstruct M = \sum_i r_i r_i' (M_hat is set NULL after Alg 1/2, so recompute)
  if (is.null(S_pre$M_hat)) {
    dotZtZ_pre <- solve(inv_pre)
    M_pre      <- dotZtZ_pre %*% S_pre$Vcr_hat %*% dotZtZ_pre
  } else {
    M_pre <- S_pre$M_hat
  }
  dimnames(M_pre) <- list(names(theta_pre), names(theta_pre))

  # Zero-extend M to (p+1) dimensions: \tilde{M}
  M_ext <- rbind(cbind(M_pre, 0), 0)
  dimnames(M_ext) <- list(nms_post, nms_post)

  # --- Prop 9: update cluster score outer product sum M^*
  # r_1 = s_1 - S_1 \hat{\theta},  \tilde{r}_1 = [r_1, 0]
  r_pre   <- unit_pre$s_i - unit_pre$S_i %*% matrix(theta_pre, ncol = 1)
  r_tilde <- rbind(r_pre, 0)

  # Post-update \dot{Z}^{*'}\dot{Z}^* and \dot{Z}^{*'}\dot{Y}^* for unit 1
  g_col      <- matrix(g_old, ncol = 1)
  S1_kappa   <- unit_pre$S_i + kappa * (g_col %*% t(g_col))   # 기존 p 차원
  s1_kappa   <- unit_pre$s_i + kappa * (g_col * doty_new)

  # Expand to (p+1) dimensions: add outer products with new time dummy
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

  e_tilde_new <- doty_new - t(g_old) %*% theta_pre   # \dot{y} - g' θ (value before computing v_new)
  g_full      <- c(g_old, 1)                          # g ∈ R^{p+1}
  r_post      <- r_tilde +
    kappa * g_full * as.numeric(e_tilde_new) -
    (S_i_tilde + kappa * outer(g_full, g_full)) %*%
    matrix(theta_post - theta_tilde, ncol = 1)

  # M^* = \tilde{M} - \tilde{r}_1 \tilde{r}_1' + r_1^* r_1^{*'}
  M_post  <- M_ext - tcrossprod(r_tilde) + tcrossprod(r_post)
  Vcr_post <- inv_post %*% M_post %*% inv_post

  # Remove unit 1's old contribution from A_N, B_N and add the new one
  S_i_old_ext <- rbind(cbind(unit_pre$S_i, 0), 0)
  s_i_old_ext <- rbind(unit_pre$s_i, 0)
  vec_old     <- matrix(vec(S_i_old_ext), ncol = 1)
  vec_new     <- matrix(vec(S_i_post),    ncol = 1)

  A_post <- A_ext - (s_i_old_ext %*% t(vec_old)) + (s_i_post %*% t(vec_new))
  B_post <- B_ext - (vec_old %*% t(vec_old))      + (vec_new %*% t(vec_new))

  # Update unit summary (to p+1 dimensions)
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
  S_post$M_hat      <- NULL   # reconstructed when entering Alg 3 next time
  S_post$A_N        <- A_post
  S_post$B_N        <- B_post

  list(S_post = S_post, unit_post = unit_post)
}


