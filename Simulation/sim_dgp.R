# =============================================================================
# sim_dgp.R
#
# Paper: "Online Updating for Linear Panel Regressions" (Hwang & Lee, 2026)
# DGP: Subsection 5.2 Simulation
#
# Model: Y_it = X_it' β + α_i + λ_t + u_it
#   β = (1.3, -0.7)'
#   α_i  ~ N(0, σ_α²),  σ_α = 1
#   λ_t  ~ N(0, σ_λ²),  σ_λ = 1
#   u_it ~ N(0, σ_u²),  σ_u = 1
#   X_it ~ N(0, I_2),   k = 2
#   T_i  ~ Uniform{T_min, ..., T},  T = 5, T_min = 2
#   S_i  ⊆ {1,...,T}: random subset of size T_i drawn without replacement
# Technical condition: every t ∈ {1,...,T} appears at least once in the sample
# =============================================================================

#' Generate panel data — Section 5.2 DGP
#'
#' @param N       number of units
#' @param T       size of the calendar time support  (default 5)
#' @param T_min   minimum number of observed time periods per unit  (default 2)
#' @param sigma_a standard deviation of individual fixed effects  (default 1)
#' @param sigma_l standard deviation of time fixed effects  (default 1)
#' @param sigma_u standard deviation of the error term  (default 1)
#' @param seed    random seed
#' @param verbose print progress messages
#' @return list:
#'   \item{df}        data.frame (id, time, y, x1, x2)
#'   \item{t_gen}     data generation time in seconds
#'   \item{N_obs}     total number of observations
#'   \item{T_i_vec}   number of observed time periods per unit (length N)
#'   \item{beta_true} true coefficient vector c(1.3, -0.7)
#'   \item{lambda_t}  realized time fixed effects (length T)
generate_panel_52 <- function(N,
                               T       = 5L,
                               T_min   = 2L,
                               sigma_a = 1,
                               sigma_l = 1,
                               sigma_u = 1,
                               seed    = 42L,
                               verbose = TRUE) {

  if (verbose)
    cat(sprintf("  Generating DGP: N = %s,  T = %d ...",
                format(N, big.mark = ","), T))
  t0 <- proc.time()

  set.seed(seed)

  beta     <- c(1.3, -0.7)                              # true coefficients
  alpha_i  <- rnorm(N,   sd = sigma_a)                  # α_i ~ N(0, σ_α²)
  lambda_t <- rnorm(T,   sd = sigma_l)                  # λ_t ~ N(0, σ_λ²)
  T_i_vec  <- sample(T_min:T, N, replace = TRUE)        # T_i ~ Uniform{T_min,...,T}
  total_obs <- sum(T_i_vec)

  # ------------------------------------------------------------------
  # Vectorize id_vec  O(n)
  # ------------------------------------------------------------------
  id_vec <- rep(seq_len(N), T_i_vec)

  # ------------------------------------------------------------------
  # time_vec: for each unit i, draw S_i ⊆ {1,...,T} of size T_i without replacement
  # Requires N calls to sample.int(T, Ti) → loop unavoidable; pre-allocate for speed
  # ------------------------------------------------------------------
  time_vec <- integer(total_obs)
  pos <- 1L
  for (i in seq_len(N)) {
    Ti  <- T_i_vec[i]
    end <- pos + Ti - 1L
    time_vec[pos:end] <- sample.int(T, Ti)   # without replacement
    pos <- end + 1L
  }

  # ------------------------------------------------------------------
  # Technical condition: ensure every t ∈ {1,...,T} appears at least once
  # (satisfied automatically for large N, but checked and fixed explicitly)
  # ------------------------------------------------------------------
  missing_t <- setdiff(seq_len(T), unique(time_vec))
  if (length(missing_t) > 0L) {
    warning(sprintf(
      "Forcing coverage of missing calendar time(s) (%s) — occurs only for very small N",
      paste(missing_t, collapse = ",")))
    cum_end <- cumsum(T_i_vec)
    cum_beg <- c(1L, cum_end[-N] + 1L)
    for (t_miss in missing_t) {
      # Replace the first observation of the first unit with T_i >= 2
      for (ii in seq_len(N)) {
        if (T_i_vec[ii] >= 2L) {
          time_vec[cum_beg[ii]] <- t_miss
          break
        }
      }
    }
  }

  # ------------------------------------------------------------------
  # Vectorize X and u  (fast)
  # ------------------------------------------------------------------
  x1_vec <- rnorm(total_obs)
  x2_vec <- rnorm(total_obs)
  u_vec  <- rnorm(total_obs, sd = sigma_u)

  # ------------------------------------------------------------------
  # Vectorize y: Y_it = α_i + λ_t + 1.3 x1 - 0.7 x2 + u
  # ------------------------------------------------------------------
  y_vec <- rep(alpha_i, T_i_vec) +
           lambda_t[time_vec]    +
           beta[1L] * x1_vec    +
           beta[2L] * x2_vec    +
           u_vec

  t_gen <- (proc.time() - t0)[["elapsed"]]

  if (verbose)
    cat(sprintf("  Done %.1f sec  (%s obs,  %.0f MB)\n",
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
