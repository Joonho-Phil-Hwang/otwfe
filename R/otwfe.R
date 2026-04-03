# =============================================================================
# Section 5: Main API
# =============================================================================

#' Online TWFE estimation
#'
#' Processes a data.frame sequentially (unit × time order) to estimate TWFE
#' coefficients and variance. Supports streaming unit by unit without loading
#' the entire dataset into memory.
#'
#' Algorithm selection (automatic):
#'   - New unit arrival          → Algorithm 1  (Props 1–3)
#'   - Existing unit, existing time → Algorithm 2  (Props 4–6)  [warm-up units only]
#'   - Existing unit, new time   → Algorithm 3  (Props 7–9)  [warm-up units only]
#'
#' @param data            data.frame containing all observations
#' @param id_col          column name for unit ID  (character)
#' @param time_col        column name for calendar time  (character, coerced to integer)
#' @param y_col           column name for the dependent variable  (character)
#' @param x_cols          character vector of covariate column names
#' @param baseline_time   time dummy to drop; NULL uses min(time) automatically
#' @param warmup_ids      unit ID vector to use as warm-up; NULL for automatic selection
#' @param warmup_n        minimum number of warm-up units (used when warmup_ids is NULL);
#'                        NULL defaults to max(3*p, 30)
#' @param track_all_units TRUE preserves per-unit state for streaming units;
#'                        FALSE (default) keeps only warm-up unit states → memory O(n_warm * p^2)
#' @param verbose         print progress messages
#' @return object of class \code{"otwfe"}
#' @export
# =============================================================================
# Vcr formula helper
# M = Σ r_i r_i',  r_i = s_i - S_i θ
# Formula: M = M_ss - A_N(θ⊗I) - [A_N(θ⊗I)]' + (θ'⊗I) B_N (θ⊗I)
# where M_ss = Σ s_i s_i',  A_N = Σ s_i vec(S_i)',  B_N = Σ vec(S_i) vec(S_i)'
# =============================================================================
.compute_vcr <- function(S, units) {
  p  <- S$p
  th <- S$theta_hat

  # M_ss_warm: sum of s_i s_i' for warm-up units
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
# Split df into chunks of chunk_size units at unit boundaries
# Assumes data is sorted by the id column
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
# otwfe_init(): create an empty handle
# =============================================================================
#' Initialize an empty Online TWFE handle
#'
#' @param x_cols        character vector of covariate column names
#' @param time_col      column name for calendar time
#' @param id_col        column name for unit ID
#' @param y_col         column name for the dependent variable
#' @param T_support     size of the calendar time support (NULL = auto-detected from first chunk)
#' @param baseline_time baseline time dummy to drop (default 1)
#' @param verbose       print progress messages
#' @return handle list
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
    time_map      = NULL,   # original time values → re-indexed mapping
    time_levels   = NULL    # original time values preserved (for output)
  )
}

# =============================================================================
# otwfe_update(): incorporate chunk_df into the handle
# =============================================================================
#' Online TWFE chunk update
#'
#' @param handle    handle returned by otwfe_init() or a previous otwfe_update()
#' @param chunk_df  data frame chunk to process (split at unit boundaries)
#' @return updated handle
#' @export
otwfe_update <- function(handle, chunk_df) {
  if (!handle$warmup_done) {
    handle <- .otwfe_first_chunk(handle, chunk_df)
  } else {
    handle <- .otwfe_stream_chunk(handle, chunk_df)
  }
  handle
}

# First chunk: warm-up selection + state initialization + Rcpp batch
.otwfe_first_chunk <- function(handle, chunk_df) {
  x_cols        <- handle$x_cols
  time_col      <- handle$time_col
  id_col        <- handle$id_col
  y_col         <- handle$y_col
  verbose       <- handle$verbose

  # 입력 검증
  for (col in c(id_col, time_col, y_col, x_cols)) {
    if (!col %in% names(chunk_df))
      stop(sprintf("Column not found: '%s'", col))
  }

  chunk_df[[time_col]] <- as.integer(chunk_df[[time_col]])
  all_times_chunk      <- sort(unique(chunk_df[[time_col]]))

  # Determine T_support (if not specified, use max time in first chunk)
  if (is.null(handle$T_support)) {
    T_support <- max(all_times_chunk)
    if (verbose) message(sprintf("T_support auto-detected from first chunk: %d", T_support))
  } else {
    T_support <- handle$T_support
  }

  # Time re-indexing (if time_map not yet set, determine it from first chunk)
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

  # Eligible units (at least 2 observations)
  unit_counts <- tapply(chunk_df[[id_col]], chunk_df[[id_col]], length)
  eligible    <- names(unit_counts)[unit_counts >= 2L]
  if (length(eligible) == 0L)
    stop("No eligible units (>= 2 observations) in first chunk.")

  # Greedy warm-up selection (ensures all calendar times are covered)
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
    cat(sprintf("Initializing state (warm-up: %d units)...\n", length(warmup_ids)))

  df_warm         <- chunk_df[id_ch_vec %in% warmup_ids, ]
  df_stream_chunk <- chunk_df[!(id_ch_vec %in% warmup_ids), ]

  # State initialization (warm-up batch OLS)
  init      <- state_init(df_warm, id_col, time_col, y_col, x_cols, baseline_time)
  S_pre     <- init$S_pre
  units     <- init$units

  handle$S             <- S_pre
  handle$units         <- units
  handle$T_support     <- T_support
  handle$baseline_time <- baseline_time
  handle$warmup_done   <- TRUE

  # Process remaining units in the first chunk via Rcpp batch
  if (nrow(df_stream_chunk) > 0L)
    handle <- .otwfe_stream_chunk(handle, df_stream_chunk)

  handle
}

# Subsequent chunks: all new units → Rcpp batch
.otwfe_stream_chunk <- function(handle, chunk_df) {
  S_pre    <- handle$S
  x_cols   <- handle$x_cols
  time_col <- handle$time_col
  id_col   <- handle$id_col
  y_col    <- handle$y_col

  # Time re-indexing (use map stored in handle)
  if (!is.null(handle$time_map))
    chunk_df[[time_col]] <- handle$time_map[as.character(chunk_df[[time_col]])]

  if (!.alg1_batch_rcpp_available)
    stop("Rcpp batch module required. Please check sourceCpp('src/alg1_batch.cpp').")

  stream_idx_list <- split(seq_len(nrow(chunk_df)),
                            as.character(chunk_df[[id_col]]))
  new_units <- names(stream_idx_list)
  if (length(new_units) == 0L) return(handle)

  new_all_idx  <- unlist(stream_idx_list, use.names = FALSE)
  t_all_new    <- as.integer(chunk_df[[time_col]])[new_all_idx]
  over_support <- t_all_new > S_pre$T_support

  if (any(over_support)) {
    warning(sprintf("%d observation(s) exceed T_support (%d) and will be dropped.",
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
# otwfe_finalize(): extract final results
# =============================================================================
#' Finalize an Online TWFE handle and return the result object
#'
#' @param handle    handle after all otwfe_update() calls are complete
#' @return object of class \code{"otwfe"}
#' @export
otwfe_finalize <- function(handle) {
  if (is.null(handle$S))
    stop("Call otwfe_update() before otwfe_finalize().")

  S       <- handle$S
  nm_list <- dimnames(S$inv_dotZtZ)

  # Compute Vcr (formula: using M_ss, A_N, B_N, and theta_final)
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

  t_start <- proc.time()   # start total elapsed time

  # --------------------------------------------------
  # chunk_size path: delegate to otwfe_init/update/finalize
  # --------------------------------------------------
  if (!is.null(chunk_size)) {
    chunk_size <- as.integer(chunk_size)
    if (chunk_size < 1L) stop("chunk_size must be a positive integer.")

    # Re-index time values across the full dataset (ensures consistency across chunks)
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
    # Time already re-indexed, so set identity map
    handle$time_map    <- setNames(seq_len(T_support_cs),
                                   as.character(seq_len(T_support_cs)))
    handle$time_levels <- time_levels_cs

    # Compute chunk boundaries as row ranges (avoids pre-copying chunks → memory savings)
    rle_cs   <- rle(data[[id_col]])
    u_lens   <- rle_cs$lengths
    n_units  <- length(u_lens)
    cum_rows <- cumsum(u_lens)
    n_chunks <- ceiling(n_units / chunk_size)

    if (verbose)
      cat(sprintf("Chunk processing: %d chunk(s) (chunk_size = %s units)\n",
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
      cat(sprintf("\nDone: N = %s, n = %s, p = %d, T_support = %d | elapsed: %.2f sec\n",
                  format(S$N, big.mark = ","), format(S$n, big.mark = ","),
                  S$p, S$T_support, t_elapsed))
    return(result)
  }

  # --------------------------------------------------
  # 1. Input validation and preprocessing (single-call path)
  # --------------------------------------------------
  stopifnot(is.data.frame(data))
  for (col in c(id_col, time_col, y_col, x_cols)) {
    if (!col %in% names(data))
      stop(sprintf("Column not found: '%s'", col))
  }

  data[[time_col]] <- as.integer(data[[time_col]])

  all_times <- sort(unique(data[[time_col]]))

  if (length(all_times) < 2L)
    stop("At least 2 distinct calendar time periods are required.")

  # Re-index calendar times to 1, 2, ..., T
  # seq_len(T_support) = 1..T must match the number of distinct observed times
  # so that the number of time dummies is determined correctly (e.g. 1935–1954 → 1–20)
  time_levels <- all_times                         # original time values preserved (for output)
  time_map    <- setNames(seq_along(time_levels),
                          as.character(time_levels))
  data[[time_col]] <- time_map[as.character(data[[time_col]])]
  all_times <- seq_along(time_levels)              # 1, 2, ..., T

  # Re-index baseline_time as well
  if (is.null(baseline_time)) {
    baseline_time <- 1L
  } else {
    baseline_time <- as.integer(time_map[as.character(as.integer(baseline_time))])
    if (is.na(baseline_time))
      stop("baseline_time value not found in data.")
  }

  # Minimum condition for within transformation: each unit needs >= 2 observations
  unit_counts <- tapply(data[[id_col]], data[[id_col]], length)
  eligible    <- names(unit_counts)[unit_counts >= 2L]

  # --------------------------------------------------
  # 2. Warm-up unit selection
  #    Key condition: warm-up data must cover all calendar times so that
  #    T_support = T_max is guaranteed after state_init().
  #    (Time dummy expansion corresponding to Algorithm 3 is handled inside the batch OLS)
  # --------------------------------------------------
  if (!is.null(warmup_ids)) {
    # User-specified warm-up IDs
    warmup_ids <- as.character(warmup_ids)
    warmup_ids <- intersect(warmup_ids, eligible)
    if (length(warmup_ids) < 1L)
      stop("None of the specified warmup_ids are eligible (require >= 2 observations).")

  } else {
    # Automatic selection: prioritize units that cover all calendar times, then meet minimum count
    p_guess  <- length(x_cols) + max(0L, length(all_times) - 1L)
    min_warm <- min(max(3L * p_guess, 30L), length(eligible))
    if (!is.null(warmup_n)) min_warm <- min(warmup_n, length(eligible))

    # Time sets observed per unit (eligible units only)
    id_ch_vec  <- as.character(data[[id_col]])
    unit_times <- tapply(as.integer(data[[time_col]]), id_ch_vec,
                         function(x) unique(x), simplify = FALSE)
    elig_times <- unit_times[eligible]

    # Greedy: select units in descending coverage order until all times are covered
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

    # Fill up to min_warm if needed
    if (length(warmup_ids) < min_warm) {
      extra      <- setdiff(elig_ord, warmup_ids)
      warmup_ids <- c(warmup_ids,
                      head(extra, min_warm - length(warmup_ids)))
    }
  }

  if (length(warmup_ids) < 1L)
    stop("No eligible warm-up units found. Each unit requires at least 2 observations.")

  df_warm   <- data[id_ch_vec %in% warmup_ids, ]
  df_stream <- data[!(id_ch_vec %in% warmup_ids), ]

  # --------------------------------------------------
  # 3. Initial state construction (batch OLS — algebraically equivalent to Alg 1/2/3)
  # --------------------------------------------------
  if (verbose) {
    cat(sprintf("Initializing state (warm-up: %d units)...\n", length(warmup_ids)))
  }

  init  <- state_init(df_warm, id_col, time_col, y_col, x_cols, baseline_time)
  S_pre <- init$S_pre
  units <- init$units   # only warm-up unit states are preserved

  # --------------------------------------------------
  # 4. Streaming: process remaining units via Algorithm 1
  #    track_all_units = FALSE (default): per-unit state not stored for streaming units
  #    → units list retains only warm-up unit states (memory savings)
  # --------------------------------------------------
  if (nrow(df_stream) == 0L) {
    t_elapsed <- (proc.time() - t_start)[["elapsed"]]
    if (verbose)
      cat(sprintf("No streaming data. Returning warm-up state. | elapsed: %.2f sec\n",
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

  # Sort by unit × time and pre-split indices by unit O(n)
  ord       <- order(df_stream[[id_col]], df_stream[[time_col]])
  df_stream <- df_stream[ord, ]

  stream_idx_list <- split(seq_len(nrow(df_stream)),
                           as.character(df_stream[[id_col]]))
  stream_unit_chs <- names(stream_idx_list)
  n_stream        <- length(stream_unit_chs)

  # Split streaming units into two groups:
  #   warm_unit_chs : warm-up units (Algorithm 2 or 3)
  #   new_unit_chs  : new units     (Algorithm 1 — Rcpp batch or R fallback)
  warm_in_stream <- stream_unit_chs[stream_unit_chs %in% names(units)]
  new_in_stream  <- stream_unit_chs[!stream_unit_chs %in% names(units)]

  if (verbose) {
    cat(sprintf("Streaming: %s units (%d warm-up / %s new), %s observations\n",
                format(n_stream, big.mark = ","),
                length(warm_in_stream),
                format(length(new_in_stream), big.mark = ","),
                format(nrow(df_stream), big.mark = ",")))
    pb <- make_progress_bar(0, n_stream)
  }

  # ------------------------------------------------------------------
  # STEP A: streaming warm-up units (Algorithm 2 / 3, R loop)
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
        # Case B: existing calendar time → Algorithm 2
        res             <- alg2_existing_unit(S_pre, units[[uid_ch]], x_obs, t_obs, y_obs)
        S_pre           <- res$S_post
        units[[uid_ch]] <- res$unit_post

      } else {
        # Case C: new calendar time T+1 → Algorithm 3
        stopifnot(t_obs == S_pre$T_support + 1L)

        res             <- alg3_new_caltime(S_pre, units[[uid_ch]], x_obs, t_obs, y_obs)
        S_pre           <- res$S_post
        units[[uid_ch]] <- res$unit_post

        # Zero-extend all other existing unit summaries
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
  # STEP B: new unit processing (Algorithm 1)
  #   track_all_units = FALSE & Rcpp available → Rcpp batch (fast)
  #   otherwise                                → pure-R loop (fallback)
  # ------------------------------------------------------------------
  use_rcpp_batch <- (!track_all_units) && .alg1_batch_rcpp_available &&
                    (length(new_in_stream) > 0L)

  if (use_rcpp_batch) {
    # ---- Rcpp batch path ----
    # Gather all new unit indices at once — vectorized without per-unit loop
    new_all_idx   <- unlist(stream_idx_list[new_in_stream], use.names = FALSE)
    t_all_new     <- as.integer(df_stream[[time_col]])[new_all_idx]
    over_support  <- t_all_new > S_pre$T_support

    if (any(over_support)) {
      warning(sprintf(
        "%d observation(s) exceed T_support (%d) and will be dropped.",
        sum(over_support), S_pre$T_support))
    }

    # Keep only observations within T_support — extract needed rows to save memory
    keep_idx   <- new_all_idx[!over_support]
    t_batch    <- t_all_new[!over_support]
    x_batch    <- as.matrix(df_stream[keep_idx, x_cols, drop = FALSE])
    y_batch    <- as.numeric(df_stream[[y_col]])[keep_idx]

    # unit_lens: actual observation count per unit after T_support filtering
    # (in the common no-filter case, lengths() is used directly for speed)
    if (any(over_support)) {
      # With filtering: compute keep flags per unit
      lens_raw  <- lengths(stream_idx_list[new_in_stream])
      # Observations are stacked contiguously per unit; use cumsum to find unit boundaries
      cum_lens  <- c(0L, cumsum(lens_raw))
      lens_vec  <- integer(length(new_in_stream))
      for (j in seq_along(new_in_stream)) {
        obs_j     <- (cum_lens[j] + 1L) : cum_lens[j + 1L]
        lens_vec[j] <- sum(!over_support[obs_j])
      }
      # Drop units with fewer than 2 observations
      valid_j   <- lens_vec >= 2L
      if (!all(valid_j)) {
        # Remove the corresponding rows from keep_idx
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

      # Restore names (Rcpp returns matrices without dimnames)
      theta_new <- rcpp_res$theta_hat
      names(theta_new) <- names(S_pre$theta_hat)
      inv_new   <- rcpp_res$inv_dotZtZ
      dimnames(inv_new) <- dimnames(S_pre$inv_dotZtZ)

      # State update: accumulate M_ss/A_N/B_N (Vcr computed via formula in otwfe_finalize)
      S_pre$inv_dotZtZ <- inv_new
      S_pre$theta_hat  <- theta_new
      S_pre$sigma2_hat <- rcpp_res$sigma2_hat
      S_pre$M_hat      <- NULL
      S_pre$N          <- S_pre$N + rcpp_res$N_add
      S_pre$n          <- S_pre$n + rcpp_res$n_add
      S_pre$M_ss       <- S_pre$M_ss + rcpp_res$M_ss_add
      S_pre$A_N        <- S_pre$A_N + rcpp_res$A_N_add
      S_pre$B_N        <- S_pre$B_N + rcpp_res$B_N_add

      # Compute Vcr (in single-call mode theta_new = theta_final, so compute immediately)
      nm_list  <- dimnames(S_pre$inv_dotZtZ)
      Vcr_new  <- .compute_vcr(S_pre, units)
      dimnames(Vcr_new) <- list(nm_list[[1]], nm_list[[2]])
      S_pre$Vcr_hat <- Vcr_new
    }

    if (verbose) {
      safe_set_progress(pb, length(warm_in_stream) + length(new_in_stream))
    }

  } else {
    # ---- pure-R loop fallback (track_all_units = TRUE or Rcpp unavailable) ----
    for (u_idx in seq_along(new_in_stream)) {
      uid_ch  <- new_in_stream[u_idx]
      df_u    <- df_stream[stream_idx_list[[uid_ch]], ]
      times_u <- df_u[[time_col]]
      y_u     <- df_u[[y_col]]
      x_mat_u <- as.matrix(df_u[, x_cols, drop = FALSE])

      if (any(times_u > S_pre$T_support)) {
        warning(sprintf(
          "Unit '%s' has observations exceeding T_support (%d); those will be dropped.",
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
    cat(sprintf("\nDone: N = %s, n = %s, p = %d, T_support = %d | elapsed: %.2f sec\n",
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
      time_levels = time_levels   # original calendar time values (for reverse remapping)
    ),
    class = "otwfe"
  )
}


