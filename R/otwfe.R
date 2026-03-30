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
      stop(sprintf("Column not found: '%s'", col))
  }

  chunk_df[[time_col]] <- as.integer(chunk_df[[time_col]])
  all_times_chunk      <- sort(unique(chunk_df[[time_col]]))

  # T_support 결정 (미지정이면 첫 청크의 max time)
  if (is.null(handle$T_support)) {
    T_support <- max(all_times_chunk)
    if (verbose) message(sprintf("T_support auto-detected from first chunk: %d", T_support))
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
    stop("No eligible units (>= 2 observations) in first chunk.")

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
    cat(sprintf("Initializing state (warm-up: %d units)...\n", length(warmup_ids)))

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
# otwfe_finalize(): 최종 결과 추출
# =============================================================================
#' Online TWFE 최종 결과 추출
#'
#' @param handle    otwfe_update() 완료 후 handle
#' @return class \code{"otwfe"} 객체
#' @export
otwfe_finalize <- function(handle) {
  if (is.null(handle$S))
    stop("Call otwfe_update() before otwfe_finalize().")

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
    if (chunk_size < 1L) stop("chunk_size must be a positive integer.")

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
  # 1. 입력 검증 및 전처리 (기존 단일 호출 경로)
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
      stop("baseline_time value not found in data.")
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
      stop("None of the specified warmup_ids are eligible (require >= 2 observations).")

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
    stop("No eligible warm-up units found. Each unit requires at least 2 observations.")

  df_warm   <- data[id_ch_vec %in% warmup_ids, ]
  df_stream <- data[!(id_ch_vec %in% warmup_ids), ]

  # --------------------------------------------------
  # 3. 초기 State 구성 (배치 OLS — Alg 1/2/3와 algebraically equivalent)
  # --------------------------------------------------
  if (verbose) {
    cat(sprintf("Initializing state (warm-up: %d units)...\n", length(warmup_ids)))
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
    cat(sprintf("Streaming: %s units (%d warm-up / %s new), %s observations\n",
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
        "%d observation(s) exceed T_support (%d) and will be dropped.",
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
      time_levels = time_levels   # 원래 calendar time 값 (재인덱싱 역변환용)
    ),
    class = "otwfe"
  )
}


