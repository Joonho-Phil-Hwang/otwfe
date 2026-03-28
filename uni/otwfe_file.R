# =============================================================================
# otwfe_file.R
#
# 대용량 CSV 패널 데이터를 메모리 초과 없이 분석하는 범용 함수
# otwfe_init / otwfe_update / otwfe_finalize 파이프라인을 내부적으로 사용
#
# 핵심 특성:
#   - 전체 데이터를 메모리에 올리지 않음
#   - chunk 단위 적응형 읽기 (unit 경계 자동 탐지)
#   - time 컬럼 전체 스캔으로 T_support 자동 탐지
#   - time 값 자동 재인덱싱 (1..T), 임의 정수/연도 형식 지원
#   - plm과 algebraically exact하게 동일한 결과 보장
#
# 전제 조건: 파일이 id 기준 비내림차순으로 정렬돼 있어야 함
#
# 작성: 2026-03-28
# =============================================================================

# --------------------------------------------------------------------------
# 의존성 확인
# --------------------------------------------------------------------------
if (!requireNamespace("data.table", quietly = TRUE))
  stop("data.table 패키지가 필요합니다: install.packages('data.table')")

# online_twfe_core.R 자동 탐지 및 로드
.uni_core_path <- local({
  self     <- tryCatch(normalizePath(sys.frame(1L)$ofile), error = function(e) NULL)
  args     <- commandArgs(trailingOnly = FALSE)
  rscript  <- tryCatch({
    f <- sub("--file=", "", args[grep("--file=", args)])
    if (length(f) == 1L && nchar(f) > 0L) normalizePath(f) else NULL
  }, error = function(e) NULL)
  cand     <- Filter(Negate(is.null), list(self, rscript))
  dir_path <- if (length(cand) > 0L) dirname(cand[[1L]]) else "."
  file.path(dir_path, "..", "online_twfe_core.R")
})
if (!exists("otwfe_finalize", mode = "function"))
  source(.uni_core_path, local = FALSE)

# =============================================================================
# 내부 헬퍼 함수
# =============================================================================

# --------------------------------------------------------------------------
# .csv_col_names(): CSV 헤더에서 컬럼명 읽기
# --------------------------------------------------------------------------
.csv_col_names <- function(path, sep) {
  names(data.table::fread(path, sep = sep, nrows = 0L, showProgress = FALSE))
}

# --------------------------------------------------------------------------
# .csv_read_chunk(): CSV에서 row_offset 이후 n_rows 행 읽기
#
# row_offset: 이미 읽은 데이터 행 수 (헤더 제외)
# n_rows    : 읽을 행 수
# 반환값    : data.frame 또는 NULL (EOF)
# --------------------------------------------------------------------------
.csv_read_chunk <- function(path, col_names, sep, row_offset, n_rows) {
  tryCatch({
    if (row_offset == 0L) {
      # 첫 번째 읽기: 헤더 포함
      chunk <- data.table::fread(path, sep = sep, nrows = n_rows,
                                  showProgress = FALSE)
    } else {
      # 이후 읽기: 헤더(1행) + 기읽은 데이터 행 건너뜀
      chunk <- data.table::fread(
        path, sep = sep,
        skip      = row_offset + 1L,   # +1은 헤더 행
        nrows     = n_rows,
        header    = FALSE,
        col.names = col_names,
        showProgress = FALSE
      )
    }
    if (nrow(chunk) == 0L) return(NULL)
    as.data.frame(chunk)
  }, error = function(e) NULL)
}

# --------------------------------------------------------------------------
# .detect_time_levels(): time 컬럼만 읽어 고유 calendar time 탐지
#
# time 컬럼 1개만 메모리에 올림 → 전체 데이터 대비 1/(k+3) 수준
# --------------------------------------------------------------------------
.detect_time_levels <- function(path, time_col, sep, verbose) {
  if (verbose) cat("  time 컬럼 스캔 중...\n")
  tmp  <- data.table::fread(path, sep = sep, select = time_col,
                              showProgress = FALSE)
  vals <- sort(unique(tmp[[1L]]))
  rm(tmp); invisible(gc())
  vals
}

# --------------------------------------------------------------------------
# .find_unit_boundary(): 마지막 완전한 unit의 마지막 행 인덱스 반환
#
# 반환값: 정수 인덱스 또는 NA (전체가 하나의 unit인 경우)
# --------------------------------------------------------------------------
.find_unit_boundary <- function(df, id_col) {
  ids     <- df[[id_col]]
  n       <- length(ids)
  if (n == 0L) return(NA_integer_)
  last_id <- ids[n]
  non_last <- which(ids != last_id)
  if (length(non_last) == 0L) return(NA_integer_)  # 전체가 같은 unit
  as.integer(max(non_last))
}

# =============================================================================
# 주요 함수: otwfe_file()
# =============================================================================
#'
#' 대용량 CSV 파일에서 Two-Way Fixed Effects 패널 회귀
#'
#' @param path       CSV 파일 경로 (id 기준 정렬 필수)
#' @param id_col     개인 식별자 컬럼명
#' @param time_col   calendar time 컬럼명
#' @param y_col      종속변수 컬럼명
#' @param x_cols     공변량 컬럼명 벡터
#' @param chunk_size 청크 당 읽을 최대 행 수 (기본 1,000,000)
#' @param sep        구분자 (기본 ",")
#' @param verbose    진행 상황 출력 여부
#' @return class "otwfe" 객체 (theta_hat, Vcr_hat 등 포함)
#'
otwfe_file <- function(path,
                       id_col,
                       time_col,
                       y_col,
                       x_cols,
                       chunk_size = 1e6L,
                       sep        = ",",
                       verbose    = TRUE) {

  chunk_size <- as.integer(chunk_size)
  t_total    <- proc.time()
  sep72      <- strrep("=", 72)
  sep_m      <- strrep("-", 72)

  # -----------------------------------------------------------------------
  # Step 0: 기본 검증
  # -----------------------------------------------------------------------
  if (!file.exists(path))
    stop(sprintf("파일을 찾을 수 없습니다: '%s'", path))

  col_names <- .csv_col_names(path, sep)
  for (col in c(id_col, time_col, y_col, x_cols))
    if (!col %in% col_names)
      stop(sprintf("컬럼을 찾을 수 없습니다: '%s'", col))

  if (verbose) {
    cat(sprintf("\n%s\n", sep72))
    cat(sprintf("  otwfe_file: %s\n", basename(path)))
    cat(sprintf("  id=%s  time=%s  y=%s  x=(%s)\n",
                id_col, time_col, y_col, paste(x_cols, collapse = ", ")))
    cat(sprintf("  chunk_size = %s 행\n", format(chunk_size, big.mark = ",")))
    cat(sprintf("%s\n", sep72))
  }

  # -----------------------------------------------------------------------
  # Step 1: T_support 자동 탐지 (time 컬럼 전체 스캔)
  # -----------------------------------------------------------------------
  if (verbose) cat("\n[Step 1] T_support 탐지\n")
  t1           <- proc.time()
  time_levels  <- .detect_time_levels(path, time_col, sep, verbose)
  T_support    <- length(time_levels)
  # 실제 time 값 → 연속 정수 인덱스 (1..T) 매핑
  time_remap   <- setNames(seq_along(time_levels), as.character(time_levels))
  baseline_idx <- 1L   # 재인덱싱 후 time=1이 기준

  if (verbose)
    cat(sprintf("  T_support = %d  |  time 범위: %s ~ %s  [%.1f초]\n",
                T_support,
                time_levels[1L], time_levels[T_support],
                (proc.time() - t1)[["elapsed"]]))

  # -----------------------------------------------------------------------
  # Step 2: 초기화 청크 구성
  #   모든 T calendar time이 포함될 때까지 청크를 누적
  #   → warm-up이 반드시 full T_support를 커버할 수 있도록 보장
  # -----------------------------------------------------------------------
  if (verbose) cat("\n[Step 2] 초기화 청크 구성 (모든 T time 포함 보장)\n")
  t2        <- proc.time()
  rows_read <- 0L
  init_df   <- NULL
  n_init_chunks <- 0L

  repeat {
    raw <- .csv_read_chunk(path, col_names, sep, rows_read, chunk_size)
    if (is.null(raw) || nrow(raw) == 0L) break

    rows_read     <- rows_read + nrow(raw)
    n_init_chunks <- n_init_chunks + 1L

    # time 값 재인덱싱 (1..T)
    raw[[time_col]] <- time_remap[as.character(raw[[time_col]])]
    init_df <- if (is.null(init_df)) raw else rbind(init_df, raw)

    times_found <- length(unique(init_df[[time_col]]))
    if (verbose)
      cat(sprintf("  초기화 청크 %d 누적: %s행,  time 커버 %d/%d\n",
                  n_init_chunks,
                  format(nrow(init_df), big.mark = ","),
                  times_found, T_support))

    if (times_found >= T_support) break   # 모든 time 포함됨
    if (nrow(raw) < chunk_size)    break   # EOF
  }
  rm(raw)

  if (is.null(init_df))
    stop("파일에서 데이터를 읽을 수 없습니다.")

  times_in_init <- length(unique(init_df[[time_col]]))
  if (times_in_init < T_support)
    warning(sprintf(
      "초기화 청크가 일부 calendar time을 포함하지 않습니다 (%d/%d). warm-up 품질이 저하될 수 있습니다.",
      times_in_init, T_support))

  # 마지막 unit 경계 탐지: 불완전 unit 행은 다음 청크로 이월
  bnd <- .find_unit_boundary(init_df, id_col)
  if (is.na(bnd)) {
    # 전체 init_df가 하나의 unit (극히 드문 경우)
    first_chunk <- init_df
    carry_over  <- NULL
  } else {
    first_chunk <- init_df[seq_len(bnd),           , drop = FALSE]
    carry_over  <- init_df[(bnd + 1L):nrow(init_df), , drop = FALSE]
  }
  rm(init_df); invisible(gc())

  if (verbose)
    cat(sprintf("  → 첫 청크: %s행  |  이월: %s행  [%.1f초]\n",
                format(nrow(first_chunk), big.mark = ","),
                format(if (!is.null(carry_over)) nrow(carry_over) else 0L,
                       big.mark = ","),
                (proc.time() - t2)[["elapsed"]]))

  # -----------------------------------------------------------------------
  # Step 3: State 초기화 (warm-up 선택 + 첫 청크 Rcpp 배치)
  # -----------------------------------------------------------------------
  if (verbose) cat("\n[Step 3] State 초기화\n")
  handle <- otwfe_init(
    x_cols        = x_cols,
    time_col      = time_col,
    id_col        = id_col,
    y_col         = y_col,
    T_support     = T_support,
    baseline_time = baseline_idx,
    verbose       = verbose
  )
  handle <- otwfe_update(handle, first_chunk)
  rm(first_chunk); invisible(gc())

  # -----------------------------------------------------------------------
  # Step 4: 이후 청크 반복 처리 (적응형 unit 경계 탐지)
  # -----------------------------------------------------------------------
  if (verbose) cat("\n[Step 4] 청크 처리\n")
  chunk_idx <- n_init_chunks  # Step 2에서 읽은 청크 수부터 시작

  repeat {
    raw <- .csv_read_chunk(path, col_names, sep, rows_read, chunk_size)
    if (is.null(raw) || nrow(raw) == 0L) break

    rows_read <- rows_read + nrow(raw)
    chunk_idx <- chunk_idx + 1L
    is_eof    <- nrow(raw) < chunk_size

    # time 재인덱싱
    raw[[time_col]] <- time_remap[as.character(raw[[time_col]])]

    # 정렬 순서 검증: carry_over 마지막 id <= raw 첫 번째 id
    if (!is.null(carry_over) && nrow(raw) > 0L) {
      if (carry_over[[id_col]][nrow(carry_over)] > raw[[id_col]][1L])
        stop(sprintf(
          "id 정렬 오류 (청크 %d): 파일이 id 기준으로 정렬되지 않았습니다.",
          chunk_idx))
    }

    # 이전 이월 행과 합치기
    combined <- if (is.null(carry_over)) raw else rbind(carry_over, raw)
    rm(raw); invisible(gc())

    if (is_eof) {
      # 파일 끝: 모든 unit 완전 처리
      carry_over <- NULL
      handle     <- otwfe_update(handle, combined)
      rm(combined); invisible(gc())
      break
    }

    # 마지막 unit 경계 탐지
    bnd <- .find_unit_boundary(combined, id_col)

    if (is.na(bnd)) {
      # 청크 전체가 하나의 unit (unit이 chunk_size보다 큰 경우)
      # → 계속 누적 (다음 청크에서 경계 탐지)
      carry_over <- combined
    } else {
      carry_over <- combined[(bnd + 1L):nrow(combined), , drop = FALSE]
      complete   <- combined[seq_len(bnd),              , drop = FALSE]
      rm(combined)
      handle     <- otwfe_update(handle, complete)
      rm(complete); invisible(gc())
    }

    if (verbose && chunk_idx %% 5L == 0L)
      cat(sprintf("  청크 %d 완료 | 누적 %s행 처리\n",
                  chunk_idx, format(rows_read, big.mark = ",")))
  }

  # 파일 마지막에 남은 carry_over 처리
  if (!is.null(carry_over) && nrow(carry_over) > 0L) {
    handle <- otwfe_update(handle, carry_over)
    rm(carry_over); invisible(gc())
  }

  # -----------------------------------------------------------------------
  # Step 5: 최종 결과 추출
  # -----------------------------------------------------------------------
  if (verbose) cat("\n[Step 5] Finalize\n")
  result <- otwfe_finalize(handle)

  # time 레이블 복원 (재인덱싱 이전 원래 값으로)
  # theta_hat의 time dummy 이름을 원래 time 값으로 복원
  orig_names <- names(result$state$theta_hat)
  for (i in seq_along(time_levels)) {
    orig_names <- gsub(
      paste0("(?<![0-9])", i, "(?![0-9])"),
      as.character(time_levels[i]),
      orig_names, perl = TRUE
    )
  }

  t_elapsed <- (proc.time() - t_total)[["elapsed"]]
  if (verbose) {
    cat(sprintf("\n%s\n", sep_m))
    cat(sprintf("  완료: N=%s units,  n=%s obs,  p=%d\n",
                format(result$state$N, big.mark = ","),
                format(result$state$n, big.mark = ","),
                result$state$p))
    cat(sprintf("  총 소요: %.1f초  (%.1f분)\n", t_elapsed, t_elapsed / 60))
    cat(sprintf("  state 크기: %.4f MB\n",
                as.numeric(object.size(result$state)) / 1e6))
    cat(sprintf("\n  추정 결과 (x 변수):\n"))
    cat(sprintf("    %-10s  %8s  %8s  %8s  %7s\n",
                "", "coef", "SE", "SE(CR)", "t"))
    cat(sprintf("    %s\n", strrep("-", 52)))
    theta <- result$state$theta_hat[x_cols]
    se_cl <- sqrt(diag(result$state$sigma2_hat *
                       result$state$inv_dotZtZ)[x_cols])
    se_cr <- sqrt(diag(result$state$Vcr_hat)[x_cols])
    for (xn in x_cols)
      cat(sprintf("    %-10s  %8.4f  %8.4f  %8.4f  %7.2f\n",
                  xn, theta[xn], se_cl[xn], se_cr[xn],
                  theta[xn] / se_cr[xn]))
    cat(sprintf("    %s\n", strrep("-", 52)))
    cat("    * SE    : classical variance 기반 표준오차\n")
    cat("    * SE(CR): HC0 cluster-robust variance 기반 표준오차 (Arellano)\n")
    cat(sprintf("%s\n", sep_m))
  }

  # 원래 time level 정보를 결과에 추가
  result$time_levels_original <- time_levels
  result$time_remap            <- time_remap
  result
}
