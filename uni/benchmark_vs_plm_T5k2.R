# =============================================================================
# benchmark_vs_plm_T5k2.R
#
# otwfe_file() vs plm 속도/정확도 비교 — T=5, k=2
# DGP: Hwang & Lee (2026), §5.2
#
# N = 1M, 5M, 10M, 20M, 50M
#   - otwfe_file(): CSV 파일 기반 청크 처리 (전체 df 불필요)
#   - plm: 단계별(pdata.frame / plm / vcov / vcovHC) 개별 타이밍 측정
#          OOM 등 에러 발생 시 에러 메시지 기록 후 계속 진행
#   - 양쪽 모두 성공한 경우 theta/V0/Vcr 정확도 비교
#
# 실행 순서:
#   1. df 생성 → CSV 저장
#   2. otwfe_file() (CSV 기반, df는 메모리에 유지)
#   3. plm 단계별 실행 (df 필요)
#   4. rm(df) + gc()
#
# 작성: 2026-03-29
# =============================================================================

`%||%` <- function(a, b) if (!is.null(a) && !is.na(a)) a else b

setwd("/Users/joonhohwang/Desktop/claude_code/uni")
source("otwfe_file.R")
source("../Simulation/sim_dgp.R")
suppressPackageStartupMessages(library(plm))
suppressPackageStartupMessages(library(data.table))

# =============================================================================
# 설정
# =============================================================================
N_VALUES   <- c(1e6L, 5e6L, 1e7L, 2e7L, 5e7L)
T_PANEL    <- 5L
X_COLS     <- c("x1", "x2")   # k = 2
SEED       <- 42L
CHUNK_SIZE <- 2e6L             # otwfe_file chunk_size

# p = k + T - 1 = 2 + 5 - 1 = 6
cat(sprintf("설정: T=%d, k=%d → p=%d\n\n", T_PANEL, length(X_COLS),
            length(X_COLS) + T_PANEL - 1L))

# =============================================================================
# plm 단계별 실행 함수
# =============================================================================
run_plm_steps <- function(df, x_cols) {
  res <- list(
    pdata_ok  = FALSE, plm_ok   = FALSE,
    vcov_ok   = FALSE, vcovhc_ok = FALSE,
    t_prep    = NA_real_, t_coef = NA_real_,
    t_v0      = NA_real_, t_vcr  = NA_real_,
    fail_step = NA_character_, fail_msg = NA_character_,
    theta = NULL, V0 = NULL, Vcr = NULL
  )

  # Step 1: pdata.frame()
  t1  <- proc.time()
  pdf <- tryCatch(
    pdata.frame(df, index = c("id", "time")),
    error = function(e) {
      res$fail_step <<- "pdata.frame()"
      res$fail_msg  <<- conditionMessage(e)
      NULL
    }
  )
  res$t_prep <- (proc.time() - t1)[["elapsed"]]
  if (is.null(pdf)) return(res)
  res$pdata_ok <- TRUE

  # Step 2: plm()
  fml <- as.formula(paste("y ~", paste(x_cols, collapse = "+"),
                          "+ factor(time)"))
  t2      <- proc.time()
  plm_fit <- tryCatch(
    plm(fml, data = pdf, model = "within", effect = "individual"),
    error = function(e) {
      res$fail_step <<- "plm()"
      res$fail_msg  <<- conditionMessage(e)
      NULL
    }
  )
  res$t_coef <- (proc.time() - t2)[["elapsed"]]
  if (is.null(plm_fit)) return(res)
  res$plm_ok <- TRUE
  res$theta  <- coef(plm_fit)[x_cols]

  # Step 3: vcov()
  t3 <- proc.time()
  V0 <- tryCatch(
    vcov(plm_fit)[x_cols, x_cols],
    error = function(e) {
      res$fail_step <<- "vcov()"
      res$fail_msg  <<- conditionMessage(e)
      NULL
    }
  )
  res$t_v0 <- (proc.time() - t3)[["elapsed"]]
  if (!is.null(V0)) {
    res$vcov_ok <- TRUE
    res$V0      <- V0
  } else return(res)

  # Step 4: vcovHC()
  t4  <- proc.time()
  Vcr <- tryCatch(
    vcovHC(plm_fit, method = "arellano", type = "HC0")[x_cols, x_cols],
    error = function(e) {
      res$fail_step <<- "vcovHC()"
      res$fail_msg  <<- conditionMessage(e)
      NULL
    }
  )
  res$t_vcr <- (proc.time() - t4)[["elapsed"]]
  if (!is.null(Vcr)) {
    res$vcovhc_ok <- TRUE
    res$Vcr       <- Vcr
  }
  res
}

# =============================================================================
# 시간 포맷 헬퍼
# =============================================================================
fmt_t <- function(x) {
  if (is.na(x)) return("    --")
  if (x < 60)   sprintf("%5.1f초", x)
  else           sprintf("%4.1f분", x / 60)
}

# =============================================================================
# 메인 벤치마크 루프
# =============================================================================
sep72 <- strrep("=", 72)
sep_m <- strrep("-", 72)

cat(sprintf("%s\n", sep72))
cat(sprintf("  Benchmark: otwfe_file vs plm  (T=%d, k=2, p=6)\n", T_PANEL))
cat(sprintf("  DGP: Hwang & Lee (2026) §5.2  |  seed=%d\n", SEED))
cat(sprintf("  N  : %s\n",
            paste(format(N_VALUES, big.mark = ","), collapse = " / ")))
cat(sprintf("  otwfe_file chunk_size: %s행\n",
            format(CHUNK_SIZE, big.mark = ",")))
cat(sprintf("  Rcpp 배치: %s\n",
            if (.alg1_batch_rcpp_available) "사용" else "미사용"))
cat(sprintf("%s\n\n", sep72))

all_res <- vector("list", length(N_VALUES))

for (idx in seq_along(N_VALUES)) {
  N <- as.integer(N_VALUES[idx])

  cat(sprintf("\n%s\n  N = %s\n%s\n",
              sep_m, format(N, big.mark = ","), sep_m))

  rec <- list(
    N = N, N_obs = NA_integer_, data_mb = NA_real_,
    t_gen = NA_real_, csv_mb = NA_real_,
    # otwfe_file
    t_otwfe = NA_real_, state_mb = NA_real_,
    theta_ot = NULL,
    # plm
    t_prep = NA_real_, t_coef = NA_real_,
    t_v0   = NA_real_, t_vcr  = NA_real_,
    # 정확도
    theta_diff = NA_real_, V0_diff = NA_real_, Vcr_diff = NA_real_,
    speedup    = NA_real_,
    # 상태
    otwfe_ok     = FALSE, plm_ok_all   = FALSE,
    otwfe_fail   = NA_character_,
    plm_fail_step = NA_character_, plm_fail_msg = NA_character_
  )

  # ── 1. 데이터 생성 ──────────────────────────────────────────────────────────
  gen <- tryCatch(
    generate_panel_52(N, T = T_PANEL, T_min = 2L, seed = SEED, verbose = TRUE),
    error = function(e) {
      cat(sprintf("  FAIL @ 데이터 생성: %s\n", conditionMessage(e)))
      rec$otwfe_fail <<- sprintf("데이터 생성: %s", conditionMessage(e))
      NULL
    }
  )
  if (is.null(gen)) { all_res[[idx]] <- rec; gc(); next }

  df           <- gen$df
  rec$N_obs    <- gen$N_obs
  rec$data_mb  <- as.numeric(object.size(df)) / 1e6
  rec$t_gen    <- gen$t_gen
  cat(sprintf("  데이터 크기: %.0f MB  (관측치 %s)\n",
              rec$data_mb, format(gen$N_obs, big.mark = ",")))

  # ── 2. CSV 저장 ──────────────────────────────────────────────────────────────
  tmp_csv <- tempfile(fileext = ".csv")
  fwrite(df, tmp_csv)
  rec$csv_mb <- file.size(tmp_csv) / 1e6
  cat(sprintf("  CSV 저장: %.0f MB\n", rec$csv_mb))

  # ── 3. otwfe_file() ──────────────────────────────────────────────────────────
  cat("  [otwfe_file] 실행 중...\n")
  t_ot   <- proc.time()
  fit_ot <- tryCatch(
    otwfe_file(
      path       = tmp_csv,
      id_col     = "id",
      time_col   = "time",
      y_col      = "y",
      x_cols     = X_COLS,
      chunk_size = CHUNK_SIZE,
      verbose    = TRUE
    ),
    error = function(e) {
      cat(sprintf("  FAIL @ otwfe_file(): %s\n", conditionMessage(e)))
      rec$otwfe_fail <<- conditionMessage(e)
      NULL
    }
  )
  rec$t_otwfe <- (proc.time() - t_ot)[["elapsed"]]
  unlink(tmp_csv)  # CSV 삭제

  if (!is.null(fit_ot)) {
    rec$otwfe_ok <- TRUE
    rec$state_mb <- as.numeric(object.size(fit_ot$state)) / 1e6
    rec$theta_ot <- fit_ot$state$theta_hat[X_COLS]
    cat(sprintf("  [otwfe_file] 소요: %.1f초  |  state: %.4f MB\n",
                rec$t_otwfe, rec$state_mb))
  } else {
    cat(sprintf("  [otwfe_file] 실패 (%.1f초)\n", rec$t_otwfe))
  }

  # ── 4. plm 단계별 ────────────────────────────────────────────────────────────
  cat("  [plm]  각 단계 실행 중...\n")
  plm_res     <- run_plm_steps(df, X_COLS)
  rec$t_prep  <- plm_res$t_prep
  rec$t_coef  <- plm_res$t_coef
  rec$t_v0    <- plm_res$t_v0
  rec$t_vcr   <- plm_res$t_vcr

  cat(sprintf("    pdata.frame() : %s%s\n",
              fmt_t(plm_res$t_prep),
              if (plm_res$pdata_ok) ""
              else sprintf("  ← FAIL: %s", plm_res$fail_msg)))
  if (plm_res$pdata_ok)
    cat(sprintf("    plm()         : %s%s\n",
                fmt_t(plm_res$t_coef),
                if (plm_res$plm_ok) ""
                else sprintf("  ← FAIL: %s", plm_res$fail_msg)))
  if (plm_res$plm_ok)
    cat(sprintf("    vcov()        : %s%s\n",
                fmt_t(plm_res$t_v0),
                if (plm_res$vcov_ok) ""
                else sprintf("  ← FAIL: %s", plm_res$fail_msg)))
  if (plm_res$vcov_ok)
    cat(sprintf("    vcovHC()      : %s%s\n",
                fmt_t(plm_res$t_vcr),
                if (plm_res$vcovhc_ok) ""
                else sprintf("  ← FAIL: %s", plm_res$fail_msg)))

  if (!is.na(plm_res$fail_step)) {
    rec$plm_fail_step <- plm_res$fail_step
    rec$plm_fail_msg  <- plm_res$fail_msg
  }
  rec$plm_ok_all <- plm_res$vcovhc_ok

  # ── 5. 정확도 비교 (양쪽 성공 시) ───────────────────────────────────────────
  if (rec$otwfe_ok && plm_res$plm_ok) {
    rec$theta_diff <- max(abs(rec$theta_ot - plm_res$theta))
  }
  if (rec$otwfe_ok && plm_res$vcov_ok) {
    V0_ot      <- fit_ot$state$sigma2_hat *
                  fit_ot$state$inv_dotZtZ[X_COLS, X_COLS]
    rec$V0_diff <- max(abs(V0_ot - plm_res$V0))
  }
  if (rec$otwfe_ok && plm_res$vcovhc_ok) {
    rec$Vcr_diff <- max(abs(fit_ot$state$Vcr_hat[X_COLS, X_COLS] -
                              plm_res$Vcr))
    t_plm_total  <- sum(c(plm_res$t_prep, plm_res$t_coef,
                           plm_res$t_v0,   plm_res$t_vcr), na.rm = TRUE)
    rec$speedup  <- t_plm_total / rec$t_otwfe
  }

  if (rec$otwfe_ok && plm_res$plm_ok)
    cat(sprintf("\n  [정확도]  theta: %.2e  V0: %.2e  Vcr: %s\n",
                rec$theta_diff,
                rec$V0_diff %||% NA,
                if (!is.na(rec$Vcr_diff)) sprintf("%.2e", rec$Vcr_diff)
                else paste0("-- (", plm_res$fail_step, ")")))

  all_res[[idx]] <- rec
  rm(df, gen); if (!is.null(fit_ot)) rm(fit_ot)
  invisible(gc())
}

# =============================================================================
# 최종 요약 표
# =============================================================================
cat(sprintf("\n\n%s\n  최종 요약  (T=%d, k=2, p=6)\n%s\n",
            sep72, T_PANEL, sep72))

# 헤더
cat(sprintf("  %-11s  %-10s  %-9s  %-8s  %-8s  %-8s  %-8s  %-8s  %-8s  %s\n",
            "N", "관측치",
            "otwfe(s)", "pdata(s)", "plm(s)", "vcov(s)", "vcovHC(s)",
            "속도비", "θ_diff", "상태"))
cat(sprintf("%s\n", strrep("-", 110)))

for (r in all_res) {
  if (is.null(r)) next

  # 시간 포맷 (초 단위, NA → "-")
  fmt <- function(x) if (is.na(x)) "       -" else sprintf("%8.1f", x)

  otwfe_str  <- if (r$otwfe_ok)    sprintf("%8.1f", r$t_otwfe) else "    FAIL"
  speed_str  <- if (!is.na(r$speedup)) sprintf("%7.1fx", r$speedup) else "       -"
  diff_str   <- if (!is.na(r$theta_diff)) sprintf("%.2e", r$theta_diff) else "       -"

  status <- if (!r$otwfe_ok) {
    "FAIL(otwfe)"
  } else if (r$plm_ok_all) {
    "OK (both)"
  } else if (!is.na(r$plm_fail_step)) {
    sprintf("plm FAIL @ %s", r$plm_fail_step)
  } else {
    "OK(otwfe only)"
  }

  cat(sprintf("  %-11s  %-10s  %s  %s  %s  %s  %s  %s  %-8s  %s\n",
              format(r$N,    big.mark = ","),
              format(r$N_obs %||% NA, big.mark = ","),
              otwfe_str,
              fmt(r$t_prep),
              fmt(r$t_coef),
              fmt(r$t_v0),
              fmt(r$t_vcr),
              speed_str,
              diff_str,
              status))

  # 에러 메시지 출력
  if (!is.na(r$plm_fail_step))
    cat(sprintf("    └─ plm FAIL (%s): %s\n",
                r$plm_fail_step,
                substr(r$plm_fail_msg %||% "", 1, 80)))
  if (!r$otwfe_ok && !is.na(r$otwfe_fail))
    cat(sprintf("    └─ otwfe FAIL: %s\n",
                substr(r$otwfe_fail, 1, 80)))
}

cat(sprintf("%s\n", sep72))
cat(sprintf("  chunk_size: %s행/청크  |  Rcpp 배치: %s\n",
            format(CHUNK_SIZE, big.mark = ","),
            if (.alg1_batch_rcpp_available) "사용됨" else "미사용"))
cat(sprintf("%s\n", sep72))
