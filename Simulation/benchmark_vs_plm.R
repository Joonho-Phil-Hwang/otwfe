# =============================================================================
# benchmark_vs_plm.R
#
# otwfe() vs plm 속도 비교 벤치마크
# DGP: Hwang & Lee (2026), Subsection 5.2
#
# 비교 항목:
#   otwfe()  : 계수 + classical variance + cluster-robust variance (1회 호출)
#   plm 측   : pdata.frame() + plm() + vcov() + vcovHC() (4단계 개별 측정)
#
# plm 실패 시: 실패 단계명과 정확한 에러 메시지를 기록
# =============================================================================

setwd("/Users/joonhohwang/Desktop/claude_code/Simulation")

source("sim_dgp.R")
source("../online_twfe_core.R")
suppressPackageStartupMessages(library(plm))

# =============================================================================
# 설정
# =============================================================================

N_VALUES    <- c(100000L, 500000L, 1000000L, 3000000L, 5000000L, 8500000L)
X_COLS      <- c("x1", "x2")
VCR_TIMEOUT <- 600L   # vcovHC 최대 허용 시간 (초)
SEED        <- 42L
T_PANEL     <- 5L

# =============================================================================
# plm 벤치마크 실행 함수
#   각 단계를 tryCatch로 감싸 실패 단계·에러 메시지·소요 시간을 정확히 기록
# =============================================================================
run_plm_steps <- function(df, x_cols, timeout) {

  res <- list(
    pdata_ok  = FALSE,  plm_ok   = FALSE,
    vcov_ok   = FALSE,  vcovhc_ok = FALSE,
    t_prep    = NA_real_, t_coef  = NA_real_,
    t_v0      = NA_real_, t_vcr   = NA_real_,
    fail_step = NA_character_,
    fail_msg  = NA_character_,
    theta = NULL, V0 = NULL, Vcr = NULL
  )

  # ── Step 1: pdata.frame() ─────────────────────────────────────────────────
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

  # ── Step 2: plm() ─────────────────────────────────────────────────────────
  fml <- as.formula(paste("y ~", paste(x_cols, collapse = "+"), "+ factor(time)"))
  t2       <- proc.time()
  plm_fit  <- tryCatch(
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

  # ── Step 3: vcov() ────────────────────────────────────────────────────────
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
  } else {
    return(res)
  }

  # ── Step 4: vcovHC()  (타임아웃 적용) ────────────────────────────────────
  t4  <- proc.time()
  Vcr <- tryCatch({
    setTimeLimit(elapsed = timeout, transient = TRUE)
    out <- vcovHC(plm_fit, method = "arellano", type = "HC0")[x_cols, x_cols]
    setTimeLimit(elapsed = Inf, transient = TRUE)
    out
  }, error = function(e) {
    setTimeLimit(elapsed = Inf, transient = TRUE)
    msg <- conditionMessage(e)
    if (grepl("elapsed time limit", msg, ignore.case = TRUE)) {
      res$fail_step <<- sprintf("vcovHC() TIMEOUT (>%ds)", timeout)
      res$fail_msg  <<- sprintf("reached elapsed time limit (%ds)", timeout)
    } else {
      res$fail_step <<- "vcovHC()"
      res$fail_msg  <<- msg
    }
    NULL
  })
  res$t_vcr <- (proc.time() - t4)[["elapsed"]]
  if (!is.null(Vcr)) {
    res$vcovhc_ok <- TRUE
    res$Vcr       <- Vcr
  }

  res
}

# =============================================================================
# 결과 출력 헬퍼
# =============================================================================
fmt_time <- function(x) {
  if (is.na(x)) return("  --  ")
  if (x < 60)   return(sprintf("%5.1f초", x))
  sprintf("%4.1f분", x / 60)
}

fmt_diff <- function(x) {
  if (is.na(x) || is.null(x)) return("   --  ")
  sprintf("%.2e", x)
}

# =============================================================================
# 메인 벤치마크 루프
# =============================================================================
cat(sprintf("\n%s\n", strrep("=", 72)))
cat("  Benchmark: otwfe vs plm\n")
cat(sprintf("  DGP: Hwang & Lee (2026) §5.2  |  T=%d, k=2, seed=%d\n",
            T_PANEL, SEED))
cat(sprintf("  vcovHC timeout: %d초\n", VCR_TIMEOUT))
cat(sprintf("%s\n\n", strrep("=", 72)))

all_res <- vector("list", length(N_VALUES))

for (idx in seq_along(N_VALUES)) {
  N <- N_VALUES[idx]
  sep <- strrep("-", 72)
  cat(sprintf("\n%s\n", sep))
  cat(sprintf("  N = %s\n", format(N, big.mark = ",")))
  cat(sprintf("%s\n", sep))

  # ── 데이터 생성 ────────────────────────────────────────────────────────────
  gen     <- generate_panel_52(N, T = T_PANEL, seed = SEED, verbose = TRUE)
  df      <- gen$df
  data_mb <- as.numeric(object.size(df)) / 1e6
  cat(sprintf("  데이터 크기: %.0f MB  (관측치 %s)\n",
              data_mb, format(gen$N_obs, big.mark = ",")))

  # ── otwfe() ────────────────────────────────────────────────────────────────
  cat("\n  [otwfe] 실행 중...\n")
  t_ot <- proc.time()
  fit  <- otwfe(data = df,
                id_col   = "id", time_col = "time",
                y_col    = "y",  x_cols   = X_COLS,
                track_all_units = FALSE,
                verbose  = TRUE)
  t_otwfe  <- (proc.time() - t_ot)[["elapsed"]]
  state_mb <- as.numeric(object.size(fit$state)) / 1e6
  cat(sprintf("  [otwfe] 소요: %.2f초  |  state 크기: %.4f MB\n",
              t_otwfe, state_mb))

  # ── plm ────────────────────────────────────────────────────────────────────
  cat("\n  [plm]  각 단계 실행 중...\n")
  plm_res <- run_plm_steps(df, X_COLS, VCR_TIMEOUT)

  # plm 단계별 결과 출력
  cat(sprintf("    pdata.frame() : %s%s\n",
              fmt_time(plm_res$t_prep),
              if (plm_res$pdata_ok) "" else
                sprintf("  ← FAIL: %s", plm_res$fail_msg)))
  if (plm_res$pdata_ok) {
    cat(sprintf("    plm()         : %s%s\n",
                fmt_time(plm_res$t_coef),
                if (plm_res$plm_ok) "" else
                  sprintf("  ← FAIL: %s", plm_res$fail_msg)))
  }
  if (plm_res$plm_ok) {
    cat(sprintf("    vcov()        : %s%s\n",
                fmt_time(plm_res$t_v0),
                if (plm_res$vcov_ok) "" else
                  sprintf("  ← FAIL: %s", plm_res$fail_msg)))
  }
  if (plm_res$vcov_ok) {
    cat(sprintf("    vcovHC()      : %s%s\n",
                fmt_time(plm_res$t_vcr),
                if (plm_res$vcovhc_ok) "" else
                  sprintf("  ← FAIL: %s", plm_res$fail_msg)))
  }

  # ── 정확도 비교 (plm 성공한 항목만) ────────────────────────────────────────
  theta_diff <- V0_diff <- Vcr_diff <- NA_real_
  if (plm_res$plm_ok) {
    theta_diff <- max(abs(fit$state$theta_hat[X_COLS] - plm_res$theta))
  }
  if (plm_res$vcov_ok) {
    V0_on  <- fit$state$sigma2_hat * fit$state$inv_dotZtZ[X_COLS, X_COLS]
    V0_diff <- max(abs(V0_on - plm_res$V0))
  }
  if (plm_res$vcovhc_ok) {
    Vcr_diff <- max(abs(fit$state$Vcr_hat[X_COLS, X_COLS] - plm_res$Vcr))
  }
  if (plm_res$plm_ok) {
    cat(sprintf("\n  [정확도]  theta: %s  |  V0: %s  |  Vcr: %s\n",
                fmt_diff(theta_diff), fmt_diff(V0_diff), fmt_diff(Vcr_diff)))
  }

  # 속도 비율 (plm 전체 vs otwfe)
  t_plm_total <- sum(c(plm_res$t_prep, plm_res$t_coef,
                        plm_res$t_v0,   plm_res$t_vcr),
                     na.rm = TRUE)
  speedup <- if (!is.na(t_plm_total) && t_plm_total > 0 && plm_res$vcovhc_ok)
               t_plm_total / t_otwfe else NA_real_

  # 결과 저장
  all_res[[idx]] <- list(
    N          = N,
    N_obs      = gen$N_obs,
    data_mb    = data_mb,
    t_gen      = gen$t_gen,
    t_otwfe    = t_otwfe,
    state_mb   = state_mb,
    t_prep     = plm_res$t_prep,
    t_coef     = plm_res$t_coef,
    t_v0       = plm_res$t_v0,
    t_vcr      = plm_res$t_vcr,
    plm_ok_all = plm_res$vcovhc_ok,
    fail_step  = plm_res$fail_step,
    fail_msg   = plm_res$fail_msg,
    speedup    = speedup,
    theta_diff = theta_diff,
    V0_diff    = V0_diff,
    Vcr_diff   = Vcr_diff
  )

  # 메모리 해제
  rm(df, fit, gen)
  if (plm_res$pdata_ok) rm(plm_res)
  invisible(gc())
}

# =============================================================================
# 최종 요약 표
# =============================================================================
cat(sprintf("\n\n%s\n", strrep("=", 72)))
cat("  최종 요약\n")
cat(sprintf("%s\n", strrep("=", 72)))
cat(sprintf("  %-12s  %-9s  %-7s  %-9s  %-8s  %-7s  %-8s  %-7s\n",
            "N", "관측치", "MB", "otwfe(초)", "plm합계", "속도비",
            "θ_diff", "상태"))
cat(sprintf("%s\n", strrep("-", 72)))

for (r in all_res) {
  if (is.null(r)) next

  t_plm_str <- if (r$plm_ok_all) {
    sprintf("%.1f초", sum(c(r$t_prep, r$t_coef, r$t_v0, r$t_vcr), na.rm=TRUE))
  } else {
    "FAIL"
  }

  speed_str <- if (!is.na(r$speedup)) sprintf("%.1fx", r$speedup) else "  --"
  diff_str  <- if (!is.na(r$theta_diff)) sprintf("%.2e", r$theta_diff) else "  --"

  status_str <- if (r$plm_ok_all) {
    "OK"
  } else if (!is.na(r$fail_step)) {
    sprintf("FAIL @ %s", r$fail_step)
  } else {
    "FAIL"
  }

  cat(sprintf("  %-12s  %-9s  %5.0f  %8.2f초  %-8s  %-7s  %-8s  %s\n",
              format(r$N, big.mark=","),
              format(r$N_obs, big.mark=","),
              r$data_mb,
              r$t_otwfe,
              t_plm_str,
              speed_str,
              diff_str,
              status_str))

  # 실패 메시지 출력
  if (!is.na(r$fail_step) && !r$plm_ok_all) {
    cat(sprintf("    └─ %s: %s\n", r$fail_step,
                substr(r$fail_msg, 1L, 65L)))
  }
}

cat(sprintf("%s\n", strrep("=", 72)))
cat(sprintf("  vcovHC timeout 설정: %d초\n", VCR_TIMEOUT))
cat(sprintf("  otwfe Rcpp 배치 모듈: %s\n",
            if (.alg1_batch_rcpp_available) "사용됨" else "미사용 (pure-R)"))
cat(sprintf("%s\n", strrep("=", 72)))
