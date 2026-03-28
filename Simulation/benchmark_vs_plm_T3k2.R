# =============================================================================
# benchmark_vs_plm_T3k2.R
#
# otwfe() vs plm 속도 비교 — 대용량 데이터 한계 탐색
# DGP: Hwang & Lee (2026), §5.2  |  T=3, k=2
#
# N = 10^6, 5×10^6, 10^7, 2×10^7, 5×10^7, 10^8
# 각 단계(데이터 생성 / otwfe / pdata.frame / plm / vcov / vcovHC)를
# 개별 tryCatch로 감싸 실패 단계와 정확한 에러 메시지를 기록
#
# 최종 수정: 2026-03-27
# =============================================================================

# NULL 병합 연산자 (base R 호환)
`%||%` <- function(a, b) if (!is.null(a) && !is.na(a)) a else b

setwd("/Users/joonhohwang/Desktop/claude_code/Simulation")
source("sim_dgp.R")
source("../online_twfe_core.R")
suppressPackageStartupMessages(library(plm))

# =============================================================================
# 설정
# =============================================================================
N_VALUES    <- c(1e6, 5e6, 1e7, 2e7, 5e7, 1e8)
T_PANEL     <- 3L      # calendar time 수
X_COLS      <- c("x1", "x2")   # k = 2
SEED        <- 42L

# p = k + T - 1 = 2 + 3 - 1 = 4
cat(sprintf("설정: T=%d, k=%d → p=%d\n", T_PANEL, length(X_COLS),
            length(X_COLS) + T_PANEL - 1L))

# =============================================================================
# plm 단계별 실행 함수
# =============================================================================
run_plm_steps <- function(df, x_cols) {
  res <- list(
    pdata_ok=FALSE, plm_ok=FALSE, vcov_ok=FALSE, vcovhc_ok=FALSE,
    t_prep=NA_real_, t_coef=NA_real_, t_v0=NA_real_, t_vcr=NA_real_,
    fail_step=NA_character_, fail_msg=NA_character_,
    theta=NULL, V0=NULL, Vcr=NULL
  )

  # Step 1: pdata.frame()
  t1  <- proc.time()
  pdf <- tryCatch(
    pdata.frame(df, index=c("id","time")),
    error=function(e) {
      res$fail_step <<- "pdata.frame()"
      res$fail_msg  <<- conditionMessage(e)
      NULL
    }
  )
  res$t_prep <- (proc.time()-t1)[["elapsed"]]
  if (is.null(pdf)) return(res)
  res$pdata_ok <- TRUE

  # Step 2: plm()
  fml <- as.formula(paste("y ~", paste(x_cols, collapse="+"), "+ factor(time)"))
  t2      <- proc.time()
  plm_fit <- tryCatch(
    plm(fml, data=pdf, model="within", effect="individual"),
    error=function(e) {
      res$fail_step <<- "plm()"
      res$fail_msg  <<- conditionMessage(e)
      NULL
    }
  )
  res$t_coef <- (proc.time()-t2)[["elapsed"]]
  if (is.null(plm_fit)) return(res)
  res$plm_ok <- TRUE
  res$theta  <- coef(plm_fit)[x_cols]

  # Step 3: vcov()
  t3 <- proc.time()
  V0 <- tryCatch(
    vcov(plm_fit)[x_cols, x_cols],
    error=function(e) {
      res$fail_step <<- "vcov()"
      res$fail_msg  <<- conditionMessage(e)
      NULL
    }
  )
  res$t_v0 <- (proc.time()-t3)[["elapsed"]]
  if (!is.null(V0)) { res$vcov_ok <- TRUE; res$V0 <- V0 } else return(res)

  # Step 4: vcovHC()
  t4  <- proc.time()
  Vcr <- tryCatch(
    vcovHC(plm_fit, method="arellano", type="HC0")[x_cols, x_cols],
    error=function(e) {
      res$fail_step <<- "vcovHC()"
      res$fail_msg  <<- conditionMessage(e)
      NULL
    }
  )
  res$t_vcr <- (proc.time()-t4)[["elapsed"]]
  if (!is.null(Vcr)) { res$vcovhc_ok <- TRUE; res$Vcr <- Vcr }
  res
}

# =============================================================================
# 메인 벤치마크 루프
# =============================================================================
sep72 <- strrep("=", 72)
sep_m <- strrep("-", 72)

cat(sprintf("\n%s\n", sep72))
cat("  Benchmark: otwfe vs plm  (대용량 한계 탐색)\n")
cat(sprintf("  DGP: Hwang & Lee (2026) §5.2  |  T=%d, k=2, p=4, seed=%d\n",
            T_PANEL, SEED))
cat(sprintf("  Rcpp 배치: %s\n",
            if (.alg1_batch_rcpp_available) "사용" else "미사용"))
cat(sprintf("%s\n\n", sep72))

all_res <- vector("list", length(N_VALUES))

for (idx in seq_along(N_VALUES)) {
  N <- as.integer(N_VALUES[idx])

  cat(sprintf("\n%s\n  N = %s\n%s\n",
              sep_m, format(N, big.mark=","), sep_m))

  rec <- list(
    N=N, N_obs=NA_integer_, data_mb=NA_real_, t_gen=NA_real_,
    t_otwfe=NA_real_, state_mb=NA_real_,
    t_prep=NA_real_, t_coef=NA_real_, t_v0=NA_real_, t_vcr=NA_real_,
    otwfe_ok=FALSE, plm_ok_all=FALSE,
    otwfe_fail=NA_character_, plm_fail_step=NA_character_, plm_fail_msg=NA_character_,
    speedup=NA_real_, theta_diff=NA_real_, V0_diff=NA_real_, Vcr_diff=NA_real_
  )

  # ── 데이터 생성 ────────────────────────────────────────────────────────────
  gen <- tryCatch(
    generate_panel_52(N, T=T_PANEL, T_min=2L, seed=SEED, verbose=TRUE),
    error=function(e) {
      cat(sprintf("  ← FAIL @ 데이터 생성: %s\n", conditionMessage(e)))
      rec$otwfe_fail <<- sprintf("데이터 생성: %s", conditionMessage(e))
      NULL
    }
  )
  if (is.null(gen)) { all_res[[idx]] <- rec; gc(); next }

  df            <- gen$df
  rec$N_obs     <- gen$N_obs
  rec$data_mb   <- as.numeric(object.size(df)) / 1e6
  rec$t_gen     <- gen$t_gen
  cat(sprintf("  데이터 크기: %.0f MB  (관측치 %s)\n",
              rec$data_mb, format(gen$N_obs, big.mark=",")))

  # ── otwfe() ────────────────────────────────────────────────────────────────
  cat("  [otwfe] 실행 중...\n")
  t_ot  <- proc.time()
  fit   <- tryCatch(
    otwfe(data=df, id_col="id", time_col="time", y_col="y",
          x_cols=X_COLS, track_all_units=FALSE, verbose=TRUE,
          chunk_size=1e6L),
    error=function(e) {
      cat(sprintf("  ← FAIL @ otwfe(): %s\n", conditionMessage(e)))
      rec$otwfe_fail <<- conditionMessage(e)
      NULL
    }
  )
  rec$t_otwfe <- (proc.time()-t_ot)[["elapsed"]]

  if (!is.null(fit)) {
    rec$otwfe_ok  <- TRUE
    rec$state_mb  <- as.numeric(object.size(fit$state)) / 1e6
    cat(sprintf("  [otwfe] 소요: %.1f초  |  state: %.4f MB\n",
                rec$t_otwfe, rec$state_mb))
  } else {
    cat(sprintf("  [otwfe] 실패 (%.1f초)\n", rec$t_otwfe))
  }

  # ── plm ────────────────────────────────────────────────────────────────────
  cat("  [plm]  각 단계 실행 중...\n")
  plm_res <- run_plm_steps(df, X_COLS)
  rec$t_prep  <- plm_res$t_prep
  rec$t_coef  <- plm_res$t_coef
  rec$t_v0    <- plm_res$t_v0
  rec$t_vcr   <- plm_res$t_vcr

  fmt_t <- function(x) if (is.na(x)) "  --" else
    if (x < 60) sprintf("%5.1f초", x) else sprintf("%4.1f분", x/60)

  cat(sprintf("    pdata.frame() : %s%s\n", fmt_t(plm_res$t_prep),
              if (plm_res$pdata_ok) "" else
                sprintf("  ← FAIL: %s", plm_res$fail_msg)))
  if (plm_res$pdata_ok)
    cat(sprintf("    plm()         : %s%s\n", fmt_t(plm_res$t_coef),
                if (plm_res$plm_ok) "" else
                  sprintf("  ← FAIL: %s", plm_res$fail_msg)))
  if (plm_res$plm_ok)
    cat(sprintf("    vcov()        : %s%s\n", fmt_t(plm_res$t_v0),
                if (plm_res$vcov_ok) "" else
                  sprintf("  ← FAIL: %s", plm_res$fail_msg)))
  if (plm_res$vcov_ok)
    cat(sprintf("    vcovHC()      : %s%s\n", fmt_t(plm_res$t_vcr),
                if (plm_res$vcovhc_ok) "" else
                  sprintf("  ← FAIL: %s", plm_res$fail_msg)))

  if (!is.na(plm_res$fail_step)) {
    rec$plm_fail_step <- plm_res$fail_step
    rec$plm_fail_msg  <- plm_res$fail_msg
  }
  rec$plm_ok_all <- plm_res$vcovhc_ok

  # ── 정확도 (둘 다 성공한 경우만) ──────────────────────────────────────────
  if (rec$otwfe_ok && plm_res$plm_ok) {
    rec$theta_diff <- max(abs(fit$state$theta_hat[X_COLS] - plm_res$theta))
  }
  if (rec$otwfe_ok && plm_res$vcov_ok) {
    V0_on <- fit$state$sigma2_hat * fit$state$inv_dotZtZ[X_COLS, X_COLS]
    rec$V0_diff <- max(abs(V0_on - plm_res$V0))
  }
  if (rec$otwfe_ok && plm_res$vcovhc_ok) {
    rec$Vcr_diff <- max(abs(fit$state$Vcr_hat[X_COLS, X_COLS] - plm_res$Vcr))
    t_plm_total  <- sum(c(plm_res$t_prep, plm_res$t_coef,
                          plm_res$t_v0,   plm_res$t_vcr), na.rm=TRUE)
    rec$speedup  <- t_plm_total / rec$t_otwfe
  }

  if (rec$otwfe_ok && plm_res$plm_ok)
    cat(sprintf("\n  [정확도]  theta: %.2e  V0: %.2e  Vcr: %s\n",
                rec$theta_diff, rec$V0_diff,
                if (!is.na(rec$Vcr_diff)) sprintf("%.2e", rec$Vcr_diff)
                else paste0("-- (", plm_res$fail_step, ")")))

  all_res[[idx]] <- rec
  rm(df, gen); if (!is.null(fit)) rm(fit)
  invisible(gc())
}

# =============================================================================
# 최종 요약 표
# =============================================================================
cat(sprintf("\n\n%s\n  최종 요약  (T=%d, k=2, p=4)\n%s\n",
            sep72, T_PANEL, sep72))
cat(sprintf("  %-11s  %-9s  %-6s  %-9s  %-10s  %-6s  %-8s  %s\n",
            "N", "관측치", "MB", "otwfe", "plm합계", "속도비", "θ_diff", "상태"))
cat(sprintf("%s\n", strrep("-", 72)))

for (r in all_res) {
  if (is.null(r)) next

  otwfe_str <- if (r$otwfe_ok)  sprintf("%.1f초", r$t_otwfe)  else "FAIL"
  t_plm_all <- sum(c(r$t_prep, r$t_coef, r$t_v0, r$t_vcr), na.rm=TRUE)
  plm_str   <- if (r$plm_ok_all) sprintf("%.1f초", t_plm_all) else "FAIL"
  speed_str <- if (!is.na(r$speedup)) sprintf("%.1fx", r$speedup) else "--"
  diff_str  <- if (!is.na(r$theta_diff)) sprintf("%.2e", r$theta_diff) else "--"

  status <- if (!r$otwfe_ok) {
    sprintf("FAIL(otwfe) @ %s", substr(r$otwfe_fail %||% "?", 1, 30))
  } else if (r$plm_ok_all) {
    "OK"
  } else if (!is.na(r$plm_fail_step)) {
    sprintf("plm FAIL @ %s", r$plm_fail_step)
  } else {
    "OK(otwfe only)"
  }

  cat(sprintf("  %-11s  %-9s  %5.0f  %-9s  %-10s  %-6s  %-8s  %s\n",
              format(r$N, big.mark=","),
              format(r$N_obs %||% NA, big.mark=","),
              r$data_mb %||% NA,
              otwfe_str, plm_str, speed_str, diff_str, status))

  if (!r$plm_ok_all && !is.na(r$plm_fail_step))
    cat(sprintf("    └─ %s: %s\n",
                r$plm_fail_step,
                substr(r$plm_fail_msg %||% "", 1, 60)))
  if (!r$otwfe_ok && !is.na(r$otwfe_fail))
    cat(sprintf("    └─ otwfe: %s\n", substr(r$otwfe_fail, 1, 60)))
}

cat(sprintf("%s\n", sep72))
cat(sprintf("  Rcpp 배치: %s  |  otwfe chunk_size: 1,000,000 units\n",
            if (.alg1_batch_rcpp_available) "사용됨" else "미사용"))
cat(sprintf("%s\n", sep72))
