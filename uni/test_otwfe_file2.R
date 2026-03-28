# =============================================================================
# test_otwfe_file2.R
#
# 대용량 검증: T=5, k=2, N = 5,000,000 / 10,000,000 / 20,000,000
# 검증 방식:
#   1. generate_panel_52()로 데이터 생성 → CSV 저장 → rm(df); gc()
#   2. otwfe_file()로 CSV 파일 분석 (전체 df 메모리 미보유 상태)
#   3. theta vs beta_true = (1.3, -0.7) 비교
#   4. 소요 시간 / 상태 크기 기록
#   5. plm은 OOM 예상되므로 생략 (N=5M+ 에서는 실행 불가)
#
# 작성: 2026-03-29
# =============================================================================

setwd("/Users/joonhohwang/Desktop/claude_code/uni")
source("otwfe_file.R")
source("../Simulation/sim_dgp.R")
suppressPackageStartupMessages(library(data.table))

# =============================================================================
# 설정
# =============================================================================
N_VALUES   <- c(5e6L, 1e7L, 2e7L)
T_PANEL    <- 5L
X_COLS     <- c("x1", "x2")
SEED       <- 42L
CHUNK_SIZE <- 2e6L       # 2M행 / 청크 (~80MB / 청크)
BETA_TRUE  <- c(x1 = 1.3, x2 = -0.7)

sep72 <- strrep("=", 72)
sep_m <- strrep("-", 72)

cat(sprintf("%s\n", sep72))
cat(sprintf("  대용량 otwfe_file() 검증: T=%d, k=2, chunk_size=%s\n",
            T_PANEL, format(CHUNK_SIZE, big.mark = ",")))
cat(sprintf("  N = %s\n",
            paste(format(N_VALUES, big.mark = ","), collapse = " / ")))
cat(sprintf("  beta_true: x1=%.1f, x2=%.1f\n", BETA_TRUE["x1"], BETA_TRUE["x2"]))
cat(sprintf("%s\n\n", sep72))

all_res <- vector("list", length(N_VALUES))

for (idx in seq_along(N_VALUES)) {
  N <- as.integer(N_VALUES[idx])
  cat(sprintf("\n%s\n  N = %s\n%s\n",
              sep_m, format(N, big.mark = ","), sep_m))

  rec <- list(N = N, N_obs = NA_integer_,
              t_gen = NA_real_, csv_mb = NA_real_,
              t_file = NA_real_, state_mb = NA_real_,
              theta_x1 = NA_real_, theta_x2 = NA_real_,
              se_cl_x1 = NA_real_, se_cl_x2 = NA_real_,
              se_cr_x1 = NA_real_, se_cr_x2 = NA_real_,
              bias_x1 = NA_real_, bias_x2 = NA_real_,
              ok = FALSE, fail_msg = NA_character_)

  # ------------------------------------------------------------------
  # 1. 데이터 생성 → CSV 저장 → df 즉시 삭제
  # ------------------------------------------------------------------
  cat("  [1/3] 데이터 생성 및 CSV 저장...\n")
  t_gen <- proc.time()
  gen <- tryCatch(
    generate_panel_52(N, T = T_PANEL, T_min = 2L, seed = SEED, verbose = TRUE),
    error = function(e) {
      cat(sprintf("  FAIL @ 데이터 생성: %s\n", conditionMessage(e)))
      rec$fail_msg <<- conditionMessage(e)
      NULL
    }
  )
  if (is.null(gen)) { all_res[[idx]] <- rec; gc(); next }

  rec$t_gen  <- (proc.time() - t_gen)[["elapsed"]]
  rec$N_obs  <- gen$N_obs

  # CSV 저장 (id 기준 이미 정렬돼 있음)
  tmp_csv <- tempfile(fileext = ".csv")
  t_write <- proc.time()
  data.table::fwrite(gen$df, tmp_csv)
  t_write_elapsed <- (proc.time() - t_write)[["elapsed"]]
  rec$csv_mb <- file.size(tmp_csv) / 1e6

  cat(sprintf("  데이터 생성: %.1f초 | CSV: %.0f MB | 쓰기: %.1f초\n",
              rec$t_gen, rec$csv_mb, t_write_elapsed))
  cat(sprintf("  관측치: %s\n", format(rec$N_obs, big.mark = ",")))

  # df 즉시 삭제 → 이후 otwfe_file은 CSV만 사용
  rm(gen); invisible(gc())

  # ------------------------------------------------------------------
  # 2. otwfe_file() 실행
  # ------------------------------------------------------------------
  cat("  [2/3] otwfe_file() 실행 중...\n")
  t_file <- proc.time()
  result <- tryCatch(
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
      rec$fail_msg <<- conditionMessage(e)
      NULL
    }
  )
  rec$t_file <- (proc.time() - t_file)[["elapsed"]]

  # CSV 삭제
  unlink(tmp_csv)

  if (is.null(result)) { all_res[[idx]] <- rec; gc(); next }

  # ------------------------------------------------------------------
  # 3. 결과 기록
  # ------------------------------------------------------------------
  S <- result$state
  rec$state_mb <- as.numeric(object.size(S)) / 1e6
  rec$theta_x1 <- S$theta_hat["x1"]
  rec$theta_x2 <- S$theta_hat["x2"]
  V0            <- S$sigma2_hat * S$inv_dotZtZ
  rec$se_cl_x1 <- sqrt(V0["x1", "x1"])
  rec$se_cl_x2 <- sqrt(V0["x2", "x2"])
  rec$se_cr_x1 <- sqrt(S$Vcr_hat["x1", "x1"])
  rec$se_cr_x2 <- sqrt(S$Vcr_hat["x2", "x2"])
  rec$bias_x1  <- rec$theta_x1 - BETA_TRUE["x1"]
  rec$bias_x2  <- rec$theta_x2 - BETA_TRUE["x2"]
  rec$ok       <- TRUE

  cat(sprintf("\n  [3/3] 결과 요약\n"))
  cat(sprintf("    otwfe_file 소요: %.1f초  |  state: %.4f MB\n",
              rec$t_file, rec$state_mb))
  cat(sprintf("    x1:  theta=%.4f (bias=%+.4f)  SE=%.4f  SE(CR)=%.4f\n",
              rec$theta_x1, rec$bias_x1, rec$se_cl_x1, rec$se_cr_x1))
  cat(sprintf("    x2:  theta=%.4f (bias=%+.4f)  SE=%.4f  SE(CR)=%.4f\n",
              rec$theta_x2, rec$bias_x2, rec$se_cl_x2, rec$se_cr_x2))

  all_res[[idx]] <- rec
  rm(result); invisible(gc())
}

# =============================================================================
# 최종 요약 표
# =============================================================================
cat(sprintf("\n\n%s\n  최종 요약  (T=%d, k=2, beta_true: x1=1.3, x2=-0.7)\n%s\n",
            sep72, T_PANEL, sep72))
cat(sprintf("  %-11s  %-10s  %-7s  %-7s  %-8s  %-8s  %-8s  %-6s  %s\n",
            "N", "관측치", "CSV(MB)", "시간(s)",
            "θ_x1", "θ_x2", "bias_x1", "state", "상태"))
cat(sprintf("%s\n", strrep("-", 90)))

for (r in all_res) {
  if (is.null(r)) next
  status <- if (r$ok) "OK" else sprintf("FAIL: %s", substr(r$fail_msg %||% "?", 1, 30))
  cat(sprintf("  %-11s  %-10s  %7.0f  %7.1f  %8.4f  %8.4f  %+8.4f  %.4f  %s\n",
              format(r$N, big.mark = ","),
              format(r$N_obs %||% NA, big.mark = ","),
              r$csv_mb %||% NA,
              r$t_file %||% NA,
              r$theta_x1 %||% NA,
              r$theta_x2 %||% NA,
              r$bias_x1  %||% NA,
              r$state_mb %||% NA,
              status))
}

cat(sprintf("%s\n", sep72))
cat(sprintf("  chunk_size: %s행 / 청크  |  Rcpp 배치: %s\n",
            format(CHUNK_SIZE, big.mark = ","),
            if (.alg1_batch_rcpp_available) "사용됨" else "미사용"))
cat(sprintf("%s\n", sep72))
