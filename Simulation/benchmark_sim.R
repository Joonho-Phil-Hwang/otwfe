# =============================================================================
# benchmark_sim.R
#
# Simulation: otwfe_file() vs plm
#   - Speed  : time each step separately
#   - Memory : process RSS (ps package) + otwfe state size
#   - Accuracy: coef / classical V / cluster-robust V vs plm (machine epsilon)
#
# DGP: Hwang & Lee (2026) §5.2  |  T=5, k=2, p=6
# =============================================================================

setwd("/Users/joonhohwang/Desktop/claude_code")

suppressPackageStartupMessages(library(devtools))
load_all(".", quiet = TRUE)
source("Simulation/sim_dgp.R")

suppressPackageStartupMessages({
  library(data.table)
  library(plm)
  library(ps)
})

# -----------------------------------------------------------------------------
# Helper: current process RSS in MB
# -----------------------------------------------------------------------------
rss_mb <- function() {
  tryCatch(ps::ps_memory_info()$rss / 1e6, error = function(e) NA_real_)
}

# -----------------------------------------------------------------------------
# Helper: run one timed block, return list(result, elapsed, err_msg)
# -----------------------------------------------------------------------------
timed <- function(expr) {
  t0      <- proc.time()
  err_msg <- NULL
  result  <- tryCatch(force(expr),
                      error = function(e) { err_msg <<- conditionMessage(e); NULL })
  el      <- (proc.time() - t0)[["elapsed"]]
  list(result = result, elapsed = el, err = err_msg)
}

# -----------------------------------------------------------------------------
# Simulation parameters
# -----------------------------------------------------------------------------
N_vec      <- c(1e6, 5e6, 1e7, 2e7, 5e7, 1e8)
T_val      <- 5L
k_val      <- 2L
chunk_size <- 2e6L
seed       <- 42L

cat(sprintf("\n%s\n", strrep("=", 72)))
cat(sprintf("  Benchmark: otwfe_file() vs plm\n"))
cat(sprintf("  T=%d, k=%d, p=%d  |  chunk_size=%s\n",
            T_val, k_val, T_val + k_val - 1L,
            format(chunk_size, big.mark = ",")))
cat(sprintf("  N: %s\n",
            paste(format(N_vec, big.mark = ",", scientific = FALSE), collapse = " / ")))
cat(sprintf("%s\n\n", strrep("=", 72)))

# -----------------------------------------------------------------------------
# Storage for results
# -----------------------------------------------------------------------------
results <- vector("list", length(N_vec))

# =============================================================================
# Main loop
# =============================================================================
for (idx in seq_along(N_vec)) {

  N <- as.integer(N_vec[idx])
  cat(sprintf("\n%s\n", strrep("-", 72)))
  cat(sprintf("  N = %s\n", format(N, big.mark = ",")))
  cat(sprintf("%s\n", strrep("-", 72)))

  # ------------------------------------------------------------------
  # Data generation
  # ------------------------------------------------------------------
  cat("  Generating data...\n")
  t_gen <- proc.time()
  gen   <- generate_panel_52(N = N, T = T_val, T_min = 2L,
                              seed = seed, verbose = FALSE)
  t_gen <- (proc.time() - t_gen)[["elapsed"]]

  tmp_csv <- file.path(getwd(), sprintf("Simulation/_bench_tmp_N%d.csv", N))
  fwrite(gen$df, tmp_csv)
  csv_mb  <- file.size(tmp_csv) / 1e6
  n_obs   <- gen$N_obs
  cat(sprintf("  Generated: %s obs  |  CSV: %.0f MB  |  %.1f sec\n",
              format(n_obs, big.mark = ","), csv_mb, t_gen))
  rm(gen); invisible(gc())

  row <- list(
    N = N, n_obs = n_obs, csv_mb = csv_mb,
    # otwfe timings
    ot_total = NA_real_, ot_state_mb = NA_real_,
    ot_rss_before = NA_real_, ot_rss_after = NA_real_,
    # plm timings
    plm_pdata  = NA_real_, plm_fit    = NA_real_,
    plm_vcov   = NA_real_, plm_vcovhc = NA_real_,
    plm_rss_before = NA_real_, plm_rss_after = NA_real_,
    # plm error info  (NULL = no failure yet)
    plm_fail_step = NULL, plm_fail_msg = NULL,
    # accuracy
    diff_coef = NA_real_, diff_V0 = NA_real_, diff_Vcr = NA_real_
  )

  # ------------------------------------------------------------------
  # [otwfe_file]
  # ------------------------------------------------------------------
  cat("  [otwfe_file] running...\n")
  row$ot_rss_before <- rss_mb()

  t_ot <- proc.time()
  fit_ot <- otwfe_file(
    path       = tmp_csv,
    id_col     = "id", time_col = "time",
    y_col      = "y", x_cols   = c("x1", "x2"),
    chunk_size = chunk_size,
    verbose    = FALSE
  )
  row$ot_total   <- (proc.time() - t_ot)[["elapsed"]]
  row$ot_rss_after  <- rss_mb()
  row$ot_state_mb   <- as.numeric(object.size(fit_ot$state)) / 1e6

  cat(sprintf("  [otwfe_file] done: %.1f sec  |  state: %.4f MB  |  RSS delta: %.0f MB\n",
              row$ot_total, row$ot_state_mb,
              row$ot_rss_after - row$ot_rss_before))

  # ------------------------------------------------------------------
  # [plm]  — each step separately, with tryCatch
  # ------------------------------------------------------------------
  cat("  [plm] running...\n")
  row$plm_rss_before <- rss_mb()

  plm_pdf <- NULL; plm_fit <- NULL; plm_coef_saved <- NULL
  plm_V0  <- NULL; plm_Vcr <- NULL

  # Step 1: pdata.frame
  cat("    pdata.frame()...")
  r1 <- timed({
    df_tmp <- fread(tmp_csv)
    pdata.frame(as.data.frame(df_tmp), index = c("id", "time"))
  })
  row$plm_pdata <- r1$elapsed
  cat(sprintf(" %.1f sec\n", r1$elapsed))
  if (!is.null(r1$err)) {
    row$plm_fail_step <- "pdata.frame"
    row$plm_fail_msg  <- r1$err
    cat(sprintf("    FAIL: %s\n", r1$err))
  } else {
    plm_pdf <- r1$result
  }

  # Step 2: plm()
  if (is.null(row$plm_fail_step)) {
    cat("    plm()...")
    r2 <- timed(plm(y ~ x1 + x2 + factor(time), data = plm_pdf,
                    model = "within", effect = "individual"))
    row$plm_fit <- r2$elapsed
    cat(sprintf(" %.1f sec\n", r2$elapsed))
    if (!is.null(r2$err)) {
      row$plm_fail_step <- "plm"
      row$plm_fail_msg  <- r2$err
      cat(sprintf("    FAIL: %s\n", r2$err))
    } else {
      plm_fit          <- r2$result
      plm_coef_saved   <- coef(r2$result)[c("x1", "x2")]
    }
    rm(plm_pdf); invisible(gc())
  }

  # Step 3: vcov() (classical)
  if (is.null(row$plm_fail_step)) {
    cat("    vcov()...")
    r3 <- timed(vcov(plm_fit)[c("x1", "x2"), c("x1", "x2")])
    row$plm_vcov <- r3$elapsed
    cat(sprintf(" %.1f sec\n", r3$elapsed))
    if (!is.null(r3$err)) {
      row$plm_fail_step <- "vcov"
      row$plm_fail_msg  <- r3$err
      cat(sprintf("    FAIL: %s\n", r3$err))
    } else {
      plm_V0 <- r3$result
    }
  }

  # Step 4: vcovHC() (cluster-robust)
  if (is.null(row$plm_fail_step)) {
    cat("    vcovHC()...")
    r4 <- timed(vcovHC(plm_fit, method = "arellano",
                       type = "HC0")[c("x1","x2"), c("x1","x2")])
    row$plm_vcovhc <- r4$elapsed
    cat(sprintf(" %.1f sec\n", r4$elapsed))
    if (!is.null(r4$err)) {
      row$plm_fail_step <- "vcovHC"
      row$plm_fail_msg  <- r4$err
      cat(sprintf("    FAIL: %s\n", r4$err))
    } else {
      plm_Vcr <- r4$result
    }
    rm(plm_fit); invisible(gc())
  }

  row$plm_rss_after <- rss_mb()

  # ------------------------------------------------------------------
  # Accuracy check (only when all plm steps succeeded)
  # ------------------------------------------------------------------
  if (is.null(row$plm_fail_step) &&
      !is.null(plm_coef_saved) && !is.null(plm_V0) && !is.null(plm_Vcr)) {
    row$diff_coef <- max(abs(coef(fit_ot)                    - plm_coef_saved))
    row$diff_V0   <- max(abs(vcov(fit_ot, type = "classical") - plm_V0))
    row$diff_Vcr  <- max(abs(vcov(fit_ot, type = "robust")   - plm_Vcr))

    cat(sprintf("  [accuracy]  coef: %.2e  V0: %.2e  Vcr: %.2e\n",
                row$diff_coef, row$diff_V0, row$diff_Vcr))
  }

  unlink(tmp_csv)
  rm(fit_ot); invisible(gc())
  results[[idx]] <- row
}

# =============================================================================
# Summary table
# =============================================================================
cat(sprintf("\n\n%s\n", strrep("=", 100)))
cat(sprintf("  SUMMARY  (T=%d, k=%d)\n", T_val, k_val))
cat(sprintf("%s\n", strrep("=", 100)))

hdr <- sprintf("  %-12s  %-12s  %8s  %10s  %8s  %8s  %8s  %8s  %s",
               "N", "obs", "ot(s)", "state(MB)", "pdata(s)", "plm(s)",
               "vcov(s)", "vcovHC(s)", "status")
cat(hdr, "\n")
cat(sprintf("  %s\n", strrep("-", 96)))

for (row in results) {
  status <- if (is.null(row$plm_fail_step)) {
    sprintf("OK | coef=%.1e  V0=%.1e  Vcr=%.1e",
            row$diff_coef, row$diff_V0, row$diff_Vcr)
  } else {
    sprintf("FAIL @ %s: %s",
            row$plm_fail_step,
            substr(row$plm_fail_msg, 1, 60))
  }

  cat(sprintf("  %-12s  %-12s  %8.1f  %10.4f  %8s  %8s  %8s  %8s  %s\n",
              format(row$N,     big.mark = ",", scientific = FALSE),
              format(row$n_obs, big.mark = ",", scientific = FALSE),
              row$ot_total,
              row$ot_state_mb,
              ifelse(is.na(row$plm_pdata),  "  -", sprintf("%5.1f", row$plm_pdata)),
              ifelse(is.na(row$plm_fit),    "  -", sprintf("%5.1f", row$plm_fit)),
              ifelse(is.na(row$plm_vcov),   "  -", sprintf("%5.1f", row$plm_vcov)),
              ifelse(is.na(row$plm_vcovhc), "  -", sprintf("%5.1f", row$plm_vcovhc)),
              status))

  if (!is.null(row$plm_fail_step))
    cat(sprintf("    └─ FAIL msg: %s\n", row$plm_fail_msg))
}

cat(sprintf("  %s\n", strrep("=", 100)))
cat(sprintf("  chunk_size: %s rows\n", format(chunk_size, big.mark = ",")))
