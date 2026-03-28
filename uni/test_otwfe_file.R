# =============================================================================
# test_otwfe_file.R
# otwfe_file() 검증: N=2000, T=3, k=2
# plm과 theta / V0 / Vcr 비교 → machine precision 수준 일치 확인
# =============================================================================

setwd("/Users/joonhohwang/Desktop/claude_code/uni")
source("otwfe_file.R")
source("../Simulation/sim_dgp.R")
suppressPackageStartupMessages(library(plm))
suppressPackageStartupMessages(library(data.table))

cat("=== otwfe_file() 검증 (N=2000, T=3, k=2) ===\n\n")

# ------------------------------------------------------------------
# 1. 합성 데이터 생성 → CSV 저장
# ------------------------------------------------------------------
set.seed(99L)
gen    <- generate_panel_52(N = 2000L, T = 3L, T_min = 2L,
                             seed = 99L, verbose = FALSE)
df     <- gen$df
x_cols <- c("x1", "x2")

# id 순 정렬 (generate_panel_52는 이미 id 순)
df <- df[order(df$id), ]

# 임시 CSV 저장
tmp_csv <- tempfile(fileext = ".csv")
data.table::fwrite(df, tmp_csv)
cat(sprintf("CSV 저장: %s  (%.1f MB, %s행)\n\n",
            basename(tmp_csv),
            file.size(tmp_csv) / 1e6,
            format(nrow(df), big.mark = ",")))

# ------------------------------------------------------------------
# 2. otwfe_file() 실행
# ------------------------------------------------------------------
cat("--- otwfe_file() ---\n")
fit_file <- otwfe_file(
  path       = tmp_csv,
  id_col     = "id",
  time_col   = "time",
  y_col      = "y",
  x_cols     = x_cols,
  chunk_size = 500L,    # 작은 chunk로 경계 처리 검증
  verbose    = TRUE
)

# ------------------------------------------------------------------
# 3. plm 실행 (offline baseline)
# ------------------------------------------------------------------
cat("\n--- plm() baseline ---\n")
pdf     <- pdata.frame(df, index = c("id", "time"))
fml     <- y ~ x1 + x2 + factor(time)
plm_fit <- plm(fml, data = pdf, model = "within", effect = "individual")
theta_plm <- coef(plm_fit)[x_cols]
V0_plm    <- vcov(plm_fit)[x_cols, x_cols]
Vcr_plm   <- vcovHC(plm_fit, method = "arellano", type = "HC0")[x_cols, x_cols]

# ------------------------------------------------------------------
# 4. 정확도 비교
# ------------------------------------------------------------------
theta_file <- fit_file$state$theta_hat[x_cols]
V0_file    <- fit_file$state$sigma2_hat * fit_file$state$inv_dotZtZ[x_cols, x_cols]
Vcr_file   <- fit_file$state$Vcr_hat[x_cols, x_cols]

theta_diff <- max(abs(theta_file - theta_plm))
V0_diff    <- max(abs(V0_file    - V0_plm))
Vcr_diff   <- max(abs(Vcr_file   - Vcr_plm))

cat("\n=== 정확도 비교 ===\n")
cat(sprintf("  theta max|diff| : %.3e  %s\n", theta_diff,
            if (theta_diff < 1e-10) "[OK]" else "[FAIL]"))
cat(sprintf("  V0    max|diff| : %.3e  %s\n", V0_diff,
            if (V0_diff    < 1e-10) "[OK]" else "[FAIL]"))
cat(sprintf("  Vcr   max|diff| : %.3e  %s\n", Vcr_diff,
            if (Vcr_diff   < 1e-10) "[OK]" else "[FAIL]"))

# ------------------------------------------------------------------
# 5. 임의의 time 값(연도 형식) 검증 — T_support 탐지, time_remap 확인
# ------------------------------------------------------------------
cat("\n=== time 값 = 연도 형식 검증 (2020, 2021, 2022) ===\n")
df2          <- df
df2$time     <- df2$time + 2019L  # {1,2,3} → {2020,2021,2022}
tmp_csv2     <- tempfile(fileext = ".csv")
data.table::fwrite(df2, tmp_csv2)

fit_yr <- otwfe_file(
  path = tmp_csv2, id_col = "id", time_col = "time",
  y_col = "y", x_cols = x_cols, chunk_size = 500L, verbose = FALSE
)
theta_yr_diff <- max(abs(fit_yr$state$theta_hat[x_cols] - theta_plm))
cat(sprintf("  theta max|diff| vs plm: %.3e  %s\n", theta_yr_diff,
            if (theta_yr_diff < 1e-10) "[OK]" else "[FAIL]"))
cat(sprintf("  탐지된 time_levels: %s\n",
            paste(fit_yr$time_levels_original, collapse = ", ")))

# 정리
unlink(tmp_csv); unlink(tmp_csv2)
cat("\n[검증 완료]\n")
