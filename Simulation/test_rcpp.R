# =============================================================================
# test_rcpp.R
# Rcpp 배치 구현 검증: N=200, T_max=5, k=3 unbalanced panel
# otwfe(Rcpp) vs plm 비교 — machine precision 수준 일치 확인
# =============================================================================

setwd("/Users/joonhohwang/Desktop/claude_code/Simulation")
source("../online_twfe_core.R")
library(plm)

set.seed(123)
N      <- 200L
T_max  <- 5L
k      <- 3L
x_cols    <- paste0("x", seq_len(k))
beta_true <- c(1.0, -0.5, 0.8)
lambda_t  <- c(0, 0.4, -0.2, 0.6, -0.1)
sigma_eps <- 0.5

T_i_vec  <- sample(2L:T_max, N, replace = TRUE)
alpha_i  <- rnorm(N)
rows <- vector("list", N)
for (i in seq_len(N)) {
  Ti    <- T_i_vec[i]
  times <- seq_len(Ti)
  X     <- matrix(rnorm(Ti * k), nrow = Ti, ncol = k,
                  dimnames = list(NULL, x_cols))
  eps   <- rnorm(Ti, sd = sigma_eps)
  y     <- alpha_i[i] + lambda_t[times] + drop(X %*% beta_true) + eps
  rows[[i]] <- as.data.frame(cbind(id = i, time = times, y = y, X))
}
df <- do.call(rbind, rows)

# ------------------------------------------------------------------
# 1. otwfe() — Rcpp 배치 사용
# ------------------------------------------------------------------
cat("=== otwfe() 실행 ===\n")
fit <- otwfe(
  data            = df,
  id_col          = "id",
  time_col        = "time",
  y_col           = "y",
  x_cols          = x_cols,
  track_all_units = FALSE,
  verbose         = TRUE
)
cat(sprintf("Rcpp 배치 모듈 사용: %s\n\n", .alg1_batch_rcpp_available))

# ------------------------------------------------------------------
# 2. plm offline benchmark
# ------------------------------------------------------------------
cat("=== plm offline refit ===\n")
pdf      <- pdata.frame(df, index = c("id", "time"))
plm_fit  <- plm(as.formula(paste("y ~", paste(x_cols, collapse = "+"),
                                 "+ factor(time)")),
                data = pdf, model = "within", effect = "individual")
V0_pl    <- vcov(plm_fit)[x_cols, x_cols]
Vcr_pl   <- vcovHC(plm_fit, method = "arellano", type = "HC0")[x_cols, x_cols]

# ------------------------------------------------------------------
# 3. 비교
# ------------------------------------------------------------------
theta_on <- fit$state$theta_hat[x_cols]
theta_pl <- coef(plm_fit)[x_cols]
V0_on    <- fit$state$sigma2_hat * fit$state$inv_dotZtZ[x_cols, x_cols]
Vcr_on   <- fit$state$Vcr_hat[x_cols, x_cols]

cat("--- Coefficients ---\n")
print(data.frame(
  true  = beta_true,
  otwfe = round(theta_on, 8),
  plm   = round(theta_pl, 8),
  diff  = theta_on - theta_pl,
  row.names = x_cols
))
cat(sprintf("theta  max|diff|: %.3e\n\n", max(abs(theta_on - theta_pl))))

cat("--- Classical SE ---\n")
print(data.frame(
  otwfe = round(sqrt(diag(V0_on)), 8),
  plm   = round(sqrt(diag(V0_pl)), 8),
  diff  = sqrt(diag(V0_on)) - sqrt(diag(V0_pl)),
  row.names = x_cols
))
cat(sprintf("V0     max|diff|: %.3e\n\n", max(abs(V0_on - V0_pl))))

cat("--- Cluster-Robust SE (Arellano HC0) ---\n")
print(data.frame(
  otwfe = round(sqrt(diag(Vcr_on)), 8),
  plm   = round(sqrt(diag(Vcr_pl)), 8),
  diff  = sqrt(diag(Vcr_on)) - sqrt(diag(Vcr_pl)),
  row.names = x_cols
))
cat(sprintf("Vcr    max|diff|: %.3e\n\n", max(abs(Vcr_on - Vcr_pl))))

# 합격 기준: machine precision 수준
tol <- 1e-10
pass_theta <- max(abs(theta_on - theta_pl))  < tol
pass_V0    <- max(abs(V0_on    - V0_pl))     < tol
pass_Vcr   <- max(abs(Vcr_on   - Vcr_pl))   < tol

cat("=== 검증 결과 ===\n")
cat(sprintf("  theta : %s (max|diff| = %.3e)\n",
            if (pass_theta) "PASS" else "FAIL", max(abs(theta_on - theta_pl))))
cat(sprintf("  V0    : %s (max|diff| = %.3e)\n",
            if (pass_V0)    "PASS" else "FAIL", max(abs(V0_on - V0_pl))))
cat(sprintf("  Vcr   : %s (max|diff| = %.3e)\n",
            if (pass_Vcr)   "PASS" else "FAIL", max(abs(Vcr_on - Vcr_pl))))

if (pass_theta && pass_V0 && pass_Vcr) {
  cat("\n[OK] 모든 검증 통과 — Rcpp 배치 구현 정확성 확인\n")
} else {
  cat("\n[FAIL] 검증 실패\n")
  quit(status = 1)
}
