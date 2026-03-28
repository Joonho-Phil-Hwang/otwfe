# =============================================================================
# test_large.R
# N=1,000,000, T_max=5, k=5 합성 패널 데이터로 otwfe() 검증 및 벤치마크
# 검증 대상: plm 패키지 offline refit과의 비교
#   - coefficients (theta_hat)
#   - classical variance (sigma2 * (dotZ'dotZ)^{-1})
#   - cluster-robust variance (Arellano HC0)
# =============================================================================

setwd("/Users/joonhohwang/Desktop/claude_code/Simulation")
source("../online_twfe_core.R")
library(plm)

# -----------------------------------------------------------------------------
# 1. 파라미터 설정
# -----------------------------------------------------------------------------
set.seed(42)
N      <- 1000000L  # 개인 수
T_max  <- 5L        # 최대 time period
k      <- 5L        # 공변량 수 (x1 ~ x5)

x_cols    <- paste0("x", seq_len(k))
beta_true <- seq(0.5, by = 0.3, length.out = k)   # 0.5, 0.8, 1.1, 1.4, 1.7
lambda_t  <- c(0, 0.5, -0.3, 0.8, -0.2)           # 시간 고정효과 (t=1..5)
sigma_eps <- 0.5

# -----------------------------------------------------------------------------
# 2. 합성 데이터 생성
# -----------------------------------------------------------------------------
cat("=== [1/4] 합성 패널 데이터 생성 중 ===\n")

T_i_vec  <- sample(2L:T_max, N, replace = TRUE)   # unbalanced (T_i ∈ {2,...,5})
total_obs <- sum(T_i_vec)
cat(sprintf("  N = %s,  총 관측치 = %s\n",
            format(N, big.mark = ","), format(total_obs, big.mark = ",")))

alpha_i <- rnorm(N)

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

data_mb <- object.size(df) / 1e6
cat(sprintf("  데이터 크기: %.1f MB\n\n", data_mb))

# -----------------------------------------------------------------------------
# 3. otwfe() 실행
# -----------------------------------------------------------------------------
cat("=== [2/4] otwfe() 실행 중 ===\n")

t_start <- proc.time()
fit <- otwfe(
  data            = df,
  id_col          = "id",
  time_col        = "time",
  y_col           = "y",
  x_cols          = x_cols,
  track_all_units = FALSE,
  verbose         = TRUE
)
t_otwfe <- (proc.time() - t_start)[["elapsed"]]

global_state_mb <- object.size(fit$state) / 1e6
total_state_mb  <- sum_object_size(fit) / 1e6

cat(sprintf("\n  실행 시간:                 %.2f 초\n",   t_otwfe))
cat(sprintf("  Global state 크기:         %.4f MB\n",    global_state_mb))
cat(sprintf("  전체 state 크기:           %.4f MB\n",    total_state_mb))
cat(sprintf("  저장된 unit 수:            %d 개 (warm-up only)\n", length(fit$units)))
cat(sprintf("  데이터 대비 전체 state:    %.1f%% 절감\n\n",
            (1 - total_state_mb / data_mb) * 100))

# -----------------------------------------------------------------------------
# 4. plm offline refit (검증 기준)
#    N=1M에서 plm 자체는 빠르지만 vcovHC는 오래 걸릴 수 있음
# -----------------------------------------------------------------------------
cat("=== [3/4] plm offline refit 실행 중 ===\n")

t_start  <- proc.time()
pdf      <- pdata.frame(df, index = c("id", "time"))
plm_fit  <- plm(as.formula(paste("y ~", paste(x_cols, collapse = "+"),
                                 "+ factor(time)")),
                data = pdf, model = "within", effect = "individual")
t_plm_coef <- (proc.time() - t_start)[["elapsed"]]
cat(sprintf("  plm 계수 추정 시간: %.2f 초\n", t_plm_coef))

t_start  <- proc.time()
V0_pl    <- vcov(plm_fit)[x_cols, x_cols]
t_plm_V0 <- (proc.time() - t_start)[["elapsed"]]
cat(sprintf("  vcov() 시간:        %.2f 초\n", t_plm_V0))

t_start   <- proc.time()
Vcr_pl    <- vcovHC(plm_fit, method = "arellano",
                    type = "HC0")[x_cols, x_cols]
t_plm_vcr <- (proc.time() - t_start)[["elapsed"]]
cat(sprintf("  vcovHC() 시간:      %.2f 초\n\n", t_plm_vcr))

# -----------------------------------------------------------------------------
# 5. 비교
# -----------------------------------------------------------------------------
cat("=== [4/4] otwfe vs plm 비교 ===\n\n")

theta_on <- fit$state$theta_hat[x_cols]
theta_pl <- coef(plm_fit)[x_cols]
V0_on    <- fit$state$sigma2_hat * fit$state$inv_dotZtZ[x_cols, x_cols]
Vcr_on   <- fit$state$Vcr_hat[x_cols, x_cols]

cat("--- Coefficients ---\n")
comp_theta <- data.frame(
  true  = beta_true,
  otwfe = round(theta_on, 8),
  plm   = round(theta_pl, 8),
  diff  = theta_on - theta_pl,
  row.names = x_cols
)
print(comp_theta)
cat(sprintf("max|diff|: %.3e\n\n", max(abs(theta_on - theta_pl))))

cat("--- Classical SE (sqrt of diagonal) ---\n")
comp_V0 <- data.frame(
  otwfe = round(sqrt(diag(V0_on)), 8),
  plm   = round(sqrt(diag(V0_pl)), 8),
  diff  = sqrt(diag(V0_on)) - sqrt(diag(V0_pl)),
  row.names = x_cols
)
print(comp_V0)
cat(sprintf("max|diff| (variance): %.3e\n\n", max(abs(V0_on - V0_pl))))

cat("--- Cluster-Robust SE (sqrt of diagonal, Arellano HC0) ---\n")
comp_Vcr <- data.frame(
  otwfe = round(sqrt(diag(Vcr_on)), 8),
  plm   = round(sqrt(diag(Vcr_pl)), 8),
  diff  = sqrt(diag(Vcr_on)) - sqrt(diag(Vcr_pl)),
  row.names = x_cols
)
print(comp_Vcr)
cat(sprintf("max|diff| (variance): %.3e\n\n", max(abs(Vcr_on - Vcr_pl))))

cat("--- 종합 ---\n")
cat(sprintf("  theta  max|diff|: %.3e\n", max(abs(theta_on - theta_pl))))
cat(sprintf("  V0     max|diff|: %.3e\n", max(abs(V0_on    - V0_pl))))
cat(sprintf("  Vcr    max|diff|: %.3e\n", max(abs(Vcr_on   - Vcr_pl))))
cat(sprintf("  otwfe %.2f초  vs  plm %.2f초 (계수+분산 합계)\n",
            t_otwfe, t_plm_coef + t_plm_V0 + t_plm_vcr))
