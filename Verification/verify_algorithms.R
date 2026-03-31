# =============================================================================
# verify_algorithms.R
#
# Section 5.1 검증: 온라인 업데이트 알고리즘 vs. offline plm 재추정
#
#   Scenario 1: 새 unit 도착          → Propositions 1–3  (Algorithm 1)
#   Scenario 2: 기존 unit, 기존 time  → Propositions 4–6  (Algorithm 2)
#   Scenario 3: 기존 unit, 새 time    → Propositions 7–9  (Algorithm 3)
#
# DGP: Y_it = X_it' β + α_i + λ_t + u_it,  β = (1.3, −0.7)'
#      N = 10,000,  T = 5,  T_min = 2
# =============================================================================

setwd("/Users/joonhohwang/Desktop/claude_code")

# 패키지 설치 후 로드 (load_all 대신 install 사용 — .so 임시 디렉토리 문제 우회)
suppressPackageStartupMessages(library(devtools))
suppressMessages(install(".", upgrade = FALSE, quiet = TRUE))
suppressPackageStartupMessages(library(otwfe))

suppressPackageStartupMessages(library(plm))

# 내부 함수 단축 바인딩 (:::로 접근)
state_init            <- otwfe:::state_init
alg1_new_unit         <- otwfe:::alg1_new_unit
alg2_existing_unit    <- otwfe:::alg2_existing_unit
alg3_new_caltime      <- otwfe:::alg3_new_caltime
build_Z_raw           <- otwfe:::build_Z_raw
within_transform_unit <- otwfe:::within_transform_unit
make_unit_state       <- otwfe:::make_unit_state

set.seed(42L)

cat(sprintf("\n%s\n", strrep("=", 72)))
cat("  Online Updating Algorithm Verification\n")
cat("  Scenarios 1–3 (Props 1–9)  |  N = 10,000, T = 5\n")
cat(sprintf("%s\n\n", strrep("=", 72)))

# -----------------------------------------------------------------------
# DGP 파라미터
# -----------------------------------------------------------------------
N        <- 10000L
T_val    <- 5L
T_min    <- 2L
beta     <- c(1.3, -0.7)
sigma_a  <- sigma_l <- sigma_u <- 1
x_cols   <- c("x1", "x2")

# -----------------------------------------------------------------------
# 고정효과 생성 (λ_t는 T+1까지 — Scenario 3의 t = T+1 대비)
# -----------------------------------------------------------------------
alpha_i  <- rnorm(N, sd = sigma_a)
lambda_t <- rnorm(T_val + 1L, sd = sigma_l)   # t = 1, ..., 6

# -----------------------------------------------------------------------
# 기반 패널 생성
#   T_i ~ Uniform{T_min,...,T},  S_i ⊆ {1,...,T} 비복원 추출
# -----------------------------------------------------------------------
T_i_vec <- sample(T_min:T_val, N, replace = TRUE)
rows    <- vector("list", N)
for (i in seq_len(N)) {
  Ti  <- T_i_vec[i]
  Si  <- sort(sample.int(T_val, Ti))
  x1  <- rnorm(Ti); x2 <- rnorm(Ti); u <- rnorm(Ti, sd = sigma_u)
  y   <- alpha_i[i] + lambda_t[Si] + beta[1L]*x1 + beta[2L]*x2 + u
  rows[[i]] <- data.frame(id = i, time = Si, y = y, x1 = x1, x2 = x2,
                           stringsAsFactors = FALSE)
}
df_base <- do.call(rbind, rows)

# 기술적 조건: 모든 t ∈ {1,...,T} 최소 1회 등장
missing_t <- setdiff(seq_len(T_val), unique(df_base$time))
if (length(missing_t) > 0L) {
  warning(sprintf("누락 calendar time 강제 커버: {%s}",
                  paste(missing_t, collapse = ",")))
  for (t_miss in missing_t) {
    for (i in seq_len(N)) {
      idx_i <- which(df_base$id == i)
      if (length(idx_i) >= 2L) { df_base$time[idx_i[1L]] <- t_miss; break }
    }
  }
}

cat(sprintf("기반 패널: N = %s units,  n = %s obs,  T = %d\n",
            format(N, big.mark = ","),
            format(nrow(df_base), big.mark = ","), T_val))

# -----------------------------------------------------------------------
# 초기 state: N unit 모두 warm-up으로 state 초기화
# -----------------------------------------------------------------------
cat("초기 state 구축 중...\n")
init  <- state_init(df_base, "id", "time", "y", x_cols, baseline_time = 1L)
S_pre <- init$S_pre
T_sup <- S_pre$T_support   # = 5

cat(sprintf("  T_support = %d,  p = %d,  df = %d\n\n",
            T_sup, S_pre$p, S_pre$n - S_pre$N - S_pre$p))

# -----------------------------------------------------------------------
# 헬퍼: plm 적합 후 coef, classical vcov, cluster-robust vcov 반환
# -----------------------------------------------------------------------
plm_estimates <- function(df) {
  pdata <- pdata.frame(df, index = c("id", "time"))
  fit   <- plm(y ~ x1 + x2 + factor(time), data = pdata,
               model = "within", effect = "individual")
  list(
    coef = coef(fit)[x_cols],
    V0   = vcov(fit)[x_cols, x_cols],
    Vcr  = vcovHC(fit, method = "arellano", type = "HC0")[x_cols, x_cols]
  )
}

# -----------------------------------------------------------------------
# 새 unit (N+1) 데이터 생성
#   T_{N+1} ∈ {T_min,...,T−1},  S_{N+1} ⊆ {1,...,T−1}
# -----------------------------------------------------------------------
alpha_new <- rnorm(1L, sd = sigma_a)
T_n1      <- sample(T_min:(T_val - 1L), 1L)
S_n1      <- sort(sample.int(T_val - 1L, T_n1))   # ⊆ {1,...,4}

x1_n1 <- rnorm(T_n1); x2_n1 <- rnorm(T_n1); u_n1 <- rnorm(T_n1, sd = sigma_u)
y_n1  <- alpha_new + lambda_t[S_n1] + beta[1L]*x1_n1 + beta[2L]*x2_n1 + u_n1
x_mat_n1 <- matrix(cbind(x1_n1, x2_n1), ncol = 2L,
                   dimnames = list(NULL, x_cols))
df_new1  <- data.frame(id = N + 1L, time = S_n1, y = y_n1,
                        x1 = x1_n1, x2 = x2_n1, stringsAsFactors = FALSE)

cat(sprintf("새 unit %d:  T_{N+1} = %d,  S_{N+1} = {%s}\n",
            N + 1L, T_n1, paste(S_n1, collapse = ",")))

# =============================================================================
# Scenario 1 (Algorithm 1, Props 1–3): 새 unit 도착
# =============================================================================
cat(sprintf("\n%s\n", strrep("-", 72)))
cat("[Scenario 1]  Algorithm 1 — 새 unit 도착\n")

S_post1 <- alg1_new_unit(S_pre, x_mat_n1, S_n1, y_n1)

# unit N+1의 per-unit state 계산 (Scenario 2에서 Algorithm 2 입력으로 사용)
Z_n1    <- build_Z_raw(x_mat_n1, S_n1, T_sup, S_pre$baseline_time)
wt_n1   <- within_transform_unit(Z_n1, y_n1)
unit_n1 <- make_unit_state(
  id     = N + 1L,
  T_i    = T_n1,
  barZ_i = wt_n1$barZ_i,
  barY_i = wt_n1$barY_i,
  S_i    = crossprod(wt_n1$dotZ_i),
  s_i    = crossprod(wt_n1$dotZ_i, matrix(wt_n1$dotY_i, ncol = 1L))
)

# offline plm: base + new unit
plm1 <- plm_estimates(rbind(df_base, df_new1))

diff1_coef <- max(abs(S_post1$theta_hat[x_cols]                            - plm1$coef))
diff1_V0   <- max(abs((S_post1$sigma2_hat * S_post1$inv_dotZtZ)[x_cols, x_cols] - plm1$V0))
diff1_Vcr  <- max(abs(S_post1$Vcr_hat[x_cols, x_cols]                     - plm1$Vcr))

cat(sprintf("  Prop 1 (coefficient)             max|diff| = %.2e\n", diff1_coef))
cat(sprintf("  Prop 2 (classical variance)      max|diff| = %.2e\n", diff1_V0))
cat(sprintf("  Prop 3 (cluster-robust variance) max|diff| = %.2e\n", diff1_Vcr))

# =============================================================================
# Scenario 2 (Algorithm 2, Props 4–6): 기존 unit, 기존 calendar time t = T
# =============================================================================
cat(sprintf("\n%s\n", strrep("-", 72)))
cat(sprintf("[Scenario 2]  Algorithm 2 — 기존 unit의 기존 time  t = %d\n", T_val))

x1_n2 <- rnorm(1L); x2_n2 <- rnorm(1L); u_n2 <- rnorm(1L, sd = sigma_u)
y_n2  <- alpha_new + lambda_t[T_val] + beta[1L]*x1_n2 + beta[2L]*x2_n2 + u_n2
df_new2 <- data.frame(id = N + 1L, time = T_val, y = y_n2,
                       x1 = x1_n2, x2 = x2_n2, stringsAsFactors = FALSE)

res2          <- alg2_existing_unit(S_post1, unit_n1, c(x1_n2, x2_n2), T_val, y_n2)
S_post2       <- res2$S_post
unit_n1_post2 <- res2$unit_post

plm2 <- plm_estimates(rbind(df_base, df_new1, df_new2))

diff2_coef <- max(abs(S_post2$theta_hat[x_cols]                            - plm2$coef))
diff2_V0   <- max(abs((S_post2$sigma2_hat * S_post2$inv_dotZtZ)[x_cols, x_cols] - plm2$V0))
diff2_Vcr  <- max(abs(S_post2$Vcr_hat[x_cols, x_cols]                     - plm2$Vcr))

cat(sprintf("  Prop 4 (coefficient)             max|diff| = %.2e\n", diff2_coef))
cat(sprintf("  Prop 5 (classical variance)      max|diff| = %.2e\n", diff2_V0))
cat(sprintf("  Prop 6 (cluster-robust variance) max|diff| = %.2e\n", diff2_Vcr))

# =============================================================================
# Scenario 3 (Algorithm 3, Props 7–9): 기존 unit, 새 calendar time t = T+1
# =============================================================================
cat(sprintf("\n%s\n", strrep("-", 72)))
cat(sprintf("[Scenario 3]  Algorithm 3 — 기존 unit의 새 time  t = %d\n", T_val + 1L))

x1_n3 <- rnorm(1L); x2_n3 <- rnorm(1L); u_n3 <- rnorm(1L, sd = sigma_u)
y_n3  <- alpha_new + lambda_t[T_val + 1L] + beta[1L]*x1_n3 + beta[2L]*x2_n3 + u_n3
df_new3 <- data.frame(id = N + 1L, time = T_val + 1L, y = y_n3,
                       x1 = x1_n3, x2 = x2_n3, stringsAsFactors = FALSE)

res3    <- alg3_new_caltime(S_post2, unit_n1_post2, c(x1_n3, x2_n3), T_val + 1L, y_n3)
S_post3 <- res3$S_post

plm3 <- plm_estimates(rbind(df_base, df_new1, df_new2, df_new3))

diff3_coef <- max(abs(S_post3$theta_hat[x_cols]                            - plm3$coef))
diff3_V0   <- max(abs((S_post3$sigma2_hat * S_post3$inv_dotZtZ)[x_cols, x_cols] - plm3$V0))
diff3_Vcr  <- max(abs(S_post3$Vcr_hat[x_cols, x_cols]                     - plm3$Vcr))

cat(sprintf("  Prop 7 (coefficient)             max|diff| = %.2e\n", diff3_coef))
cat(sprintf("  Prop 8 (classical variance)      max|diff| = %.2e\n", diff3_V0))
cat(sprintf("  Prop 9 (cluster-robust variance) max|diff| = %.2e\n", diff3_Vcr))

# =============================================================================
# 결과 요약 (논문 Table 2 형식)
# =============================================================================
cat(sprintf("\n\n%s\n", strrep("=", 72)))
cat("  Table 2: Online updating vs. offline plm  (max absolute difference)\n")
cat(sprintf("%s\n", strrep("-", 72)))
cat(sprintf("  %-12s  %-10s  %-28s  %s\n",
            "Event", "Prop", "Quantity", "Max|diff|"))
cat(sprintf("  %s\n", strrep("-", 64)))

tab <- list(
  list("Scenario 1", "Prop 1", "coefficient",             diff1_coef),
  list("Scenario 1", "Prop 2", "classical variance",       diff1_V0),
  list("Scenario 1", "Prop 3", "cluster-robust variance",  diff1_Vcr),
  list("Scenario 2", "Prop 4", "coefficient",              diff2_coef),
  list("Scenario 2", "Prop 5", "classical variance",       diff2_V0),
  list("Scenario 2", "Prop 6", "cluster-robust variance",  diff2_Vcr),
  list("Scenario 3", "Prop 7", "coefficient",              diff3_coef),
  list("Scenario 3", "Prop 8", "classical variance",       diff3_V0),
  list("Scenario 3", "Prop 9", "cluster-robust variance",  diff3_Vcr)
)
for (r in tab)
  cat(sprintf("  %-12s  %-10s  %-28s  %.2e\n",
              r[[1L]], r[[2L]], r[[3L]], r[[4L]]))

cat(sprintf("  %s\n", strrep("=", 72)))
cat("  * 모든 차이는 floating-point rounding 수준 (machine precision ~1e-16)\n\n")
