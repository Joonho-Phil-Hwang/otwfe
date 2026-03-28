# =============================================================================
# test_bench.R
# Rcpp 배치 vs pure-R fallback 속도 비교 (N=10,000, T_max=5, k=5)
# =============================================================================

setwd("/Users/joonhohwang/Desktop/claude_code/Simulation")
source("../online_twfe_core.R")

set.seed(42)
N      <- 10000L
T_max  <- 5L
k      <- 5L
x_cols    <- paste0("x", seq_len(k))
beta_true <- seq(0.5, by = 0.3, length.out = k)
lambda_t  <- c(0, 0.5, -0.3, 0.8, -0.2)
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
cat(sprintf("N=%s, 총 관측치=%s\n\n",
            format(N, big.mark=","), format(nrow(df), big.mark=",")))

# --- Rcpp 배치 ---
cat("=== Rcpp 배치 ===\n")
t1 <- proc.time()
fit_rcpp <- otwfe(data=df, id_col="id", time_col="time",
                  y_col="y", x_cols=x_cols,
                  track_all_units=FALSE, verbose=FALSE)
t_rcpp <- (proc.time()-t1)[["elapsed"]]
cat(sprintf("  시간: %.2f 초\n", t_rcpp))

# --- pure-R fallback (Rcpp 비활성화) ---
cat("=== pure-R fallback ===\n")
saved_flag <- .alg1_batch_rcpp_available
.alg1_batch_rcpp_available <<- FALSE
t2 <- proc.time()
fit_r <- otwfe(data=df, id_col="id", time_col="time",
               y_col="y", x_cols=x_cols,
               track_all_units=FALSE, verbose=FALSE)
t_r <- (proc.time()-t2)[["elapsed"]]
cat(sprintf("  시간: %.2f 초\n", t_r))
.alg1_batch_rcpp_available <<- saved_flag

cat(sprintf("\n속도 향상: %.1fx\n", t_r / t_rcpp))

# 결과 일치 확인
theta_rcpp <- fit_rcpp$state$theta_hat[x_cols]
theta_r    <- fit_r$state$theta_hat[x_cols]
cat(sprintf("theta max|diff| (Rcpp vs R): %.3e\n", max(abs(theta_rcpp - theta_r))))
