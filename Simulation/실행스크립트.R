source("/Users/joonhohwang/Downloads/online_twfe.R")

# N_grid <- c(1000000L, 5000000L, 10000000L, 20000000L, 50000000L)
N_grid <- c(500L)
res <- run_benchmark_grid(
  N_grid = N_grid,
  N_warm = 5L,
  do_offline = TRUE,
  offline_maxN_attempt = Inf,
  verify_smallN = TRUE,
  verify_maxN = 1000000L,
  verbose = TRUE,
  alg1_batch_size = 50000L
)


