# =============================================================================
# practice.R
#
# Hands-on practice script for the otwfe package
# Run this line by line in RStudio (or source it all at once)
#
# Sections:
#   0. Setup
#   1. In-memory estimation: otwfe()
#   2. S3 methods: coef / vcov / nobs / confint / print / summary
#   3. Comparison with plm (accuracy check)
#   4. File-based estimation: otwfe_file()
#   5. Large-scale demo (N = 500,000)
# =============================================================================

# =============================================================================
# 0. Setup
# =============================================================================
# Install and load the package
#   - First run or outdated version (missing otwfe() or Rcpp functions): install from GitHub
#   - Up-to-date version already installed: load directly
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
needs_install <- !requireNamespace("otwfe", quietly = TRUE) ||
                 !"otwfe" %in% getNamespaceExports("otwfe") ||
                 !exists("alg1_batch_cpp", where = asNamespace("otwfe"), inherits = FALSE)
if (needs_install) {
  # upgrade = "never": suppress interactive prompt asking whether to update dependencies
  remotes::install_github("Joonho-Phil-Hwang/otwfe", upgrade = "never")
}
suppressPackageStartupMessages(library(otwfe))

# Load DGP helper: use local file if available, otherwise source from GitHub
if (file.exists("Simulation/sim_dgp.R")) {
  source("Simulation/sim_dgp.R")
} else {
  source("https://raw.githubusercontent.com/Joonho-Phil-Hwang/otwfe/main/Simulation/sim_dgp.R")
}

suppressPackageStartupMessages({
  library(plm)         # offline benchmark
  library(data.table)  # fast CSV I/O
})

cat("Setup complete.\n")


# =============================================================================
# 1. In-memory estimation: otwfe()
# =============================================================================
# Generate a small synthetic panel: N=1000 units, T_max=5, k=2 covariates
set.seed(42)
gen <- generate_panel_52(N = 1000L, T = 5L, T_min = 2L, seed = 42L, verbose = FALSE)
df  <- gen$df

cat("\n--- Dataset summary ---\n")
cat(sprintf("Units (N): %s\n", format(gen$N, big.mark = ",")))
cat(sprintf("Obs   (n): %s\n", format(gen$N_obs, big.mark = ",")))
cat(sprintf("Columns  : %s\n", paste(names(df), collapse = ", ")))
cat(sprintf("True beta: x1 = 1.3,  x2 = -0.7\n"))

# Fit the model
fit <- otwfe(
  data     = df,
  id_col   = "id",
  time_col = "time",
  y_col    = "y",
  x_cols   = c("x1", "x2"),
  verbose  = TRUE
)


# =============================================================================
# 2. S3 methods
# =============================================================================

# --- coef() ------------------------------------------------------------------
cat("\n\n--- coef() ---\n")
coef(fit)                        # x1, x2 only (default)
coef(fit, which = "all")         # includes time FE dummies


# --- vcov() ------------------------------------------------------------------
cat("\n--- vcov() ---\n")
vcov(fit)                        # cluster-robust (default)
vcov(fit, type = "classical")    # homoskedastic


# --- nobs() ------------------------------------------------------------------
cat("\n--- nobs() ---\n")
nobs(fit)                        # number of observations


# --- confint() ---------------------------------------------------------------
cat("\n--- confint() ---\n")
confint(fit)                              # 95% CI, cluster-robust (default)
confint(fit, level = 0.99)               # 99% CI
confint(fit, type = "classical")          # classical SE
confint(fit, parm = "x1")                # single coefficient


# --- print() -----------------------------------------------------------------
cat("\n--- print() ---\n")
print(fit)                               # compact output with significance stars
# or just: fit


# --- summary() ---------------------------------------------------------------
cat("\n--- summary() ---\n")
sm <- summary(fit)                        # returns a summary.otwfe object
print(sm)                                 # prints detailed table

# You can also access the coefficient table directly:
sm$coefficients


# =============================================================================
# 3. Comparison with plm
# =============================================================================
cat("\n\n--- plm comparison ---\n")

pdf     <- pdata.frame(df, index = c("id", "time"))
fml     <- y ~ x1 + x2 + factor(time)
plm_fit <- plm(fml, data = pdf, model = "within", effect = "individual")

# Coefficients
cat("coef diff (max abs):",
    max(abs(coef(fit) - coef(plm_fit)[c("x1", "x2")])), "\n")

# Classical VCV
V0_plm  <- vcov(plm_fit)[c("x1", "x2"), c("x1", "x2")]
V0_ot   <- vcov(fit, type = "classical")
cat("V0 diff   (max abs):", max(abs(V0_ot - V0_plm)), "\n")

# Cluster-robust VCV
Vcr_plm <- vcovHC(plm_fit, method = "arellano", type = "HC0")[c("x1","x2"), c("x1","x2")]
Vcr_ot  <- vcov(fit, type = "robust")
cat("Vcr diff  (max abs):", max(abs(Vcr_ot - Vcr_plm)), "\n")
# All three should be < 1e-10 (machine precision)


# =============================================================================
# 4. File-based estimation: otwfe_file()
# =============================================================================
cat("\n\n=== otwfe_file() demo ===\n")

# Save dataset to a temporary CSV
tmp_csv <- tempfile(fileext = ".csv")
fwrite(df, tmp_csv)
cat(sprintf("CSV saved: %.1f MB\n", file.size(tmp_csv) / 1e6))

# Run otwfe_file() — reads in chunks, never loads full file
fit_file <- otwfe_file(
  path       = tmp_csv,
  id_col     = "id",
  time_col   = "time",
  y_col      = "y",
  x_cols     = c("x1", "x2"),
  chunk_size = 500L,     # small chunk to demonstrate boundary handling
  verbose    = TRUE
)
unlink(tmp_csv)

# All S3 methods work on the otwfe_file() result too
coef(fit_file)
summary(fit_file)

# Accuracy vs in-memory otwfe()
cat("\notwfe vs otwfe_file (coef diff):",
    max(abs(coef(fit) - coef(fit_file))), "\n")


# =============================================================================
# 4b. otwfe_file() with year-format time column
# =============================================================================
cat("\n--- otwfe_file: year-format time (2020-2024) ---\n")

df_yr       <- df
df_yr$time  <- df_yr$time + 2019L   # {1,...,5} -> {2020,...,2024}
tmp_yr      <- tempfile(fileext = ".csv")
fwrite(df_yr, tmp_yr)

fit_yr <- otwfe_file(
  path     = tmp_yr,
  id_col   = "id",  time_col = "time",
  y_col    = "y",   x_cols   = c("x1", "x2"),
  verbose  = FALSE
)
unlink(tmp_yr)

cat("Detected time levels:", paste(fit_yr$time_levels_original, collapse=", "), "\n")
cat("coef diff vs original:", max(abs(coef(fit) - coef(fit_yr))), "\n")


# =============================================================================
# 5. Large-scale demo: N = 500,000
# =============================================================================
cat("\n\n=== Large-scale demo: N = 500,000 ===\n")

gen_big <- generate_panel_52(N = 500000L, T = 5L, T_min = 2L,
                              seed = 42L, verbose = TRUE)
tmp_big <- tempfile(fileext = ".csv")
fwrite(gen_big$df, tmp_big)
cat(sprintf("CSV: %.0f MB  (%s obs)\n",
            file.size(tmp_big) / 1e6,
            format(gen_big$N_obs, big.mark = ",")))
rm(gen_big); invisible(gc())

t0 <- proc.time()
fit_big <- otwfe_file(
  path       = tmp_big,
  id_col     = "id",
  time_col   = "time",
  y_col      = "y",
  x_cols     = c("x1", "x2"),
  chunk_size = 1e6L,
  verbose    = TRUE
)
unlink(tmp_big)
cat(sprintf("Elapsed: %.1f sec\n", (proc.time() - t0)[["elapsed"]]))

# Results
print(fit_big)
summary(fit_big)

# Bias check (true beta = 1.3, -0.7)
cat(sprintf("\nBias: x1 = %+.4f,  x2 = %+.4f\n",
            coef(fit_big)["x1"] - 1.3,
            coef(fit_big)["x2"] - (-0.7)))
