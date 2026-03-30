library(otwfe)

# DGP는 Simulation/sim_dgp.R 에서 로드
source(system.file("../Simulation/sim_dgp.R", package = "otwfe",
                   mustWork = FALSE) |>
       (\(p) if (nchar(p) > 0) p else "../../Simulation/sim_dgp.R")())

set.seed(42)
gen <- generate_panel_52(N = 500L, T = 5L, T_min = 2L, seed = 42L, verbose = FALSE)
df  <- gen$df

fit <- otwfe(
  data     = df,
  id_col   = "id",
  time_col = "time",
  y_col    = "y",
  x_cols   = c("x1", "x2"),
  verbose  = FALSE
)

test_that("otwfe returns otwfe object", {
  expect_s3_class(fit, "otwfe")
})

test_that("coef returns named numeric vector", {
  b <- coef(fit)
  expect_named(b, c("x1", "x2"))
  expect_true(is.numeric(b))
})

test_that("vcov returns symmetric matrix", {
  V <- vcov(fit)
  expect_equal(dim(V), c(2L, 2L))
  expect_equal(V, t(V))
  # diagonal entries must be positive
  expect_true(all(diag(V) > 0))
})

test_that("nobs returns total observation count", {
  expect_equal(nobs(fit), gen$N_obs)
})

test_that("confint has correct dimensions and column names", {
  ci <- confint(fit)
  expect_equal(dim(ci), c(2L, 2L))
  expect_equal(rownames(ci), c("x1", "x2"))
})

test_that("otwfe matches plm to machine precision", {
  skip_if_not_installed("plm")
  library(plm)
  pdf     <- pdata.frame(df, index = c("id", "time"))
  fml     <- y ~ x1 + x2 + factor(time)
  plm_fit <- plm(fml, data = pdf, model = "within", effect = "individual")

  coef_diff <- max(abs(coef(fit) - coef(plm_fit)[c("x1", "x2")]))
  expect_lt(coef_diff, 1e-10)

  V0_plm <- vcov(plm_fit)[c("x1", "x2"), c("x1", "x2")]
  V0_ot  <- vcov(fit, type = "classical")
  expect_lt(max(abs(V0_ot - V0_plm)), 1e-10)
})
