# otwfe

**Online Updating for Two-Way Fixed Effects Panel Regression**

An R package implementing algebraically exact online updating of two-way fixed effects (TWFE) panel regression estimators. Coefficients, classical variance, and cluster-robust variance (HC0, Arellano 1987) are updated sequentially without revisiting historical data.

Based on Hwang & Lee (2026).

---

## Installation

```r
# Install from GitHub
remotes::install_github("Joonho-Phil-Hwang/otwfe")
```

Or load directly from source for development:

```r
library(devtools)
load_all(".")   # run from the package root directory
```

---

## Usage

### In-memory estimation: `otwfe()`

```r
library(otwfe)

fit <- otwfe(
  data     = df,
  id_col   = "id",
  time_col = "time",
  y_col    = "y",
  x_cols   = c("x1", "x2")
)

coef(fit)                      # coefficients (x covariates only)
coef(fit, which = "all")       # includes time FE dummies
vcov(fit)                      # cluster-robust VCV (default)
vcov(fit, type = "classical")  # homoskedastic VCV
confint(fit)                   # 95% CI, cluster-robust
nobs(fit)                      # number of observations
print(fit)
summary(fit)
```

### File-based estimation: `otwfe_file()`

For datasets too large to load into memory. The file must be sorted by `id_col`.

```r
fit <- otwfe_file(
  path       = "panel_data.csv",
  id_col     = "id",
  time_col   = "time",
  y_col      = "y",
  x_cols     = c("x1", "x2"),
  chunk_size = 1e6L        # rows per chunk (default: 1,000,000)
)
```

All S3 methods (`coef`, `vcov`, `nobs`, `confint`, `print`, `summary`) work identically on results from both functions.

### Low-level streaming API

For fine-grained control over chunked processing:

```r
handle <- otwfe_init(x_cols = c("x1", "x2"), time_col = "time",
                     id_col = "id", y_col = "y", T_support = 5L)
handle <- otwfe_update(handle, chunk1)
handle <- otwfe_update(handle, chunk2)
fit    <- otwfe_finalize(handle)
```

---

## How It Works

Standard TWFE eliminates individual fixed effects via within transformation and includes time dummies as regressors. When new data arrives, a full re-estimation requires rebuilding the entire design matrix — infeasible when historical data cannot be revisited.

This package derives closed-form update formulas for three scenarios:

| Algorithm | Scenario |
|-----------|----------|
| **Algorithm 1** | New individual unit arrives |
| **Algorithm 2** | New observation for an existing unit, within known time support |
| **Algorithm 3** | New calendar time period (expands parameter dimension by 1) |

Each update modifies only a small set of stored summary objects. Results are **algebraically identical** to a full offline re-estimation (verified to machine precision against `plm`).

### What is stored

The state object contains:

| Object | Description |
|--------|-------------|
| `inv_dotZtZ` | Inverse of the within-demeaned design matrix cross-product |
| `theta_hat` | Coefficient vector |
| `sigma2_hat` | Residual variance |
| `Vcr_hat` | Cluster-robust variance matrix (HC0, Arellano) |
| `A_N`, `B_N`, `M_ss` | Aggregates for cluster-robust updating |

State size is **O(p²)** regardless of sample size — approximately 0.02 MB for k=2, T=5.

### Key property of Algorithm 3

When a new calendar time arrives, the parameter dimension increases (p → p+1). This is valid because the new time dummy equals zero for all pre-update observations, so the update reduces to a rank-1 extension without accessing historical data.

---

## Requirements

- R ≥ 4.0.0
- `Rcpp`, `RcppArmadillo` (compiled at installation)
- `data.table` (for `otwfe_file()`)

---

## Reference

Hwang, J. & Lee, S. (2026). *Online Updating for Linear Panel Regressions.*
