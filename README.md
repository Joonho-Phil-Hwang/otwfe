# otwfe: Online Updating for Two-Way Fixed Effects Panel Regression

<!-- badges -->
![Status](https://img.shields.io/badge/status-research--prototype-yellow)
![Language](https://img.shields.io/badge/language-R-blue)
![License](https://img.shields.io/badge/license-MIT-green)

> **Hwang & Lee (2026)** — *Online Updating for Linear Panel Regressions*

---

## Overview

`otwfe` implements **algebraically exact online updating** for Two-Way Fixed Effects (TWFE) panel regression. When new data arrives sequentially, `otwfe` updates the TWFE estimates **without ever revisiting the full historical dataset**.

The key guarantee is not an approximation:

> The online estimate after processing *N* units is **identical** to the offline full-sample OLS estimate up to floating-point precision.

This makes `otwfe` suitable for settings where:

- The full panel is too large to hold in memory at once
- Raw historical data cannot be re-accessed for privacy or security reasons
- Estimates must be updated in real time as new observations stream in

---

## The Problem

Consider a linear panel regression model with individual fixed effects $\alpha_i$ and time fixed effects $\lambda_t$:

$$Y_{it} = \alpha_i + \lambda_t + X_{it}'\beta + \varepsilon_{it}$$

Standard estimation requires storing and processing the entire dataset. For panels with $N = 10^6$ units and $T = 5$ periods, this means holding tens of millions of rows in memory and re-running OLS every time new data arrives, which is time-consuming and (sometimes) infeasible.

**`otwfe` breaks this bottleneck.** It maintains a compact *state object* of size $O(p^2)$ (where $p = k + T - 1$ is the number of regressors), and updates the state incrementally as each new unit or observation arrives. No raw data needs to be stored after processing.

---

## What Gets Updated

For each streaming update, `otwfe` produces **exact** updates to:

| Quantity | Description |
|---|---|
| $\hat{\beta}$ | TWFE coefficient vector |
| $(\dot{Z}'\dot{Z})^{-1}$ | Inverse of the within-demeaned Gram matrix |
| $\hat{\sigma}^2$ | Classical variance estimate |
| $\hat{V}_{\text{classical}} = \hat{\sigma}^2 (\dot{Z}'\dot{Z})^{-1}$ | Classical (homoskedastic) variance–covariance matrix |
| $\hat{V}_{\text{CR}} = (\dot{Z}'\dot{Z})^{-1} \hat{M} (\dot{Z}'\dot{Z})^{-1}$ | Cluster-robust variance (Arellano HC0) |

Here $\dot{Z}_{it} = Z_{it} - \bar{Z}_i$ denotes the within-demeaned design matrix (individual FE removed by within transformation; time FE included as dummies).

---

## Algorithm Structure

The streaming updates are governed by three algorithms, each handling a different type of new data arrival.

### Algorithm 1 — New Unit Arrives

A unit $i = N+1$ with $T_i$ observations $(Y_{N+1,t}, X_{N+1,t})$ arrives for the first time.

**Key idea:** The within-demeaned contribution $S_i = \dot{Z}_i'\dot{Z}_i$ of the new unit is independent of all previous units. The Woodbury matrix identity gives a rank-$T_i$ update to $(\dot{Z}'\dot{Z})^{-1}$ without re-inverting the full matrix.

**State updates (Props 1–3):**

$$(\dot{Z}^{*\prime}\dot{Z}^*)^{-1} = (\dot{Z}'\dot{Z})^{-1} - (\dot{Z}'\dot{Z})^{-1}\dot{Z}_{N+1}'H_{N+1}^{-1}\dot{Z}_{N+1}(\dot{Z}'\dot{Z})^{-1}$$

$$\hat{\beta}^* = \hat{\beta} + (\dot{Z}'\dot{Z})^{-1}\dot{Z}_{N+1}'H_{N+1}^{-1}\tilde{e}_{N+1}$$

where $H_{N+1} = I_{T_i} + \dot{Z}_{N+1}(\dot{Z}'\dot{Z})^{-1}\dot{Z}_{N+1}'$ and $\tilde{e}_{N+1} = \dot{Y}_{N+1} - \dot{Z}_{N+1}\hat{\beta}$.

**Stored per warm-up unit:** $S_i$, $s_i = \dot{Z}_i'\dot{Y}_i$, $\bar{Z}_i$, $\bar{Y}_i$, $T_i$

### Algorithm 2 — New Observation for Existing Unit (Existing Calendar Time)

Unit $i$ already in the sample receives a new observation at a calendar time $t \leq T_{\max}$ already seen.

**Key idea:** The individual mean $\bar{Z}_i$ shifts when $T_i \to T_i + 1$. The within deviation of the new observation depends on the current $\bar{Z}_i$. A scalar Sherman–Morrison update suffices.

**State updates (Props 4–6):**

$$(\dot{Z}^{*\prime}\dot{Z}^*)^{-1} = \left(I - s_{i,t^*}(\dot{Z}'\dot{Z})^{-1}\dot{z}_{i,t^*}\dot{z}_{i,t^*}'\right)(\dot{Z}'\dot{Z})^{-1}$$

where $s_{i,t^*} = \kappa_i / (1 + \kappa_i h_{i,t^*})$, $\kappa_i = T_i/(T_i+1)$, $h_{i,t^*} = \dot{z}_{i,t^*}'(\dot{Z}'\dot{Z})^{-1}\dot{z}_{i,t^*}$.

> **Note:** Only the within-transformed quantities of unit $i$ change. All other units' contributions remain unchanged.

### Algorithm 3 — New Calendar Time (New Time Dummy Required)

A new calendar time $T+1$ not previously observed appears. The parameter dimension increases: $p \to p+1$.

**Key idea:** The new time dummy equals 0 for all pre-update observations. This zero-extension property allows the existing $(\dot{Z}'\dot{Z})^{-1}$ to be embedded exactly into the larger $(p+1) \times (p+1)$ inverse matrix without any costly recomputation.

**State updates (Props 7–9):**

$$(\dot{Z}^{*\prime}\dot{Z}^*)^{-1} = \begin{pmatrix} (\dot{Z}'\dot{Z})^{-1} + \frac{q q'}{g_2^2} & -q/g_2 \\ -q'/g_2 & (1+\kappa h)/(\kappa g_2^2) \end{pmatrix}$$

where $q = (\dot{Z}'\dot{Z})^{-1}g$, $h = g'q$, $g_2 = 1$ (within deviation of new dummy = 1).

> **Zero-extension:** All stored per-unit objects ($S_i$, $s_i$, $\bar{Z}_i$) are extended to $p+1$ dimensions by appending zeros—valid because the new time dummy was always 0 for pre-update observations.

---

## Memory Architecture

`otwfe` uses a two-tier memory design:

```
┌─────────────────────────────────────────────────────────────┐
│  Global State  (always stored, O(p²))                       │
│  ─────────────────────────────────────────────────────────  │
│  inv_dotZtZ   : (dot{Z}'dot{Z})^{-1}         [p × p]       │
│  theta_hat    : β-hat                          [p]          │
│  sigma2_hat   : σ²-hat                         [1]          │
│  Vcr_hat      : cluster-robust VCV             [p × p]      │
│  A_N, B_N     : aggregate matrices for Vcr     [p × p²]    │
│                                                [p² × p²]    │
└─────────────────────────────────────────────────────────────┘
┌─────────────────────────────────────────────────────────────┐
│  Per-Unit State  (warm-up units only, O(n_warm × p²))       │
│  ─────────────────────────────────────────────────────────  │
│  S_i  : dot{Z}_i' dot{Z}_i                    [p × p]      │
│  s_i  : dot{Z}_i' dot{Y}_i                    [p]          │
│  barZ_i, barY_i : within means                [p], [1]     │
│  T_i  : number of observations                [1]          │
└─────────────────────────────────────────────────────────────┘
         Streaming units (Algorithm 1 only):
         Per-unit state is NOT stored by default.
         Only the global state is updated.
```

### Warm-up Strategy

Before streaming begins, a small *warm-up sample* is selected. This sample must collectively span **all calendar time periods** in the data. The warm-up batch OLS initializes the global state $(\dot{Z}'\dot{Z})^{-1}$, $\hat{\beta}$, $\hat{\sigma}^2$, $\hat{V}_{\text{CR}}$, ensuring that `T_support = T_max` from the start. After warm-up, all streaming units are new and processed via Algorithm 1 only.

The greedy warm-up selector chooses the minimum number of units (typically $\max(3p, 30)$) that cover all $T$ calendar times.

---

## Performance

### Memory

| Setting | Full Dataset | otwfe State | Reduction |
|---|---|---|---|
| N=100K, T=5, k=5 | ~100 MB | ~0.14 MB | **>99%** |
| N=1M, T=5, k=5 | ~1 GB | ~0.14 MB | **>99.9%** |

The global state size is $O(p^2)$, independent of $N$.

### Speed (Apple M-series, arm64)

| N | plm offline | otwfe (Rcpp) | otwfe (pure R) |
|---|---|---|---|
| 10,000 | < 1s | **0.6s** | 3.1s |
| 100,000 | ~5s | **5.6s** | 57s |
| 1,000,000 | ~61s* | **~55s** | ~82min |

\* plm does not compute cluster-robust SE at N=1M in the table above; `vcovHC` would take additional time.

Speed comes from Rcpp (C++ via RcppArmadillo): streaming new units are processed in a two-pass C++ loop, replacing $N$ sequential R-level Woodbury updates with a single matrix inversion.

---

## Installation

```r
# Development version (from source)
# Requires: R >= 4.0, Rcpp, RcppArmadillo, plm (for benchmarking)
install.packages(c("Rcpp", "RcppArmadillo"))

# Then in the project directory:
setwd("path/to/claude_code")
source("online_twfe_core.R")
```

> **Note:** A full CRAN-installable package (`otwfe`) is in development.

---

## Usage

```r
library(plm)    # for comparison
source("online_twfe_core.R")

# Synthetic unbalanced panel: N=1,000, T_max=5, k=3
set.seed(42)
N <- 1000L;  T_max <- 5L;  k <- 3L
x_cols <- paste0("x", seq_len(k))

T_i_vec <- sample(2L:T_max, N, replace = TRUE)
df <- do.call(rbind, lapply(seq_len(N), function(i) {
  Ti  <- T_i_vec[i]
  X   <- matrix(rnorm(Ti * k), Ti, k, dimnames = list(NULL, x_cols))
  y   <- rnorm(1) + c(0, .5, -.3, .8, -.2)[seq_len(Ti)] + X %*% c(1, -.5, .8) + rnorm(Ti, sd=.5)
  data.frame(id = i, time = seq_len(Ti), y = y, X)
}))

# --- Online estimation ---
fit <- otwfe(
  data            = df,
  id_col          = "id",
  time_col        = "time",
  y_col           = "y",
  x_cols          = x_cols,
  track_all_units = FALSE,   # memory-efficient mode (default)
  verbose         = TRUE
)

# Coefficients
fit$state$theta_hat[x_cols]

# Classical standard errors
sqrt(diag(fit$state$sigma2_hat * fit$state$inv_dotZtZ[x_cols, x_cols]))

# Cluster-robust standard errors (Arellano HC0)
sqrt(diag(fit$state$Vcr_hat[x_cols, x_cols]))
```

### Key Parameters

| Parameter | Description | Default |
|---|---|---|
| `track_all_units` | Store per-unit state for every streaming unit | `FALSE` |
| `warmup_ids` | Manually specify warm-up unit IDs | `NULL` (auto) |
| `warmup_n` | Minimum number of warm-up units | `max(3p, 30)` |
| `baseline_time` | Omitted time dummy for identification | `min(time)` |
| `verbose` | Print progress | `TRUE` |

### `track_all_units`

- `FALSE` (default): Only warm-up unit states are stored. Streaming units processed via Rcpp batch Algorithm 1. Memory usage is $O(p^2)$ regardless of $N$. Algorithms 2 and 3 are unavailable for streaming units.
- `TRUE`: Per-unit state stored for all units. Supports Algorithms 2 and 3 for future updates. Memory usage is $O(N \times p^2)$.

---

## Accuracy

All three estimates match the offline `plm` benchmark at machine precision:

```
N = 1,000,000 | T_max = 5 | k = 5 | unbalanced T_i ∈ {2,...,5}

  theta  max|diff| vs plm : 6.2e-14   (< 1e-10 ✓)
  V0     max|diff| vs plm : 5.5e-21   (< 1e-10 ✓)
  Vcr    max|diff| vs plm : 8.9e-21   (< 1e-10 ✓)
```

---

## File Structure

```
claude_code/
├── online_twfe_core.R       # Main implementation (otwfe + Algorithms 1–3)
├── src/
│   └── alg1_batch.cpp       # Rcpp batch Algorithm 1 (two-pass C++)
├── online_twfe_sim.R        # Simulation utilities
├── test_large.R             # Benchmark: N=1M vs plm
├── test_rcpp.R              # Accuracy verification (N=200)
├── test_bench.R             # Speed comparison: Rcpp vs pure-R
├── online_algorithms.pdf    # Working paper (Hwang & Lee, 2026)
└── README.md
```

---

## Technical Notes

### Why Individual FE via Within Transformation, Not Dummies?

With $N = 10^6$ units, including individual dummies would require inverting a $(k + N + T - 2) \times (k + N + T - 2)$ matrix. `otwfe` instead applies the **within (Frisch–Waugh) transformation** to eliminate individual effects analytically, reducing the problem to a $p = k + T - 1$ dimensional system that stays small regardless of $N$.

### Why Exact, Not Approximate?

Methods like stochastic gradient descent or incremental SVD produce approximations that converge over time. `otwfe` provides **algebraically exact** updates via the Woodbury and Sherman–Morrison identities—the online estimate equals the offline full-sample estimate to machine precision at every step.

### Cluster-Robust Variance Without Raw Data

The Arellano (1987) HC0 cluster-robust estimator is:

$$\hat{V}_{\text{CR}} = (\dot{Z}'\dot{Z})^{-1} \hat{M} (\dot{Z}'\dot{Z})^{-1}, \qquad \hat{M} = \sum_{i=1}^N r_i r_i', \quad r_i = \dot{Z}_i'\hat{e}_i$$

Updating $\hat{M}$ requires $r_i = s_i - S_i\hat{\beta}$, which depends on the **final** $\hat{\beta}$. In the Rcpp batch implementation, a two-pass approach computes all $S_i$, $s_i$ in Pass 1 (to obtain $\hat{\beta}_{\text{new}}$), then recomputes $r_i$ in Pass 2—without storing raw data.

---

## References

- **Hwang, J. & Lee, S. (2026).** *Online Updating for Linear Panel Regressions.* Working paper.
- Arellano, M. (1987). Computing robust standard errors for within-groups estimators. *Oxford Bulletin of Economics and Statistics*, 49(4), 431–434.
- Woodbury, M. A. (1950). Inverting modified matrices. *Memorandum Report 42, Statistical Research Group, Princeton University.*
- Sherman, J. & Morrison, W. J. (1950). Adjustment of an inverse matrix corresponding to a change in one element of a given matrix. *Annals of Mathematical Statistics*, 21(1), 124–127.
