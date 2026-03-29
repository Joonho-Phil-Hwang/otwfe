# otwfe: Online Updating for Two-Way Fixed Effects Panel Regression

<!-- badges -->
![Status](https://img.shields.io/badge/status-research--prototype-yellow)
![Language](https://img.shields.io/badge/language-R-blue)
![License](https://img.shields.io/badge/license-MIT-green)

> **Hwang & Lee (2026)** — *Online Updating for Linear Panel Regressions*

---

## Overview

`otwfe` implements **algebraically exact online updating** for Two-Way Fixed Effects (TWFE) panel regression. When new data arrives sequentially—unit by unit or observation by observation—`otwfe` updates the TWFE estimates **without ever revisiting the full historical dataset**.

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

Standard estimation requires storing and processing the entire dataset. For panels with $N = 10^6$ units and $T = 5$ periods, this means holding tens of millions of rows in memory and re-running OLS every time new data arrives.

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

> **Note:** Only the within-transformed quantities of unit $i$ change. All other units' contributions remain unchanged.

### Algorithm 3 — New Calendar Time (New Time Dummy Required)

A new calendar time $T+1$ not previously observed appears. The parameter dimension increases: $p \to p+1$.

**Key idea:** The new time dummy equals 0 for all pre-update observations. This zero-extension property allows the existing $(\dot{Z}'\dot{Z})^{-1}$ to be embedded exactly into the larger $(p+1) \times (p+1)$ inverse matrix without any costly recomputation.

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

## Installation

```r
# Requires: R >= 4.0, Rcpp, RcppArmadillo, data.table
install.packages(c("Rcpp", "RcppArmadillo", "data.table"))

# Source the core file from the project root:
setwd("path/to/claude_code")
source("online_twfe_core.R")

# For large-scale CSV analysis:
source("uni/otwfe_file.R")
```

> **Note:** `online_twfe_core.R` automatically compiles `src/alg1_batch.cpp` via Rcpp on first load. If compilation fails, it falls back to a pure-R implementation.

---

## Usage

### In-memory estimation (`otwfe`)

For moderate-sized panels that fit in memory:

```r
source("online_twfe_core.R")

fit <- otwfe(
  data     = df,
  id_col   = "id",
  time_col = "time",
  y_col    = "y",
  x_cols   = c("x1", "x2"),
  verbose  = TRUE
)

# Coefficients
fit$state$theta_hat[c("x1", "x2")]

# Classical standard errors
sqrt(diag(fit$state$sigma2_hat * fit$state$inv_dotZtZ)[c("x1", "x2")])

# Cluster-robust standard errors (Arellano HC0)
sqrt(diag(fit$state$Vcr_hat)[c("x1", "x2")])
```

### Large-scale CSV analysis (`otwfe_file`)

For datasets too large to load into memory, `otwfe_file()` reads directly from a CSV file in chunks. The full dataset is **never materialized in R**—only one chunk (~80–100 MB) and a compact state object (~0.02 MB) reside in memory at any time.

```r
source("uni/otwfe_file.R")

result <- otwfe_file(
  path       = "/path/to/large_panel.csv",
  id_col     = "id",
  time_col   = "time",
  y_col      = "y",
  x_cols     = c("x1", "x2"),
  chunk_size = 2e6L,   # rows per chunk (default 1e6)
  verbose    = TRUE
)
```

**Prerequisites for `otwfe_file()`:**

- File must be **sorted by `id_col`** (non-decreasing). Each unit's observations must be contiguous.
- File format: **CSV** with a header row.
- `time_col` can contain arbitrary integer values (e.g., years `{2020, 2021, 2022}`); they are automatically remapped to `{1, ..., T}` internally.

**Output:**

```
  추정 결과 (x 변수):
                    coef        SE    SE(CR)        t
    ----------------------------------------------------
    x1            1.3002    0.0002    0.0002  6498.97
    x2           -0.7002    0.0002    0.0002  -3500.48
    ----------------------------------------------------
    * SE    : classical variance 기반 표준오차
    * SE(CR): HC0 cluster-robust variance 기반 표준오차 (Arellano)
```

---

## How `otwfe_file()` Works: 5-Step Pipeline

```
CSV File (sorted by id)
        │
        ▼
┌───────────────────────────────────────────────────────────┐
│  Step 1: T_support Detection                              │
│  ─────────────────────────────────────────────────────    │
│  - File connection 유지 → 순차 읽기 (재스캔 없음)          │
│  - 100K행 청크 단위로 time 컬럼만 추출                     │
│  - 새 time 값이 없는 청크가 5회 연속이면 조기 종료          │
│  - 정렬된 패널에서 보통 수십만 행 이내에서 완료             │
│  → T_support, time_remap({year} → {1,...,T}) 확정         │
└───────────────────────────────────────────────────────────┘
        │
        ▼
┌───────────────────────────────────────────────────────────┐
│  Step 2: Initialization Chunk Assembly                    │
│  ─────────────────────────────────────────────────────    │
│  - chunk_size 행씩 읽어 누적                               │
│  - 모든 T calendar time이 포함될 때까지 누적               │
│  - 마지막 unit 경계(unit boundary) 탐지                    │
│  - 불완전한 마지막 unit은 다음 단계로 이월(carry-over)      │
└───────────────────────────────────────────────────────────┘
        │
        ▼
┌───────────────────────────────────────────────────────────┐
│  Step 3: State Initialization (Warm-up)                   │
│  ─────────────────────────────────────────────────────    │
│  - 초기 청크에서 greedy 방식으로 warm-up units 선택         │
│  - warm-up은 모든 T time period를 반드시 커버              │
│  - otwfe_init() + otwfe_update()로 초기 state 구성        │
└───────────────────────────────────────────────────────────┘
        │
        ▼
┌───────────────────────────────────────────────────────────┐
│  Step 4: Adaptive Chunked Processing                      │
│  ─────────────────────────────────────────────────────    │
│  - chunk_size 행씩 반복 읽기                               │
│  - 각 청크: carry-over 병합 → unit boundary 탐지           │
│  - 완전한 units만 otwfe_update()로 전달                    │
│  - 불완전한 마지막 unit → 다음 청크로 carry-over            │
│  - Rcpp 배치 처리: C++ 2-pass로 다수 unit 일괄 업데이트     │
└───────────────────────────────────────────────────────────┘
        │
        ▼
┌───────────────────────────────────────────────────────────┐
│  Step 5: Finalize                                         │
│  ─────────────────────────────────────────────────────    │
│  - otwfe_finalize()로 최종 theta_hat, Vcr_hat 계산        │
│  - 결과 출력: coef / SE / SE(CR) / t                      │
└───────────────────────────────────────────────────────────┘
        │
        ▼
   State (~0.02 MB)
   theta_hat, Vcr_hat, inv_dotZtZ, sigma2_hat
```

### Unit Atomicity

`otwfe_file()` ensures that all observations for each unit are processed together. The `chunk_size` parameter controls the number of **rows** per read, but each chunk is trimmed at a unit boundary before being passed to `otwfe_update()`. Observations for a unit that straddles a chunk boundary are carried over to the next chunk.

### Time Remapping

Any monotone integer sequence can serve as `time_col` values (e.g., `{2020, 2021, 2022}` or `{1, 3, 5}`). `otwfe_file()` automatically detects all distinct time values in Step 1 and remaps them to `{1, ..., T}` before passing to the core algorithm.

---

## Performance

### Memory

The global state is $O(p^2)$, fixed regardless of $N$:

| Setting | CSV Size | otwfe State | plm |
|---|---|---|---|
| N=1M, T=5, k=2 | 221 MB | **0.02 MB** | 112 MB (in-memory) |
| N=5M, T=5, k=2 | 1.1 GB | **0.02 MB** | 560 MB (in-memory) |
| N=10M, T=5, k=2 | 2.2 GB | **0.02 MB** | 1.1 GB (in-memory) |
| N=20M, T=5, k=2 | 4.5 GB | **0.02 MB** | OOM |
| N=50M, T=5, k=2 | 11 GB | **0.02 MB** | OOM |

### Speed: `otwfe_file` vs `plm` (T=5, k=2, Apple M-series)

| N | otwfe_file | pdata.frame | plm() | vcov() | vcovHC() | θ diff vs plm |
|---|---|---|---|---|---|---|
| 1M | 15초 | 1.4초 | 6.7초 | 0.0초 | 16.2초 | **8.9e-15** |
| 5M | 63초 | 8.8초 | 43초 | 0.0초 | 80초 | **3.0e-14** |
| 10M | 131초 | 13초 | 151초 | 0.0초 | **FAIL** (OOM) | 3.6e-14 |
| 20M | 308초 | 49초 | **FAIL** (OOM) | — | — | — |
| 50M | 2218초 | **FAIL** (OOM) | — | — | — | — |

- plm OOM 한계: `vcovHC()` → N≥10M, `plm()` → N≥20M, `pdata.frame()` → N≥50M
- `otwfe_file`은 N=50M (관측치 1.75억, CSV 11GB)에서도 **state 0.02MB**만으로 완료

### Step 1 (T_support Detection) Optimization

파일 커넥션 기반 순차 읽기 + 조기 종료로 time 탐지 시간이 대폭 단축됩니다:

| N (CSV) | 이전 (`fread` 전체 스캔) | 이후 (조기 종료) | 스캔 행 수 |
|---|---|---|---|
| 5M (1.1 GB) | 2.1초 | 1.0초 | 600K행 |
| 10M (2.2 GB) | 1.4초 | 2.2초 | 600K행 |
| 20M (4.5 GB) | **1,097초** | **19초** | 600K행 |

---

## Accuracy

`otwfe_file()` matches offline `plm` at machine precision (where plm is computable):

```
N = 1,000,000 | T = 5 | k = 2 | unbalanced T_i ∈ {2,...,5}

  theta  max|diff| vs plm : 8.9e-15   (< 1e-10 ✓)
  V0     max|diff| vs plm : 8.8e-19   (< 1e-10 ✓)
  Vcr    max|diff| vs plm : 5.5e-20   (< 1e-10 ✓)

N = 5,000,000 | T = 5 | k = 2

  theta  max|diff| vs plm : 3.0e-14   (< 1e-10 ✓)
  V0     max|diff| vs plm : 2.0e-19   (< 1e-10 ✓)
  Vcr    max|diff| vs plm : 2.1e-20   (< 1e-10 ✓)
```

---

## Key Parameters

### `otwfe()`

| Parameter | Description | Default |
|---|---|---|
| `data` | `data.frame` (id 기준 정렬 권장) | — |
| `id_col` | 개인 식별자 컬럼명 | — |
| `time_col` | calendar time 컬럼명 | — |
| `y_col` | 종속변수 컬럼명 | — |
| `x_cols` | 공변량 컬럼명 벡터 | — |
| `track_all_units` | 모든 unit의 per-unit state 저장 | `FALSE` |
| `chunk_size` | 청크 당 unit 수 (Rcpp 배치) | `1e6L` |
| `verbose` | 진행 상황 출력 | `TRUE` |

### `otwfe_file()`

| Parameter | Description | Default |
|---|---|---|
| `path` | CSV 파일 경로 (id 기준 정렬 필수) | — |
| `id_col` | 개인 식별자 컬럼명 | — |
| `time_col` | calendar time 컬럼명 | — |
| `y_col` | 종속변수 컬럼명 | — |
| `x_cols` | 공변량 컬럼명 벡터 | — |
| `chunk_size` | 청크 당 최대 행 수 | `1e6L` |
| `sep` | CSV 구분자 | `","` |
| `verbose` | 진행 상황 출력 | `TRUE` |

---

## File Structure

```
claude_code/
├── online_twfe_core.R          # 핵심 알고리즘 (otwfe_init/update/finalize)
├── online_twfe_core_Ronly.R    # 순수 R 구현 (Rcpp 없이)
├── online_twfe.R               # 래퍼 함수 (otwfe)
├── src/
│   └── alg1_batch.cpp          # Rcpp 배치 Algorithm 1 (C++ 2-pass)
├── uni/
│   ├── otwfe_file.R            # 대용량 CSV 분석 함수 (otwfe_file)
│   ├── test_otwfe_file.R       # 정확도 검증: N=2000, plm 비교
│   ├── test_otwfe_file2.R      # 대용량 검증: T=5, N=5M/10M/20M
│   └── benchmark_vs_plm_T5k2.R # otwfe_file vs plm 벤치마크
├── Simulation/
│   ├── sim_dgp.R               # DGP: Hwang & Lee (2026) §5.2
│   ├── benchmark_vs_plm.R      # 벤치마크 (T=3)
│   ├── benchmark_vs_plm_T3k2.R # 벤치마크 (T=3, k=2, 대용량)
│   └── ...                     # 기타 시뮬레이션/검증 스크립트
├── Application/
│   ├── application_design.pdf  # 설계 문서
│   └── application_design.tex
├── online_algorithms.pdf        # 워킹 페이퍼 (Hwang & Lee, 2026)
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

### `fread(skip=n)` vs File Connection

`data.table::fread(skip=n)` re-reads the file from the beginning each time, causing $O(n^2)$ total I/O for $n$ chunks. `otwfe_file()` uses:

- **Step 1**: `file()` connection + `readLines()` — position maintained across calls, no re-scanning
- **Step 2–4**: `fread(skip=n)` — acceptable because chunk count is small (typically 10–100 chunks for large files)

---

## References

- **Hwang, J. & Lee, S. (2026).** *Online Updating for Linear Panel Regressions.* Working paper.
- Arellano, M. (1987). Computing robust standard errors for within-groups estimators. *Oxford Bulletin of Economics and Statistics*, 49(4), 431–434.
- Woodbury, M. A. (1950). Inverting modified matrices. *Memorandum Report 42, Statistical Research Group, Princeton University.*
- Sherman, J. & Morrison, W. J. (1950). Adjustment of an inverse matrix corresponding to a change in one element of a given matrix. *Annals of Mathematical Statistics*, 21(1), 124–131.

---

## Author

**Joonho Hwang** — PhD student, Econometrics
Research focus: Online inference for large-scale panel data
