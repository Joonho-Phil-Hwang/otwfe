setwd("/Users/joonhohwang/Desktop/claude_code/Simulation")
source("sim_dgp.R")
source("../online_twfe_core.R")
suppressPackageStartupMessages(library(plm))

# sim_dgp 구조 확인
cat("=== sim_dgp.R 검증 ===\n")
gen <- generate_panel_52(N=1000L, T=5L, seed=42L, verbose=TRUE)
cat(sprintf("beta_true: %s\n", paste(gen$beta_true, collapse=", ")))
cat(sprintf("T_i 분포: min=%d, max=%d, mean=%.2f\n",
            min(gen$T_i_vec), max(gen$T_i_vec), mean(gen$T_i_vec)))
cat(sprintf("등장 calendar time: %s\n",
            paste(sort(unique(gen$df$time)), collapse=",")))
cat(sprintf("총 관측치: %d\n\n", gen$N_obs))

# run_plm_steps 함수 (benchmark_vs_plm.R에서 동일하게 사용)
run_plm_steps <- function(df, x_cols, timeout) {
  res <- list(
    pdata_ok=FALSE, plm_ok=FALSE, vcov_ok=FALSE, vcovhc_ok=FALSE,
    t_prep=NA_real_, t_coef=NA_real_, t_v0=NA_real_, t_vcr=NA_real_,
    fail_step=NA_character_, fail_msg=NA_character_,
    theta=NULL, V0=NULL, Vcr=NULL
  )
  t1  <- proc.time()
  pdf <- tryCatch(
    pdata.frame(df, index=c("id","time")),
    error=function(e){ res$fail_step<<-"pdata.frame()"; res$fail_msg<<-conditionMessage(e); NULL }
  )
  res$t_prep <- (proc.time()-t1)[["elapsed"]]
  if (is.null(pdf)) return(res)
  res$pdata_ok <- TRUE

  fml <- as.formula(paste("y ~", paste(x_cols,collapse="+"), "+ factor(time)"))
  t2  <- proc.time()
  plm_fit <- tryCatch(
    plm(fml, data=pdf, model="within", effect="individual"),
    error=function(e){ res$fail_step<<-"plm()"; res$fail_msg<<-conditionMessage(e); NULL }
  )
  res$t_coef <- (proc.time()-t2)[["elapsed"]]
  if (is.null(plm_fit)) return(res)
  res$plm_ok <- TRUE
  res$theta  <- coef(plm_fit)[x_cols]

  t3 <- proc.time()
  V0 <- tryCatch(
    vcov(plm_fit)[x_cols,x_cols],
    error=function(e){ res$fail_step<<-"vcov()"; res$fail_msg<<-conditionMessage(e); NULL }
  )
  res$t_v0 <- (proc.time()-t3)[["elapsed"]]
  if (!is.null(V0)) { res$vcov_ok <- TRUE; res$V0 <- V0 } else return(res)

  t4  <- proc.time()
  Vcr <- tryCatch({
    setTimeLimit(elapsed=timeout, transient=TRUE)
    out <- vcovHC(plm_fit, method="arellano", type="HC0")[x_cols,x_cols]
    setTimeLimit(elapsed=Inf, transient=TRUE)
    out
  }, error=function(e){
    setTimeLimit(elapsed=Inf, transient=TRUE)
    msg <- conditionMessage(e)
    if (grepl("elapsed time limit", msg, ignore.case=TRUE)){
      res$fail_step <<- sprintf("vcovHC() TIMEOUT (>%ds)", timeout)
      res$fail_msg  <<- msg
    } else {
      res$fail_step <<- "vcovHC()"
      res$fail_msg  <<- msg
    }
    NULL
  })
  res$t_vcr <- (proc.time()-t4)[["elapsed"]]
  if (!is.null(Vcr)) { res$vcovhc_ok <- TRUE; res$Vcr <- Vcr }
  res
}

# N=1000 정확도 검증
cat("=== N=1000 otwfe vs plm 정확도 검증 ===\n")
gen2 <- generate_panel_52(N=1000L, T=5L, seed=1L, verbose=FALSE)
df   <- gen2$df
x_cols <- c("x1","x2")

fit     <- otwfe(data=df, id_col="id", time_col="time", y_col="y",
                 x_cols=x_cols, track_all_units=FALSE, verbose=FALSE)
plm_res <- run_plm_steps(df, x_cols, timeout=60L)

cat(sprintf("plm 단계별 성공: pdata=%s  plm=%s  vcov=%s  vcovhc=%s\n",
            plm_res$pdata_ok, plm_res$plm_ok, plm_res$vcov_ok, plm_res$vcovhc_ok))
cat(sprintf("theta max|diff| : %.3e\n",
            max(abs(fit$state$theta_hat[x_cols] - plm_res$theta))))
V0_on <- fit$state$sigma2_hat * fit$state$inv_dotZtZ[x_cols, x_cols]
cat(sprintf("V0    max|diff| : %.3e\n", max(abs(V0_on - plm_res$V0))))
cat(sprintf("Vcr   max|diff| : %.3e\n",
            max(abs(fit$state$Vcr_hat[x_cols,x_cols] - plm_res$Vcr))))
cat("[OK] 검증 완료\n")
