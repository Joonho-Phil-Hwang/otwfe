# =============================================================================
# Section 6: S3 Methods — coef / vcov / nobs / confint / tidy / print / summary
# =============================================================================

# --------------------------------------------------------------------------
# Internal helper: build coefficient table for x_cols
# --------------------------------------------------------------------------
.otwfe_coef_table <- function(object, type = c("robust", "classical"),
                               conf_level = 0.95) {
  type  <- match.arg(type)
  S     <- object$state
  theta <- S$theta_hat[object$x_cols]
  V     <- if (type == "classical") {
    (S$sigma2_hat * S$inv_dotZtZ)[object$x_cols, object$x_cols]
  } else {
    S$Vcr_hat[object$x_cols, object$x_cols]
  }
  se   <- sqrt(pmax(diag(V), 0))
  df   <- S$n - S$N - S$p
  tval <- theta / se
  pval <- 2 * pt(abs(tval), df = df, lower.tail = FALSE)
  q    <- qt((1 + conf_level) / 2, df = df)
  data.frame(
    term      = names(theta),
    estimate  = unname(theta),
    std.error = unname(se),
    statistic = unname(tval),
    p.value   = unname(pval),
    conf.low  = unname(theta - q * se),
    conf.high = unname(theta + q * se),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

# Significance stars helper
.sig_stars <- function(p) {
  ifelse(p < 0.001, "***",
  ifelse(p < 0.01,  "** ",
  ifelse(p < 0.05,  "*  ",
  ifelse(p < 0.1,   ".  ", "   "))))
}

# --------------------------------------------------------------------------
# coef.otwfe
# --------------------------------------------------------------------------
#' Extract coefficients from an otwfe object
#'
#' @param object An \code{otwfe} object
#' @param which  \code{"x"} (default) returns only x covariates;
#'               \code{"all"} also includes time FE dummies
#' @export
coef.otwfe <- function(object, which = c("x", "all"), ...) {
  which <- match.arg(which)
  theta <- object$state$theta_hat
  if (which == "x") theta[object$x_cols] else theta
}

# --------------------------------------------------------------------------
# vcov.otwfe
# --------------------------------------------------------------------------
#' Extract variance-covariance matrix from an otwfe object
#'
#' @param object An \code{otwfe} object
#' @param type   \code{"robust"} (default) for cluster-robust (HC0, Arellano);
#'               \code{"classical"} for homoskedastic
#' @export
vcov.otwfe <- function(object, type = c("robust", "classical"), ...) {
  type <- match.arg(type)
  S    <- object$state
  V    <- if (type == "classical") S$sigma2_hat * S$inv_dotZtZ else S$Vcr_hat
  V[object$x_cols, object$x_cols]
}

# --------------------------------------------------------------------------
# nobs.otwfe
# --------------------------------------------------------------------------
#' Number of observations in an otwfe object
#' @export
nobs.otwfe <- function(object, ...) object$state$n

# --------------------------------------------------------------------------
# confint.otwfe
# --------------------------------------------------------------------------
#' Confidence intervals for otwfe coefficients
#'
#' @param object     An \code{otwfe} object
#' @param parm       Parameter names to include (default: all x covariates)
#' @param level      Confidence level (default 0.95)
#' @param type       \code{"robust"} (default) or \code{"classical"}
#' @export
confint.otwfe <- function(object, parm, level = 0.95,
                           type = c("robust", "classical"), ...) {
  type  <- match.arg(type)
  tab   <- .otwfe_coef_table(object, type = type, conf_level = level)
  ci    <- as.matrix(tab[, c("conf.low", "conf.high")])
  rownames(ci) <- tab$term
  alpha        <- 1 - level
  colnames(ci) <- c(sprintf("%.1f %%", 100 * alpha / 2),
                    sprintf("%.1f %%", 100 * (1 - alpha / 2)))
  if (!missing(parm)) ci[parm, , drop = FALSE] else ci
}

# --------------------------------------------------------------------------
# tidy.otwfe  (broom-compatible)
# --------------------------------------------------------------------------
#' Tidy method for otwfe objects
#'
#' @param x              An \code{otwfe} object
#' @param vcov           \code{"robust"} (default) or \code{"classical"}
#' @param conf_level     Confidence level (default 0.95)
#' @param include_time_fe Include time FE dummies in output (default FALSE)
#' @export
tidy.otwfe <- function(x, vcov = c("robust", "classical"),
                        conf_level = 0.95, include_time_fe = FALSE, ...) {
  vcov <- match.arg(vcov)
  tab  <- .otwfe_coef_table(x, type = vcov, conf_level = conf_level)

  if (include_time_fe) {
    S     <- x$state
    theta <- S$theta_hat
    fe_nm <- setdiff(names(theta), x$x_cols)
    if (length(fe_nm) > 0L) {
      V  <- if (vcov == "classical") S$sigma2_hat * S$inv_dotZtZ else S$Vcr_hat
      df <- S$n - S$N - S$p
      q  <- qt((1 + conf_level) / 2, df = df)
      se_fe <- sqrt(pmax(diag(V)[fe_nm], 0))
      tv_fe <- theta[fe_nm] / se_fe
      pv_fe <- 2 * pt(abs(tv_fe), df = df, lower.tail = FALSE)
      fe_tab <- data.frame(
        term      = fe_nm,
        estimate  = unname(theta[fe_nm]),
        std.error = unname(se_fe),
        statistic = unname(tv_fe),
        p.value   = unname(pv_fe),
        conf.low  = unname(theta[fe_nm] - q * se_fe),
        conf.high = unname(theta[fe_nm] + q * se_fe),
        row.names = NULL, stringsAsFactors = FALSE
      )
      tab <- rbind(tab, fe_tab)
    }
  }
  tab
}

# --------------------------------------------------------------------------
# print.otwfe
# --------------------------------------------------------------------------
#' Print method for otwfe objects
#' @export
print.otwfe <- function(x, digits = 4, ...) {
  S   <- x$state
  tab <- .otwfe_coef_table(x, type = "robust")

  cat("Online TWFE Estimator  (Hwang & Lee, 2026)\n")
  cat(sprintf("  N = %-12s  n = %-12s  T = %d  k = %d  p = %d\n",
              format(S$N, big.mark = ","),
              format(S$n, big.mark = ","),
              S$T_support, length(x$x_cols), S$p))
  cat("\nCoefficients (cluster-robust SE):\n")

  fmt <- paste0("%.", digits, "f")
  coef_df <- data.frame(
    Estimate   = sprintf(fmt, tab$estimate),
    `Std. Err` = sprintf(fmt, tab$std.error),
    `t value`  = sprintf("%.2f",  tab$statistic),
    `Pr(>|t|)` = format.pval(tab$p.value, digits = 3),
    ` `        = .sig_stars(tab$p.value),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  rownames(coef_df) <- tab$term
  print(coef_df)

  cat("---\nSignif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  cat(sprintf("\nResidual sigma^2: %.4f   df: %s\n",
              S$sigma2_hat, format(S$n - S$N - S$p, big.mark = ",")))
  invisible(x)
}

# --------------------------------------------------------------------------
# summary.otwfe  — returns a summary.otwfe object
# --------------------------------------------------------------------------
#' Summary method for otwfe objects
#'
#' @param object     An \code{otwfe} object
#' @param vcov_type  \code{"robust"} (default) or \code{"classical"}
#' @param conf_level Confidence level for intervals (default 0.95)
#' @export
summary.otwfe <- function(object, vcov_type = c("robust", "classical"),
                           conf_level = 0.95, ...) {
  vcov_type <- match.arg(vcov_type)
  tab       <- .otwfe_coef_table(object, type = vcov_type, conf_level = conf_level)
  structure(
    list(
      vcov_type    = vcov_type,
      conf_level   = conf_level,
      coefficients = tab,
      state        = object$state,
      x_cols       = object$x_cols,
      y_col        = object$y_col,
      id_col       = object$id_col,
      time_col     = object$time_col,
      time_levels  = object$time_levels
    ),
    class = "summary.otwfe"
  )
}

# --------------------------------------------------------------------------
# print.summary.otwfe
# --------------------------------------------------------------------------
#' @export
print.summary.otwfe <- function(x, digits = 4, ...) {
  S   <- x$state
  tab <- x$coefficients
  df  <- S$n - S$N - S$p

  sep56 <- strrep("=", 56)
  cat(sep56, "\n")
  cat("Online TWFE Estimator  (Hwang & Lee, 2026)\n")
  cat(sep56, "\n")
  cat(sprintf("Dependent variable : %s\n", x$y_col))
  cat(sprintf("Covariates         : %s\n", paste(x$x_cols, collapse = ", ")))
  cat(sprintf("ID column          : %s\n", x$id_col))
  cat(sprintf("Time column        : %s\n", x$time_col))
  if (!is.null(x$time_levels) && length(x$time_levels) > 1L)
    cat(sprintf("Time periods       : %s - %s  (T = %d)\n",
                x$time_levels[1L], x$time_levels[length(x$time_levels)],
                S$T_support))
  cat(sprintf("N (units)          : %s\n",  format(S$N, big.mark = ",")))
  cat(sprintf("n (observations)   : %s\n",  format(S$n, big.mark = ",")))
  cat(sprintf("Parameters (p)     : %d\n",  S$p))
  cat(sprintf("Residual df        : %s\n",  format(df, big.mark = ",")))
  cat(sprintf("sigma^2            : %.6f\n", S$sigma2_hat))
  cat(sprintf("Variance type      : %s\n",
              if (x$vcov_type == "robust")
                "Cluster-robust (HC0, Arellano)"
              else "Classical (homoskedastic)"))
  cat(strrep("-", 56), "\n")

  # CI column labels
  alpha <- 1 - x$conf_level
  ci_lo <- sprintf("%.1f %%", 100 * alpha / 2)
  ci_hi <- sprintf("%.1f %%", 100 * (1 - alpha / 2))

  fmt <- paste0("%.", digits, "f")
  coef_df <- data.frame(
    Estimate   = sprintf(fmt, tab$estimate),
    `Std. Err` = sprintf(fmt, tab$std.error),
    `t value`  = sprintf("%.2f", tab$statistic),
    `Pr(>|t|)` = format.pval(tab$p.value, digits = 3),
    ` `        = .sig_stars(tab$p.value),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  coef_df[[ci_lo]] <- sprintf(fmt, tab$conf.low)
  coef_df[[ci_hi]] <- sprintf(fmt, tab$conf.high)
  rownames(coef_df) <- tab$term
  print(coef_df)

  cat("---\n")
  cat("Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  cat(sep56, "\n")
  invisible(x)
}
