#' @noRd
#' @param obj A list object. The analysis results and settings.
#' @param isSurvival A logical. If TRUE, survival end point.
#' 
#' @include dtrProcedure.R utils.R
#' @keywords internal
.bootstrap <- function(obj, isSurvival) {
  
  n_cases <- nrow(obj$data)
  
  # default # of samples in each bootstrap iteration is the number of participants
  if (is.null(obj$boot.controls$M) || obj$boot.controls$M == 0L) {
    obj$boot.controls$M <- n_cases
  }
  
  # note time for progress reports
  ptm <- proc.time()
  psi_boot <- list()
  
  # continue indicator
  cont <- "y"
  
  bobj <- obj

  for (i in seq_len(obj$boot.controls$B)) {
    boot_sample <- sample(seq_len(n_cases), obj$boot.controls$M, replace = TRUE)
    if (!isSurvival) bobj$outcome <- obj$outcome[boot_sample]
    bobj$data <- obj$data[boot_sample, , drop = FALSE]
    bobj$var.estim <- "none"
    bobj$boot.controls <- list()
    
    if (is.list(obj$tx.wgt.man)) {
      bobj$tx.wgt.man <- lapply(obj$tx.wgt.man, "[", boot_sample)
    }
    
    psi_boot[[i]] <- .dtrProcedure(obj = bobj, quiet = TRUE, isSurvival)$psi

    # if verbose, display ETA and give option to abort
    if (obj$boot.controls$verbose & i >= 10) {
      # only display if projected > 30s, then only display every 30s
      eta <- {obj$boot.controls$B - i}*{{proc.time() - ptm}[3L] / i}
      # if very long (> 10 mins) give option to abort but if ignored continue
      if (i > 50 & eta > 600 & cont != "y" & obj$boot.controls$interrupt) {
        cont <- readline(paste("Estimated run time", round(eta / 60),
                               "minutes. Continue? y/n: "))
        if (cont %in% c("n", "no", "N", "NO")) {
          message("bootstrap aborted")
          return(list("covmat" = "aborted", "psi.boot" = "aborted"))
        } else {
          cont <- "y"
        }
      }
      if (i == 10L) last <- eta + 31.0
      if (eta > 30.0 & eta < {last - 30.0}) {
        message("Approximately ", round(eta), " seconds remaining.")
        last <- eta
      }
    }
  }
  
  psi_boot <- do.call(function(...) mapply(rbind, ..., SIMPLIFY = FALSE), 
                      psi_boot)

  # truncation function for truncate option
  .trunc <- function(x,p) {
    l <- quantile(x, p)
    u <- quantile(x, 1.0 - p)
    x[x < l] <- l
    x[x > u] <- u
    x
  }
  
  # if requested, truncate estimates based on specified percentile
  if (obj$boot.controls$truncate > 0) {
    psi_boot <- lapply(psi_boot,
                       FUN = function(x) {
                               apply(x, 2L, .trunc, obj$boot.controls$truncate)
                       })
  }
  covmat <- lapply(psi_boot, var)
  list("covmat" = covmat, "psi.boot" = psi_boot)
}

#' Non-parametric Bootstrap Algorithm for Survival Data
#' 
#' @noRd
#' @param step.obj A list. The step specific models and parameter estimates
#'   as defined at the end of each stage analysis in .step1Surv(). Must contain
#'   psi, beta, A, n, D.hat, models, optimization, K, time.vars, 
#'   status.vars, treat.vars, data, type, drop)
#' @param obj A list. The analysis controls.
#' @param residuals A numeric vector. The residuals of the current model.
#' @param status.var A character object. The column header of `data` containing
#'   the 0/1 status indicator.
#' @param time.var A character object. The column header of `data` containing 
#'   the current time variable.
#' @param tx.wgt.man NULL or a vector object. User specified treatment weights
#'   for complete cases.
#' @param data A data.frame. The full data for complete cases.
#' @param isSurvival A logical. If TRUE, survival analysis.
#' @param d.hat A numeric vector. The estimated event probability.
#' 
#' @importFrom graphics hist
#' @importFrom stats rbinom rnorm runif sd var
.bootstrap_nonpar <- function(step.obj, obj, residuals, status.var, tx.var, time.var,
                              tx.wgt.man, data, isSurvival, d.hat) {
  
  # Matrix to hold psi estimates for each bootstrap sample
  psi_boot <- matrix(NA, ncol = obj$boot.controls$B, nrow = length(step.obj$psi))

  # To avoid repeated calculations of the outcome
  # Ensure that data contains the observed treatment
  data[[step.obj$tx.var]] <- step.obj$A

  if (step.obj$cts == "cts.q") {
    data[["l__txvar2__l"]] <- step.obj$A^2
  }
  
  hold_Y <- do.call(ifelse(isSurvival, exp, c),
                    list(predict(step.obj$outcome.fit, data)))
  
  # mean and sd of residuals; used to sample N(mean, sd) if bootstrap option is 
  # "normal"
  mean.res <- mean(residuals)
  sd.res <- stats::sd(residuals)
  
  # binned residuals; used to sample U() if bootstrap option is "empirical"
  reshist <- graphics::hist(residuals, plot = FALSE)
  
  # create single step version of main object
  single_stage_obj <- obj
  single_stage_obj$K <- 1L
  single_stage_obj$models <- list(step.obj$models)
  single_stage_obj$var.estim <- "none"
  single_stage_obj$dependent.vars$treat <- tx.var
  if (isSurvival) {
    single_stage_obj$dependent.vars$status <- status.var
    single_stage_obj$dependent.vars$time <- time.var
  }
  single_stage_obj$data <- data
  if (is.null(tx.wgt.man)) {
    single_stage_obj["tx.wgt.man"] <- list(NULL)
  } else {
    single_stage_obj$tx.wgt.man <- tx.wgt.man
  }
  single_stage_obj$boot.controls <- list()
  
  for (k in seq_len(obj$boot.controls$B)) {
    # generate new status variable based on fitted probabilities
    if (isSurvival) {
      single_stage_obj$data[, status.var] <- stats::rbinom(step.obj$n, 1, d.hat)
    }
      
    if (obj$boot.controls$type == "empirical") {
      bins <- with(reshist, 
                   sample(length(mids), 
                          sum(single_stage_obj$data[, status.var]), 
                          prob = density, replace = TRUE))
      newres <- stats::runif(step.obj$n, 
                             reshist$breaks[bins], reshist$breaks[bins + 1])
    } else if (obj$boot.controls$type == "normal") {
      newres <- stats::rnorm(step.obj$n, mean = mean.res, sd = sd.res)
    }
    
    if (isSurvival) {
      single_stage_obj$data[, time.var] <- hold_Y*exp(newres)
    } else {
      single_stage_obj$outcome <- hold_Y + newres
    }

    psi_boot[, k] <- .dtrProcedure(single_stage_obj, quiet = TRUE, isSurvival)$psi[[1L]]
  }
  
  psi_boot <- t(psi_boot)
  list("covmat" = stats::var(psi_boot), "psi.boot" = psi_boot)
}