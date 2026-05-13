#' @noRd
#' @param obj A list object. The analysis results and settings.
#' @param isSurvival A logical. If TRUE, survival end point.
#' 
#' @include dtrProcedure.R utils.R
#' @keywords internal
.resampleObj <- function(obj, boot_sample, isSurvival) {
  bobj <- obj
  if (!isSurvival) bobj$outcome <- obj$outcome[boot_sample]
  bobj$data <- obj$data[boot_sample, , drop = FALSE]
  bobj$var.estim <- "none"
  bobj$boot.controls <- list()
  if (is.list(obj$tx.wgt.man)) {
    bobj$tx.wgt.man <- lapply(obj$tx.wgt.man, "[", boot_sample)
  }
  bobj
}

#' @noRd
#' @keywords internal
.trunc <- function(x,p) {
  l <- quantile(x, p)
  u <- quantile(x, 1.0 - p)
  x[x < l] <- l
  x[x > u] <- u
  x
}

#' @noRd
#' @keywords internal
.stackPsiBoot <- function(psi_boot) {
  do.call(function(...) mapply(rbind, ..., SIMPLIFY = FALSE), 
          psi_boot)
}

.truncPsiBoot <- function(truncate, psi_boot) {
  if (truncate > 0) {
    lapply(psi_boot, FUN = function(x) {
      apply(x, 2L, .trunc, truncate)
    })
  } else {
    psi_boot
  }
}

#' @noRd
#' @keywords internal
.computeCovmat <- function(psi_boot) {
  lapply(psi_boot, var)
}

.reportProgress <- function(ptm, i, B, cont, verbose, interrupt) {
  if (!verbose || i < 10L) return(cont)
  
  eta <- {B - i} * {{proc.time() - ptm}[3L] / i}
  
  if (i > 50L && eta > 600 && cont == "y" && interrupt) {
    last_val <- attr(cont, "last")
    response <- readline(paste("Estimated run time", round(eta / 60L),
                               "minutes. Continue? y/n: "))
    if (response %in% c("n", "no", "N", "NO")) return("abort")
    cont <- "asked" 
    attr(cont, "last") <- last_val
  }
  
  if (i == 10L) attr(cont, "last") <- eta + 31.0
  last <- attr(cont, "last") %||% Inf
  
  if (eta > 30.0 && eta < {last - 30.0}) {
    message("Approximately ", round(eta), " seconds remaining.")
    attr(cont, "last") <- eta
  }
  
  cont
}

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
  last <- Inf
  cont <- "y"
  
  bobj <- obj

  for (i in seq_len(obj$boot.controls$B)) {
    boot_sample <- sample(seq_len(n_cases), obj$boot.controls$M, replace = TRUE)
    bobj <- .resampleObj(obj, boot_sample, isSurvival)

    result <- .dtrProcedure(obj = bobj, quiet = TRUE, isSurvival)
    psi_boot[[i]] <- lapply(result$stages, "[[", "psi")

    cont <- .reportProgress(ptm, i, obj$boot.controls$B, cont,
                            obj$boot.controls$verbose, obj$boot.controls$interrupt)
    if (identical(cont, "abort")) {
      message("bootstrap aborted")
      return(list("covmat" = "aborted", "psi.boot" = "aborted"))
    }
  }
  
  psi_boot <- .stackPsiBoot(psi_boot)
  
  if (any(sapply(psi_boot, function(x) any(is.na(x))))) {
    stop("psi contains NA values", call. = FALSE)
  }
  
  # truncation function for truncate option
  # if requested, truncate estimates based on specified percentile
  psi_boot <- .truncPsiBoot(obj$boot.controls$truncate, psi_boot)
  
  covmat <- .computeCovmat(psi_boot)
  
  list("covmat" = covmat, "psi.boot" = psi_boot)
}

#' @noRd
#' @keywords internal
.computeHoldY <- function(step.obj, data, isSurvival) {
  do.call(ifelse(isSurvival, exp, c),
          list(predict(step.obj$outcome.fit, data)))
}

#' @noRd
#' @keywords internal
.buildSingleStageObj <- function(obj, step.obj, tx.var, status.var, time.var, 
                                 data, tx.wgt.man, isSurvival) {
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
  single_stage_obj
}

#' @noRd
#' @keywords internal
.sampleResiduals <- function(type, reshist, mean.res, sd.res, n, n_sample) {
  if (type == "empirical") {
    bins <- with(reshist,
                  sample(length(mids), n_sample, prob = density, replace = TRUE))
    newres <- stats::runif(n, reshist$breaks[bins], reshist$breaks[bins + 1])
  } else if (type == "normal") {
    newres <- stats::rnorm(n, mean = mean.res, sd = sd.res)
  }
  newres
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

  data <- step.obj$cts.obj$prep(data, step.obj$A)
  
  hold_Y <- .computeHoldY(step.obj, data, isSurvival)
  
  # mean and sd of residuals; used to sample N(mean, sd) if bootstrap option is 
  # "normal"
  mean.res <- mean(residuals)
  sd.res <- stats::sd(residuals)
  
  # binned residuals; used to sample U() if bootstrap option is "empirical"
  reshist <- graphics::hist(residuals, plot = FALSE)
  
  # create single step version of main object
  single_stage_obj <- .buildSingleStageObj(obj, step.obj, tx.var, status.var,
                                           time.var, data, tx.wgt.man, isSurvival)
  
  for (k in seq_len(obj$boot.controls$B)) {
    # generate new status variable based on fitted probabilities
    if (isSurvival) {
      single_stage_obj$data[, status.var] <- stats::rbinom(step.obj$n, 1, d.hat)
    }
    
    n_sample <- if (isSurvival) sum(single_stage_obj$data[, status.var]) else step.obj$n
    
    newres <- .sampleResiduals(obj$boot.controls$type, reshist,
                               mean.res, sd.res, step.obj$n, n_sample)
    
    if (isSurvival) {
      single_stage_obj$data[, time.var] <- hold_Y*exp(newres)
    } else {
      single_stage_obj$outcome <- hold_Y + newres
    }

    psi_boot[, k] <- .dtrProcedure(single_stage_obj, quiet = TRUE, isSurvival)$stages[[1L]]$psi
  }
  
  psi_boot <- t(psi_boot)
  list("covmat" = stats::var(psi_boot), "psi.boot" = psi_boot)
}