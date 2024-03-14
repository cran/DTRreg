#' Censoring or Complete Case Probabilities
#' 
#' @noRd
#' @param model A formula object. The censoring or complete case model.
#' @param delta An integer object. Indicator of censored or complete case.
#' @param data A data.frame object. All covariate and treatments.
#' @param censoring.modeled A logical object. If TRUE, censoring is modeled.
#' @param dp An integer. The stage under consideration.
#' @param quiet A logical. If TRUE messaging is suppressed.
#' 
#' @return A list containing the fitted model `$cens.model.fitted` (if estimated),
#'   censoring weight `$wt.cen`, and estimated probabilities `$D.hat`.
#'   
#' @importFrom stats fitted glm update.formula
#' @keywords internal
.processCensoring <- function(model, delta, data, censoring.modeled, dp, quiet) {
  res <- list()
  
  res$cens.model.fitted <- NA

  # censoring model: binary logistic
  if (any(delta == 0L) && censoring.modeled) {
    model <- stats::update.formula(l__delta__l ~ ., model)
    res$cens.model.fitted <- stats::glm(model, cbind(data, "l__delta__l" = delta), 
                                        family = "binomial")
    d.hat <- stats::fitted(res$cens.model.fitted)
  } else if (all(delta == 1L) && censoring.modeled) {
    
    d.hat <- delta
    if (!quiet) {
      message("no censoring in stage", dp, ", censoring model will be ignored.")
    }
    
  } else {
    d.hat <- delta
  }
  res$D.hat <- d.hat
  res$wt.cen <- 1.0 / {delta * d.hat + {1.0 - delta} * {1.0 - d.hat}}

  res
  
}

#' Calculate weights for censoring or complete cases for all decision points.
#'
#' @noRd
#' @param obj A list object. Analysis settings. Must include
#'   K: An integer object. The number of decision points (dp); 
#'   outcome: A numeric vector object. The outcome of interest;
#'   models: A list of length K, each element a list of dp specific models; 
#'   data: A data.frame object. The full covariate and treatment data; 
#'   dependent.vars: A list of dependent variables;
#'   censoring.modeled: A logical. If TRUE, element `$models` must include
#'     censoring models
#'   manual.censor.weight: A logical. If TRUE, censoring weights were provided
#'     by the user through treat.wgt.man.
#' @param quiet A logical. If TRUE messaging is suppressed.
#'  
#' @return A list object containing `$last.stage`, an integer vector of length
#'   n indicating the last stage for which each participant had complete data;
#'   `$prob.complete.case`, an n x K matrix of inverse weights; and
#'   `$cens.mod.fitted` a list of `glm` objects and/or NA_character_ (indicating
#'   that a fit was not performed for the stage).
#'   
#' @importFrom stats complete.cases
#' @include utils.R
.completeCaseProbability <- function(obj, quiet, isSurvival) {
  rqrd_elements <- c("K", "models", "data", "censoring.modeled",
                     "dependent.vars")
  if (isSurvival) rqrd_elements <- c(rqrd_elements, "manual.censor.weight")
  
  stopifnot(
    "`obj` does not contain the required information" = 
      is.list(obj) && all(rqrd_elements %in% names(obj))
  )

  n_cases <- nrow(obj$data)
  
  # this vector will indicate the last stage at which each participant
  # has complete data
  last_stage <- integer(n_cases)
  
  # this matrix will contain the inverse probabilities
  # for each decision point. "missing" cases will be NA
  prob_complete <- matrix(NA_real_, nrow = n_cases, ncol = obj$K)
  d.hat <- matrix(NA_real_, nrow = n_cases, ncol = obj$K)
  
  # temporary vector to track cases that are still in the analysis at the
  # stage under consideration
  still_in <- rep(TRUE, n_cases)
  
  # if complete cases or censoring is not modeled, we need to remove all data 
  # for cases with missing data. Do not want to require complete data when 
  # censoring is provided through the manual treatments
  if ({!obj$censoring.modeled && !isSurvival} ||
      {!obj$censoring.modeled && isSurvival && !obj$manual.censor.weight}) {
    still_in <- complete.cases(obj$data)
  }

  fitted_censor_models <- list()
  
  # NOTE this is in the FORWARD direction
  # This ensures that once someone has missing data that are removed from 
  # all subsequent stages
  for (k in 1L:obj$K) {
    # Pull treatment and possibly status and time variables
    dependent.vars <- NULL
    for (dv in obj$dependent.vars) {
      dependent.vars <- cbind(dependent.vars,
                              tryCatch(obj$data[, dv[min(k, length(dv))]],
                                       error = function(e) {
                                         stop("unable to extract variable ", 
                                              dv[k], "\n\t", e$message, 
                                              call. = FALSE)
                                       }))
    }

    model_matrices <- lapply(obj$models[[k]], .getModelMatrix, data = obj$data)

    # indicator is 1 if all stage covariates are present and the participant is
    # still in the cohort
    cc <- stats::complete.cases(dependent.vars, model_matrices) & still_in

    # for the participants that remain in the cohort (no missing data yet)
    # change last_stage to indicate their inclusion in this stage
    last_stage[cc] <- k

    if (!isSurvival) {
      delta <- rep(1L, n_cases)
      delta[!cc & still_in] <- 0L
      delta <- delta[still_in]
    } else {
      still_in <- still_in & cc
      delta <- obj$data[still_in, obj$dependent.vars$status]
    }

    cen_fit <- .processCensoring(obj$models[[k]]$cens, 
                                 delta = delta, 
                                 data = obj$data[still_in, , drop = FALSE], 
                                 censoring.modeled = obj$censoring.modeled, 
                                 dp = k,
                                 quiet = quiet)

    prob_complete[still_in, k] <- cen_fit$wt.cen
    d.hat[still_in, k] <- cen_fit$D.hat
    fitted_censor_models[[k]] <- cen_fit$cens.model.fitted

    # for those that are lost at or before this stage, set as NA
    prob_complete[!cc, k] <- NA_real_
    
    # adjust cases that are still under consideration
    still_in <- still_in & cc
  }
  
  list("last.stage" = last_stage, 
       "prob.complete.case" = prob_complete,
       "cens.mod.fitted" = fitted_censor_models,
       "d.hat" = d.hat)
}

#' Calculate treatment weights
#' 
#' @noRd
#' @param A A numeric or factor vector object. The observed treatments.
#' @param Ahat A numeric vector or matrix object of length/rows equivalent to `A`. 
#'   The estimated treatment probabilities.
#' @param tx.weight A character object. If "iptw", the inverse probability of
#'   treatment is used; others, abs(A-Ahat).
#' @param cts.obj An R6 object. The treatment-type specific definitions.
#' @param n.bins An integer object. The number of bins (levels) to be used for 
#'   categorizing continuous doses. This input is required only when
#'   `treat.type` = "cont" and `weight` = "wo" or `weight` = "qpom".
#' 
#' @return A numeric vector object of length equivalent to input `A`
#'
#' @import R6
#'
#' @include treatmentClasses.R
#'
#' @keywords internal
.treatmentWeights <- function(A, A.hat, tx.weight, tx.mod.fitted, cts.obj, n.bins, data) {

  stopifnot(
    "`A` must be a numeric or factor vector" = !missing(A) &&
      {{is.numeric(A) && is.vector(A)} || is.factor(A)} && length(A) > 1L,
    "`A.hat` must be a numeric vector or matrix with same length/nrow as `A`" = 
      !missing(A.hat) &&
      is.numeric(A.hat) && {{is.vector(A.hat) && length(A.hat) == length(A)} ||
          {is.matrix(A.hat) && nrow(A.hat) == length(A)}},
    "`tx.weight` must be a character object" = !missing(tx.weight) &&
      is.character(tx.weight) && is.vector(tx.weight) && length(tx.weight) == 1L,
    "`tx.mod.fitted` must be provided" = !missing(tx.mod.fitted),
    "`cts.obj` must be an R6 object" = !missing(cts.obj) && R6::is.R6(cts.obj),
    "`n.bins` must be an integer" = !missing(n.bins) && is.integer(n.bins)
  )    
  if ("MultiNom" %in% class(cts.obj)) {
    n.bins <- length(levels(A))
  }

  if (tx.weight == "ipw") {
    tx_wgt <- cts.obj$ipw(A, A.hat = A.hat, tx.mod.fitted = tx.mod.fitted)
  } else if (tx.weight == "cipw") {
    weights <- cts.obj$ipw(A, A.hat = A.hat, tx.mod.fitted = tx.mod.fitted)
    cap <- quantile(weights, 0.99)
    tx_wgt <- pmin(weights, cap)
  } else if (tx.weight == "qpom") {
    tx_wgt <- q.pom(A, tx.mod.fitted, data, n.bins, cts.obj)
  } else if (tx.weight == "wo") {
    tx_wgt <- bin.o.fcn(A, tx.mod.fitted, data, n.bins, cts.obj)
  } else if (tx.weight == "abs") {
    if ("MuliNom" %in% class(cts.obj)) {
      stop("`abs` weight is not available for multinomial treatments",
           call. = FALSE)
    }
    tx_wgt <- abs(A - A.hat)
  } else {
    tx_wgt <- rep(1.0, length(A))
  }

  tx_wgt
}




#########################################################
### Different weight functions for two stage function ###
#########################################################

# different weight functions
# @noRd
# quantile binning based on POM
q.pom <- function(A, tx.mod.fitted, data, m, cts.obj) {
  pom <- cts.obj$.pom(A, tx.mod.fitted, data, m)
  {1.0 / m} / pom$prob
}

# @noRd
# overlap weights
bin.o.fcn <- function(A, tx.mod.fitted, data, m, cts.obj) {
  pom <- cts.obj$.pom(A, tx.mod.fitted, data, m)
  {1.0 / pom$prob} / pom$sum.prob
  
}