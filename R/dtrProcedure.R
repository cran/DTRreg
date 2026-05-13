.getY <- function(obj, isSurvival) {
  if (isSurvival) {
    # numeric vector of 0's
    Y <- obj$data |> nrow() |> numeric()
  } else {
    # provided outcome
    Y <- obj$outcome
  }
  Y
}

.getCompleteCaseProbability <- function(obj, quiet, isSurvival) {
  ### Complete case weights
  # object will be a list containing `$last.stage`, `$prob.complete.case`, and
  # `$fitted.censor.models`.
  # `$last.stage` is an integer vector object of length n containing the 
  #   final dp for which each participant has complete data; 
  # `$prob.complete.case` is a numeric matrix object containing the  
  #   weight for each participant at each stage. Note that 
  #   prob.complete.case = NA_real_ indicates that the participant has been lost 
  #   to follow-up.
  # `$cens.mod.fitted` is a list containing NA_character_ values if censoring
  #    was not fit and/or `glm` objects.
  complete_case_info <- .completeCaseProbability(obj, quiet, isSurvival)
  
  if (!isSurvival) {
    # using cumulative probabilities when not survival
    complete_case_info$prob.complete.case <-
      apply(complete_case_info$prob.complete.case, 1L, cumprod) |> t() |>
      matrix(ncol = obj$K)
  }
  complete_case_info
}

.getTxInfo <- function(obj) {
  # Need to ensure that treatments align with user specification and
  # are in expected form (0/1, factor, numeric)

  # process stage treatment information
  # elements include A, A.hat, cts, cts.obj, and treat.mod.fitted
  # Note that we keep A here because we modify `data` to be the optimal
  # treatment as we go along; i.e., A is always the observed treatment
  .processTreatment(stage.models = obj$models, 
                    data = obj$data, 
                    tx.var = obj$dependent.vars$treat, 
                    tx.range = obj$tx.range, 
                    tx.type = obj$tx.type,
                    tx.family = obj$tx.family)
  
}

.addYdelta <- function(obj, step.obj, Y, isSurvival) {
  
  if (isSurvival) {
    # survival outcome is log(t_k + sum_{j = k+1}^{K} t_hat)
    step.obj$Y <- obj$data[, obj$dependent.vars$time]
    step.obj$Y <- log(step.obj$Y + Y)
    
    # status will be NA if interactive and no censoring indicated
    # if not interactive, the user has to specify a status that is
    # expected to contain only 1s
    if (!is.na(obj$dependent.vars$status)) {
      step.obj$delta <- obj$data[, obj$dependent.vars$status]
    } else {
      step.obj$delta <- rep(1.0, step.obj$n)
    }
  } else {
    step.obj$Y <- Y
    step.obj$delta <- rep(1.0, step.obj$n)
  }
  step.obj
}


.addTxWgt <- function(obj, step.obj) {
  
  if (obj$tx.weight == "none" || obj$method == "gest") {
    step.obj$tx.wgt <- rep(1.0, step.obj$n)
  } else if (obj$tx.weight %in% c("manual")) {
    step.obj$tx.wgt <- obj$tx.wgt.man
  } else if (obj$tx.weight %in% c("manual.with.censor")) {
    step.obj$tx.cens.wgt <- obj$tx.wgt.man
    step.obj$cens.wgt <- NULL
  } else {
    step.obj$tx.wgt <- .treatmentWeights(A = step.obj$A, 
                                         A.hat = step.obj$A.hat,
                                         tx.weight = obj$tx.weight,
                                         tx.mod.fitted = step.obj$tx.mod.fitted,
                                         cts.obj = step.obj$cts.obj,
                                         n.bins = obj$n.bins,
                                         data = obj$data)
  }
  step.obj
}

.addWgts <- function(obj, step.obj, complete.case.info) {
  
  step.obj$cens.wgt <- complete.case.info$prob.complete.case
  wgts <- step.obj$cens.wgt
  
  step.obj <- .addTxWgt(obj, step.obj)
  
  if (obj$tx.weight %in% c("manual.with.censor")) {
    wgts <- step.obj$tx.cens.wgt
  } else if (obj$tx.weight != "none" && obj$method != "gest") {
    wgts <- wgts * step.obj$tx.wgt
  }
  wgts <- wgts * step.obj$delta
  step.obj$wgts <- wgts
  step.obj
}

.performMethod <- function(obj, step.obj) {
  
  if (obj$method == "dwols") {
    delta_is_1 <- step.obj$delta > 0.5
    .dwols(Y = step.obj$Y[delta_is_1],
           A = step.obj$A[delta_is_1], 
           data = obj$data[delta_is_1, , drop = FALSE],
           wgts = step.obj$wgts[delta_is_1],
           cts.obj = step.obj$cts.obj, 
           tx.var = step.obj$tx.var)
  } else if (obj$method == "gest") {
    .gest(A = step.obj$A, 
          wgts = step.obj$wgts,
          data = obj$data,
          Y = step.obj$Y, 
          A.hat = step.obj$A.hat,
          tx.var = step.obj$tx.var,
          cts.obj = step.obj$cts.obj, 
          treat.mod.fitted = step.obj$tx.mod.fitted)
  }
}

.getStage <- function(obj, k, quiet, complete.case.info, isSurvival, Y) {
  
  if (!quiet) message("Stage ", k, " analysis")
  
  stage_cases <- complete.case.info$last.stage >= k
  
  # Subset all vectors/matrices to stage_cases and k
  obj$data <- obj$data[stage_cases, , drop = FALSE]
  obj$tx.wgt.man <- if (is.list(obj$tx.wgt.man)) obj$tx.wgt.man[[k]][stage_cases] else obj$tx.wgt.man
  obj$dependent.vars$time <- obj$dependent.vars$time[k]
  obj$dependent.vars$treat <- obj$dependent.vars$treat[k]
  obj$models <- obj$models[[k]]
  obj$tx.range <- obj$tx.range[[k]]
  
  Y <- Y[stage_cases]
  
  complete.case.info$last.stage <- complete.case.info$last.stage[stage_cases]
  complete.case.info$prob.complete.case <- complete.case.info$prob.complete.case[stage_cases, k]
  complete.case.info$d.hat <- complete.case.info$d.hat[stage_cases, k]
  complete.case.info$cens.mod.fitted <- complete.case.info$cens.mod.fitted[[k]]
  
  tx_info <- .getTxInfo(obj)
  obj$data[, obj$dependent.vars$treat] <- tx_info$A
  
  # initialize step specific result list
  step_obj <- c(
    list(dp = k,
         tx.var = obj$dependent.vars$treat,
         models = obj$models,
         n = nrow(obj$data)),
    tx_info
  )
  
  step_obj <- .addYdelta(obj, step_obj, Y, isSurvival)

  step_obj <- .addWgts(obj, step_obj, complete.case.info)
  
  step_obj <- c(step_obj, .performMethod(obj, step_obj))
  
  if (!quiet) {
    message("beta: ", paste(step_obj$beta, collapse = ", "))
    message("psi: ", paste(step_obj$psi, collapse = ", "))
  }
  
  step_obj$fitted.values <- fitted(step_obj$outcome.fit)
  step_obj$residuals <- step_obj$Y[step_obj$delta > 0.5] - step_obj$fitted.values
  step_obj$blip.data <- obj$data[step_obj$delta > 0.5, , drop = FALSE]
  
  step_obj  
}

#' Main procedure for estimating parameters across all decision points for
#'   DTRreg
#'
#' @noRd
#' @param obj A list object. Analysis settings. Expected to include
#'   K: An integer object. The number of decision points (dp); 
#'   outcome: A numeric vector object. The outcome of interest; only
#'     present if not survival
#'   models: A list of length K, each element a list of dp specific models; 
#'   data: A data.frame object. The full covariate and treatment data; 
#'   type: A character object. If "DTR", DTR estimation using regret; otherwise
#'     effect estimation using blip; 
#'   dependent.vars: A list of dependent variables. If not survival, this 
#'     includes only 
#'     treat.vars: A character vector object of length K. The treatment variable 
#'       for each dp; 
#'     if survival, also inclust status.vars and time.vars
#'   tx.range: NULL or a list of numeric vector objects of length 2. The range of
#'     allowed treatment values when tx is continuous.
#'   method: A character object. Must be one of "dwols" or "gest"; 
#'   var.estim: A character object. The covariance estimation method; must be
#'     one of "bootstrap", "sandwich", or "none".
#'   tx.weight: A character object. The type of treatment weighting.
#'   tx.wgt.man: A list of user specified treatment weights. 
#'   tx.type: A character indicating if treatment is binary, multi-nomial, or cont
#'   n.bins: An integer required when tx is continuous and weight is qpom or pom
#'   tx.family: A family object
#'   censoring.modeled: A logical object. If TRUE, element `$models` must 
#'     include censoring models.
#'   boot.controls: A list of control variables when bootstrap variance
#'     is requested.
#'   manual.censor.weight: A logical object. If TRUE, censoring weights are
#'     included in the user provided treatment weights. Provided only if
#'     procedure being called from DWsurv.
#' @param quiet A logical. If TRUE, screen prints are suppressed. Used to quiet
#'   call made by calls from bootstrap
#' @param isSurvival A logical. If TRUE analysis is for survival data.
#'   
#' @return A list object containing
#' 
#' @importFrom stats complete.cases model.matrix model.response
#' @include Ahat.R bootstrapStep.R dwols.R gest.R sandwich.R utils.R weights.R
#' @keywords internal
.dtrProcedure <- function(obj, quiet, isSurvival = FALSE) {

  # Confirm that everything is provided as expected
  obj_musts <- c("K", "dependent.vars", "models", "tx.range", 
                 "var.estim", "method", "data", "type", "tx.weight",
                 "censoring.modeled", "tx.wgt.man",
                 "tx.type", "n.bins",
                 "tx.family", "boot.controls", "full.cov")
  if (!isSurvival) {
    obj_musts <- c("outcome", obj_musts)
  } else {
    obj_musts <- c("manual.censor.weight", obj_musts)
  }

  stopifnot(
    "`obj` must be a list" = is.list(obj) && all(obj_musts %in% names(obj)),
    "`quiet` must be a logical" = is.logical(quiet) && is.vector(quiet) &&
      length(quiet) == 1L,
    "`isSurvival` must be a logical" = is.logical(isSurvival) && 
      is.vector(isSurvival) && length(isSurvival) == 1L
  )
  
  obj <- obj[obj_musts]

  # make sure that we aren't calling gest for multinomial treatment
  if (obj$method == "gest" && obj$tx.type == "multi") {
    stop("G-estimation for multinomial treatments is not yet implemented",
         call. = FALSE)
  }
  
  # get initial Y value
  Y <- .getY(obj, isSurvival)

  # Complete case weights
  complete_case_info <- .getCompleteCaseProbability(obj, quiet, isSurvival)
  
  # stages will contain the analysis results for each decision point
  # in the order of stages (NOT steps)
  obj$stages <- vector("list", length = obj$K)

  for (k in obj$K:1L) {
    
    obj$stages[[k]] <- .getStage(obj, k, quiet, complete_case_info, isSurvival, Y)
    
    # There is the possibility that the optimal treatment identified in this
    # stage will change the optimal treatment in later stages that have already
    # been evaluated. Need to re-evaluate the optimal treatment in each later 
    # stage and include regret/blip adjustment at the (potentially) new optimal

    obj$stages <- .stageTx(obj$stages, complete_case_info, obj$data, quiet)

    if (isSurvival) {
      Y <- .stageYSurvival(obj$dependent.vars$time,
                           obj$stages[k:obj$K],
                           complete.case.info = complete_case_info, 
                           data = obj$data, 
                           type = obj$type)
    } else {
      Y <- .stageY(obj$outcome, obj$stages[k:obj$K], 
                   complete.case.info = complete_case_info, 
                   data = obj$data, 
                   type = obj$type)
    }
  }
  
  obj <- c(obj, complete_case_info)
  
  obj$regret <- .stageRegret(obj$stages, complete_case_info, obj$data, obj$type)

  # Y if all optimal treatments followed
  obj$opt.Y <- Y
  
  obj$treat.vars <- obj$dependent.vars$treat
  
  obj
}
