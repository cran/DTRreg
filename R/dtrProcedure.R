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

  obj_musts <- c("K", "dependent.vars", "models", "tx.range", 
                 "var.estim", "method", "data", "type", "tx.weight",
                 "censoring.modeled", "tx.wgt.man",
                 "tx.type", "n.bins",
                 "tx.family", "boot.controls", "full.cov")
  if (!isSurvival) obj_musts <- c("outcome", obj_musts)

  stopifnot(
    "`obj` must be a list" = is.list(obj) && all(obj_musts %in% names(obj)),
    "`quiet` must be a logical" = is.logical(quiet) && is.vector(quiet) &&
      length(quiet) == 1L,
    "`isSurvival` must be a logical" = is.logical(isSurvival) && 
      is.vector(isSurvival) && length(isSurvival) == 1L
  )
  
  obj <- obj[obj_musts]

  if (obj$method == "gest" && obj$tx.type == "multi") {
    stop("G-estimation for multinomial treatments is not yet implemented",
         call. = FALSE)
  }
  
  if (isSurvival) {
    Y <- obj$data |> nrow() |> numeric()
  } else {
    Y <- obj$outcome
  }

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
      apply(complete_case_info$prob.complete.case, 1L, cumprod) |> t()
  }
  
  # Need to ensure that treatments align with user specification and
  # are in expected form (0/1, factor, numeric)
  tx_info <- list()
  for (k in 1L:obj$K) {

    stage_cases <- complete_case_info$last.stage >= k
    stage_data <- obj$data[stage_cases, , drop = FALSE]

    # process stage treatment information
    # elements include A, A.hat, cts, cts.obj, and treat.mod.fitted
    # Note that we keep A here because we modify `data` to be the optimal
    # treatment as we go along; i.e., A is always the observed treatment
    tx_tmp <- .processTreatment(stage.models = obj$models[[k]], 
                                data = stage_data, 
                                tx.var = obj$dependent.vars$treat[k], 
                                tx.range = obj$tx.range[[k]], 
                                tx.type = obj$tx.type,
                                tx.family = obj$tx.family)

    new_A <- NULL
    if (tx_tmp$cts == "bin") {
      new_A <- rep(NA_integer_, nrow(obj$data))
    } else if (tx_tmp$cts == "multinom") {
      new_A <- factor(rep(NA, nrow(obj$data)), levels = levels(tx_tmp$A))
    } else {
      new_A <- rep(NA, nrow(obj$data))
    }
    new_A[stage_cases] <- tx_tmp$A
    obj$data[, obj$dependent.vars$treat[k]] <- new_A

    tx_info[[k]] <- tx_tmp
  }
  
  
  # stages will contain the analysis results for each decision point
  # in the order of stages (NOT steps)
  obj$stages <- vector("list", length = obj$K)

  for (k in obj$K:1L) {
    if (!quiet) message("Stage ", k, " analysis")
 
    step_obj <- list()
    
    # some info that will be handy to carry around for the backward analysis
    step_obj$dp <- k
    step_obj$tx.var <- obj$dependent.vars$treat[k]

    step_obj$models <- obj$models[[k]]

    stage_cases <- complete_case_info$last.stage >= k
    step_obj$n <- sum(stage_cases)
    
    if (isSurvival) {
      # survival outcome is log(t_k + sum_{j = k+1}^{K} t_hat)
      step_obj$Y <- obj$data[stage_cases, obj$dependent.vars$time[k]]
      step_obj$Y <- log(step_obj$Y + Y[stage_cases])
      
      # status will be NA if interactive and no censoring indicated
      # if not interactive, the user has to specify a status that is
      # expected to contain only 1s
      if (!is.na(obj$dependent.vars$status)) {
        step_obj$delta <- obj$data[stage_cases, obj$dependent.vars$status]
      } else {
        step_obj$delta <- rep(1.0, step_obj$n)
      }
    } else {
      step_obj$Y <- Y[stage_cases]
      step_obj$delta <- rep(1.0, step_obj$n)
    }
    delta_is_1 <- step_obj$delta > 0.5
    
    stage_data <- obj$data[stage_cases, , drop = FALSE]
    
    # process stage treatment information
    # elements added to step_obj include A, A.hat, cts, cts.obj, and
    # treat.mod.fitted
    step_obj <- c(step_obj, tx_info[[k]])

    step_obj$cens.wgt <- complete_case_info$prob.complete.case[stage_cases, k]
    wgts <- step_obj$cens.wgt
    
    if (obj$tx.weight == "none" || obj$method == "gest") {
      step_obj$tx.wgt <- NA
    } else if (obj$tx.weight == "manual") {
      step_obj$tx.wgt <- obj$tx.wgt.man[[k]]
      wgts <- wgts * step_obj$tx.wgt[stage_cases]
    } else {
      step_obj$tx.wgt <- rep(NA_real_, length(stage_cases))
      step_obj$tx.wgt[stage_cases] <- .treatmentWeights(A = step_obj$A, 
                                                        A.hat = step_obj$A.hat,
                                                        tx.weight = obj$tx.weight,
                                                        tx.mod.fitted = step_obj$tx.mod.fitted,
                                                        cts.obj = step_obj$cts.obj,
                                                        n.bins = obj$n.bins,
                                                        data = stage_data)
        wgts <- wgts * step_obj$tx.wgt[stage_cases]
    }
    wgts <- wgts * step_obj$delta
    step_obj$wgts <- wgts

    if (obj$method == "dwols") {
      step_obj <- c(step_obj, 
                    .dwols(Y = step_obj$Y[delta_is_1],
                           A = step_obj$A[delta_is_1], 
                           data = stage_data[delta_is_1, ],
                           wgts = wgts[delta_is_1],
                           cts.obj = step_obj$cts.obj, 
                           tx.var = step_obj$tx.var))
    } else if (obj$method == "gest") {
      step_obj <- c(step_obj,
                    .gest(A = step_obj$A, 
                          wgts = wgts,
                          data = stage_data,
                          Y = step_obj$Y, 
                          A.hat = step_obj$A.hat,
                          tx.var = step_obj$tx.var,
                          cts.obj = step_obj$cts.obj, 
                          treat.mod.fitted = step_obj$tx.mod.fitted))
    }
    
    if (!quiet) {
      message("beta: ", paste(step_obj$beta, collapse = ", "))
      message("psi: ", paste(step_obj$psi, collapse = ", "))
    }
    
    # leave procedure if parametric bootstrap requested
    # adds covmat and psi.boot to step_obj
    # quiet = TRUE only when already in bootstrap mode
    if (!quiet && obj$var.estim == "bootstrap" &&
        obj$boot.controls$type %in% c("empirical", "normal")) {
      message("starting ", obj$boot.controls$type, " nonparametric bootstrap procedure")

      if (step_obj$cts == "cts.q") {
        stage_data[["l__txvar2__l"]] <- step_obj$A^2
      }

      residuals <- {step_obj$Y - 
          predict(step_obj$outcome.fit, stage_data)}[delta_is_1]

      step_obj <- c(step_obj, 
                    .bootstrap_nonpar(step.obj = step_obj, 
                                      obj = obj, 
                                      residuals = residuals,
                                      status.var = obj$dependent.vars$status,
                                      tx.var = obj$dependent.vars$treat[k],
                                      time.var = obj$dependent.vars$time[k],
                                      tx.wgt.man = step_obj$tx.wgt.man,
                                      data = stage_data,
                                      isSurvival = isSurvival,
                                      d.hat = complete_case_info$d.hat[stage_cases,k]))
    }
    
    # Sandwich variance estimator
    if (!quiet && obj$var.estim == "sandwich") {
      step_obj$covmat <- .sandwich(cts.obj = step_obj$cts.obj, 
                                   outcome.fit = step_obj$outcome.fit, 
                                   Y = step_obj$Y, 
                                   A = step_obj$A, 
                                   A.hat = step_obj$A.hat,
                                   wgts = wgts, 
                                   tx.var = step_obj$tx.var, 
                                   data = stage_data,
                                   treat.mod.fitted = step_obj$tx.mod.fitted,
                                   method = obj$method)

      blip_vars <- step_obj$cts.obj$blip_params(step_obj$outcome.fit$coefficients)
      colnames(step_obj$covmat)[blip_vars] <- names(blip_vars)
      rownames(step_obj$covmat)[blip_vars] <- names(blip_vars)
      
      if (!obj$full.cov) {  
        step_obj$covmat <- step_obj$covmat[blip_vars, blip_vars]
      }
    }
    
    obj$stages[[k]] <- step_obj
    
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

  # To conform to original return object, shift list of individual results to 
  # individual results as lists
  # Don't want to keep 
  #   dp, tx.var, models, or cts.obj
  obj$stages <- lapply(obj$stages,
                       FUN = function(stg) {
                         stg[c("dp", "tx.var", "models", "cts.obj")] <- NULL
                         stg
                       })
  obj <- c(obj, do.call(mapply, c(obj$stages, FUN = "list", SIMPLIFY = FALSE)))
  obj$stages <- NULL

  # Y if all optimal treatments followed
  obj$opt.Y <- Y
  
  # bootstrap standard errors
  if (!quiet && 
      obj$var.estim == "bootstrap" && 
      obj$boot.controls$type == "standard") {
    obj <- c(obj,  .bootstrap(obj = obj, isSurvival = isSurvival))
  }
  
  # non-regularity
  if (!quiet & obj$cts[[1L]] == "bin") {
    obj$nonreg <- .varest(obj)
  }
  
  obj$treat.vars <- obj$dependent.vars$treat
  
  obj
}