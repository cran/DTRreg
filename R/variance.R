.finalizeCovmat <- function(covmat, blip_vars, full.cov) {
  colnames(covmat)[blip_vars] <- names(blip_vars)
  rownames(covmat)[blip_vars] <- names(blip_vars)
  if (!full.cov) covmat <- covmat[blip_vars, blip_vars, drop = FALSE]
  covmat
}

.getVariance <- function(obj, step.obj, complete.case.info, isSurvival, k) {
  
  stage_cases <- complete.case.info$last.stage >= k
  
  # Subset all vectors/matrices to stage_cases and k
  stage_data <- obj$data[stage_cases, , drop = FALSE]
  stage_data[, step.obj$tx.var] <- step.obj$A
  step_data <- step.obj$cts.obj$prep(stage_data, step.obj$A)
  
  complete.case.info$d.hat <- complete.case.info$d.hat[stage_cases, k]
  
  if (obj$var.estim == "bootstrap" &&
      obj$boot.controls$type %in% c("empirical", "normal")) {
    
    message("starting ", obj$boot.controls$type, " nonparametric bootstrap procedure")
    
    residuals <- {step.obj$Y - 
        predict(step.obj$outcome.fit, step_data)}[step.obj$delta > 0.5]
    
   step.obj <- c(step.obj, 
                  .bootstrap_nonpar(step.obj = step.obj, 
                                    obj = obj, 
                                    residuals = residuals,
                                    status.var = obj$dependent.vars$status,
                                    tx.var = obj$dependent.vars$treat[k],
                                    time.var = obj$dependent.vars$time[k],
                                    tx.wgt.man = if (is.list(obj$tx.wgt.man)) obj$tx.wgt.man[[k]] else NULL,    
                                    data = step_data,
                                    isSurvival = isSurvival,
                                    d.hat = complete.case.info$d.hat))
  }
  
  # Sandwich variance estimator
  if (obj$var.estim == "sandwich") {
    step.obj$covmat <- .sandwich(cts.obj = step.obj$cts.obj, 
                                 outcome.fit = step.obj$outcome.fit, 
                                 Y = step.obj$Y, 
                                 A = step.obj$A, 
                                 A.hat = step.obj$A.hat,
                                 wgts = step.obj$wgts, 
                                 tx.var = step.obj$tx.var, 
                                 data = step_data,
                                 treat.mod.fitted = step.obj$tx.mod.fitted,
                                 method = obj$method)
    
    blip_vars <- step.obj$cts.obj$blip_params(step.obj$outcome.fit$coefficients)
    step.obj$covmat <- .finalizeCovmat(step.obj$covmat, blip_vars, obj$full.cov)    
  }
  step.obj
}

#' @noRd
#' @keywords internal
#' @include bootstrapStep.R sandwich.R varest.R
.computeVariance <- function(obj, isSurvival = FALSE) {
  
  # Complete case weights
  complete_case_info <- obj[c("last.stage", "prob.complete.case", "d.hat", "cens.mod.fitted")]
  
  for (k in obj$K:1L) {
    obj$stages[[k]] <- .getVariance(obj, obj$stages[[k]], complete_case_info, isSurvival, k)
  }
  
  # bootstrap standard errors
  if (obj$var.estim == "bootstrap" && 
      obj$boot.controls$type == "standard") {
    obj <- c(obj, .bootstrap(obj = obj, isSurvival = isSurvival))
  }

  obj
}