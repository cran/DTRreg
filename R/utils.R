#' Extract Model Matrix
#' 
#' @noRd
#' @param model A formula object.
#' @param data A data.frame object.
#'
#' @return A matrix. The design matrix
#' 
#' @keywords internal
.getModelMatrix <- function(model, data) {
  tryCatch(stats::model.frame(model, data, na.action = "na.pass") |>
             stats::model.matrix(object = model),
           error = function(e) {
             stop("unable to extract model frame\n\t", e$message, call. = FALSE)
           })
}

#' Recalculate optimal treatment for all analyzed stages
#'
#' @noRd
#' @param steps A list object. All evaluated stage results.
#' @param complete.case.info A list object. The complete cases information,
#'   `$last.stage` and `$prob.complete.case`
#' @param data A data.frame object. The full covariate and treatment dataset.
#'   Must contain all data for all participants.
#' 
#' @return A numeric vector object.
#' @keywords internal
.stageTx <- function(steps, complete.case.info, data, quiet) {
  
  stopifnot(
    "`steps` must be a list" = is.list(steps),
    "`complete.case.info` must be a list" = is.list(complete.case.info) &&
      all(c("last.stage", "prob.complete.case") %in% names(complete.case.info)),
    "`data` must be a data.frame" = is.data.frame(data)
  )
  
  # NOTE this moves in the FORWARD direction
  for (i in seq_len(length(steps))) {
    if (is.null(steps[[i]])) next
    
    # identify the complete cases for the stage
    stage_cases <- complete.case.info$last.stage >= steps[[i]]$dp

    # estimate optimal treatment
    opt.tx <- steps[[i]]$cts.obj$opt(steps[[i]]$outcome.fit,
                                     data[stage_cases, , drop = FALSE], quiet) |> 
      drop()
    if (steps[[i]]$cts == "multi") {
      tmp <- factor(rep(NA, nrow(data)), levels = levels(opt.tx))
    } else {
      tmp <- rep(NA_real_, nrow(data))
    }
    tmp[stage_cases] <- opt.tx
    
    steps[[i]]$opt.treat <- tmp
    data[, steps[[i]]$tx.var] <- steps[[i]]$opt.treat
  }
  steps
}

#' Recalculate Y based using all model estimates
#'
#' @noRd
#' @param Y A numeric vector object of length N. The observed outcome of 
#'   interest for all participants.
#' @param steps A list object. All evaluated stage results.
#' @param complete.case.info A list object. The complete cases information,
#'   `$last.stage` and `$prob.complete.case`
#' @param data A data.frame object. The full covariate and treatment dataset.
#'   Must contain all data for all participants.
#' @param type A character object. The shift (effect/optimal).
#' @param isSurvival A logical object. TRUE if survival outcome.
#' 
#' @return A numeric vector object.
#' @keywords internal
.stageY <- function(Y, steps, complete.case.info, data, type) {
  
  stopifnot(
    "`Y` must be a numeric vector" = is.numeric(Y) && is.vector(Y),
    "`steps` must be a list" = is.list(steps),
    "`complete.case.info` must be a list" = is.list(complete.case.info) &&
      all(c("last.stage", "prob.complete.case") %in% names(complete.case.info)),
    "`data` must be a data.frame" = is.data.frame(data) && nrow(data) == length(Y),
    "`type` must be a character" = is.character(type) && is.vector(type) &&
      length(type) == 1L
  )

  # NOTE this moves in the FORWARD direction
  for (stp in steps) {
    # identify the complete cases for the stage
    stage_cases <- complete.case.info$last.stage >= stp$dp
    
    shift <- stp$cts.obj$shiftY(type = type, 
                                outcome.fit = stp$outcome.fit,
                                data = data[stage_cases, , drop = FALSE],
                                opt = stp$opt.treat[stage_cases], 
                                A = stp$A)
    Y[stage_cases] <- Y[stage_cases] + shift
  }
  Y
}

#' Recalculate Y based using all model estimates
#'
#' @noRd
#' @param times A character vector. The time variables
#' @param steps A list object. All evaluated stage results.
#' @param complete.case.info A list object. The complete cases information,
#'   `$last.stage` and `$prob.complete.case`
#' @param data A data.frame object. The full covariate and treatment dataset.
#'   Must contain all data for all participants.
#' @param type A character object. The shift (effect/optimal).
#' @param isSurvival A logical object. TRUE if survival outcome.
#' 
#' @return A numeric vector object.
#' @keywords internal
.stageYSurvival <- function(times, steps, complete.case.info, data, type) {
  
  stopifnot(
    "`times` must be a character vector" = is.character(times) && is.vector(times),
    "`steps` must be a list" = is.list(steps),
    "`complete.case.info` must be a list" = is.list(complete.case.info) &&
      all(c("last.stage", "prob.complete.case") %in% names(complete.case.info)),
    "`data` must be a data.frame" = is.data.frame(data),
    "`type` must be a character" = is.character(type) && is.vector(type) &&
      length(type) == 1L
  )
  
  Y <- numeric(nrow(data))
  
  # NOTE this moves in the FORWARD direction
  for (stp in steps) {
    # identify the complete cases for the stage
    stage_cases <- complete.case.info$last.stage >= stp$dp
    
    shift <- stp$cts.obj$shiftY(type = type, 
                                outcome.fit = stp$outcome.fit,
                                data = data[stage_cases, , drop = FALSE],
                                opt = stp$opt.treat[stage_cases], 
                                A = stp$A)
    
    Y[stage_cases] <- Y[stage_cases] + data[stage_cases, times[stp$dp]] * exp(shift)
  }
  Y
}

#' Recalculate regret using all model estimates
#'
#' @noRd
#' @param Y A numeric vector object of length N. The observed outcome of 
#'   interest for all participants.
#' @param steps A list object. All evaluated stage results.
#' @param complete.case.info A list object. The complete cases information,
#'   `$last.stage` and `$prob.complete.case`
#' @param data A data.frame object. The full covariate and treatment dataset.
#'   Must contain all data for all participants.
#' @param type A character object. The shift (effect/optimal).
#' @param isSurvival A logical object. TRUE if survival outcome.
#' 
#' @return A numeric vector object.
#' @keywords internal
.stageRegret <- function(steps, complete.case.info, data, type) {
  
  stopifnot(
    "`steps` must be a list" = is.list(steps),
    "`complete.case.info` must be a list" = is.list(complete.case.info) &&
      all(c("last.stage", "prob.complete.case") %in% names(complete.case.info)),
    "`data` must be a data.frame" = is.data.frame(data),
    "`type` must be a character" = is.character(type) && is.vector(type) &&
      length(type) == 1L
  )
  
  regret <- list()
  
  # NOTE this moves in the FORWARD direction
  for (i in 1L:length(steps)) {
    stp <- steps[[i]]
    # identify the complete cases for the stage
    stage_cases <- complete.case.info$last.stage >= stp$dp
    
    regret[[i]] <- stp$cts.obj$shiftY(type = type, 
                                      outcome.fit = stp$outcome.fit,
                                      data = data[stage_cases, , drop = FALSE],
                                      opt = stp$opt.treat[stage_cases], 
                                      A = stp$A)
    
  }
  regret
}