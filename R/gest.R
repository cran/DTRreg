#' Estimate parameters using g-estimation
#'   
#' @noRd
#' @param A A numeric vector object of length n. Must contain only 0/1 values.
#' @param wgts A numeric vector object of length n. The complete case + tx weights.
#' @param Y A numeric vector object of length n. The outcome of interest.
#' @param A.hat A numeric vector object of length n. The estimated treatment.
#' @param tx.var A character object. The treatment variable name.
#' @param cts.obj An R6 object.
#' @param treat.mod.fitted A fitted regression object.
#' @param data A data.frame object with n rows. The covariate and treatment data.
#' 
#' @return A list object. Element `$outcome.fit` a list of class GEST providing 
#'   the full outcome model and the estimated parameters; element
#'   `$beta` is the parameter estimates of the treatment-free component of
#'   the model; and `$psi` is the parameter estimates of the blip component.
#'
#' @importFrom R6 is.R6
#' @importFrom stats terms
#'   
#' @include treatmentClasses.R
#' @keywords internal
.gest <- function(A, wgts, data, Y, A.hat, tx.var, cts.obj, treat.mod.fitted, ...) {
  
  stopifnot(
    "`A` must be a numeric vector object" = !missing(A) &&
      is.numeric(A) && is.vector(A) && length(A) > 1L,
    "`wgts` must be a numeric vector object with length equivalent to that of `A`" = 
      !missing(wgts) && is.numeric(wgts) && is.vector(wgts) && 
      length(wgts) == length(A),
    "`Y` must be a numeric vector object" = !missing(Y) &&
      is.numeric(Y) && is.vector(Y) && length(Y) == length(A),
    "`data` must be a data.frame object" = !missing(data) &&
      is.data.frame(data) && nrow(data) == length(A),
    "`A.hat` must be a numeric vector object" = !missing(A.hat) &&
      is.numeric(A.hat) && is.vector(A.hat) && length(A.hat) == length(A),
    "`tx.var` must be a character" = !missing(tx.var) &&
      is.character(tx.var) && is.vector(tx.var) && length(tx.var) == 1L,
    "`cts.obj` must be an R6 object" = !missing(cts.obj) && R6::is.R6(cts.obj),
    "`treat.mod.fitted` must be NULL a fitted lm regression object" = 
      !missing(treat.mod.fitted) && 
      {is.list(treat.mod.fitted) || inherits(treat.mod.fitted, "glm")}
  )
  
  Hd <- cts.obj$Hd(data, A)
  
  Hw <- cts.obj$Hw(data = data, A = A, A.hat = A.hat, wgt = wgts, 
                   treat.mod.fitted = treat.mod.fitted)

  # get estimates
  est <- tryCatch(solve(crossprod(Hw, Hd), crossprod(Hw, Y)),
                  error = function(e) {
                    stop("unable to invert matrix\n\t", e$message,
                         call. = FALSE)
                    })
  names(est) <- colnames(Hd)
  
  # creating a fake regression object to allow for a new S3 prediction method
  outcome.fit <- list("formula" = cts.obj$full.model,
                      "fitted.values" = drop(Hd %*% est),
                      "coefficients" = est)
  class(outcome.fit) <- c("GEST", class(outcome.fit))
  
  blip_vars <- cts.obj$blip_params(est)

  psi <- est[blip_vars]
  names(psi) <- names(blip_vars)

  list("outcome.fit" = outcome.fit,
       "beta" = est[-blip_vars],
       "psi" = psi)
}

#' Internal prediction method for parameters obtained by G-estimation
#'
#' @noRd
#' @param object An S3 object. The full outcome model and estimated parameters.
#' @param newdata A data.frame object. The data on which to make predictions#
#' 
#' @return A numeric vector. The predicted outcome
#' 
#' @importFrom stats model.matrix predict
#' @keywords internal
predict.GEST <- function(object, newdata, ...) {
  mm <- stats::model.matrix(object$formula, newdata)
  mm %*% object$coefficients |> drop()
}