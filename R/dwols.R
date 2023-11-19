#' Estimate parameters using dynamic weighted ordinary least squares
#'   (only appropriate for binary treatments)
#'   
#' It is expected that the dimensions of the inputs are limited to only the
#'   cases included in the stage analysis.
#'   
#' @noRd
#' @param Y A numeric vector object of length n. The outcome of interest.
#' @param A A vector object of length n. The stage treatment
#' @param data A data.frame object. The full covariate set.
#' @param wgts A numeric vector object of length n. The complete case + tx weights.
#' @param cts.obj An R6 object defining treatment + blip specific functions.
#' @param tx.var A character object. The treatment variable name.
#' 
#' @return A list object. Element `$outcome.fit` is the lm object; element
#'   `$beta` is the parameter estimates of the treatment-free component of
#'   the model; and `$psi` is the parameter estimates of the blip component.
#'   
#' @importFrom R6 is.R6
#' @importFrom stats coef lm terms
#' 
#' @include treatmentClasses.R
#' 
#' @keywords internal
.dwols <- function(Y, A, data, wgts, cts.obj, tx.var, ...) {
  
  stopifnot(
    "`Y` must be a numeric vector object" = !missing(Y) &&
      is.numeric(Y) && is.vector(Y) && length(Y) > 1L,
    "`A` must be a vector object" = !missing(A) && {is.vector(A) || is.factor(A)} && 
      length(A) == length(Y),
    "`data` must be a data.frame" = !missing(data) && is.data.frame(data) &&
      nrow(data) == length(Y),
    "`wgts` must be a numeric vector object with length equivalent to that of `Y`" = 
      !missing(wgts) &&
      is.numeric(wgts) && is.vector(wgts) && length(wgts) == length(Y),
    "`cts.obj` must be an R6 object" = !missing(cts.obj) && R6::is.R6(cts.obj),
    "`tx.var` must be a character" = !missing(tx.var) && is.character(tx.var)
  )
  
  # update the outcome model to include the internally defined outcome
  full_model <- stats::update.formula(cts.obj$full.model, Y_internal_Y ~ .)
  
  # add internally defined outcome and stage treatment to analysis data.frame
  data[["Y_internal_Y"]] <- Y
  data[[tx.var]] <- A

  # if a quadratic blip function, update the internal
  # I(a^2) variable with the square of the current treatment variable value
  if ("ContQuadraticBlip" %in% class(cts.obj)) {
    data[["l__txvar2__l"]] <- A^2
  }
  l___wgts___l <- wgts

  fit <- tryCatch(stats::lm(full_model, cbind(data, "l___wgts___l" = wgts), 
                            weights = l___wgts___l),
                  error = function(e) {
                    stop("unable to complete DWOLS\n\t",
                         e$message, call. = FALSE)
                    })

  # need blip terms to extract appropriate elements of the parameter estimates
  est <- stats::coef(fit)

  blip_vars <- cts.obj$blip_params(est)

  psi <- est[blip_vars]
  names(psi) <- names(blip_vars)
  
  list("outcome.fit" = fit,
       "beta" = est[-blip_vars],
       "psi" = psi)
}