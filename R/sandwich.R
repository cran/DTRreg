#' Sandwich estimator for UNADJUSTED estimating equation
#' 
#' @noRd
#' @param cts.obj An R6 object. The treatment type specific tools.
#' @param outcome.fit A regression object or S3 GEST used to make predictions.
#' @param Y A numeric vector object of length n. The outcome of interest.
#' @param A A numeric vector object of length n. The observed treatment. 
#'   (may be recast as 0/1)
#' @param wgts A numeric vector object of length n. The stage weights.
#' @param tx.var A character object. The treatment variable name.
#' @param data A data.frame object. The complete data for the participants
#'   currently under analysis.
#' 
#' @return A numeric matrix object {n_beta + n_psi} x {n_beta + n_psi}. The 
#'   sandwich estimator of variance.
#' 
#' @keywords internal
.sandwich <- function(cts.obj, outcome.fit, Y, A, A.hat, wgts, tx.var, data,
                      treat.mod.fitted, method) {

  stopifnot(
    "`cts.obj` must be an R6 object" = !missing(cts.obj) && R6::is.R6(cts.obj),
    "`outcome.fit` must be provided" = !missing(outcome.fit),
    "`Y` must be a numeric vector" = is.numeric(Y) && is.vector(Y),
    "`A` must be a numeric vector" = {is.vector(A) || is.factor(A)} &&
      length(A) == length(Y),
    "`wgts` must be a numeric vector" = is.numeric(wgts) && is.vector(wgts) &&
      length(wgts) == length(Y),
    "`tx.var` must be a character object" = !missing(tx.var) && 
      is.character(tx.var) && length(tx.var) == 1L,
    "`data` must be a data.frame" = !missing(data) && is.data.frame(data) &&
      nrow(data) == length(Y),
    "`method` must be a character" = !missing(method) && is.character(method)
  )

  # n x n_theta
  X_theta <- cts.obj$Hd(data, A)

  # n x n_theta
  if (method == "dwols") {
    HW <- wgts * X_theta
  } else {
    HW <- cts.obj$Hw(data, A, A.hat, wgts, treat.mod.fitted = treat.mod.fitted)
  }
  
  n <- nrow(HW)
  
  # n x n_theta
  # (Y - Yhat) wgt (A-Ahat) H_psi
  data[, tx.var] <- A
  if ("ContQuadraticBlip" %in% class(cts.obj)) {
    data[, "l__txvar2__l"] <- A^2
  }
  U <- drop(Y - predict(outcome.fit, data)) * HW

  # n_theta x n_theta
  dtheta <- {-1.0 / n} * crossprod(HW, X_theta)
  
  # derivative with respect to parameter of interest
  inv_mat <- tryCatch(solve(dtheta, t(U)),
                      error = function(e) {
                        stop("unable to invert matrix\n\t", e$message,
                             call. = FALSE)
                      })

  var <- {1.0 / n} * var(t(inv_mat))

  var
}
