# STH SHOULD THIS USE OPTIMAL TREATMENT OR OBSERVED?

#' Estimate the confidence intervals
#' 
#' @noRd
#' @param obj A list object. Must contain `$K`, `$models`, `$data`, `$covmat`, 
#'   `$psi`, `$tx.type`
#'   
#' @return A numeric vector object of length `obj$K`
#' 
#' @keywords internal
.varest <- function(obj) {
  
  if (is.null(obj$covmat)) return(NA)
  if (obj$tx.type != "bin") return(NA)
  
  .each_stage <- function(models, data, covmat, psi) {
    
    H_psi <- stats::model.matrix(models$blip, data)
    H_psi <- H_psi[stats::complete.cases(H_psi), , drop = FALSE]

    # lower/upper limits on psi
    tmp <- covmat |> diag() |> sqrt()
    psi_l <- psi - 1.96 * tmp
    psi_u <- psi + 1.96 * tmp

    # look at max/min value of blip based on sign of covariates and max/min of 
    # parameter CIs
    H_psi_p <- pmax(H_psi, 0.0)
    H_psi_n <- pmin(H_psi, 0.0)

    # lower blip is positive values * lower CI + negative values * upper CI
    blip_l <- H_psi_p %*% psi_l + H_psi_n %*% psi_u
    blip_u <- H_psi_p %*% psi_u + H_psi_n %*% psi_l
    
    sum(blip_l < 0 & blip_u > 0) / nrow(H_psi)
  }
    
  res <- numeric(obj$K)
  
  for (j in obj$K:1L) {
    res[j] <- .each_stage(obj$models[[j]], obj$data, obj$covmat[[j]], obj$psi[[j]])
  }
  res
}