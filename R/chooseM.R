#' Adaptive Choice of the Bootstrap Resample Size M
#' 
#' 
#' Implementation of a double-bootstrap algorithm for choosing the bootstrap 
#'   resample size \eqn{m}{m} in a data-adaptive manner. The function returns the 
#'   resample size to be used to apply the m-out-of-n bootstrap with \link{DTRreg}.
#'
#' The m-out-of-n bootstrap is an adequate tool for constructing valid 
#'   confidence intervals for the first stage parameters in \link{DTRreg}. The 
#'   resample size \eqn{m}{m} is: 
#'   \eqn{m = n^{\frac{1 + \alpha(1-\hat{p})}{1+\alpha}}}{%
#'   m = n^[1+alpha(1-pHat)/(1+alpha)]}. The estimated non-regularity level is 
#'   computed by \link{DTRreg}. The double-bootstrap algorithm is a cross-validation 
#'   tool for choosing the tuning parameter \eqn{\alpha}{alpha} in a data-driven way.
#'  
#' The current implementation is valid for a two-stage DTR. Moreover, the 
#'   current implementation may be unstable when there are many missing data.
#'   
#' @inheritParams DTRreg
#' @param B1 Number of first-level bootstrap resamples.
#' @param B2 Number of second-level bootstrap resamples.
#' 
#' @return A list with a single element 
#'   \item{m }{Resample size for using in the m-out-of-n bootstrap.}
#'   
#' @references
#'   Chakraborty, B., Moodie, E. E. M. (2013) \emph{Statistical Methods for 
#'   Dynamic Treatment Regimes}. New York: Springer.
#'   
#'   Efron B., Tibshirani R. J. (1994) An Introduction to the Bootstrap. 
#'   \emph{CRC press}.
#'   
#'   Wallace, M. P., Moodie, E. M. (2015) Doubly-Robust Dynamic Treatment 
#'   Regimen Estimation Via Weighted Least Squares. \emph{Biometrics} 
#'   \bold{71}(3), 636--644 (doi:10.1111/biom.12306.)
#' % Eventually update with reference to my paper
#' @author Gabrielle Simoneau
#' 
#' @examples 
#' data(twoStageCont)
#' @template model_definitions
#' @examples
#' 
#' # perform dWOLS without calculating confidence intervals
#' mod1 <- DTRreg(twoStageCont$Y, blip.mod, treat.mod, tf.mod, 
#'                data = twoStageCont, method = "dwols")
#'   
#' # choose m adaptively for that model
#' \dontrun{
#'   m <- chooseM(twoStageCont$Y, blip.mod, treat.mod, tf.mod, 
#'                data = twoStageCont, method = "dwols",
#'                B1 = 200, B2 = 200)$m
#' }
#' m <- 94
#'   
#' # dWOLS with confidence intervals from the m-out-of-n bootstrap
#' mod2 <- DTRreg(twoStageCont$Y, blip.mod, treat.mod, tf.mod, 
#'                data = twoStageCont, method = "dwols",
#'                var.estim = "bootstrap", 
#'                bootstrap.controls = list(M = m))
#' @concept dynamic treatment regimens
#' @concept adaptive treatment strategies
#' @concept personalized medicine
#' @concept g-estimation
#' @concept dynamic weighted ordinary least squares
#' @concept m-out-of-n bootstrap
#' @concept double bootstrap
#' @export
chooseM <- function(outcome, 
                    blip.mod, treat.mod, tf.mod, 
                    data = NULL, 
                    method = c("gest", "dwols", "qlearn"), 
                    treat.type = c("bin", "multi", "cont"),
                    treat.fam = gaussian(link = "identity"),
                    weight = c("abs", "ipw", "cipw", "qpom", "wo", "none", "manual"), 
                    n.bins = 3L, 
                    treat.wgt.man = NULL,
                    treat.range = NULL, 
                    missing = c("drop", "ipw"), 
                    missing.mod = NULL,
                    B1 = 500, B2 = 500) {
  
  stopifnot(
    "`outcome`, `blip.mod`, `treat.mod`, and `tf.mod` must be provided" =
      !missing(outcome) && !missing(blip.mod) && !missing(treat.mod) &&
      !missing(tf.mod)
  )
  
  Y <- outcome
  
  # fit initial DTRreg model to the data to get original estimates + phat
  mod1 <- DTRreg(outcome = outcome, 
                 blip.mod = blip.mod, 
                 treat.mod = treat.mod, 
                 tf.mod = tf.mod, 
                 data = data, 
                 method = method, 
                 weight = weight, 
                 n.bins = n.bins,
                 missing = missing, 
                 missing.mod = missing.mod,
                 treat.type = treat.type,
                 treat.fam = treat.fam,
                 treat.wgt.man = treat.wgt.man, 
                 treat.range = treat.range,
                 var.estim = "bootstrap", 
                 bootstrap.controls = list(M = 0, truncate = 0.0, B = 50))
  
  blip.psi.1 <- mod1$psi[[1L]]
  p <- mod1$nonreg[2L]
  
  # deal with missing data here
  
  # initialize
  alpha <- 0.025
  n <- nrow(mod1$training_data$data)
  
  # loop until alpha reaches 0.5, with break when the coverage reaches the 
  # nominal 95% rate
  while (alpha <= 0.5) { 
    # reset coverage
    coverage <- 0
    
    # reset matrix to save estimates
    est <- matrix(NA, ncol = B2 + 3L, nrow = 1L)
    
    # loop over B1 first stage bootstrap samples
    for (j in seq_len(B1)) {
      
      # draw a n-out-of-n bootstrap sample
      index <- sample(seq_len(n), n, replace = TRUE)
      boot1 <- mod1$training_data$data[index, , drop = FALSE]
      Y <- outcome[index]
      
      if (is.list(treat.wgt.man)) {
        tx.wgt.man <- lapply(treat.wgt.man, "[", index)
      } else {
        tx.wgt.man <- NULL
      }
      
      # fit the model to b1-th bootstrap sample
      res1 <- tryCatch(DTRreg(outcome = Y, 
                              blip.mod = blip.mod, 
                              treat.mod = treat.mod, 
                              tf.mod = tf.mod, 
                              data = boot1, 
                              method = method, 
                              weight = weight, 
                              n.bins = n.bins,
                              missing = missing, 
                              missing.mod = missing.mod,
                              treat.type = treat.type,
                              treat.fam = treat.fam,
                              treat.wgt.man = tx.wgt.man, 
                              treat.range = treat.range,
                              var.estim = "bootstrap"),
                       error = function(e) {
                         stop("unable to complete DTRreg\n\t",
                              e$message, call. = FALSE)
                       })
      
      esb1 <- res1$psi[[1L]][1L] # only consider main effect of treatment
      est[1L] <- esb1
      
      # estimate m for each b1 bootstrap sample
      phat <- res1$nonreg[2L]
      est[2L] <- phat
      
      # resampling size
      m <- n^({1.0 + alpha * {1.0 - phat}}/{1.0 + alpha})
      est[3L] <- m
      
      # loop over B2 second stage bootstrap samples
      for (k in seq_len(B2) ) {
        # resample with replacement
        index <- sample(seq_len(n), floor(m), replace = TRUE)

        if (is.list(treat.wgt.man)) {
          tx_wgt_man <- lapply(tx.wgt.man, "[", index)
        } else {
          tx_wgt_man <- NULL
        }
        
        # fit the model to bootstrap sample
        res2 <- tryCatch(DTRreg(outcome = Y[index], 
                                blip.mod = blip.mod, 
                                treat.mod = treat.mod, 
                                tf.mod = tf.mod, 
                                data = boot1[index, ], 
                                method = method, 
                                weight = weight, 
                                n.bins = n.bins,
                                missing = missing, 
                                missing.mod = missing.mod,
                                treat.type = treat.type,
                                treat.fam = treat.fam,
                                treat.wgt.man = tx_wgt_man, 
                                treat.range = treat.range,
                                var.estim = "none"),
                         error = function(e) {
                           stop("unable to complete DTRreg\n\t",
                                e$message, call. = FALSE)
                         })
        
        # save the (b1,b2) bootstrap estimates in the (k+3) column
        est[k + 3L] <- res2$psi[[1L]][1L]
      }
      quan <- quantile(sqrt(m) * {est[-c(1:3)] - esb1}, 
                       probs = c(0.025, 0.975))
      coverage <- coverage + {{esb1 - quan[2L] / sqrt(m) <= blip.psi.1[1L]} & 
                              {esb1 - quan[1L] / sqrt(m) >= blip.psi.1[1L]}}
    }
    
    coverage <- coverage / B1
    
    # print coverage and alpha
    message("alpha= ", format(alpha, digits = 4L), 
            ", coverage = ", format(coverage, digits = 4L))
    
    if ( coverage >= 0.95 ) {
      message("final value of alpha = ", format(alpha, digits = 4L))
      m <- n^({1.0 + alpha * {1.0 - p}} / {1.0 + alpha})
      message("selected subsample size m = ", floor(m))
      break
    } else {
      alpha <- alpha + 0.025
    }
  }
  list("m" = as.numeric(floor(m)))
}