#' Toy Two-Stage Trial Datasets
#' 
#' These datasets are provided only to facilitate examples. They are not based 
#'   on or representative of any real-world applications.
#' 
#' @name data
#' @rdname data
#' 
#' @usage data(twoStageCont)
#'
#' @format twoStageCont is a dataset generated to mimic a simple two-stage
#'   trial. The data.frame contains 1000 observations with 5 columns: 
#' \describe{
#' \item{X1}{The first stage covariate. A normally distributed continuous variable.}
#' \item{A1}{The first stage treatment. A binary variable.}
#' \item{X2}{The second stage covariate. A normally distributed continuous variable.}
#' \item{A2}{The second stage treatment. A binary variable.}
#' \item{Y}{The outcome. A continuous variable.}
#' }
#'
#' @keywords datasets
"twoStageCont"



#' @noRd
#' @keywords internal
.generateTwoStageCont <- function(seed = 1234L, n = 1000L) {
  set.seed(seed)
 
  # expit function
  expit <- function(x) { 1.0 / (1.0 + exp(-x)) }

  # variables (X = patient information, A = treatment)
  X1 <- rnorm(n)
  A1 <- rbinom(n, 1, expit(X1))
  X2 <- rnorm(n)
  A2 <- rbinom(n, 1, expit(X2))
   
  # blip functions
  gamma1 <- A1 * (1 + X1)
  gamma2 <- A2 * (1 + X2)

  # observed outcome: treatment-free outcome plus blip functions
  Y <- exp(X1) + exp(X2) + gamma1 + gamma2 + rnorm(n)
  
  data.frame("X1" = X1, "A1" = A1,
             "X2" = X2, "A2" = A2, "Y" = Y)
}

#' @rdname data
#' 
#' @usage data(twoStageCens)
#'
#' @format twoStageCens is a dataset generated to mimic a simple two-stage
#'   trial with right-censoring. The data.frame contains 1000 observations 
#'   with 9 columns: 
#' \describe{
#' \item{X11}{A first stage covariate. A normally distributed continuous variable.}
#' \item{X12}{A first stage covariate. A continuous variable = X11^4.}
#' \item{A1}{The first stage treatment. A binary variable.}
#' \item{T1}{The time from the beginning of the first stage to the event or to 
#'           stage 2 entry, whichever comes first.}
#' \item{X21}{A second stage covariate. A normally distributed continuous variable.}
#' \item{X22}{A second stage covariate. A continuous variable = X21^3.}
#' \item{A2}{The second stage treatment. A binary variable.}
#' \item{T2}{The time from the beginning of the second stage to the event 
#'           defined only for subjects who enter the second stage.}
#' \item{delta}{Event indicator.}
#' }
#' 
#' Note: For participants who experienced the event during stage 1, i.e., did
#'   not continue to stage 2, the "survival time" is T1. For participants that
#'   entered stage 2, the "survival time" is T1 + T2.
#' 
"twoStageCens"

#' Generate toy data for two-stage survival with censoring
#' @noRd
#' @importFrom stats rbinom rnorm runif
#' 
#' @keywords internal
.generateTwoStageCensored <- function(seed = 1234L, n = 1000L) {
  
  set.seed(seed)
  
  # expit function
  expit <- function(x) {1.0 / (1.0 + exp(-x))}
  
  theta1 <- c(6.3, 1.5, -0.8, 0.1, 0.1)
  theta2 <- c(4, 1.1, -0.2, -0.9, 0.6, -0.1)
  lambda <- 1/300
  p <- 0.9
  beta <- 2
  # covariates and treatment (X = patient information, A = treatment)
  X1 <- runif(n, 0.1, 1.29)
  X14 <- X1^4
  A1 <- rbinom(n, size = 1, prob = expit(2*X1 - 1))
  X2 <- runif(n, 0.9, 2)
  X23 <- X2^3
  A2 <- rbinom(n, size = 1, prob = expit(-2*X2 + 2.8))
  delta <- rbinom(n, size = 1, prob = expit(2*X1 - 0.4))
  eta2 <- rbinom(n, 1, prob = 0.8)
  delta2 <- delta[eta2 == 1]
  # survival time
  logY2 <- logT2 <- theta2[1] + theta2[2]*X2[eta2 == 1] + theta2[3]*X23[eta2 == 1] 
  + theta2[4]*A2[eta2 == 1] + theta2[5]*A2[eta2 == 1]*X2[eta2 == 1] 
  + theta2[6]*X1[eta2 == 1] + rnorm(sum(eta2), sd = 0.3)
  trueA2opt <- ifelse(theta2[4]*A2[eta2 == 1] 
                      + theta2[5]*A2[eta2 == 1]*X2[eta2 == 1] > 0, 1, 0)
  logT2opt <- logT2 + (trueA2opt - A2[eta2 == 1])*(theta2[4]*A2[eta2 == 1] 
                                                   + theta2[5]*A2[eta2 == 1]*X2[eta2 == 1])
  logT <- theta1[1] + theta1[2]*X1 + theta1[3]*X14 + theta1[4]*A1 
  + theta1[5]*A1*X1 + rnorm(n, sd = 0.3)
  T1 <- exp(logT[eta2 == 1 & delta == 1]) - exp(logT2opt[delta2 == 1])
  logT[eta2 == 1 & delta == 1] <- log(T1 + exp(logT2[delta2 == 1]))
  # censoring time
  C <- (- log(runif(n - sum(delta), 0, 1))/(lambda * exp(beta * X1[delta == 0])))^(1/p)
  eta2d0 <- eta2[delta == 0]
  C1 <- rep(NA, length(C))
  C2 <- rep(NA, length(C))
  for(i in 1:length(C))
  {
    if(eta2d0[i] == 0){
      C1[i] <-  C[i]
      C2[i] <- 0
    }else{
      C1[i] <- runif(1, 0, C[i])
      C2[i] <- C[i] - C1[i]
    }
  }
  # observed survival time
  Y2 <- rep(NA, n)
  Y1 <- rep(NA, n)
  Y2[delta == 0] <- C2
  Y1[delta == 0] <- C1
  Y1[delta == 1 & eta2 == 1] <- T1
  Y1[delta == 1 & eta2 == 0] <- exp(logT[delta == 1 & eta2 == 0])
  Y2[delta == 1 & eta2 == 0] <- 0
  Y2[delta == 1 & eta2 == 1] <- exp(logT2[delta2 == 1])
  logY <- log(Y1 + Y2)
  logY2 <- log(Y2[eta2 == 1])

  data.frame("X11" = X1, "X12" = X14, "A1" = A1, 
             "X21" = X2, "X22" = X23, "A2" = A2, 
             "delta" = delta, "T1" = Y1, "T2" = Y2)
  
}

#' @rdname data
#' 
#' @usage data(twoStageSurv)
#'
#' @format twoStageSurv is a dataset generated to mimic a simple two-stage
#'   trial without censoring. The data.frame contains 1000 observations 
#'   with 9 columns: 
#' \describe{
#' \item{X11}{A first stage covariate. A normally distributed continuous variable.}
#' \item{X12}{A first stage covariate. A continuous variable = X11^4.}
#' \item{A1}{The first stage treatment. A binary variable.}
#' \item{T1}{The time from the beginning of the first stage to the event or to 
#'           stage 2 entry, whichever comes first.}
#' \item{X21}{A second stage covariate. A normally distributed continuous variable.}
#' \item{X22}{A second stage covariate. A continuous variable = X21^3.}
#' \item{A2}{The second stage treatment. A binary variable.}
#' \item{T2}{The time from the beginning of the second stage to the event 
#'           non-zero only for subjects who did not have an event in Stage I.}
#' }
#' 
#' Note: The "survival time" is T1 + T2.
#' 
"twoStageSurv"

#' Generate toy data for two-stage survival with censoring
#' @noRd
#' @importFrom stats rbinom rnorm runif
#' @keywords internal
.generateTwoStageSurv <- function(seed = 1234L, n = 1000L) {
  
  set.seed(seed)
  
  # expit function
  expit <- function(x) {1.0 / (1.0 + exp(-x))}
  
  theta1 <- c(6.3, 1.5, -0.8, 0.1, 0.1)
  theta2 <- c(4, 1.1, -0.2, -0.9, 0.6, -0.1)
  lambda <- 1/300
  p <- 0.9
  beta <- 2
  # covariates and treatment (X = patient information, A = treatment)
  X1 <- runif(n, 0.1, 1.29)
  X14 <- X1^4
  A1 <- rbinom(n, size = 1, prob = expit(2*X1 - 1))
  X2 <- runif(n, 0.9, 2)
  X23 <- X2^3
  A2 <- rbinom(n, size = 1, prob = expit(-2*X2 + 2.8))
  delta <- rbinom(n, size = 1, prob = expit(2*X1 - 0.4))
  eta2 <- rbinom(n, 1, prob = 0.8)
  delta2 <- delta[eta2 == 1]
  # survival time
  logY2 <- logT2 <- theta2[1] + theta2[2]*X2[eta2 == 1] + theta2[3]*X23[eta2 == 1] 
  + theta2[4]*A2[eta2 == 1] + theta2[5]*A2[eta2 == 1]*X2[eta2 == 1] 
  + theta2[6]*X1[eta2 == 1] + rnorm(sum(eta2), sd = 0.3)
  trueA2opt <- ifelse(theta2[4]*A2[eta2 == 1] 
                      + theta2[5]*A2[eta2 == 1]*X2[eta2 == 1] > 0, 1, 0)
  logT2opt <- logT2 + (trueA2opt - A2[eta2 == 1])*(theta2[4]*A2[eta2 == 1] 
                                                   + theta2[5]*A2[eta2 == 1]*X2[eta2 == 1])
  logT <- theta1[1] + theta1[2]*X1 + theta1[3]*X14 + theta1[4]*A1 
  + theta1[5]*A1*X1 + rnorm(n, sd = 0.3)

  delta <- rep(1,n)
  delta2 <- delta[eta2 == 1]
  T1 <- exp(logT[eta2 == 1 & delta == 1]) - exp(logT2opt[delta2 == 1])
  logT[eta2 == 1 & delta == 1] <- log(T1 + exp(logT2[delta2 == 1]))
  # observed survival time
  Y2 <- rep(NA, n)
  Y1 <- rep(NA, n)
  Y1[delta == 1 & eta2 == 1] <- T1
  Y1[delta == 1 & eta2 == 0] <- exp(logT[delta == 1 & eta2 == 0])
  Y2[delta == 1 & eta2 == 0] <- 0
  Y2[delta == 1 & eta2 == 1] <- exp(logT2[delta2 == 1])
  logY <- log(Y1 + Y2)
  logY2 <- log(Y2[eta2 == 1])
  
  data.frame("X11" = X1, "X12" = X14, "A1" = A1, 
             "X21" = X2, "X22" = X23, "A2" = A2, 
             "delta" = rep(1L, n), "T1" = Y1, "T2" = Y1)
}