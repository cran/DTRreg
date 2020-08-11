\name{DWSurv}
\alias{DWSurv}
\title{DTR estimation and inference for time-to-event data using DWSurv}
\description{
Dynamic treatment regimen estimation and inference via dynamic weighted survival modeling (DWSurv).  Inference for the blip estimators with single- and multi-stage data.
}

\usage{
DWSurv(time, blip.mod, treat.mod, tf.mod, cens.mod, data = NULL, weight = "default",
        var.estim = "none", asymp.opt = "adjusted", boot.opt = "standard", B = 500,
        optimization = "max", quiet = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{time}{A list of formula specifying the survival time variable for each stage in order. The time variable should be specified on the right hand side of the formula. No dependent variable should be specified. The list should be as long as the maximum number of stages.}
  \item{blip.mod}{A list of formula objects specifying covariates of a (linear) blip function for each stage in order. No dependent variable should be specified.}
  \item{treat.mod}{A list of formula objects specifying the treatment model for each stage in order. The treatment variable should be binary and included as the dependent variable. Logistic regression models are used.}
  \item{tf.mod}{A list of formula objects specifying covariates of a (linear) treatment-free model for each stage in order. No dependent variable should be specified.}
  \item{cens.mod}{A list of formula objects specifying the censoring model for each stage in order. The event indicator, which takes value 1 if an event was observed and 0 otherwise, should be included as the dependent variable and should be the same across stages. In the absence of censoring, one still needs to specify an event indicator with 1s on the right-hand side of the formula and leave the left-hand side empty (see example below).}
  \item{data}{A data frame containing all necessary covariates contained in the above models.}
  \item{weight}{A user-supplied function for the weights to be used in DWSurv. The function must have the four following arguments: treatment received A, probability of receiving treatment A=1, status, probability of being observed status = 1. Default is the inverse probability of censoring weights combined with |A - E[A|...]|.}
  \item{var.estim}{Covariance matrix estimation method, either "asymptotic", "bootstrap" or "none" (default).}
 \item{asymp.opt}{If the asymptotic variance estimation is used, specify either the "adjusted" (default) or "naive" version.}
  \item{boot.opt}{If bootstrap is used for variance estimation, specify either the "standard" (default), "empirical" or "normal". The last two are parametric bootstraps.}
  \item{B}{Number of bootstrap resamples, if applicable.}
  \item{optimization}{If "max" (default), it is assumed that larger values/longer survival times are preferred. Set to "min" if the sequence of optimal decision rules should minimize survival time.}
  \item{quiet}{To suppress warnings when bootstrapping.}
}
\details{
The function \code{DWSurv()} allows estimating an optimal dynamic treatment regime from multi-stage trials or observational data when the outcome of interest is survival time subject to right-censoring. The dynamic weighted survival modeling (DWSurv) algorithm is implemented.  The method focuses on estimating the parameters of the blip: a model of the difference in expected outcome under the observed treatment and some reference treatment (usually a control) at a given stage, assuming identical histories and optimal treatment thereafter.

The method requires the specification of four models for each stage of the analysis: a treatment model (conditional mean of the treatment variable), a censoring model, a treatment-free model (conditional mean of outcome assuming only reference treatments are used), and a blip model.  Only the blip model must be correctly specified (or over-specified), with consistent parameter estimates obtainable if at least one of the treatment-free or the treatment and censoring models are correctly specified.  Note that all of these must be specified as lists of formula objects, even if only one stage of treatment is considered.

Note that as is conventional, it is assumed a larger survival time is preferred (which can be easily achieved via transformation of your data if necessary).

%When treatment is binary, if confidence intervals are computed (via specification of var.estim other than 'none'), then DTRreg will calculate the proportion of subjects at each stage for whom optimal treatment is non-unique.  If this proportion exceeds 0.05 a non-regularity warning will be displayed, along with the proportion of subjects for whom this is the case.  Note that this warning is only displayed if a variance estimation option is selected.
}
\value{An object of class \code{DTR}, a list including elements
 \item{psi}{Blip parameter estimates for each stage of treatment.}
 \item{opt.treat}{Optimal treatment decisions for each subject at each stage of treatment.}
 \item{covmat}{Covariance matrix of blip parameter estimates.}
 \item{log.regret}{Estimates of the log-transformed regret for each subject based on observed treatment and blip parameter estimates.}
 \item{beta}{Treatment-free model parameter estimates (note that these may not be consistent).}
 \item{opt.Y}{Predicted optimal survival time under recommended regimen.}
 \item{nonreg}{Non-regularity estimates.}
 \item{psi.boot}{If applicable, the B bootstrap estimates of the blip parameters across stages.}
The functions coef and confint may be used with such model objects.  The first has specific help files for their implementation, while confint is used in the same way as the standard confint command, with an additional \code{type} options which can be set to "percentile" when bootstrap is used to derive confidence intervals. The parm option is not available.
}
\references{
Simoneau, G., Moodie, E. E. M., Wallace, M.P., Platt, R. W. (2020) Optimal Dynamic Treatment Regimes with Survival Endpoints: Introducing DWSurv in the R package DTRreg. \emph{Journal of Statistical Computation and Simulation (in press)}. 

Simoneau, G., Moodie, E. E. M., Nijjar, J. S., Platt, R. W. (2019) Estimating Optimal Dynamic Treatment with Survival Outcomes. \emph{Journal of the American Statistical Association}, pp.1-9 (doi:10.1080/01621459.2019.1629939).

Wallace, M. P., Moodie, E. E. M., Stephens, D. A. (2017) Dynamic Treatment Regimen Estimation via Regression-Based Techniques: Introducing {R} Package {DTRreg}. \emph{Journal of Statistical Software} \bold{80}(2), 1--20 (doi:10.18637/jss.v080.i02).
}
\author{Gabrielle Simoneau}

\examples{
##################
#### example single run of a 2-stage DWSurv analysis
set.seed(1)
# expit function
expit <- function(x) {1 / (1 + exp(-x))}
# sample size and parameters
n <- 1000
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
# data and run DWSurv
mydata <- data.frame(X1,X14,A1,X2,X23,A2,delta,Y1,Y2)
mod <- DWSurv(time = list(~Y1, ~Y2), blip.mod = list(~X1, ~X2), 
  treat.mod = list(A1~X1, A2~X2), tf.mod = list(~X1 + X14, ~X2 + X23 + X1), 
  cens.mod = list(delta~X1, delta~X1), var.estim = "asymptotic", data = mydata)
mod

#### example in the absence of censoring
# create an event indicator
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
# data and run DWSurv
mydata <- data.frame(X1,X14,A1,X2,X23,A2,delta,Y1,Y2)
mod_nocensoring <- DWSurv(time = list(~Y1, ~Y2), blip.mod = list(~X1, ~X2), 
  treat.mod = list(A1~X1, A2~X2), tf.mod = list(~X1 + X14, ~X2 + X23 + X1), 
  cens.mod = list(~delta, ~delta), var.estim = "asymptotic", data = mydata)
mod_nocensoring
##################
}

\concept{dynamic treatment regimens}
\concept{adaptive treatment strategies}
\concept{personalized medicine}
\concept{dynamic weighted survival modeling}
\concept{accelerated failure time}
