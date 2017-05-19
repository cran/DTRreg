\name{chooseM}
\alias{chooseM}
\title{Adaptive Choice of the Bootstrap Resample Size M for the m-out-of-n Bootstrap with for DTR Estimation}
\description{
Implementation of a double-bootstrap alogrithm for choosing the bootstrap resample size m in a data-adaptive manner. The function returns a resample size m to be used to apply the m-out-of-n bootstrap with DTRreg. }
\usage{
chooseM(outcome, blip.mod, treat.mod, tf.mod, data = NULL,
        method = "gest", weight = "default", missing = "default", 
                        treat.mod.man = NULL, B1 = 500, B2 = 500)
}
\arguments{
  \item{outcome}{The outcome variable.}
  \item{blip.mod}{A list of formula objects specifying covariates of a (linear) blip function for each stage in order. No dependent variable should be specified.}
  \item{treat.mod}{A list of formula objects specifying the treatment model for each stage in order. Treatment variable should be included as the dependent variable. If treatment is binary a logistic regression model will be used, otherwise a linear regression model will be used.}
  \item{tf.mod}{A list of formula objects specifying covariates of a (linear) treatment-free model for each stage in order. No dependent variable should be specified.}
  \item{data}{A data frame containing all necessary covariates contained in the above models.}
  \item{method}{The DTR method to be used, choose "dwols" for dynamic WOLS, "gest" for G-estimation, or "qlearn" for Q-learning.}
  \item{weight}{If using dynamic WOLS the option for the weights used. Default is the form |A - E[A|...]|, "iptw" gives inverse probability of treatment style weights.}
  \item{missing}{If set to "ipcw" and data are missing then inverse probability of censored weights is used with the probability of censoring estimated via logistic regression on the full covariate history up to that point.}
  \item{treat.mod.man}{A list of vectors of known treatment weights can be specified to be used instead of those estimated by the routine.}
  \item{B1}{Number of first-level boostrap resamples.}
  \item{B2}{Number of second-level boostrap resamples.}
}
\details{
The m-out-of-n bootstrap is an adequate tool for constructing valid confidence intervals for the first stage parameters in DTRreg. The resample size m is: \eqn{m = n^{\frac{1+alpha(1-pHat)}{1+alpha}}}{% 
      m = n^[1+alpha(1-pHat)/(1+alpha)]}. The estimated non-regularity level is computed by \code{DTRreg}. The double-bootstrap algorithm is a cross-validation tool for choosing the tuning parameter alpha in a data-driven way. 
      
The current implementation is valid for a two-stage DTR. Moreover, the current implementation may be unstable when there are many missing data.
}
\value{
\item{m }{Resample size for using in the m-out-of-n bootstrap.}
}
\references{
Chakraborty, B., Moodie, E. E. M. (2013) \emph{Statistical Methods for Dynamic Treatment Regimes}. New York: Springer.

Efron B., Tibshirani R. J. (1994) An Introduction to the Bootstrap. \emph{CRC press}.

Wallace, M. P., Moodie, E. M. (2015) Doubly-Robust Dynamic Treatment Regimen Estimation Via Weighted Least Squares. \emph{Biometrics} \bold{71}(3), 636--644 (doi:10.1111/biom.12306.)
} % Eventually update with reference to my paper
\author{Gabrielle Simoneau}

\examples{
##################
# example single run of a 2-stage g-estimation analysis
set.seed(1)
# expit function
expit <- function(x) {1 / (1 + exp(-x))}
# sample size
n <- 100
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
# models to be passed to DTRreg
# blip model
blip.mod <- list(~X1, ~X2)
# treatment model (correctly specified)
treat.mod <- list(A1~X1, A2~X2)
# treatment-free model (incorrectly specified)
tf.mod <- list(~X1, ~X2)

# perform dWOLS without calculating confidence intervals
mod1 <- DTRreg(Y, blip.mod, treat.mod, tf.mod, method = "dwols")

# choose m adaptively for that model
# m <- chooseM(Y, blip.mod, treat.mod, tf.mod, method = "dwols", B1 = 200, B2 = 200)$m
m <- 94

# perform dWOLS and calculate confidence intervals with the m-out-of-n bootstrap
mod2 <- DTRreg(Y, blip.mod, treat.mod, tf.mod, method = "dwols", var.estim = "bootstrap", M = m)

##################
}

\keyword{dynamic treatment regimens,adaptive treatment strategies,personalized medicine,g-estimation,dynamic weighted ordinary least squares,m-out-of-n bootstrap,double bootstrap}
