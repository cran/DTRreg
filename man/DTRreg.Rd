% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DTRreg.R, R/s3-methods.R
\name{DTRreg}
\alias{DTRreg}
\alias{print.DTRreg}
\alias{summary.DTRreg}
\alias{coef.DTRreg}
\title{DTR Estimation and Inference via G-estimation, Dynamic WOLS, or Q-learning}
\usage{
DTRreg(
  outcome,
  blip.mod,
  treat.mod,
  tf.mod,
  data = NULL,
  method = c("gest", "dwols", "qlearn"),
  interactive = FALSE,
  treat.type = c("bin", "multi", "cont"),
  treat.fam = gaussian(link = "identity"),
  weight = c("abs", "ipw", "cipw", "qpom", "wo", "none", "manual"),
  n.bins = 3L,
  treat.range = NULL,
  treat.wgt.man = NULL,
  var.estim = c("none", "bootstrap", "sandwich"),
  full.cov = FALSE,
  bootstrap.controls = list(B = 100L, M = nrow(data), type = "standard", truncate = 0,
    verbose = FALSE, interrupt = FALSE),
  missing = c("drop", "ipw"),
  missing.mod = NULL,
  dtr = TRUE
)

\method{print}{DTRreg}(x, ...)

\method{summary}{DTRreg}(object, ...)

\method{coef}{DTRreg}(object, ...)
}
\arguments{
\item{outcome}{The outcome variable. Missing data will result in a stopping
error.}

\item{blip.mod}{A list of formula objects specifying covariates of the
blip function for each stage in order. No dependent variable should be 
specified. Note that this input should include the treatment variable 
ONLY if the blip model is quadratic in treatment. See Details for further 
clarification.}

\item{treat.mod}{A list of formula objects specifying the treatment model for
each stage in order. Treatment variable should be included as the dependent
variable. If treatment is binary \code{\link[stats]{glm}(family = binomial)}
will be 
used to obtain parameter estimates; if treatment is multi-nomial, 
\code{\link[nnet]{multinom}()} will be used to obtain parameter estimates; and if 
treatment is continuous, \code{\link[stats]{lm}()} will be used.}

\item{tf.mod}{A list of formula objects specifying covariates of the
treatment-free model for each stage in order. No dependent variable should 
be specified.}

\item{data}{A data frame containing all necessary covariates and treatments
contained in the models. Missing data should be coded as \code{NA}.}

\item{method}{The DTR method to be used, choose "dwols" for dynamic WOLS, 
"gest" for G-estimation, or "qlearn" for Q-learning.}

\item{interactive}{If \code{TRUE} on-screen prompts will guide the user through 
the specification of blip, treatment, and treatment-free models.}

\item{treat.type}{A character object. Must be one of \{"bin", "multi", "cont"\}
indicating that the treatments at each stage are binary, multinomial,
or continuous, respectively. Each stage must have the same treatment type.}

\item{treat.fam}{A character or family object. 
The description of the dose distribution along with the link 
function to be used in the treatment model for computing weights; should be 
specified in a similar format as that used in \code{\link[stats]{glm}()}. 
If character object, must be one of \{"gaussian", "Gamma"\}, for which
\code{\link[stats]{gaussian}(link = "identity")} or 
\code{\link[stats]{Gamma}(link = "log")} 
will be used,
respectively. Input is ignored for \code{treat.type = "bin"} and 
\code{treat.type = "multi"}.}

\item{weight}{The form of the treatment weight. See details.}

\item{n.bins}{An integer object. The number of bins (levels) to be used for 
categorizing continuous doses. This input is required only when
\code{treat.type = "cont"} and \code{weight = "wo"} or \code{weight = "qpom"}.}

\item{treat.range}{For continuous treatments. Specify the maximum/minimum 
value that treatments can be take. If unspecified then the minimum/maximum 
value of observed treatments is used. If you wish to have unrestricted 
treatments set this option to \code{c(-Inf, +Inf)}. If each stage has its own
range, provide as a list, the ith element providing the min and max
for the ith stage treatment.}

\item{treat.wgt.man}{NULL or a list of vectors of known treatment weights can be 
specified to be used instead of hard-coded treatment weight options.
The \eqn{i^{th}}{ith} element of the list contains the multiplicative weights 
for the \eqn{i^{th}}{ith} stage. Each vector must be of length \eqn{n}{n}, 
the number of participants. Used only for \code{method = "dwols"}.}

\item{var.estim}{Covariance matrix estimation method, either "bootstrap"
or "sandwich" for sandwich estimation.}

\item{full.cov}{A logical. If \code{TRUE}, the full covariance matrix will be
returned. If \code{FALSE}, only the terms pertaining to the blip parameters
are returned.}

\item{bootstrap.controls}{A named list specifying control parameters of the
bootstrap if \code{var.estim = "bootstrap"}. Available controls are:
\describe{
  \item{B: }{The number of bootstrap samples.} 
  \item{M: }{The subsample size for m out of n bootstrap.}
  \item{type: }{The type of bootstrap. Must be one of \{"standard", 
  "empirical", "normal"\}. The last two are parametric bootstraps.}
  \item{truncate: }{A number between 0 and 0.5. The lowest and highest 
  specified proportion of parameter estimates will be replaced by the
  relevant quantiles affording some robustness to extreme values 
  when estimating covariance.}
  \item{verbose: }{If TRUE, estimated time to completion will be printed 
  to the console every ~30 seconds.}
  \item{interrupt: }{If TRUE then user will be given the option 
  to abort the bootstrap without error if estimated time to completion
  exceeds 10 minutes.}
}}

\item{missing}{A character object. Must be one of \{"drop", "ipw"\}.
If set to "ipw" and covariate or treatment data are missing then inverse 
probability 
weights are used. The complete case probability is estimated 
via logistic regression. If set to "drop" and data are missing, participants
with missing data are dropped for all stage analyses.}

\item{missing.mod}{An optional list of formula objects specifying the model
for the inverse probability of weights for each stage in order.
No dependent variable should be specified. If \code{missing = "ipw"} and 
\code{missing.mod = NULL}, then the models are assumed to be linear comprising
the full covariate history derived from all of the previous stage models.}

\item{dtr}{A logical object. If \code{TRUE}, use the DTR estimation approach, which
estimates the stage pseudo-outcome by adding a regret function. If \code{FALSE},
use an 'effect estimation' approach, which treats the observed outcome
as being equal to an outcome assuming no treatment is received at any
stage, plus a blip component at each stage; each stage pseudo-outcome is 
generated by subtracting a blip function. Note that most of the 
DTR-specific output will either be suppressed or irrelevant.}

\item{x}{An object of class `DTRreg`.}

\item{...}{Ignored.}

\item{object}{An object of class `DTRreg`.}
}
\value{
An object of class \code{DTRreg}, a list including elements
    \item{K: }{The number of decision points.}
    \item{beta: }{A list. The ith element contains the parameter estimates of
      the ith stage treatment-free model.}
    \item{psi: }{A list. The ith element contains the parameter estimates of
      the ith stage blip model.}
    \item{covmat: }{A list. The ith element contains covariance matrix of 
      the ith stage blip parameter estimates.}
    \item{nonreg: }{Non-regularity estimates.}
    \item{setup: }{A list detailing the input parameter settings used for the
      analysis}
      \describe{
        \item{models: }{A list of the models used for the analysis.}
        \item{method: }{The parameter estimation method.}
        \item{var.estim: }{The variance esetimation method.}
        \item{cc.modeled: }{If TRUE, missing data was modeled. If FALSE, cases
           with missing data were removed from the analysis.}
        \item{tx.weight: }{The treatment weighting used for the analysis.}
        \item{tx.type: }{Treatment was binary, multinomial, or continuous.}
        \item{n.bins: }{The number of bins (levels) used for categorizing
          continuous doses when \code{tx.weight = "wo"} or 
          \code{tx.weight = "qpom"}.}
        \item{tx.wgt.man: }{Any user provided treatment weights.}
        \item{tx.range: }{For continuous treatments, the range of allowed
          treatment values.}
        \item{tx.family: }{The description of the dose distribution along 
          with the link function used in the continuous treatment model.}
        \item{boot.controls: }{A list of the bootstrap controls.}
        \item{type: }{The type of effect. Dynamic treatment regime or treatment
          effect.}
      }
    \item{training_data: }{A list containing the training data.}
      \describe{
        \item{data: }{The covariates and treatment data.}
        \item{outcome: }{The outcome of interest.}
        \item{A: }{The treatment variables, possibly recoded to adhere to internal
          code requirements.}
      }
    \item{analysis: }{A list containing the primary results of each stage analysis.}
      \describe{
        \item{n: }{The number of participants included in the stage analysis.}
        \item{last.stage: }{The last stage each participant was included in 
          the analysis.}
        \item{prob.cc: }{The complete case probabilities.}
        \item{cc.mod.fitted: }{The regression objects returned for estimating
          the complete case probabilities.}
        \item{cc.wgt: }{The complete case weights.}
        \item{cts: }{The treatment type at each stage.}
        \item{tx.mod.fitted: }{The regression objects returned for estimating
          the treatment probabilities.}
        \item{A.hat: }{The estimated or provided treatment probabilities.}
        \item{tx.wgt: }{The treatment weights.}
        \item{outcome.fit: }{The regression objects returned for each stage outcome
          regression.}
        \item{Y: }{The pseudo-outcomes.}
        \item{regret: }{Estimates of the regret for each subject based on observed 
          treatment and blip parameter estimates.}
        \item{opt.treat: }{Optimal treatment decisions for each subject at each 
          stage of treatment.}
        \item{opt.Y: }{Predicted optimal outcome under recommended regimen.}
     }
    \item{call: }{The original function call.}
  The functions \code{coef()}, \code{predict()} and 
  \code{confint()} may be used with such 
  model objects. The first two have specific help files for their 
  implementation, while \code{confint()} is used in the same way as 
  the standard 
  \code{\link[stats]{confint}()} command, with the exception of the \code{parm} 
  option, which is not available.
}
\description{
Dynamic treatment regimen estimation and inference via G-estimation and 
  dynamic WOLS.  Estimation of blip model parameters for multi-stage data.
}
\details{
\code{DTRreg()} allows the estimation of optimal dynamic treatment regimens 
  (DTRs, also known as adaptive treatment strategies) from multi-stage 
  trials using G-estimation, dynamic weighted ordinary least squares 
  (dWOLS), and genearlized dWOLS. All methods focus on estimating the 
  parameters of the blip: a 
  model of the difference in expected outcome under the observed treatment 
  and some reference treatment (usually a control) at a given stage, assuming
  identical histories and optimal treatment thereafter. The reader is 
  referred to Chakraborty and Moodie (2013) for a thorough introduction and 
  review of DTR methods. The dWOLS method may be used to obtain parameter 
  estimates identical to those from Q-learning (by setting \code{weight = "none"}).
  This option is intended primarily for exploratory purposes; the authors 
  note that there is a dedicated R package for Q-learning (qLearn), although 
  it is limited to the 2-stage setting; multi-stage settings are available
  in R package DynTxRegime.
  
  This implementation assumes an outcome regression model of the form
  E(Y|X=x,A=a) = tf.mod + a blip.mod. That is -- the input \code{blip.mod} 
  formula should include the treatment variable \emph{ONLY} if it is quadratic. 
  For example, if the full blip model is linear in the treatment variable 
  \deqn{\sim a \psi_0 + a x \psi_1,}{~ a psi_0 + a x psi_1,} then the input
  should model should be \code{blip.mod = ~ x}. 
  If the full blip model is quadratic in the treatment variable
  \deqn{\sim a \psi_0 + a^2 \psi_1 + a x \psi_2 + a^2 x \psi_3,}{
  ~ a psi_0 + a^2 psi_1 + a x psi_2 + a^2 x psi_3,} \code{blip.mod =
  ~ a*x}. For continuous treatments, only quadratic blip
  functions are supported.
  
  All methods require the specification of three models for each 
  stage of the analysis: a treatment model (conditional mean of the treatment 
  variable), a treatment-free model (conditional mean of outcome assuming 
  only reference treatments are used), and a blip model. Only the blip model
  must be correctly specified (or over-specified), with consistent parameter 
  estimates obtainable if at least one of the other two models is correctly 
  specified. Note that all of these must be specified as lists of formula 
  objects, even if only one stage of treatment is considered.
  
  Note that as is conventional, it is assumed a larger value of the outcome 
  is preferred (which can be easily achieved via transformation of your data 
  if necessary).
  
  When treatment is binary, if confidence intervals are computed (via 
  specification of \code{var.estim} other than "none"), then DTRreg will calculate 
  the proportion of subjects at each stage for whom optimal treatment is 
  non-unique. If this proportion exceeds 0.05 a non-regularity warning will 
  be displayed, along with the proportion of subjects for whom this is the 
  case. Note that this warning is only displayed if a variance estimation 
  option is selected.
  
  Several treatment weight function options have been implemented within the 
  package:
  \itemize{ 
    \item "none": No treatment weights applied. If \code{method = "dWOLS"}, this 
      selection results in the implementation of Q-learning, modified 
      slightly to use the G-estimation or dWOLS style pseudo-outcome 
      (computed using the observed outcome modified by the estimated 
       treatment effect) rather than the traditional Q-learning outcome 
      (predicted based on model only, rather than observed outcome with 
      treatment effect).   
    \item "ipw": weights based on the inverse probability of 
      treatment. For binary treatments, a logistic regression is used.
      For multinomial, a multinomial log-linear model is fit using 
      \code{\link[nnet]{multinom}}. For continuous treatments, a GLM with the specified 
      family and link function provided in the \code{treat.fam} argument is used.
    \item "cipw": inverse probability of treatment weights as described for
      "ipw" and capped at the 99th percentile of the observed weights.
    \item "qpom": weights based on the stabilized inverse 
      probability of treatment applied to the categorized (into n.bins bins) 
      continuous doses or multinomial treatments; probabilities are 
      calculated using a proportional odds model.
      This weight is appropriate only for continuous and multinomial treatments.
    \item "wo": overlap weights for the categorized continuous doses 
      or multinomial treatments (Li and Li, 2019).
      This weight is appropriate only for continuous treatments.
    \item "abs": Absolute difference \eqn{|A - E[A|...]|}{|A - E[A|...]|}. 
      appropriate only for binary treatments.
    \item "manual": User provides weights through input \code{treat.wgt.man}.
      Manual treatments are only used in dwols.
   }
}
\examples{
data(twoStageCont)

# models to be passed to DTRreg
# blip model
blip.mod <- list(~ X1, ~ X2)
# treatment model (correctly specified)
treat.mod <- list(A1 ~ X1, A2 ~ 1)
# treatment-free model (incorrectly specified)
tf.mod <- list(~ X1, ~ X2)

# perform G-estimation
mod1 <- DTRreg(twoStageCont$Y, blip.mod, treat.mod, tf.mod, 
               data = twoStageCont, method = "gest")
mod1
}
\references{
Chakraborty, B., Moodie, E. E. M. (2013) \emph{Statistical Methods for 
  Dynamic Treatment Regimes}. New York: Springer.
  
  Robins, J. M. (2004) \emph{Optimal structural nested models for optimal 
  sequential decisions}. In Proceedings of the Second Seattle Symposium on 
  Biostatistics, D. Y. Lin and P. J. Heagerty (eds), 189--326. New York: 
  Springer.
  
  Wallace, M. P., Moodie, E. E. M. (2015) Doubly-Robust Dynamic Treatment 
  Regimen Estimation Via Weighted Least Squares. \emph{Biometrics} 
  \bold{71}(3), 636--644 (doi:10.1111/biom.12306.)
  
  Simoneau, G., Moodie, E. E. M., Nijjar, J. S., and Platt, R. W. (2020)
  Finite Sample Variance Estimation for Optimal Dynamic Treatment
  Regimes of Survival Outcomes. \emph{Statistics in Medicine} \bold{39},
  4466-4479.
  
  Efron, B., and Tibshirani, R. (1986)
  Bootstrap Methods for Standard Errors, Confidence Intervals, and Other 
  Measures of Statistical Accuracy \emph{Source: Statistical Science} \bold{1}
  54-75.
}
\author{
Michael Wallace

Shannon T. Holloway
}
\concept{adaptive treatment strategies}
\concept{dynamic treatment regimens}
\concept{dynamic weighted ordinary least squares}
\concept{g-estimation}
\concept{personalized medicine}
