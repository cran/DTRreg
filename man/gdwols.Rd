\name{gdwols}
\alias{gdwols}
\title{DTR Estimation and Inference via Generalized dWOLS}
\description{
Dynamic treatment regimen estimation via generalized dynamic weighted ordinary least squares for continuous treatment. Estimation of blip model parameters for multi-stage data.
}

\usage{
gdwols(outcome, Xpsi1, Xpsi2, treat.mod, treat.fam, tf.mod,
        weight.fcn, data = NULL, m, k, treat.range)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{outcome}{The outcome variable.}
  \item{Xpsi1, Xpsi2}{A list of formula objects specifying the tailoring variables of a (quadratic) blip function for each stage in order. Xpsi1 specifies the covariates interacting with the treatment (A); Xpsi2 specifies the covariates interacting with squared-treatment (A^2). No dependent variable should be specified.}
  \item{treat.mod}{A list of formula objects specifying the treatment model for each stage in order. The treatment variable should be included as the dependent variable.}
   \item{treat.fam}{A description of the dose distribution along with the link function to be used in the treatment model for computing weights; should be specified in a similar format as that used in the standard glm function. Options include the gaussian and gamma distributions; ex: gaussian(link = "identity"), Gamma(link = "log"). If unspecified, default will be gaussian(link = "identity").}
  \item{tf.mod}{A list of formula objects specifying covariates of a (linear) treatment-free model for each stage in order. No dependent variable should be specified.}
  \item{weight.fcn}{The weight function to be used in calculating the balancing weights. Options include "ipw", "cipw", "qpom", "wo"; see details below.}
  \item{data}{The data (across all stages), including the outcome, stage-specific covariates and doses.}
  \item{m}{The number of bins (levels) to be used for categorizing the continuous doses. (This argument is only pertinent in the case where the weight function used is based on a categorization of the continuous doses.)}
  \item{k}{The number of treatment stages.}
  \item{treat.range}{The range of permissible dose values for continuous treatments; the optimal treatment will be restricted to lie withinin the specified treatment range. If unspecified then the minimum/maximum value of the observed treatments is used. For unrestricted treatments, set this option to c(-Inf,+Inf). In the case of an unrestricted treatment range, a warning message will be shown whenever the optimal treatment is in fact a minimum for at least one observation.}
}
\details{
GdWOLS (generalized dynamic weighted ordinaly least squares) allows for the estimation of optimal dynamic treatment regimes (DTRs, also known as adaptive treatment strategies) for multi-stage trials in the case where treatment is measured on a continuous scale. This method allows for the estimation of the blip model parameters. The GdWOLS approach requires the specification of three models at each stage: a treatment model (conditional mean of the treatment variable), a treatment-free model (conditional mean of the outcome assuming the reference treatment levels are used, in this case, 0), and a blip model (the difference in expected outcome under the observed treatment and the reference treatment at a given stage, assuming identical histories and optimal treatment thereafter). In the GdWOLS framework, only the blip model must be correctly specified (or over-specified); consistent parameter estimates are obtained as long as at least one of the other two models is correctly specified.  

Several weight function options have been implemented within the package, including the IPW, capped-IPW, Q-POM and W-O weights:

"ipw": (IPW) weights based on the inverse probability of treatment, probabilities are calculated using a GLM with the specified family and link funciton provided in the treat.fam argument. 

"cipw": (capped-IPW) IPW weights capped at the 99th percentile of the observed weights

"qpom": (Q-POM) weights based on the stabilized inverse probability of treatment applied to the categorized (into m bins) doses, probabilities are calculated using a proportional odds model

"wo": (W-O) overlap weights for the categorized continuous doses (Li and Li, 2019)

Note:
- All the model components (Xpsi1, Xpsi2, treatment and treatment-free models) must be specified as lists of formula objects, even if only one stage of treatment is considered.
- The implementation in the gdwols function only considers quadratic blip functions of the form psi1*Xpsi1*A + psi2*Xpsi2*(A^2). This ensures that the optimal treatment is given a.opt=-(psi1*X1)/2(psi2*X2), provided that psi2*X2<0, see paper for more details.
- As is conventional, it is assumed that a larger value of the outcome is preferred (which can be easily achieved via transformation of the data if necessary).

}
\value{An object of class \code{DTR}, a list including elements
 \item{psi}{Blip parameter estimates for each stage of treatment.}
 \item{beta}{Treatment-free model parameter estimates for each stage of treatment.}
 \item{opt.treat}{Optimal treatment value for each stage of treatment.}
}
\references{

Schulz, J., and Moodie, E. E. M. (2020) Doubly Robust Estimation of Optimal Dosing Strategies. Journal of the American Statistical Association. (DOI: 10.1080/01621459.2020.1753521)

Li, F., and Li, F. (2019) "Propensity Score Weighting for Causal Inference With Multiple Treatments," The Annals of Applied Statistics, 13, 2389-2415.

}
\author{Juliana Schulz}

\examples{

##################
# example single run of a 2-stage G-dWOLS analysis

set.seed(1)
# number of stages
k<-2
# number of bins for categorized weight functions
m<-5
# sample size
n<-1000

# parameter values
psi.mat<-matrix(rep(c(1,1,-1),2),byrow=TRUE,nrow=2)
alpha.mat<-matrix(rep(c(-1,1),2),byrow=TRUE,nrow=2)

### generate data
# stage 1
x1<-abs(rnorm(n,10,1))
a1<-rnorm(n,alpha.mat[1,1]+alpha.mat[1,2]*x1,1)
# stage 2
x2<-abs(rnorm(n,10,1))
a2<-rnorm(n,alpha.mat[2,1]+alpha.mat[2,2]*x2,1)
# blips
gamma1<-as.matrix(cbind(a1,a1*x1,a1^2)) %*% psi.mat[1,]
gamma2<-as.matrix(cbind(a2,a2*x2,a2^2)) %*% psi.mat[2,]
# y: outcome
# y <- trmt free + blip
y<-log(x1)+sin(x1)+log(x2)+sin(x2)+gamma1+gamma2 + rnorm(n,0,1)
# convert to a vector for formatting into data frame
y <- as.vector(y)
# data
data<-data.frame(cbind(y,x1,x2,a1,a2))

# models to be passed to gdwols function
# tailoring variables that interact with a1
Xpsi1<-list(~x1,~x2)
# tailoring variables that interact with a2
Xpsi2<-list(~1,~1)
# treatment model at each stage
treat.mod<-list(a1~x1,a2~x2)
# treatment-free model at each stage (misspecified)
tf.mod<-list(~x1,~log(x2)+sin(x2))


out1 <- gdwols(y, Xpsi1, Xpsi2, treat.mod, treat.fam = gaussian(link = "identity"),
               tf.mod, weight.fcn="ipw", data, m, k)
out1$psi

out2 <- gdwols(y, Xpsi1, Xpsi2, treat.mod, treat.fam = gaussian(link = "identity"),
               tf.mod, weight.fcn="cipw", data, m, k)
out2$psi

out3 <- gdwols(y, Xpsi1, Xpsi2, treat.mod, treat.fam = gaussian(link = "identity"),
               tf.mod, weight.fcn="qpom", data, m, k)
out3$psi

out4 <- gdwols(y, Xpsi1, Xpsi2, treat.mod, treat.fam = gaussian(link = "identity"),
               tf.mod, weight.fcn="wo", data, m, k)
out4$psi


##################
}

\concept{dynamic treatment regimens}
\concept{adaptive treatment strategies}
\concept{personalized medicine}
\concept{continuous treatment}
\concept{generalized dynamic weighted ordinary least squares}
