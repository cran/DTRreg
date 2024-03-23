#' DTR estimation and inference for time-to-event data using DWSurv
#' 
#' Dynamic treatment regimen estimation and inference via dynamic weighted 
#'   survival modeling (DWSurv).  Inference for the blip estimators with single- 
#'   and multi-stage data.
#'   
#' The function \code{DWSurv()} allows estimating an optimal dynamic treatment regime
#'   from multi-stage trials or observational data when the outcome of interest
#'   is survival time subject to right-censoringg. The dynamic weighted survival 
#'   modeling (DWSurv) algorithm is implemented.  The method focuses on 
#'   estimating the parameters of the blip: a model of the difference in 
#'   expected outcome under the observed treatment and some reference treatment 
#'   (usually a control) at a given stage, assuming identical histories and 
#'   optimal treatment thereafter.
#'   
#' The method requires the specification of four models for each stage of the 
#'   analysis: a treatment model (conditional mean of the treatment variable), 
#'   a censoring model, a treatment-free model (conditional mean of outcome 
#'   assuming only reference treatments are used), and a blip model.  Only the 
#'   blip model must be correctly specified (or over-specified), with consistent 
#'   parameter estimates obtainable if at least one of the treatment-free or the 
#'   treatment and censoring models are correctly specified.  Note that all of 
#'   these must be specified as lists of formula objects, even if only one stage
#'   of treatment is considered.
#'   
#' Note that as is conventional, it is assumed a larger survival time is 
#'   preferred (which can be easily achieved via transformation of your data if 
#'   necessary).
#'   
#' Several treatment weight function options have been implemented within the 
#'   package:
#'   \itemize{ 
#'     \item "none": No treatment weights applied. If \code{method = "dWOLS"}, this 
#'       selection results in the implementation of Q-learning, modified 
#'       slightly to use the dWOLS style pseudo-outcome 
#'       (computed using the observed outcome modified by the estimated 
#'        treatment effect) rather than the traditional Q-learning outcome 
#'       (predicted based on model only, rather than observed outcome with 
#'       treatment effect).  
#'     \item "ipw": weights based on the inverse probability of 
#'       treatment. For binary treatments, a logistic regression is used.
#'       For multinomial, a multinomial log-linear model is fit using 
#'       \code{\link[nnet]{multinom}}. For continuous treatments, a GLM with the specified 
#'       family and link function provided in the \code{treat.fam} argument is used.
#'     \item "cipw": inverse probability of treatment weights as described for
#'       "ipw" and capped at the 99th percentile of the observed weights.
#'     \item "qpom": weights based on the stabilized inverse 
#'       probability of treatment applied to the categorized (into n.bins bins) 
#'       continuous doses or multinomial treatments; probabilities are 
#'       calculated using a proportional odds model.
#'       This weight is appropriate only for continuous and multinomial treatments.
#'     \item "wo": overlap weights for the categorized continuous doses 
#'       or multinomial treatments (Li and Li, 2019).
#'       This weight is appropriate only for continuous treatments.
#'     \item "abs": Absolute difference \eqn{|A - E[A|...]|}{|A - E[A|...]|}. 
#"       This weight is 
#'       appropriate only for binary treatments.
#'     \item "manual": User provides treatment weights through input 
#'     \code{treat.wgt.man}.
#'     \item "manual.with.censor": User provides combined treatment * censoring
#'       weights through input \code{treat.wgt.man}. Note that `cens.mod` should
#'       be specified with the event indicator on the right-hand side of the 
#'       formula (e.g., \code{~ status}).
#'    }
#'   
#' @inheritParams DTRreg
#' @param method The DTR method to be used, choose "dwols" for dynamic WOLS, 
#'   or "qlearn" for Q-learning.
#' @param time A list of formula specifying the survival time variable for each 
#'   stage in order. The time variable should be specified on the right hand 
#'   side of the formula. No dependent variable should be specified. The list 
#'   should be as long as the number of stages.
#' @param cens.mod A list of formula objects specifying the censoring
#'   model for each stage in order. The event indicator, which takes value 1 if
#'   an event was observed and 0 otherwise, should be included as the dependent
#'   variable and should be the same across stages. In the absence of censoring
#'   or if censoring weights are provided by the user through `treat.wgt.man`, 
#'   (i.e., \code{weight = 'manual.with.censor'})
#'   one still needs to specify an event indicator on the right-hand
#'   side of the formula and leave the left-hand side empty (see example below).
#' @param treat.wgt.man NULL or a list of vectors of known treatment 
#'   (or treatment * censoring) weights can be 
#'   specified to be used instead of hard-coded treatment weight options.
#'   The \eqn{i^{th}}{ith} element of the list contains the multiplicative weights 
#'   for the \eqn{i^{th}}{ith} stage. Each vector must be of length \eqn{n}{n}, 
#'   the number of participants. Used only for \code{method = "dwols"}. If
#'   providing the treatment * censoring weights, \code{cens.mod = NA} must
#'   be used.
#'    
#' @return An object of class \code{DWSurv}, a list including elements
#'     \item{K: }{The number of decision points.}
#'     \item{beta: }{A list. The ith element contains the parameter estimates of
#'       the ith stage treatment-free model.}
#'     \item{psi: }{A list. The ith element contains the parameter estimates of
#'       the ith stage blip model.}
#'     \item{covmat: }{A list. The ith element contains covariance matrix of 
#'       the ith stage blip parameter estimates.}
#'     \item{nonreg: }{Non-regularity estimates.}
#'     \item{setup: }{A list detailing the input parameter settings used for the
#'       analysis}
#'       \itemize{
#'         \item models: A list of the models used for the analysis.
#'         \item method: The parameter estimation method.
#'         \item var.estim: The variance estimation method.
#'         \item cc.modeled: If TRUE, missing data was modeled. If FALSE, cases
#'            with missing data were removed from the analysis.
#'         \item tx.weight: The treatment weighting used for the analysis.
#'         \item tx.type: Treatment was binary, multinomial, or continuous.
#'         \item n.bins: The number of bins (levels) used for categorizing
#'           continuous doses when \code{tx.weight = "wo"} or 
#'           \code{tx.weight = "qpom"}.
#'         \item tx.wgt.man: Any user provided treatment weights.
#'         \item tx.range: For continuous treatments, the range of allowed
#'           treatment values.
#'         \item tx.family: The description of the dose distribution along 
#'           with the link function used in the continuous treatment model.
#'         \item boot.controls: A list of the bootstrap controls.
#'         \item type: The type of effect. Dynamic treatment regime or treatment
#'           effect.
#'       }
#'     \item{training_data: }{A list containing the training data.}
#'       \itemize{
#'         \item data: The covariates and treatment data.
#'         \item outcome: The outcome of interest.
#'         \item A: The treatment variables, possibly recoded to adhere to internal
#'           code requirements.
#'       }
#'     \item{analysis: }{A list containing the primary results of each stage analysis.}
#'       \itemize{
#'         \item n: The number of participants included in the stage analysis.
#'         \item last.stage: The last stage each participant was included in 
#'           the analysis.
#'         \item prob.cens: The complete case probabilities.
#'         \item cens.mod.fitted: The regression objects returned for estimating
#'           the complete case probabilities.
#'         \item cens.wgt: The complete case weights.
#'         \item cts: The treatment type at each stage.
#'         \item tx.mod.fitted: The regression objects returned for estimating
#'           the treatment probabilities.
#'         \item A.hat: The estimated or provided treatment probabilities.
#'         \item tx.wgt: The treatment weights.
#'         \item outcome.fit: The regression objects returned for each stage outcome
#'           regression.
#'         \item Y: The pseudo-outcomes.
#'         \item regret: Estimates of the regret for each subject based on observed 
#'           treatment and blip parameter estimates.
#'         \item opt.treat: Optimal treatment decisions for each subject at each 
#'           stage of treatment.
#'         \item opt.Y: Predicted optimal outcome under recommended regimen.
#'      }
#'     \item{call: }{The original function call.}
#'     
#'   The functions \code{coef()}, \code{predict()} and 
#'   \code{confint()} may be used with such 
#'   model objects. The first two have specific help files for their 
#'   implementation, while \code{confint()} is used in the same way as 
#'   the standard 
#'   \code{\link[stats]{confint}()} command, with the exception of the \code{parm} 
#'   option, which is not available.
#'
#' @references
#' Simoneau, G., Moodie, E. E. M., Wallace, M.P., Platt, R. W. (2020) Optimal 
#'   Dynamic Treatment Regimes with Survival Endpoints: Introducing DWSurv in the 
#'   R package DTRreg. \emph{Journal of Statistical Computation and Simulation}.
#'   \bold{90}, 2991-3008. (doi:10.1080/00949655.2020.1793341)
#'   
#' Simoneau, G., Moodie, E. E. M., Nijjar, J. S., Platt, R. W. (2019) Estimating 
#'   Optimal Dynamic Treatment with Survival Outcomes. \emph{Journal of the 
#'   American Statistical Association}, \bold{115}, 1531-1539 
#'   (doi:10.1080/01621459.2019.1629939).
#'   
#' Wallace, M. P., Moodie, E. E. M., Stephens, D. A. (2017) Dynamic Treatment 
#'   Regimen Estimation via Regression-Based Techniques: Introducing {R} Package 
#'   {DTRreg}. \emph{Journal of Statistical Software} \bold{80}(2), 1--20 
#'   (doi:10.18637/jss.v080.i02).
#'   
#' Simoneau, G., Moodie, E. E. M., Nijjar, J. S., and Platt, R. W. (2020)
#'   Finite Sample Variance Estimation for Optimal Dynamic Treatment
#'   Regimes of Survival Outcomes. \emph{Statistics in Medicine} \bold{39},
#'   4466-4479.
#'   
#' Efron, B., and Tibshirani, R. (1986)
#'   Bootstrap Methods for Standard Errors, Confidence Intervals, and Other 
#'   Measures of Statistical Accuracy \emph{Source: Statistical Science} \bold{1}
#'   54-75.
#'   
#' @examples
#' #### example single run of a 2-stage DWSurv analysis
#' data(twoStageCens)
#' mod <- DWSurv(time = list(~ T1, ~ T2), 
#'               blip.mod = list(~ X11, ~ X21), 
#'               treat.mod = list(A1 ~ X11, A2 ~ 1), 
#'               tf.mod = list(~ X11 + X12, ~ X21 + X22 + X11), 
#'               cens.mod = list(delta ~ 1, delta ~ X11), 
#'               var.estim = "sandwich", 
#'               data = twoStageCens)
#' mod
#'   
#' #### example in the absence of censoring
#' data(twoStageSurv)
#' mod_nocensoring <- DWSurv(time = list(~ T1, ~ T2), 
#'                           blip.mod = list(~ X11, ~ X21), 
#'                           treat.mod = list(A1 ~ X11, A2 ~ 1), 
#'                           tf.mod = list(~ X11 + X12, ~ X21 + X22 + X11), 
#'                           cens.mod = list(~ delta, ~ delta), 
#'                           var.estim = "sandwich", 
#'                           data = twoStageSurv)
#' mod_nocensoring
#' 
#' @concept dynamic treatment regimens
#' @concept adaptive treatment strategies
#' @concept personalized medicine
#' @concept dynamic weighted survival modeling
#' @concept accelerated failure time
#'
#' @importFrom stats as.formula family Gamma gaussian get_all_vars model.frame terms
#' @importFrom stats na.pass
#' @include dtrProcedure.R inputProcessing.R interactive.R
#' @export
DWSurv <- function(time, 
                   blip.mod, treat.mod, tf.mod, cens.mod,
                   data = NULL, 
                   method = c("dwols", "qlearn"), 
                   interactive = FALSE,
                   treat.type = c("bin", "multi", "cont"), 
                   treat.fam = gaussian(link = "identity"), 
                   weight = c("abs", "ipw", "cipw", "qpom", "wo", "none", 
                              "manual", "manual.with.censor"),
                   n.bins = 3L, 
                   treat.range = NULL, 
                   treat.wgt.man = NULL,
                   var.estim = c("none", "bootstrap", "sandwich"), 
                   bootstrap.controls = list(B = 100L, 
                                             M = 0L, 
                                             type = "standard", 
                                             truncate = 0.0, 
                                             verbose = FALSE, 
                                             interrupt = FALSE), 
                   dtr = TRUE,
                   full.cov = FALSE) {
  
  # if interactive mode chosen, build standardized input from user input
  if (interactive) {
    # returned list contains 
    # $models A list. The models (blip, tf, treat) grouped according to 
    #   decision point
    # $method A character. The DTR method (dwols).
    # $var.estim A character. The covariance matrix estimation method 
    #   ("none", "bootstrap", "sandwich")
    # $time A list of formulae. RHS is the time variable for each stage.
    # $K An integer. The number of decision points.
    obj <- .interactive(isSurvival = TRUE, prompt.cens = TRUE, data = data)
    
    time <- obj$time
    obj$time <- NULL
    
  } else {
    if (inherits(time, "formula")) time <- list(time)
    if (inherits(cens.mod, "formula")) cens.mod <- list(cens.mod)
    stopifnot(
      "`time`, `blip.mod`, `treat.mod`, `tf.mod`, and `cens.mod` must be provided" =
        !missing(time) && !missing(blip.mod) & !missing(treat.mod) & 
        !missing(tf.mod) && !missing(cens.mod),
      "`time` must be a list of formulae of the form ~ RHS" = 
        is.list(time) && all(sapply(time, inherits, what = "formula")) &&
        all(lengths(time) == 2L),
      "`cens.mod` must be a list of formulae" = !is.null(cens.mod) &&
        {is.list(cens.mod) && all(sapply(cens.mod, inherits, what = "formula")) &&
         (all(lengths(cens.mod) == 2L) || all(lengths(cens.mod) == 3L))},
      "`data` must be NULL or a data.frame" = is.null(data) || is.data.frame(data)
    )
    
    obj <- list()
    
    # group models by stage
    # this is the model testing for non-survival, so it is missing tests of 
    # time and censoring
    obj$models <- .modelsTest(blip.mod, treat.mod, tf.mod)

    # tests for time and censoring    
    if (!is.list(time)) time <- list(time)
    if (!is.list(cens.mod)) cens.mod <- list(cens.mod)
    
    if (length(obj$models) != length(time) || 
        length(obj$models) != length(cens.mod)) {
      stop("`treat.mod`, `blip.mod`, `tf.mod`, `cens.mod`, and `time` must be of ",
           "the same length", call. = FALSE)
    }
    
    # add appropriate censoring model to each set of models
    obj$models <- mapply(function(m, cens) {
                           m[["cens"]] <- cens
                           m
                         }, obj$models, cens.mod, SIMPLIFY = FALSE)
    
    obj$method <- match.arg(method)
    obj$var.estim <- match.arg(var.estim)
    
    obj$K <- length(obj$models)
  }
  
  if (obj$method == "qlearn") {
    method <- "qlearn"
    weight <- "none"
    obj$method <- "dwols"
  } else {
    method <- obj$method
  }
  
  # time formula should only contain the survival time as the RHS
  if (any(lengths(lapply(time, all.vars)) != 1L)) {
    stop("the survival time must be supplied as a formula of the form ~ RHS, ",
         "where the RHS contains only the survival time variable", 
         call. = FALSE)
  }
  
  # If all censoring models are of length 2, no dependent variable given and
  # RHS provides the status variable. If of length 3, censoring is to be modeled
  # and LHS is status variable and RHS is model covariates
  # NOTE this is all or nothing. Can't model some stages but not others
  if (!is.null(cens.mod)) {
    obj$censoring.modeled <- length(obj$models[[1L]]$cens) == 3L
  } else {
    obj$censoring.modeled <- FALSE
  }
  
  # get all vars from the models
  data_vars <- obj$models |> 
    unlist() |> 
    lapply(FUN = all.vars) |>
    unlist()
  data_vars <- c(data_vars, 
                 time |> 
                   unlist() |> 
                   lapply(FUN = all.vars) |>
                   unlist())
  data_vars <- data_vars |> unique() |>
    paste(collapse = " + ")
  data_vars <- stats::as.formula(paste("~", data_vars))
  
  if (is.null(data)) {
    # if no data set, create one from environment using the models
    
    # get all vars from the models
    obj$data <- tryCatch(stats::get_all_vars(data_vars),
                         error = function(e) {
                           stop("unable to retrieve all covariates from the ",
                                "environment\n\t", e$message, call. = FALSE)
                         })
  } else if (is.data.frame(data)) {
    obj$data <- tryCatch(stats::model.frame(data_vars, data, 
                                            na.action = stats::na.pass),
                         error = function(e) {
                           stop("unable to retrieve all covariates from the ",
                                "`data`\n\t", e$message, call. = FALSE)
                         })
    rm(data, data_vars)
  } else {
    stop("if provided, `data` must be a data.frame", call. = FALSE)
  }
  
  rownames(obj$data) <- seq_len(nrow(obj$data))

    
  ### Treatment related inputs
  
  obj$tx.weight <- match.arg(weight)
  obj$tx.type <- match.arg(treat.type)
  
  if(obj$tx.weight == "abs" && obj$tx.type != "bin") {
    stop("cannot select `weight = 'abs'` for multinomial or continuous treatments",
         call. = FALSE)
  }
  
  obj$manual.censor.weight <- obj$tx.weight == "manual.with.censor"
  # can't have weight = "manual.with.censor" and specify a full censoring model
  if (obj$censoring.modeled && obj$manual.censor.weight) {
    stop("censoring models should be of the form ~ status when ",
         "`weight` = 'manual.with.censor'",
         call. = FALSE)
  }
  
  obj <- c(obj,
           .treatmentTest(weight = obj$tx.weight, 
                          treat.type = obj$tx.type, 
                          n.bins = n.bins, 
                          treat.wgt.man = treat.wgt.man, 
                          treat.range = treat.range, 
                          treat.fam = treat.fam, 
                          K = obj$K, 
                          data = obj$data,
                          models = obj$models))
  
  if (obj$var.estim != "none" && is.list(obj$tx.mod.man)) {
    warningMsg <- paste("Treatment weights were supplied as an argument",
                        "(generated outside the DWSurv function);",
                        "standard errors do not account for estimation of the",
                        "treatment score model.\n")
    warning(warningMsg, call. = FALSE)
    addWarning <- TRUE
  } else {
    addWarning <- FALSE
  }
  
  ### Status and Time variables  
  obj$dependent.vars$time <- lapply(time, 
                                    function(x) {
                                      rownames(attr(stats::terms(x), "factors"))
                                    }) |> unlist()
  
  # Note that if user opted for interactive mode, the censoring model will
  # be NULL if data is not censored. 
  obj$dependent.vars$status <- lapply(obj$models, 
                                      function(x) {
                                        if (is.null(x$cens)) return(NA)
                                        as.character(x$cens)[2L]
                                      }) |> unlist() |> unique()

  # Ensure that status variable is 0/1 if data is censored
  obj$data <- .statusTest(obj$dependent.vars$status, obj$data)
  
  ### Variance related inputs
  
  obj$boot.controls <- .varestimTest(var.estim = obj$var.estim, 
                                     method = obj$method, 
                                     bootstrap.controls = bootstrap.controls, 
                                     n = nrow(obj$data))
  
  stopifnot("`full.cov` must be logical" = is.logical(full.cov))
  obj$full.cov <- full.cov
  
  ### Analysis details
  
  if (is.logical(dtr)) {
    obj$type <- ifelse(dtr, "DTR", "effect")
  } else {
    stop("`dtr` must be a logical", call. = FALSE)
  }
  
  obj <- .dtrProcedure(obj = obj, quiet = FALSE, isSurvival = TRUE)
  
  obj$call <- match.call()
  
  analysis_settings <- c("models" = "models", 
                         "method" = "method", 
                         "var.estim" = "var.estim", 
                         "boot.controls" = "boot.controls", 
                         "cens.modeled" = "censoring.modeled",
                         "tx.weight" = "tx.weight", 
                         "tx.type" = "tx.type", 
                         "n.bins" = "n.bins", 
                         "tx.wgt.man" = "tx.wgt.man",
                         "tx.range" = "tx.range", 
                         "tx.family" = "tx.family", 
                         "type" = "type",
                         "tx.vars" = "treat.vars")
  setup <- obj[analysis_settings]
  obj[analysis_settings] <- NULL
  names(setup) <- names(analysis_settings)
  obj$setup <- setup
  
  training_data <- c("data", "outcome", "A")
  training <- obj[training_data]
  obj[training_data] <- NULL
  obj$training_data <- training
  
  analysis_names <- c("n" = "n", 
                      "last.stage" = "last.stage", 
                      "prob.cens" = "prob.complete.case", 
                      "cens.mod.fitted" = "cens.mod.fitted",
                      "cens.wgt" = "cens.wgt", 
                      "cts" = "cts", 
                      "tx.mod.fitted" = "tx.mod.fitted", 
                      "A.hat" = "A.hat", 
                      "tx.wgt" = "tx.wgt",
                      "outcome.fit" = "outcome.fit", 
                      "Y" = "Y",
                      "regret" = "regret",
                      "opt.treat" = "opt.treat", 
                      "opt.Y" = "opt.Y")
  analysis <- obj[analysis_names]
  obj[analysis_names] <- NULL
  names(analysis) <- names(analysis_names)
  obj$analysis <- analysis
  
  if (is.null(obj$covmat)) obj$covmat <- NA
  if (is.null(obj$nonreg)) obj$nonreg <- NA

  for_return <- c("K", "beta", "psi", "covmat", "psi.boot", "nonreg", 
                  "setup", "training_data", "analysis", "call")
  obj <- obj[for_return[for_return %in% names(obj)]]
  
  if (addWarning) obj$warn <- warningMsg
  
  obj$setup$method <- method
  
  # return
  class(obj) <- c("DWSurv", "DTRreg", class(obj))
  obj
}