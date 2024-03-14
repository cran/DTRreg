#' Process treatment and predicted treatment; setup Binary variable
#' 
#' @noRd
#' @param A A vector. The treatment coded as 0/1.
#' @param data A data.frame object. The complete data for the participants
#'   under analysis in the stage.
#' @param tx.var A character. The treatment variable.
#' @param models A list object. The models.
#' 
#' @return A list object containing 
#'   `$A`, the observed stage treatments, 
#'   `$A.hat`, the predicted treatment probabilities,
#'   `$cts`, a character object indicating binary treatment
#'   `$cts.obj`, an R6 object of class Binary,
#'   `$treat.mod.fitted`, a fitted glm object if user did not provide 
#'     treatment probabilities.
#'   
#' @importFrom stats glm predict
#' @include treatmentClasses.R
#' @keywords internal
.Ahat_binary <- function(A, data, tx.var, models, ...) {
  
  stopifnot(
    "`A` must be a vector" = !missing(A) && is.vector(A) && length(A) > 1L,
    "`data` must be a data.frame" = !missing(data) && is.data.frame(data) &&
      nrow(data) == length(A),
    "`tx.var` must be a character" = !missing(tx.var) && is.character(tx.var),
    "`models` must be a list" = !missing(models) && is.list(models) &&
      "treat" %in% names(models)
  )
  
  # Ensure that A is 0/1; convert if needed
  if (!is.factor(A))  A <- as.factor(A)
  A <- {unclass(A) - 1L} |> as.integer()
  attributes(A) <- NULL
  
  cts_obj <- Binary$new(tf.model = models$tf, 
                        blip.model = models$blip, 
                        tx.var = tx.var)

  data[[tx.var]] <- A
  fit <- tryCatch(stats::glm(models$treat, data, family = "binomial"),
                  error = function(e) {
                    stop("unable to fit treatment model\n\t",
                         e$message, call. = FALSE)
                  })
    
  Ahat <- stats::predict(fit, type = "response")
  if (length(Ahat) != length(A)) {
    stop("A.hat predictions not returned for all cases; NAs likely present\n",
         "contact developer", call. = FALSE)
  }
  
  list("A" = A, 
       "A.hat" = Ahat,
       "cts" = "bin", 
       "cts.obj" = cts_obj, 
       "tx.mod.fitted" = fit)
}

#' Process treatment and predicted treatment; setup MulitNom variable
#' 
#' @noRd
#' @param A A vector. The treatment as given in the original data.frame.
#' @param data A data.frame object. The complete data for the participants
#'   under analysis in the stage.
#' @param tx.var A character. The treatment variable.
#' @param models A list object. The models.
#' 
#' @return A list object containing 
#'   `$A`, the observed stage treatments, 
#'   `$A.hat`, a matrix, the predicted treatment probabilities for all treatments,
#'   `$cts.obj`, an R6 object,
#'   `$cts`, a character object indicating multi-nomial treatment
#'   `$treat.mod.fitted`, a fitted multinom object if user did not provide 
#'     treatment probabilities.
#'   
#' @importFrom nnet multinom
#' @importFrom stats predict
#' @importFrom utils capture.output
#'
#' @include treatmentClasses.R
#'
#' @keywords internal
.Ahat_multinom <- function(A, data, tx.var, models, ...) {

  stopifnot(
    "`A` must be a vector" = !missing(A) && {is.vector(A) || is.factor(A)} && 
      length(A) > 1L,
    "`data` must be a data.frame" = !missing(data) && is.data.frame(data) &&
      nrow(data) == length(A),
    "`tx.var` must be a character" = !missing(tx.var) && is.character(tx.var),
    "`models` must be a list" = !missing(models) && 
      all(c("tf", "blip", "treat") %in% names(models))
    )
  
  # Ensure that A is a factor variable
  if (!is.factor(A)) A <- factor(A)
  
  levs <- levels(A)
  n_levels <- length(levs)
  if (n_levels < 3L) {
    stop("provided treatment variable is not multinomial", call. = FALSE)
  }
  
  cts_obj <- MultiNom$new(tf.model = models$tf,
                          blip.model = models$blip,
                          tx.var = tx.var, 
                          tx.levels = levs)
  
  data[[tx.var]] <- A
  utils::capture.output(fit <- tryCatch(
      suppressMessages(nnet::multinom(models$treat, data)),
      error = function(e) {
        stop("unable to fit treatment model\n\t",
             e$message, call. = FALSE)
        }))
  Ahat <- stats::predict(fit, type = "probs")
    
  if (nrow(Ahat) != length(A)) {
    stop("A.hat predictions not returned for all cases; NAs likely present\n",
          "contact developer", call. = FALSE)
  }
  if (ncol(Ahat) == {n_levels - 1L}) {
    Ahat <- cbind(1.0 - rowSums(Ahat), Ahat)
  } else if (ncol(Ahat) != length(levels(A))) {
    stop("too few A.hat predictions returned; contact developer",
          call. = FALSE)
  }

  list("A" = A, 
       "A.hat" = Ahat,
       "cts" = "multinom", 
       "cts.obj" = cts_obj, 
       "tx.mod.fitted" = fit)
}

#' Process treatment and predicted treatment; setup Linear/Quadratic variable
#' 
#' @noRd
#' @param A A numeric vector object. The treatment as given in the original 
#'   data.frame. It is expected that participants that are not under analysis
#'   are removed.
#' @param data A data.frame object. The complete data for participants under
#'   analysis.
#' @param tx.var A character object. The treatment variable.
#' @param models A list object. The models.
#' @param tx.range A numeric vector. The min/max of the allowed treatment values.
#' @param tx.family A character object. Name of the family function passed to
#'   `glm` to fit treatment model.
#' 
#' @return A list object containing 
#'   `$A`, the observed stage treatments, 
#'   `$A.hat`, the predicted treatments,
#'   `$cts`, a character object indicating linear or quadratic blip
#'   `$cts.obj`, an R6 object,
#'   `$treat.mod.fitted`, a fitted lm object of user did not provide 
#'     treatment predictions.
#'
#' @importFrom stats glm predict
#'
#' @include treatmentClasses.R
#'
#' @keywords internal
.Ahat_cont <- function(A, data, tx.var, models, tx.range, tx.family, ...) {
  stopifnot(
    "`A` must be a vector" = !missing(A) && is.vector(A) && length(A) > 1L,
    "`data` must be a data.frame" = !missing(data) && is.data.frame(data) &&
      nrow(data) == length(A),
    "`tx.var` must be a character" = !missing(tx.var) && is.character(tx.var),
    "`models` must be a list" = !missing(models) && is.list(models) &&
      all(c("tf", "blip", "treat") %in% names(models)),
    "`tx.range` must be a numeric vector of length 2" = 
      is.numeric(tx.range) && is.vector(tx.range) && length(tx.range) == 2L,
    "`tx.family` must be a character in {'gaussian', 'Gamma'}" = 
      !missing(tx.family) && 
      {{is.character(tx.family) && tx.family %in% c("gaussian", "Gamma")} ||
      {inherits(tx.family, "family") && tx.family$family %in% c("gaussian", "Gamma")}}
  )
  
  data[[tx.var]] <- A
  fit <- tryCatch(stats::glm(models$treat, data, family = tx.family),
                  error = function(e) {
                    stop("unable to fit treatment model\n\t",
                         e$message, call. = FALSE)
                  })

  Ahat <- stats::predict(fit, type = "response")

  if (any(Ahat > tx.range[2L] | Ahat < tx.range[1L])) {
    warning("one or more predicted treatments are outside of allowed ",
            "treatment range; predicted treatment range: ", 
            format(range(Ahat), digits = 4L), "\n", 
            call. = FALSE)
  }
  
  psi_vars <- all.vars(models$blip)
  
  if (tx.var %in% psi_vars) {
    # continuous with quadratic blip
    cts_obj <- ContQuadraticBlip$new(tf.model = models$tf,
                                     blip.model = models$blip,
                                     treat.range = tx.range,
                                     tx.var = tx.var)
    cts <- "cts.q"
  } else {
    # continuous with linear blip
    cts_obj <- ContLinearBlip$new(tf.model = models$tf,
                                  blip.model = models$blip,
                                  treat.range = tx.range,
                                  tx.var = tx.var)
    cts <- "cts.l"
  }

  list("A" = A,
       "A.hat" = Ahat,
       "cts" = cts,
       "cts.obj" = cts_obj,
       "tx.mod.fitted" = fit)
}

#' Process stage treatment information
#' 
#' @noRd
#' @param stage.models A list object. Contains the 3 stage model named
#'   `$blip`, `$tf`, and `$treat`.
#' @param data A data.frame object. The complete covariate and treatment dataset
#'   for all participants under analysis in current step.
#' @param tx.var A character object. The column header of `data` pertaining to
#'   the treatment variable.
#' @param tx.range A numeric vector object of length 2. For continuous tx,
#'   the min and max of the allowed treatments. Ignored if not continuous tx.
#' @param tx.type A character object. Must be one of {"bin", "multi", "cont"}
#'   indicating if the treatments are binary, multinomial, or continuous,
#'   respectively.
#' 
#' @return A list object containing 
#'   `$A`, the observed tx variable (recoded if required); 
#'   `$A.hat`, the estimated or user specified tx probability; 
#'   `$treat.mod.fitted`, the regression object, 
#'   `$cts.obj` an R6 object
#'   `$cts`, a character object indicating the type of tx (bin/lin/quad); and
#'   Note that all vectors are of the length sum(stage.cases).
#'   
#' @keywords internal
.processTreatment <- function(stage.models, data, tx.var, tx.range, 
                              tx.type, tx.family) {
  
  stopifnot(
    "`stage.models` must be a list of formula" = !missing(stage.models) &&
      is.list(stage.models) && "treat" %in% names(stage.models) && 
      all(sapply(stage.models, inherits, what = "formula")),
    "`data` must be a data.frame" = !missing(data) && is.data.frame(data),
    "`tx.var` must be a character object" = !missing(tx.var) &&
      is.character(tx.var) && is.vector(tx.var) && length(tx.var) == 1L,
    "`tx.range` must be a numeric object or NULL" = !missing(tx.range) &&
      any({is.null(tx.range) | is.na(tx.range)}) ||
      (is.numeric(tx.range) && is.vector(tx.range) && length(tx.range) == 2L),
    "`tx.type` must be a character" = !missing(tx.type) && is.character(tx.type) &&
      tx.type %in% c("bin", "multi", "cont")  )
  
  # extract stage treatment variable
  A <- tryCatch(data[, tx.var],
                error = function(e) {
                  stop("unable to extract treatment variable ", 
                       tx.var, "\n\t", e$message, call. = FALSE)
                })
  
  if (tx.type == "bin") {
    
    # ensure that A is in fact binary
    # note that only the complete cases are provided, so NAs should not be
    # an issue here
    if (length(unique(A)) != 2L) {
      stop("treatment does not appear to be binary", call. = FALSE)
    }
    
    tx_info <- .Ahat_binary(A = A, 
                            data = data,
                            tx.var = tx.var,
                            models = stage.models)
    
  } else if (tx.type == "multi") {
    
    # ensure that 3 or more categories are present in the data
    # note that only the complete cases are provided, so NAs should not be
    # an issue here
    if (length(unique(A)) < 3L) {
      stop("treatment does not appear to be multinomial", call. = FALSE)
    }
    
    if (length(unique(A)) > 10L) {
      warning("a large number of treatment options identified; ",
              "verify that tx is multinomial\n", call. = FALSE)
    }
    
    tx_info <- .Ahat_multinom(A = A, 
                              data = data,
                              tx.var = tx.var,
                              models = stage.models)
    
  } else {
    
    tx_info <- .Ahat_cont(A = A, 
                          data = data,
                          tx.var = tx.var,
                          models = stage.models,
                          tx.range = tx.range,
                          tx.family = tx.family)
  } 
  tx_info
}