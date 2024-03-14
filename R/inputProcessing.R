#' Verify and group models
#'
#' @noRd
#' @param blip.mod A list of formulae.
#' @param treat.mod A list of formulae.
#' @param tf.mod A list of formulae
#'
#' @return A list of the input formula grouped by stage
#' 
#' @keywords internal
.modelsTest <- function(blip.mod, treat.mod, tf.mod) {
  
  # allow for the possibility of providing single decision point analysis
  # as single formula rather than one element lists
  if (!is.list(treat.mod)) treat.mod <- list(treat.mod)
  if (!is.list(blip.mod)) blip.mod <- list(blip.mod)
  if (!is.list(tf.mod)) tf.mod <- list(tf.mod)

  stopifnot(
    "`treat.mod` must be a list of formulae of the form LHS ~ RHS" = 
      is.list(treat.mod) && all(sapply(treat.mod, inherits, what = "formula")) &&
      all(lengths(treat.mod) == 3L),
    "`blip.mod` must be a list of formulae of the form ~ RHS" = 
      is.list(blip.mod) && all(sapply(blip.mod, inherits, what = "formula")) &&
      all(lengths(blip.mod) == 2L),
    "`tf.mod` must be a list of formulae of the form ~ RHS" = 
      is.list(tf.mod) && all(sapply(tf.mod, inherits, what = "formula")) &&
      all(lengths(blip.mod) == 2L)
  )
  
  if (length(treat.mod) != length(blip.mod) || 
      length(treat.mod) != length(tf.mod)) {
    stop("`treat.mod`, `blip.mod`, and `tf.mod` must be of ",
         "the same length", call. = FALSE)
  }
  
  models <- mapply(list, 
                   "blip" = blip.mod,
                   "treat" = treat.mod,
                   "tf" = tf.mod,
                   SIMPLIFY = FALSE)
  models
  
}


#' Create, verify, or incorporate missingness models
#'
#' Note: Call this function only if user specifies that missingness should
#'   be modeled
#' 
#' @noRd
#' @param models A list. One element for each stage. Each element is itself
#'   a list containing "blip", "treat", and "tf". If interactive session
#'   was used to construct this list, it may also contain "cens".
#' @param missing.mod NULL or a list of formula provided as input to dtr
#' 
#' @return Input `models` possibly augmented with missingness model
#' 
#' @importFrom stats as.formula
#'
#' @keywords internal
.constructMissingness <- function(models, missing.mod) {
  
  K <- length(models)
  
  if (is.null(models[[1L]]$cens) && is.null(missing.mod)) {
    # No missingness model provided as input or specified in interactive
    # session, create model from full covariate history
    # NOTE this logic does not maintain interaction terms of models
    
    missing.mod <- vector(mode = "list", length = K)
    missing.mod[[1L]] <- ~ 1
    
    vars <- lapply(models[[1L]], all.vars) |> unlist() |> unique()
    for (i in 2L:K) {
      # missing model defined by all previous stage model covariates
      missing.mod[[i]] <- stats::as.formula(paste("~", paste(vars, collapse = "+")))
      message("Stage ", i, " missingness model defined as ",
              deparse(missing.mod[[i]]))
      
      # include ith stage model covariates in vector of covariates
      vars <- c(vars, lapply(models[[i]], all.vars) |> unlist()) |> unique()
    }
  }
  
  # allow for the possibility of a single stage scenario with formula provided
  if (inherits(missing.mod, "formula")) missing.mod <- list(missing.mod)
  
  if (is.list(missing.mod)) {
    # `missing.mod` either constructed above or provided by user
    if (is.null(models[[1L]]$cens)) {
      # if no censoring/missingness model stored in models list add `missing.mod`
      if (length(missing.mod) == K) {
        # ensure that all elements of the list are formulae
        tst <- sapply(missing.mod, inherits, what = "formula") |> all()
        # add to models list
        if (tst) {
          models <- mapply(function(x, y) {x[["cens"]] <- y; x},
                           models,
                           missing.mod, SIMPLIFY = FALSE)
        } else {
          stop("`missing.mod` must be a list of formulae" , call. = FALSE)
        }
      } else {
        stop("`missing.mod` is not of appropriate length, K = ", K, 
             "length(missing.mod) = ", length(missing.mod), call. = FALSE)
      }
    } else {
      warning("missingness models specified in interactive mode;",
              "input `missing.mod` will be ignored", call. = FALSE)
    }
  } else {
    stop("`missing.mod` must be a list of models, 1 for each decision pt",
         call. = FALSE)
  }
  models
}

#' Verify treatment related inputs
#'
#' @noRd
#' @param weight A character. The type of treatment weights
#' @param treat.type A character. The type of treatment variable (binary,
#'   multinomial, continuous).
#' @param n.bins An integer. If continuous treatment at weight function is pom
#'   the number of discrete bins to discretize the treatment.
#' @param treat.wgt.man NULL or a list of numeric vectors. User provided treatment
#'   weights.
#' @param treat.range NULL or a two element vector. If treatment is continuous,
#'   the allowed range of treatment values.
#' @param treat.fam NULL, character, or family object. If treatment is 
#'   continuous, the family object for obtaining parameter estimates
#' @param K An integer. The number of stages
#' @param models A list. The models.
#' @param data A data.frame.
#'
#' @return A list containing all verified treatment related information
#' 
#' @importFrom stats family Gamma gaussian terms
#' 
#' @keywords internal
.treatmentTest <- function(weight, treat.type, n.bins, treat.wgt.man,
                           treat.range, treat.fam, K, models, data) {
  local_obj <- list()
  
  local_obj$n.bins <- NA_integer_
  if (weight %in% c("qpom", "wo")) {
    if (treat.type == "cont") {
      if (is.numeric(n.bins) && is.vector(n.bins) && length(n.bins) == 1L &&
          isTRUE(all.equal(n.bins, round(n.bins)))) {
        local_obj$n.bins <- n.bins
      } else {
        stop("`n.bins` must be an integer", call. = FALSE)
      }
    } else if (treat.type != "multi") {
      stop("binary treatments are not appropriate for `weights` 'qpom' or 'wo'", 
           call. = FALSE)
    }
  }
  
  if (is.null(treat.wgt.man)) {
    if (weight %in% c("manual", "manual.with.censor")) {
      stop("manual weights must be provided if `weight` = '", weight, "'",
           call. = FALSE)
    }
    local_obj$tx.wgt.man <- NA
  } else if (is.list(treat.wgt.man)) {
    if (!{weight %in% c("manual", "manual.with.censor")}) {
      stop("manual weights cannot be provided if `weight` = '", weight, "'",
           call. = FALSE)
    }
    # must be a list of length K
    if (length(treat.wgt.man) != K || 
        any(!sapply(treat.wgt.man, is.numeric)) ||
        any(!sapply(treat.wgt.man, is.vector, mode = "numeric"))) {
      stop("if provided, `treat.wgt.man` must be a list of K numeric vectors",
           call. = FALSE)
    }
    
    # each element must be a vector of length n
    if (any(lengths(treat.wgt.man) != nrow(data))) {
      stop("all participants must be included in each element of `treat.wgt.man`",
           call. = FALSE)
    }
    
    # all weights must be positive (though they can be NA)
    if (lapply(treat.wgt.man, "<", 1e-8) |> unlist() |> any(na.rm = TRUE)) {
      stop("treatment weights must be positive", call. = FALSE)
    }
    
    local_obj$tx.wgt.man <- treat.wgt.man
  } else {
    stop("`treat.wgt.man` must be NULL or a list", call. = FALSE)
  }

  # this element facilitates blending of survival
  local_obj$dependent.vars <- list()
  local_obj$dependent.vars$treat <- lapply(models,
                                           function(x) {
                                             as.character(x$treat)[2L]
                                             }) |> unlist()
  
  if (is.null(treat.range) && treat.type == "cont") {
    As <- tryCatch(data[, local_obj$dependent.vars$treat],
                   error = function(e) {
                     stop("unable to retrieve treatment variables from `data`\n\t",
                          e$message, call. = FALSE)
                   })
    local_obj$tx.range <- lapply(As, range, na.rm = TRUE)
    message("no treatment range specified, [min, max] of observed ",
            "treatments used; \n")
    for (i in 1L:length(local_obj$tx.range)) {
      message("Stage ", i, " [", 
              paste(format(local_obj$tx.range[[i]], digits = 4), collapse = ", "), "]\n")
    }
    
  } else if (treat.type == "cont") {
    if (is.numeric(treat.range) && 
        is.vector(treat.range) &&
        length(treat.range) == 2L) {
      local_obj$tx.range <- vector(mode = "list", length = K)
      local_obj$tx.range <- lapply(local_obj$tx.range, function(x) { treat.range })
    } else if (is.list(treat.range) && length(treat.range) == K && 
               all(lengths(treat.range) == 2L)) {
      local_obj$tx.range <- treat.range
    } else {
      stop("if provided, `treat.range` must be a list of length K each element ",
           "a numeric vector of length 2",
           call. = FALSE)
    }
  } else {
    local_obj$tx.range <- as.list(rep(NA, K))
  }
  
  if (treat.type == "cont") {
    if (is.null(treat.fam)) {
      local_obj$tx.family <- stats::gaussian(link = "identity")
    } else if (is.character(treat.fam) && length(treat.fam) == 1L) {
      if (grepl("gaussian", treat.fam)) {
        local_obj$tx.family <- stats::gaussian(link = "identity")
      } else if (grepl("Gamma", treat.fam)) {
        local_obj$tx.family <- stats:: Gamma(link = "log")
      } else {
        stop("if specificed as a character, `treat.fam` must be one of ",
             "'gaussian' or 'Gamma'", call. = FALSE)
      }
    } else if (inherits(treat.fam, "family")) {
      if (treat.fam$family == "gaussian" && treat.fam$link != "identity") {
        stop("unsupported family specified for treatment model", call. = FALSE)
      } else if (treat.fam$family == "Gamma" && treat.fam$link != "log") {
        stop("unsupported family specified for treatment model", call. = FALSE)
      }
      local_obj$tx.family <- treat.fam
    } else {
      stop("`treat.fam` must be NULL, a character object, or a family object",
           call. = FALSE)
    }
  } else {
    local_obj$tx.family <- NA
  }
  local_obj
}

#' Verify status variable
#' 
#' @noRd
#' @param status A character vector. The status variable name repeated for each
#'   stage or NA if data is not censored
#' @param data A data.frame. The full data.
#' 
#' @return Input `data` possibly updated to ensure that status variable is
#'   binary 0/1.
#'   
#' @keywords internal
.statusTest <- function(status, data) {
  
  if (!all(status == status[1L])) {
    stop("more than 1 status variable has been provided", call. = FALSE)
  }
  status <- status[1L]
  
  if (is.na(status)) {
    message("data is not censored")
  } else {
    if (any(is.na(data[, status]))) { 
      stop("missing status values are not allowed", call. = FALSE)
    }
    
    d.bin <- {data[, status] |> unique() |> length()} <= 2L
    if (!d.bin) {
      stop("status variables must be binary", call. = FALSE)
    }
      
    if (!all(data[, status] %in% c(0L, 1L))) {
      # recode as 0/1
      data[, status] <- {factor(data[, status]) |> unclass()} - 1L 
    }
  }
  data
}

#' Verify and accept user specified bootstrap controls if provided
#' 
#' @noRd
#' @param bootstrap.controls A list. Optional bootstrap contraols
#' @param n An integer. The number of participants
#'
#' @returns The full bootstrap.controls list
#' 
#' @keywords internal
.bootstrapControlsTest <- function(bootstrap.controls, n) {

  if (is.list(bootstrap.controls)) {
    
    # ensure that list contains only expected elements
    allowed_names <- c("B", "M", "type", "truncate", "verbose", "interrupt")
    if (!all(names(bootstrap.controls) %in% allowed_names)) {
      stop("`bootstrap.control` has unrecognized element names", call. = FALSE)
    }
    
    # define default values for those not provided as input
    default_values <- list("B" = 100L, 
                           "M" = n, 
                           "type" = "standard", 
                           "truncate" = 0.0, 
                           "verbose" = FALSE, 
                           "interrupt" = FALSE)
    
    # replace defaults with provided values
    default_values[names(bootstrap.controls)] <- bootstrap.controls
    
    boot.controls <- default_values
  } else {
    stop("if provided, `bootstrap.controls` must be a named list", 
         call. = FALSE)
  }

  if (is.null(boot.controls$M) || boot.controls$M < 1e-8) boot.controls$M <- n

  stopifnot(
    "`bootstrap.controls$B` must be a positive integer" = 
      is.numeric(boot.controls$B) && 
      isTRUE(all.equal(boot.controls$B, round(boot.controls$B))) &&
      boot.controls$B > 0L,
    "`bootstrap.controls$M` must be a positive integer" = 
      is.numeric(boot.controls$M) && 
      isTRUE(all.equal(boot.controls$M, round(boot.controls$M))) &&
      boot.controls$M > 0L,
    "`bootstrap.controls$type must be one of {'standard', 'empirical', 'normal'}" =
      is.character(boot.controls$type) && length(boot.controls$type) == 1L &&
      boot.controls$type %in% c("standard", "empirical", "normal"),
    "`bootstrap.controls$truncate must be [0.0, 0.5]" = 
      is.numeric(boot.controls$truncate) && 
      length(boot.controls$truncate) == 1L &&
      boot.controls$truncate > -1e-8 && boot.controls$truncate <= 0.5,
    "`bootstrap.controls$verbose` must be logical" = 
      is.logical(boot.controls$verbose),
    "`bootstrap.controls$interrupt` must be logical" = 
      is.logical(boot.controls$interrupt)
  )
  boot.controls
}

#' Verify variance related inputs
#' 
#' @noRd
#' @param var.estim A character. One of "none", "sandwich", "bootstrap"
#' @param method A character. One of "dwols" "gest"
#' @param bootstrap.controls NULL or a list
#' @param n An integer. The number of participants
#'
#' @returns NA if bootstrap not selected; otherwise a list of bootstrap
#'   controls
#'   
#' @keywords internal
.varestimTest <- function(var.estim, method, bootstrap.controls, n) {
  
  if (var.estim == "bootstrap") {
    boot.controls <- .bootstrapControlsTest(bootstrap.controls, n)
  } else {
    boot.controls <- NA
  }
  
  boot.controls
}