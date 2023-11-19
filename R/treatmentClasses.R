#' Specialized function for Binary treatments
#' 
#' @noRd
#' @import R6
#' @keywords internal
Binary <- R6::R6Class(
  "Binary",
  public = list(
    blip.model = NULL,
    # the full outcome model
    full.model = NULL,
    
    initialize = function(tf.model, blip.model, tx.var, ...) {
        private$tx.var <- tx.var
      
        # add the treatment variable dependence to the blip formula
        self$blip.model <- paste("~ - 1 +", tx.var, "+", 
                                 paste(tx.var, attr(stats::terms(blip.model), "term.labels"), 
                                       sep = ":", collapse = "+")) |> as.formula()
        
        tf_vars <- paste(attr(stats::terms(tf.model), "term.labels"), collapse = "+")
        blip_vars <- paste(attr(stats::terms(self$blip.model), "term.labels"), 
                           collapse = "+")
        
        self$full.model <- paste("~", tf_vars, "+", blip_vars) |> as.formula()
    },
    
    blip_params = function(coefs) {
      possible_names <- private$tx.var
      nms <- names(coefs) |> strsplit(":")
      nms_list <- lapply(nms,
                         function(x) {
                           any(x %in% possible_names)
                         }) |> unlist()
      blip_vars <- which(nms_list)
      names(blip_vars) <- names(coefs)[blip_vars]
      
      blip_vars
    },
    
    ipw = function(A, A.hat, ...) {
      1.0 / {A * A.hat + {1.0 - A} * {1.0 - A.hat}}
    }, 
    
    ipw.cap = function(A, tx.mod.fitted, ...) {
      weights <- self$ipw(A, tx.mod.fitted)
      cap <- quantile(weights, 0.99)
      pmin(weights, cap)
    },
    
    .pom = function(...) {
      stop(".pom() is not appropriate for binary treatments", call. = FALSE)
    },
    
    Hd = function(data, A, ...) {
      data[, private$tx.var] <- A
      model.matrix(self$full.model, data)
    },
    
    Hpsi = function(data, A, ...) {
      data[, private$tx.var] <- A
      model.matrix(self$blip.model, data)
    },
    
    Hw = function(data, A, A.hat, wgt, ...) {
        data[, private$tx.var] <- {A - A.hat} * wgt
        model.matrix(self$full.model, data)
      },
    
    # If DTR, add regret function; otherwise, subtract stage blip at observed
    # treatment
    shiftY = function(type, outcome.fit, data, opt, A, ...) {
        if (type == "DTR") {
          self$regret(outcome.fit, data, opt, A)
        } else {
          # setting optimal treatment to 0 is equivalent to stage blip at
          # observed treatment
          self$regret(outcome.fit, data, 0L, A)
        }
      },
    
    regret = function(outcome.fit, data, opt, A) {
      # predicted outcome at optimal treatment
      data[, private$tx.var] <- opt
      at_opt <- predict(outcome.fit, data)
      
      # predicted outcome at observed treatment
      data[, private$tx.var] <- A
      at_A <- predict(outcome.fit, data)
      
      at_opt - at_A
    },
    
    # return optimal treatment
    opt = function(outcome.fit, data, ...) {
        data[, private$tx.var] <- 1L
        y_1 <- predict(outcome.fit, data)
        data[, private$tx.var] <- 0L
        y_0 <- predict(outcome.fit, data)
        as.integer({y_1 - y_0} > 1e-8) |> drop()
      }
  ),
  private = list(
    tx.var = NULL
  )
)

#' Specialized function for Multinomial treatments
#' 
#' @noRd
#' @keywords internal
MultiNom <- R6::R6Class(
  "MultiNom",
  public = list(
    blip.model = NULL,
    # the full outcome model
    full.model = NULL,

    initialize = function(tf.model, blip.model, tx.var, tx.levels, ...) {
      private$tx.var <- tx.var
      private$tx.levels <- tx.levels
      
      # add the treatment variable dependence to the blip formula
      self$blip.model <- paste("~ - 1 +", tx.var, "+", 
                               paste(tx.var, attr(stats::terms(blip.model), "term.labels"), 
                                     sep = ":", collapse = "+")) |> as.formula()
      
      tf_vars <- paste(attr(stats::terms(tf.model), "term.labels"), collapse = "+")
      blip_vars <- paste(attr(stats::terms(self$blip.model), "term.labels"), 
                         collapse = "+")
      
      self$full.model <- paste("~", tf_vars, "+", blip_vars) |> as.formula()
    },
    
    blip_params = function(coefs) {
      possible_names <- paste0(private$tx.var, private$tx.levels[-1L])
      nms <- names(coefs) |> strsplit(":")
      nms_list <- lapply(nms,
                         function(x) {
                           any(x %in% possible_names)
                         }) |> unlist()
      blip_vars <- which(nms_list)
      names(blip_vars) <- names(coefs)[blip_vars]

      blip_vars
    },
    
    ipw = function(A, A.hat, ...) {
      idx <- match(A, private$tx.levels)
      deno <- A.hat[cbind(seq_len(length(A)), idx)]
      1.0 / drop(deno)
    }, 
    
    ipw.cap = function(A, tx.mod.fitted, ...) {
      weights <- self$ipw(A, tx.mod.fitted)
      cap <- quantile(weights, 0.99)
      pmin(weights, cap)
    },
    
    .pom = function(A, tx.mod.fitted, data, m) {
      
      if (all(is.null(tx.mod.fitted)) || all(is.na(tx.mod.fitted))) {
        stop("`tx.weight` cannot be 'qpom' if treatment is not modelled",
             call. = FALSE)
      }
      X_alpha <- stats::model.matrix(tx.mod.fitted, 
                                     stats::model.frame(tx.mod.fitted$terms, data))

      a_cat <- A

      if (ncol(X_alpha) == 1L && attr(terms(tx.mod.fitted), "intercept") == 1L) {
        pom <- MASS::polr(a_cat ~ 1)
      } else{
        pom <- MASS::polr(a_cat ~ X_alpha[, -1L])
      }
      A <- unclass(A)
      
      prob <- pom$fitted.values[cbind(seq_along(A), A)]
      sum_prob <- rowSums(1.0 / pom$fitted.values)
      
      list("prob" = prob, "sum.prob" = sum_prob)
    },
    
    Hd = function(data, A, ...) {
      data[, private$tx.var] <- A
      model.matrix(self$full.model, data)
    },
    
    Hpsi = function(data, A, ...) {
      stop("Multinomial with g-estimation is not yet supported",
           call. = FALSE)
    },
    
    Hw = function(data, A, A.hat, wgt, ...) {
      stop("Multinomial with g-estimation is not yet supported",
           call. = FALSE)
    },
    
    # If DTR, add regret function; otherwise, subtract stage blip at observed
    # treatment
    shiftY = function(type, outcome.fit, data, opt, A, ...) {
      if (type == "DTR") {
        self$regret(outcome.fit, data, opt, A)
      } else {
        # setting optimal treatment to base level is equivalent to stage blip at
        # observed treatment
        self$regret(outcome.fit, data, private$tx.levels[1L], A)
      }
    },
    
    regret = function(outcome.fit, data, opt, A, ...) {
      # predict outcome at optimal level
      data[, private$tx.var] <- factor(opt, levels = private$tx.levels)
      at_opt <- predict(outcome.fit, data)
      
      # predict outcome at observed level
      data[, private$tx.var] <- factor(A, levels = private$tx.levels)
      at_A <- predict(outcome.fit, data)
      
      at_opt - at_A
    },
    
    # estimate optimal treatment
    opt = function(outcome.fit, data, ...) {
      est_Y <- matrix(0.0, nrow = nrow(data), ncol = length(private$tx.levels))
      for (i in seq_along(private$tx.levels)) {
        data[, private$tx.var] <- private$tx.levels[i]
        est_Y[, i] <- predict(outcome.fit, data)
      }
      max_tx <- apply(est_Y, 1L, which.max)
      private$tx.levels[max_tx]
    }
  ),
  private = list(
    tx.var = NULL,
    tx.levels = NULL
  )
)

#' Specialized function for Continuous treatment with linear blip model
#' 
#' @noRd
#' @keywords internal
ContLinearBlip <- R6::R6Class(
  "ContLinearBlip",
  public = list(
    blip.model = NULL,
    # the full outcome model
    full.model = NULL,
    
    initialize = function(tf.model, blip.model, tx.var, treat.range) {
      
      stop("linear blip functions are not supported", call. = FALSE)
      
      private$max.tx = max(treat.range)
      private$min.tx = min(treat.range)
      private$tx.var = tx.var
        
      # add the treatment variable dependence to the blip formula
      self$blip.model <- paste("~ - 1 +", tx.var, "+", 
                               paste(tx.var, attr(stats::terms(blip.model), "term.labels"), 
                                     sep = ":", collapse = "+")) |> as.formula()
        
      tf_vars <- paste(attr(stats::terms(tf.model), "term.labels"), collapse = "+")
      blip_vars <- paste(attr(stats::terms(self$blip.model), "term.labels"), 
                           collapse = "+")
        
      self$full.model <- paste("~", tf_vars, "+", blip_vars) |> as.formula()
    },
    
    blip_params = function(coefs) {
      possible_names <- paste0(private$tx.var, private$tx.levels[-1L])
      nms <- names(coefs) |> strsplit(":")
      nms_list <- lapply(nms,
                         function(x) {
                           any(x %in% possible_names)
                         }) |> unlist()
      blip_vars <- which(nms_list)
      names(blip_vars) <- names(coefs)[blip_vars]
      
      blip_vars
    },
    
    ipw = function(A, tx.mod.fitted, ...) {
      res_working <- (A - tx.mod.fitted$fitted.values) / tx.mod.fitted$fitted.values
      dispersion <- sum(res_working^2) / (length(A) - tx.mod.fitted$df)
      
      if (tx.mod.fitted$family$family == "gaussian") {
        1.0 / dnorm(A, mean = tx.mod.fitted$fitted.values, 
                    sd = sqrt(dispersion))
      } else if (tx.mod.fitted$family$family == "Gamma") {
        1.0 / dgamma(A, shape = 1.0 / dispersion,
                     scale = tx.mod.fitted$fitted.values * dispersion)
      } else {
        stop("unsupported family", call. = FALSE)
      }
    }, 
    
    ipw.cap = function(A, tx.mod.fitted, ...) {
      weights <- self$ipw(A, tx.mod.fitted)
      cap <- quantile(weights, 0.99)
      pmin(weights, cap)
    },
    
    Hd = function(data, A, ...) {
      data[, private$tx.var] <- A
      model.matrix(self$full.model, data)
    },
    
    Hpsi = function(data, A, ...) {
      data[, private$tx.var] <- A
      model.matrix(self$blip.model, data)
    },
    
    Hw = function(data, A, A.hat, wgt, ...) {
        data[, private$tx.var] <- {A - A.hat} * wgt
        model.matrix(self$full.model, data)
      },
    
    # If DTR, add regret function; otherwise, subtract stage blip at observed
    # treatment
    shiftY = function(type, outcome.fit, data, opt, A, ...) {
        if (type == "DTR") {
          self$regret(outcome.fit, data, opt, A)
        } else {
          self$regret(outcome.fit, data, 0.0, A)
        }
      },
    
    regret = function(outcome.fit, data, opt, A) {
        # predict treatment at optimal treatment
        data[, private$tx.var] <- opt
        at_opt <- predict(outcome.fit, data)
        
        # predict treatment at observed treatment
        data[, private$tx.var] <- A
        at_A <- predict(outcome.fit, data)

        at_opt - at_A
    },
    
    # estimate optimal treatment
    opt = function(outcome.fit, data, ...) {
        # use treatment values to extract only the blip component
        data[, private$tx.var] <- 1.0
        y_full <- predict(outcome.fit, data)
        data[, private$tx.var] <- 0.0
        y_me <- predict(outcome.fit, data)
        blip <- y_full - y_me
        
        tst <- as.numeric(blip > 0) |> drop() |> as.numeric()
        res <- tst * private$max.tx + (1.0 - tst) * private$min.tx
        res
      }
    
  ),
  private = list(
    min.tx = NULL,
    max.tx = NULL,
    tx.var = NULL
  )
) 

#' Specialized function for Continuous treatment with quadratic blip model
#' 
#' @noRd
#' @keywords internal
ContQuadraticBlip <- R6::R6Class(
  "ContQuadraticBlip",
  public = list(
    blip.model = NULL,
    # the full outcome model
    full.model = NULL,
    
    initialize = function(tf.model, blip.model, tx.var, treat.range) {
        private$max.tx <- max(treat.range)
        private$min.tx <- min(treat.range)
        private$tx.var <- tx.var
        private$tx.var.b <- "l__txvar2__l"
        
        # need to identify quad and linear
        # this should never be an intercept only model
        orig_vars <- attr(stats::terms(blip.model), "term.labels")
        
        if (length(orig_vars) == 0L) {
          stop("cannot provide an intercept only model for quadratic blip",
               call. = FALSE)
        }
        
        mod_vars <- sapply(orig_vars,
                           function(x) {
                             elems <- strsplit(x, ":")[[1L]]
                             if (any(elems == tx.var)) {
                               # any terms that includes tx.var is a quadratic
                               # convert I(a^2) to a new variable to facilitate
                               # later methods
                               elems[elems == tx.var] <- "l__txvar2__l"
                             } else {
                               # add treatment as interaction term
                               elems <- c(tx.var, elems)
                             }
                             paste(elems, collapse = ":")
                           })

        self$blip.model <- paste("~ -1 +", tx.var, "+",
                                 paste(mod_vars, collapse = "+")) |> 
          as.formula()
        
        tf_vars <- paste(attr(stats::terms(tf.model), "term.labels"), collapse = "+")
        blip_vars <- paste(attr(stats::terms(self$blip.model), "term.labels"), 
                           collapse = "+")
        self$full.model <- paste("~", tf_vars, "+", blip_vars) |> as.formula()
    },
    
    blip_params = function(coefs) {
      possible_names <- paste0(private$tx.var, private$tx.levels[-1L])
      nms <- names(coefs) |> strsplit(":")
      nms_list <- lapply(nms,
                         function(x) {
                           any(x %in% possible_names)
                         }) |> unlist()
      blip_vars <- which(nms_list)
      names(blip_vars) <- names(coefs)[blip_vars]
      
      blip_vars
    },

    ipw = function(A, tx.mod.fitted, ...) {
      dispersion <- summary(tx.mod.fitted)$dispersion

      if (tx.mod.fitted$family$family == "gaussian") {
        1.0 / dnorm(A, mean = tx.mod.fitted$fitted.values, 
                    sd = sqrt(dispersion))
      } else if (tx.mod.fitted$family$family == "Gamma") {
        1.0 / dgamma(A, shape = 1.0 / dispersion,
                     scale = tx.mod.fitted$fitted.values * dispersion)
      } else {
        stop("unsupported family", call. = FALSE)
      }
    }, 
    
    ipw.cap = function(A, tx.mod.fitted, ...) {
      weights <- self$ipw(A, tx.mod.fitted)
      cap <- quantile(weights, 0.99)
      pmin(weights, cap)
    },
    
    
    .pom = function(A, tx.mod.fitted, data, m) {
      if (all(is.null(tx.mod.fitted)) || all(is.na(tx.mod.fitted))) {
        stop("cannot use `weight = 'qpom' or 'wo' if treatment not modeled",
             call. = FALSE)
      }
      
      if (!inherits(tx.mod.fitted, "glm")) {
        stop("`tx.weight` cannot be 'qpom' if treatment is not modelled",
             call. = FALSE)
      }
      
      X_alpha <- stats::model.matrix(tx.mod.fitted, 
                                     stats::model.frame(tx.mod.fitted, data))
      
      a_binned <- dplyr::ntile(A, m)
      a_cat <- as.factor(a_binned)
      
      if (ncol(X_alpha) == 1L && attr(terms(tx.mod.fitted), "intercept") == 1L) {
        pom <- MASS::polr(a_cat ~ 1)
      } else{
        pom <- MASS::polr(a_cat ~ X_alpha[, -1L])
      }
      
      prob <- pom$fitted.values[cbind(seq_along(a_binned), a_binned)]
      sum_prob <- rowSums(1.0 / pom$fitted.values)
      
      list("prob" = prob, "sum.prob" = sum_prob)
    },
    
    Hd = function(data, A, ...) {
      data[, private$tx.var] <- A
      data[, private$tx.var.b] <- A^2
      model.matrix(self$full.model, data)
    },
    
    Hw = function(data, A, A.hat, wgt, treat.mod.fitted, ...) {
      A2hat <- summary(treat.mod.fitted)$dispersion + A.hat * A.hat

      data[, private$tx.var] <- {A - A.hat} * wgt
      data[, private$tx.var.b] <- {A * A - A2hat} * wgt

      model.matrix(self$full.model, data)
      
     },
    
    shiftY = function(type, outcome.fit, data, opt, A, ...) {
        if (type == "DTR") {
          self$regret(outcome.fit, data, opt, A)
        } else {
          self$regret(outcome.fit, data, 0.0, A)
        }
      },
    
    regret = function(outcome.fit, data, opt, A) {
        data[, private$tx.var] <- opt
        data[, private$tx.var.b] <- opt^2
        at_opt <- predict(outcome.fit, data)

        data[, private$tx.var] <- A
        data[, private$tx.var.b] <- A^2
        at_A <- predict(outcome.fit, data)
      
        at_opt - at_A
      },
    
    opt = function(outcome.fit, data, quiet = FALSE, ...) {
        # both treatment variables to 0 to get the non-blip term
        data[, private$tx.var] <- 0.0
        data[, private$tx.var.b] <- 0.0
        y_no_tx <- predict(outcome.fit, data)
        
        # only linear terms of blip model
        data[, private$tx.var] <- 1.0
        data[, private$tx.var.b] <- 0.0
        y_linear <- predict(outcome.fit, data) - y_no_tx

        # only quadratic terms of blip model
        data[, private$tx.var] <- 0.0
        data[, private$tx.var.b] <- 1.0
        y_quad <- predict(outcome.fit, data) - y_no_tx

        if (!is.infinite(private$min.tx) && !is.infinite(private$max.tx)) {
          data[, private$tx.var] <- private$min.tx
          data[, private$tx.var.b] <- private$min.tx^2
          predict_min <- predict(outcome.fit, data)
          
          data[, private$tx.var] <- private$max.tx
          data[, private$tx.var.b] <- private$max.tx^2
          predict_max <- predict(outcome.fit, data)
          
          opt <- -0.5 * y_linear / y_quad
          opt[y_quad >= 0.0 & predict_min >= predict_max] <- private$min.tx
          opt[y_quad >= 0.0 & predict_min < predict_max] <- private$max.tx
        } else {
          opt <- -0.5 * y_linear / y_quad
          if (any(y_quad > 0.0) && !quiet) {
            warning(" optimal treatment may be a minimum for some observations.\n",
                    call. = FALSE)
          }
          
        }
        pmin(pmax(opt, private$min.tx), private$max.tx)
      }
    
  ),
  private = list(
    max.tx = NULL,
    min.tx = NULL,
    tx.var = NULL,
    tx.var.b = NULL
  )
)