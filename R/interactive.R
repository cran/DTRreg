#' Test to see if user requested to exit interactive mode
#'
#' @noRd
#' @param x A character object.
#' @returns The unmodified input
#' @keywords internal
.test_quit <- function(x) {
  if (tolower(x) == "quit") stop("interactive mode cancelled", call. = FALSE)
  x
}

.getVariables <- function(x, data) {
  if (is.null(data)) {
    x <- tryCatch(get(x),
                  error = function(e) {
                    stop("unable to retreive ", x, " from environment; ",
                         "`data` must be provided as input\n\t",
                         e$message, call. = FALSE)
                  })
  } else {
    x <- tryCatch(data[, x],
                  error = function(e) {
                    stop("unable to retreive ", x, " from `data`\n\t",
                         e$message, call. = FALSE)
                  })
  }
  x
}

#' Function obtains required inputs from user through prompted interactions
#'
#' @noRd
#' @param isSurvival A logical. If TRUE messaging changes a bit.
#' @param prompt.cens A logical. If TRUE prompt for missingness/censoring model.
#'
#' @return A list containing the following input information
#' \itemize{
#' \item{models }{A list object. The ith element contains the models for the ith
#'   decision point. Models are "blip", "tf", and "treat"}
#' \item{method }{A character object. The DTR method to be used, one of 
#'   "dwols", "gest", or "qlearn"}
#' \item{var.estim }{A character object. Covariance matrix estimation method, 
#'   one of "none", "bootstrap", or "sandwich".}
#' \item{outcome }{A character object. The outcome of interest.}
#' }
#' @keywords internal
.interactive <- function(isSurvival = FALSE, prompt.cens = TRUE, data = NULL) {
  
  cat("DTR estimation interactive mode!\n\nEnter 'quit' at anytime to cancel.\n\n")
  cat("Note: Do not use quotation marks.\n")
  
  if (isSurvival) {
    has_censored_data <- readline("Is the data censored? (y/n) ") |> .test_quit()
    prompt.cens <- tolower(has_censored_data) == "y"
    if (prompt.cens) {
      status_var <- readline("What is the status variable? ")
    } else {
      status_var <- NA
    }
  } else {
    # prompt.cens is set by calling function
    # if `missing` == "ipw", prompt.cens = TRUE and user will be
    # prompted for a model
    outcome <- readline("Enter outcome variable: ") |> .test_quit()
    outcome <- .getVariables(outcome, data)
    status_var <- NA
    
    if (prompt.cens) {
      # if `missing` == "ipw", allow user to take default models
      message("You have set `missing` = 'ipw'.")
      use_default_missingness <- readline("Use default missingness models? (y/n) ")
      prompt.cens <- tolower(use_default_missingness) == "n"
    }
  }
  
  while (TRUE) {
    K <- readline("Enter number of stages: ") |> .test_quit() |> as.integer()
    if (K > 0L) break
    cat("K must be a positive integer\n")
  }
  
  model_type <- c("blip model", "treatment-free model", "treatment model",
                  ifelse(isSurvival, "censoring model", "missingness model"))
  models <- vector(mode = "list", length = K)
  model_names <- c("blip", "tf", "treat", "cens")

  time <- list()
  tx.var <- character(K)
  
  for (j in 1:K) {

    models[[j]] <- list()

    if (isSurvival) {
      out_time <- readline(paste("Enter stage", j, "time variable: ")) |> 
        .test_quit()
      time[[j]] <- paste("~", out_time) |> as.formula()
    }
    
    out_treat <- readline(paste("Enter stage", j, "treatment variable: ")) |> 
      .test_quit()
    tx.var[j] <- out_treat
    
    for (m in seq_along(model_names)) {

      if (model_names[m] == "treat") {
        # only the treatment model requires a dependent variable
        out <- tx.var[j]
      } else if (model_names[m] == "cens") {
        # if survival, `cens.mod` model has dependent variable only if 
        # censored is to be modeled. if not survival, and missingness is
        # to be modeled, there is no dependent variable
        if (!prompt.cens) next
        out <- ifelse(isSurvival, status_var, "")
      } else {
        out <- ""
      }
      
      while (TRUE) {
        new.var <- readline(prompt = paste("Enter all stage", j, model_type[m],
                                           "covariates separated by a comma;",
                                           "press RETURN when finished: "))
        .test_quit(new.var)
        new.var <- strsplit(new.var, ",") |> unlist() |> trimws()
        new.var <- new.var[nchar(new.var) > 0L]
        if (length(new.var) == 0L) new.var <- "1"
        form <-  paste(out, " ~ ", paste(new.var, collapse = " + "))
        cat(model_type[m], ": ", form, "\n")
        move_on <- readline(prompt = "Is this correct (y/n)? ")
        if (tolower(move_on) == "y") break
      }
      
      models[[j]][[model_names[m]]] <- form |> trimws() |> as.formula()
    }
  }
  
  while (TRUE) {
    method <- readline(paste("Dynamic WOLS (dwols), G-estimation (gest), or Q-learning (qlearn)? ")) |> 
      .test_quit()
    if (method %in% c("dwols", "gest", "qlearn")) break
    cat("method must be 'dwols', 'gest', or 'qlearn' (do not use quotation marks)\n")
  }
  
  while (TRUE) {
    var_estim <- readline("Variance estimation: none, bootstrap or sandwich? ") |> 
      .test_quit()
    if (var_estim %in% c("none", "bootstrap", "sandwich")) break
    cat("variance must be 'none', 'bootstrap', or 'sandwich' (do not use quotation marks)\n")
  }
  
  res <- list("models" = models, 
              "method" = method, 
              "var.estim" = var_estim,
              "K" = K)
  
  if (isSurvival) {
    res$time <- time
  } else {
    res$outcome <- outcome
  }
  res
}