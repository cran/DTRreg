set.seed(42L)
n <- 20L

# Complete data - two stage, no missingness
dat_complete <- data.frame(
  X1 = rnorm(n),
  A1 = rbinom(n, 1L, 0.5),
  X2 = rnorm(n),
  A2 = rbinom(n, 1L, 0.5),
  Y = rnorm(n)
)

# Data with one NA in stage 2 covariate
dat_na <- dat_complete
dat_na$X2[3L] <- NA_real_

models_2stage <- list(
  list(blip = ~X1, treat = A1 ~ X1, tf = ~X1, cens = ~1),
  list(blip = ~X2, treat = A2 ~ X2, tf = ~X2, cens = ~X1)
)

# Minimal obj for non-survival
make_obj <- function(data = dat_complete,
                     censoring.modeled = FALSE,
                     isSurvival = FALSE) {
  obj <- list(
    K = 2L,
    data = data,
    models = models_2stage,
    dependent.vars = list(treat = c("A1", "A2")),
    censoring.modeled = censoring.modeled
  )
  if (isSurvival) obj$manual.censor.weight <- FALSE
  obj
}


test_that(".initStillIn: all TRUE when censoring.modeled=TRUE, non-survival", {
  obj <- make_obj(censoring.modeled = TRUE)
  res <- .initStillIn(obj, isSurvival = FALSE)
  expect_true(all(res))
  expect_length(res, n)
})

test_that(".initStillIn: all TRUE when censoring.modeled=TRUE, survival", {
  obj <- make_obj(censoring.modeled = TRUE, isSurvival = TRUE)
  res <- .initStillIn(obj, isSurvival = TRUE)
  expect_true(all(res))
})

test_that(".initStillIn: uses complete.cases when censoring.modeled=FALSE, non-survival", {
  obj <- make_obj(data = dat_na, censoring.modeled = FALSE)
  res <- .initStillIn(obj, isSurvival = FALSE)
  expect_equal(res, complete.cases(dat_na))
  expect_false(res[3L])
})

test_that(".initStillIn: uses complete.cases when censoring.modeled=FALSE, survival, no manual weight", {
  obj <- make_obj(data = dat_na, censoring.modeled = FALSE, isSurvival = TRUE)
  obj$manual.censor.weight <- FALSE
  res <- .initStillIn(obj, isSurvival = TRUE)
  expect_equal(res, complete.cases(dat_na))
})

test_that(".initStillIn: all TRUE when censoring.modeled=FALSE, survival, manual weight provided", {
  obj <- make_obj(data = dat_na, censoring.modeled = FALSE, isSurvival = TRUE)
  obj$manual.censor.weight <- TRUE
  res <- .initStillIn(obj, isSurvival = TRUE)
  expect_true(all(res))
})


test_that(".dependentVarsCompleteCases: returns logical vector of length n", {
  obj <- make_obj()
  res <- .dependentVarsCompleteCases(obj, K = 1L)
  expect_true(is.logical(res))
  expect_length(res, n)
})

test_that(".dependentVarsCompleteCases: all TRUE when no missing treatment data", {
  obj <- make_obj()
  res <- .dependentVarsCompleteCases(obj, K = 1L)
  expect_true(all(res))
})

test_that(".dependentVarsCompleteCases: FALSE for rows with missing treatment", {
  dat_na_tx <- dat_complete
  dat_na_tx$A1[5L] <- NA_integer_
  obj <- make_obj(data = dat_na_tx)
  res <- .dependentVarsCompleteCases(obj, K = 1L)
  expect_false(res[5L])
  expect_true(all(res[-5L]))
})

test_that(".dependentVarsCompleteCases: stops when treatment column absent from data", {
  obj <- make_obj()
  obj$dependent.vars$treat <- c("A1", "A99")
  expect_error(.dependentVarsCompleteCases(obj, K = 2L), "unable to extract variable")
})


test_that(".modelVarsCompleteCases: returns logical vector of length n", {
  obj <- make_obj()
  res <- .modelVarsCompleteCases(obj, K = 1L)
  expect_true(is.logical(res))
  expect_length(res, n)
})

test_that(".modelVarsCompleteCases: all TRUE when no missing covariate data", {
  obj <- make_obj()
  res <- .modelVarsCompleteCases(obj, K = 1L)
  expect_true(all(res))
})

test_that(".modelVarsCompleteCases: FALSE for rows with missing model covariate", {
  obj <- make_obj(data = dat_na)
  # X2 is NA for row 3 - stage 2 model uses X2
  res <- .modelVarsCompleteCases(obj, K = 2L)
  expect_false(res[3L])
})


test_that(".completeCases: oracle - conjunction of dependent and model complete cases", {
  obj <- make_obj(data = dat_na)
  dep <- .dependentVarsCompleteCases(obj, K = 2L)
  mod <- .modelVarsCompleteCases(obj, K = 2L)
  expect_equal(.completeCases(obj, K = 2L), dep & mod)
})

test_that(".completeCases: all TRUE for complete data", {
  obj <- make_obj()
  expect_true(all(.completeCases(obj, K = 1L)))
  expect_true(all(.completeCases(obj, K = 2L)))
})


test_that(".getDelta: non-survival returns vector of length sum(still.in)", {
  still_in <- c(TRUE, TRUE, FALSE, TRUE, TRUE)
  cc <- c(TRUE, FALSE, FALSE, TRUE, TRUE)
  res <- .getDelta(still_in, cc, status = NULL)
  expect_length(res, sum(still_in))
})

test_that(".getDelta: non-survival oracle - 0 for dropouts, 1 for complete cases", {
  still_in <- c(TRUE, TRUE, TRUE, TRUE, TRUE)
  cc <- c(TRUE, FALSE, TRUE, FALSE, TRUE)
  res <- .getDelta(still_in, cc, status = NULL)
  expected <- c(1L, 0L, 1L, 0L, 1L)
  expect_equal(res, expected)
})

test_that(".getDelta: non-survival excludes participants not still_in", {
  still_in <- c(TRUE, TRUE, FALSE, TRUE, TRUE)
  cc <- c(TRUE, FALSE, FALSE, TRUE, TRUE)
  res <- .getDelta(still_in, cc, status = NULL)
  # row 3 is not still_in so excluded; row 2 is still_in but not cc so delta=0
  expected <- c(1L, 0L, 1L, 1L)
  expect_equal(res, expected)
})

test_that(".getDelta: survival returns status for still_in & cc subset", {
  still_in <- c(TRUE, TRUE, TRUE, TRUE, TRUE)
  cc <- c(TRUE, FALSE, TRUE, TRUE, FALSE)
  status <- c(1L, 0L, 1L, 0L, 1L)
  res <- .getDelta(still_in, cc, status = status)
  # still_in & cc = rows 1, 3, 4
  expect_equal(res, status[still_in & cc])
})


local({
  n_upd <- 6L
  K_upd <- 2L
  
  ml <- list(
    last.stage = integer(n_upd),
    prob.complete.case = matrix(NA_real_, nrow = n_upd, ncol = K_upd),
    d.hat = matrix(NA_real_, nrow = n_upd, ncol = K_upd),
    cens.mod.fitted = list()
  )
  
  still_in <- c(TRUE,  TRUE,  TRUE,  TRUE,  FALSE, FALSE)
  cc <- c(FALSE,  FALSE,  FALSE, FALSE, FALSE, FALSE)
  k <- 1L
  
  cen_fit <- list(
    wt.cen = c(1.0, 1.0, 1.0, 1.0),
    D.hat = c(0.9, 0.8, 0.7, 0.6),
    cens.model.fitted = NA
  ) 
  
  res <- .updateCensorMatrices(ml, still_in, cc, k, cen_fit)

  test_that(".updateCensorMatrices: last.stage updated for cc rows", {
    expect_equal(res$last.stage[cc], rep(k, sum(cc)))
    expect_equal(res$last.stage[!cc], rep(0L, sum(!cc)))
  })
  
  test_that(".updateCensorMatrices: prob.complete.case updated", {
    expect_equal(res$prob.complete.case[, k], rep(NA_real_, length(still_in)))
  })
  
  test_that(".updateCensorMatrices: prob.complete.case is NA for !cc rows", {
    expect_true(all(is.na(res$prob.complete.case[!cc, k])))
  })
  
  test_that(".updateCensorMatrices: d.hat updated for still_in rows", {
    expect_equal(res$d.hat[still_in, k], cen_fit$D.hat)
  })
  
  test_that(".updateCensorMatrices: cens.mod.fitted stored at correct stage index", {
    expect_identical(res$cens.mod.fitted[[k]], cen_fit$cens.model.fitted)
  })
})

local({
  n_upd <- 6L
  K_upd <- 2L
  
  ml <- list(
    last.stage = integer(n_upd),
    prob.complete.case = matrix(NA_real_, nrow = n_upd, ncol = K_upd),
    d.hat = matrix(NA_real_, nrow = n_upd, ncol = K_upd),
    cens.mod.fitted = list()
  )
  
  still_in <- c(TRUE,  TRUE,  TRUE,  TRUE,  FALSE, FALSE)
  cc <- c(TRUE,  TRUE,  TRUE, TRUE, FALSE, FALSE)
  k <- 1L
  
  cen_fit <- list(
    wt.cen = c(1.0, 1.0, 1.0, 1.0),
    D.hat = c(0.9, 0.8, 0.7, 0.6),
    cens.model.fitted = NA
  ) 
  
  res <- .updateCensorMatrices(ml, still_in, cc, k, cen_fit)
  
  test_that(".updateCensorMatrices: last.stage updated for cc rows", {
    expect_equal(res$last.stage[cc], rep(k, sum(cc)))
    expect_equal(res$last.stage[!cc], rep(0L, sum(!cc)))
  })
  
  test_that(".updateCensorMatrices: prob.complete.case updated for still_in rows", {
    expect_equal(res$prob.complete.case[still_in, k], cen_fit$wt.cen)
  })
  
  test_that(".updateCensorMatrices: prob.complete.case is NA for !cc rows", {
    expect_true(all(is.na(res$prob.complete.case[!cc, k])))
  })
  
  test_that(".updateCensorMatrices: d.hat updated for still_in rows", {
    expect_equal(res$d.hat[still_in, k], cen_fit$D.hat)
  })
  
  test_that(".updateCensorMatrices: cens.mod.fitted stored at correct stage index", {
    expect_identical(res$cens.mod.fitted[[k]], cen_fit$cens.model.fitted)
  })
})

# Minimal data for censoring model
dat_cens <- data.frame(X = rnorm(n), A = rbinom(n, 1L, 0.5))

test_that(".processCensoring: returns list with required element names", {
  delta <- rep(1L, n)
  res <- .processCensoring(~1, delta, dat_cens,
                           censoring.modeled = FALSE, dp = 1L, quiet = TRUE)
  expect_named(res, c("cens.model.fitted", "D.hat", "wt.cen"), ignore.order = TRUE)
})

test_that(".processCensoring: path 3 - censoring not modeled, d.hat equals delta", {
  delta <- rep(1L, n)
  res <- .processCensoring(~1, delta, dat_cens,
                           censoring.modeled = FALSE, dp = 1L, quiet = TRUE)
  expect_equal(as.numeric(res$D.hat), as.numeric(delta))
})

test_that(".processCensoring: path 2 - all delta=1, censoring modeled, messages and d.hat=delta", {
  delta <- rep(1L, n)
  expect_message(
    res <- .processCensoring(~X, delta, dat_cens,
                             censoring.modeled = TRUE, dp = 1L, quiet = FALSE),
    "censoring model will be ignored"
  )
  expect_equal(as.numeric(res$D.hat), as.numeric(delta))
  expect_true(is.na(res$cens.model.fitted[[1L]]))
})

test_that(".processCensoring: path 1 - some delta=0, censoring modeled, fits glm", {
  delta <- c(rep(1L, n - 5L), rep(0L, 5L))
  res <- .processCensoring(~X, delta, dat_cens,
                           censoring.modeled = TRUE, dp = 1L, quiet = TRUE)
  expect_true(inherits(res$cens.model.fitted, "glm"))
  expect_length(res$D.hat, n)
})

test_that(".processCensoring: wt.cen oracle - 1/(delta*d.hat + (1-delta)*(1-d.hat))", {
  delta <- rep(1L, n)
  res <- .processCensoring(~1, delta, dat_cens,
                           censoring.modeled = FALSE, dp = 1L, quiet = TRUE)
  d.hat <- as.numeric(res$D.hat)
  delta_n <- as.numeric(delta)
  expected <- 1.0 / (delta_n * d.hat + (1.0 - delta_n) * (1.0 - d.hat))
  expect_equal(res$wt.cen, expected)
})

test_that(".processCensoring: wt.cen oracle - mixed delta with modeled censoring", {
  delta <- c(rep(1L, n - 5L), rep(0L, 5L))
  res <- .processCensoring(~X, delta, dat_cens,
                           censoring.modeled = TRUE, dp = 1L, quiet = TRUE)
  d.hat <- as.numeric(res$D.hat)
  delta_n <- as.numeric(delta)
  expected <- 1.0 / (delta_n * d.hat + (1.0 - delta_n) * (1.0 - d.hat))
  expect_equal(unname(res$wt.cen), expected)
})


A_bin <- as.numeric(rbinom(n, 1L, 0.5))
Ahat_v <- runif(n, 0.3, 0.7)

cts_bin <- Binary$new(tf.model = ~X1,
                      blip.model = ~X1,
                      tx.var = "A1")

fit_tx <- glm(A1 ~ X1, dat_complete, family = "binomial")

# Convenience wrapper
tw <- function(A = A_bin,
               A.hat = Ahat_v,
               tx.weight = "ipw",
               tx.mod.fitted = fit_tx,
               cts.obj = cts_bin,
               n.bins = NA_integer_,
               data = dat_complete) {
  .treatmentWeights(A = A, A.hat = A.hat, tx.weight = tx.weight,
                    tx.mod.fitted = tx.mod.fitted, cts.obj = cts.obj,
                    n.bins = n.bins, data = data)
}

test_that(".treatmentWeights: ipw matches direct cts.obj$ipw call", {
  res <- tw(tx.weight = "ipw")
  expected <- cts_bin$ipw(A = A_bin, A.hat = Ahat_v,
                          tx.mod.fitted = fit_tx)
  expect_equal(res, expected)
})

test_that(".treatmentWeights: cipw oracle - pmin(ipw, quantile(ipw, 0.99))", {
  res <- tw(tx.weight = "cipw")
  weights <- cts_bin$ipw(A = A_bin, A.hat = Ahat_v, tx.mod.fitted = fit_tx)
  cap <- quantile(weights, 0.99)
  expected <- pmin(weights, cap)
  expect_equal(res, expected)
})

test_that(".treatmentWeights: cipw values are <= ipw values", {
  res_ipw <- tw(tx.weight = "ipw")
  res_cipw <- tw(tx.weight = "cipw")
  expect_true(all(res_cipw <= res_ipw))
})

test_that(".treatmentWeights: overlap binary oracle - abs(A - A.hat)", {
  res <- tw(tx.weight = "overlap")
  expected <- abs(A_bin - Ahat_v)
  expect_equal(res, expected)
})

test_that(".treatmentWeights: unrecognised weight returns vector of 1s", {
  res <- tw(tx.weight = "none")
  expect_equal(res, rep(1.0, n))
})

test_that(".treatmentWeights: stops when A has length 1", {
  expect_error(tw(A = 1.0, A.hat = 1.0), "`A`")
  expect_error(tw(A.hat = Ahat_v[1L:5L]), "`A.hat`")
  expect_error(tw(tx.weight = 1L), "`tx.weight`")
  expect_error(tw(cts.obj = list()), "`cts.obj`")
  expect_error(tw(n.bins = 5.0), "`n.bins`")
})

local({
  # ContQuadraticBlip with .pom() implemented
  dat_cont <- data.frame(
    X = cos(seq(0.0, pi, length.out = n)),
    A = runif(n, 0.1, 0.9)
  )
  cts_cont <- ContQuadraticBlip$new(tf.model = ~X,
                                    blip.model = ~A,
                                    tx.var = "A",
                                    treat.range = c(0.1, 0.9))
  fit_cont <- glm(A ~ X, dat_cont, family = gaussian())
  m <- 5L
  
  test_that("q.pom: oracle - (1/m) / pom$prob with named arguments to .pom()", {
    pom <- cts_cont$.pom(A = dat_cont$A,
                         tx.mod.fitted = fit_cont,
                         data = dat_cont,
                         m = m)
    expected <- (1.0 / m) / pom$prob
    result <- q.pom(dat_cont$A, fit_cont, dat_cont, m, cts_cont)
    expect_equal(result, expected)
  })
  
  test_that("q.pom: returns numeric vector of length n", {
    result <- q.pom(dat_cont$A, fit_cont, dat_cont, m, cts_cont)
    expect_true(is.numeric(result))
    expect_length(result, n)
  })
  
  test_that("bin.o.fcn: oracle - (1/pom$prob) / pom$sum.prob with named arguments", {
    pom <- cts_cont$.pom(A = dat_cont$A,
                         tx.mod.fitted = fit_cont,
                         data = dat_cont,
                         m = m)
    expected <- (1.0 / pom$prob) / pom$sum.prob
    result <- bin.o.fcn(dat_cont$A, fit_cont, dat_cont, m, cts_cont)
    expect_equal(result, expected)
  })
  
  test_that("bin.o.fcn: returns numeric vector of length n", {
    result <- bin.o.fcn(dat_cont$A, fit_cont, dat_cont, m, cts_cont)
    expect_true(is.numeric(result))
    expect_length(result, n)
  })
})


test_that(".completeCaseProbability: stops when required element missing from obj", {
  obj_bad <- make_obj()
  obj_bad$K <- NULL
  expect_error(.completeCaseProbability(obj_bad, quiet = TRUE, isSurvival = FALSE),
               "required information")
})

test_that(".completeCaseProbability: returns list with required element names", {
  obj <- make_obj()
  obj$models[[1L]]$cens <- ~1
  obj$models[[2L]]$cens <- ~X1
  res <- .completeCaseProbability(obj, quiet = TRUE, isSurvival = FALSE)
  expect_true(all(c("last.stage", "prob.complete.case",
                    "d.hat", "cens.mod.fitted") %in% names(res)))
})

test_that(".completeCaseProbability: last.stage is integer vector of length n", {
  obj <- make_obj()
  obj$models[[1L]]$cens <- ~1
  obj$models[[2L]]$cens <- ~X1
  res <- .completeCaseProbability(obj, quiet = TRUE, isSurvival = FALSE)
  expect_true(is.integer(res$last.stage))
  expect_length(res$last.stage, n)
})

test_that(".completeCaseProbability: prob.complete.case is n x K matrix", {
  obj <- make_obj()
  obj$models[[1L]]$cens <- ~1
  obj$models[[2L]]$cens <- ~X1
  res <- .completeCaseProbability(obj, quiet = TRUE, isSurvival = FALSE)
  expect_true(is.matrix(res$prob.complete.case))
  expect_equal(dim(res$prob.complete.case), c(n, 2L))
})

test_that(".completeCaseProbability: all last.stage = K for complete data", {
  obj <- make_obj()
  obj$models[[1L]]$cens <- ~1
  obj$models[[2L]]$cens <- ~X1
  res <- .completeCaseProbability(obj, quiet = TRUE, isSurvival = FALSE)
  expect_true(all(res$last.stage == obj$K))
})

test_that(".completeCaseProbability: last.stage < K for participants with missing stage 2 data", {
  obj <- make_obj(data = dat_na)
  obj$models[[1L]]$cens <- ~1
  obj$models[[2L]]$cens <- ~X1
  res <- suppressMessages(
    .completeCaseProbability(obj, quiet = TRUE, isSurvival = FALSE)
  )
  # row 3 has NA in X2 - cannot complete stage 2
  expect_lt(res$last.stage[3L], obj$K)
})

test_that(".treatmentWeights: overlap continuous oracle — matches bin.o.fcn", {
  dat_cont <- data.frame(X = cos(seq(0.0, pi, length.out = n)),
                         A = runif(n, 0.1, 0.9))
  cts_cont <- ContQuadraticBlip$new(tf.model    = ~X,
                                    blip.model  = ~A,
                                    tx.var      = "A",
                                    treat.range = c(0.1, 0.9))
  fit_cont <- glm(A ~ X, dat_cont, family = gaussian())
  m        <- 5L
  res      <- .treatmentWeights(A = dat_cont$A, A.hat = fit_cont$fitted.values,
                                tx.weight = "overlap", tx.mod.fitted = fit_cont,
                                cts.obj = cts_cont, n.bins = m, data = dat_cont)
  expected <- bin.o.fcn(dat_cont$A, fit_cont, dat_cont, m, cts_cont)
  expect_equal(res, expected)
})

test_that(".treatmentWeights: overlap multinomial oracle — matches bin.o.fcn", {
  dat_multi <- data.frame(
    X = rnorm(n),
    A = factor(sample(c("low", "med", "high"), n, replace = TRUE),
               levels = c("low", "med", "high"))
  )
  cts_multi  <- MultiNom$new(tf.model = ~X, blip.model = ~X, tx.var = "A",
                             tx.levels = c("low", "med", "high"))
  fit_multi  <- suppressMessages(nnet::multinom(A ~ X, dat_multi))
  Ahat_mat   <- predict(fit_multi, type = "probs")
  m          <- length(levels(dat_multi$A))
  res        <- .treatmentWeights(A = dat_multi$A, A.hat = Ahat_mat,
                                  tx.weight = "overlap", tx.mod.fitted = fit_multi,
                                  cts.obj = cts_multi, n.bins = m, data = dat_multi)
  expected   <- bin.o.fcn(dat_multi$A, fit_multi, dat_multi, m, cts_multi)
  expect_equal(res, expected)
})