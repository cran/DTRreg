data(twoStageCont)
n <- nrow(twoStageCont)

# Minimal valid obj for non-survival, two-stage, binary treatment, DWOLS
make_obj <- function(method = "dwols",
                     var.estim = "none",
                     tx.weight = "ipw",
                     type = "DTR") {
  list(
    K = 2L,
    outcome = twoStageCont$Y,
    data = twoStageCont,
    models = list(
      list(blip = ~X1, treat = A1 ~ X1, tf = ~X1, cens = ~1),
      list(blip = ~X2, treat = A2 ~ X2, tf = ~X2, cens = ~X1)
    ),
    dependent.vars = list(treat = c("A1", "A2")),
    tx.range = list(NA_real_, NA_real_),
    method = method,
    var.estim = var.estim,
    type = type,
    tx.weight = tx.weight,
    censoring.modeled = FALSE,
    tx.wgt.man = NA,
    tx.type = "bin",
    n.bins = NA_integer_,
    tx.family = NA,
    boot.controls = list(),
    full.cov = FALSE
  )
}

test_that(".getY: non-survival returns obj$outcome", {
  obj <- make_obj()
  expect_equal(.getY(obj, isSurvival = FALSE), obj$outcome)
})

test_that(".getY: survival returns numeric zero vector of length n", {
  obj <- make_obj()
  res <- .getY(obj, isSurvival = TRUE)
  expect_true(is.numeric(res))
  expect_length(res, n)
  expect_true(all(res == 0.0))
})

local({
  obj <- make_obj()
  res <- .getCompleteCaseProbability(obj, quiet = TRUE, isSurvival = FALSE)
  
  test_that(".getCompleteCaseProbability: returns list with required names", {
    expect_true(all(c("last.stage", "prob.complete.case",
                      "d.hat", "cens.mod.fitted") %in% names(res)))
  })
  
  test_that(".getCompleteCaseProbability: prob.complete.case is n x K matrix", {
    expect_true(is.matrix(res$prob.complete.case))
    expect_equal(dim(res$prob.complete.case), c(n, 2L))
  })
  
  test_that(".getCompleteCaseProbability: non-survival applies cumprod across stages", {
    # Oracle: for complete data, cumprod of all-1 weights stays 1
    # For stage 1: prob[,1] should equal the raw weight (no previous stages)
    # For stage 2: prob[,2] should be product of stage 1 and stage 2 weights
    raw <- .completeCaseProbability(obj, quiet = TRUE, isSurvival = FALSE)
    expected <- apply(raw$prob.complete.case, 1L, cumprod) |> t() |> matrix(ncol = 2L)
    expect_equal(res$prob.complete.case, expected)
  })
  
  test_that(".getCompleteCaseProbability: last.stage is integer vector of length n", {
    expect_true(is.integer(res$last.stage))
    expect_length(res$last.stage, n)
  })
})

local({
  obj <- make_obj()
  # Subset to stage 1 as .getStage would
  cc <- .getCompleteCaseProbability(obj, quiet = TRUE, isSurvival = FALSE)
  sc <- cc$last.stage >= 1L
  obj$data <- obj$data[sc, , drop = FALSE]
  obj$dependent.vars$treat <- obj$dependent.vars$treat[1L]
  obj$models <- obj$models[[1L]]
  obj$tx.range <- obj$tx.range[[1L]]
  
  res <- .getTxInfo(obj)
  
  test_that(".getTxInfo: returns list with required element names", {
    expect_named(res, c("A", "A.hat", "cts", "cts.obj", "tx.mod.fitted"),
                 ignore.order = TRUE)
  })
  
  test_that(".getTxInfo: cts is 'bin' for binary treatment", {
    expect_equal(res$cts, "bin")
  })
  
  test_that(".getTxInfo: cts.obj is Binary R6 object", {
    expect_true(inherits(res$cts.obj, "Binary"))
  })
  
  test_that(".getTxInfo: A is 0/1 integer vector", {
    expect_true(all(res$A %in% c(0L, 1L)))
  })
})

local({
  n_s <- 50L
  obj_ns <- list(
    data = data.frame(X = rnorm(n_s)),
    dependent.vars = list(status = NA_character_, time = NA_character_)
  )
  step_obj <- list(n = n_s)
  Y <- rnorm(n_s)
  
  test_that(".addYdelta: non-survival sets Y to provided Y", {
    res <- .addYdelta(obj_ns, step_obj, Y, isSurvival = FALSE)
    expect_equal(res$Y, Y)
  })
  
  test_that(".addYdelta: non-survival sets delta to all 1s", {
    res <- .addYdelta(obj_ns, step_obj, Y, isSurvival = FALSE)
    expect_equal(res$delta, rep(1.0, n_s))
  })
  
  test_that(".addYdelta: survival with NA status sets delta to all 1s", {
    obj_surv <- obj_ns
    obj_surv$data$T1 <- abs(rnorm(n_s)) + 0.1
    obj_surv$dependent.vars$time <- "T1"
    obj_surv$dependent.vars$status <- NA_character_
    Y_surv <- abs(rnorm(n_s))
    res <- .addYdelta(obj_surv, step_obj, Y_surv, isSurvival = TRUE)
    expect_equal(res$delta, rep(1.0, n_s))
  })
  
  test_that(".addYdelta: survival with status column sets delta from data", {
    obj_surv <- obj_ns
    delta_true <- rbinom(n_s, 1L, 0.7)
    obj_surv$data$T1 <- abs(rnorm(n_s)) + 0.1
    obj_surv$data$status <- delta_true
    obj_surv$dependent.vars$time <- "T1"
    obj_surv$dependent.vars$status <- "status"
    Y_surv <- abs(rnorm(n_s))
    res <- .addYdelta(obj_surv, step_obj, Y_surv, isSurvival = TRUE)
    expect_equal(res$delta, delta_true)
  })
  
  test_that(".addYdelta: survival Y is log(time + Y)", {
    obj_surv <- obj_ns
    T1 <- abs(rnorm(n_s)) + 0.1
    obj_surv$data$T1 <- T1
    obj_surv$dependent.vars$time <- "T1"
    obj_surv$dependent.vars$status <- NA_character_
    Y_surv <- abs(rnorm(n_s))
    res <- .addYdelta(obj_surv, step_obj, Y_surv, isSurvival = TRUE)
    expect_equal(res$Y, log(T1 + Y_surv))
  })
})

local({
  n_s <- 30L
  set.seed(42L)
  dat_s <- data.frame(X = rnorm(n_s), A = rbinom(n_s, 1L, 0.5))
  cts_b <- Binary$new(tf.model = ~X, blip.model = ~X, tx.var = "A")
  fit_tx <- glm(A ~ X, dat_s, family = "binomial")
  Ahat <- predict(fit_tx, type = "response")
  
  step_s <- list(
    n = n_s,
    A = as.numeric(dat_s$A),
    A.hat = Ahat,
    tx.mod.fitted = fit_tx,
    cts.obj = cts_b
  )
  
  make_obj_s <- function(tx.weight, tx.wgt.man = NA, method = "dwols") {
    list(tx.weight = tx.weight, method = method,
         tx.wgt.man = tx.wgt.man, n.bins = NA_integer_,
         data = dat_s)
  }
  
  test_that(".addTxWgt: 'none' sets tx.wgt to all 1s", {
    res <- .addTxWgt(make_obj_s("none"), step_s)
    expect_equal(res$tx.wgt, rep(1.0, n_s))
  })
  
  test_that(".addTxWgt: 'gest' method sets tx.wgt to all 1s", {
    res <- .addTxWgt(make_obj_s("ipw", method = "gest"), step_s)
    expect_equal(res$tx.wgt, rep(1.0, n_s))
  })
  
  test_that(".addTxWgt: 'manual' sets tx.wgt to tx.wgt.man", {
    wgt <- runif(n_s, 0.5, 2.0)
    res <- .addTxWgt(make_obj_s("manual", tx.wgt.man = wgt), step_s)
    expect_equal(res$tx.wgt, wgt)
  })
  
  test_that(".addTxWgt: 'manual.with.censor' sets tx.cens.wgt and nulls cens.wgt", {
    wgt <- runif(n_s, 0.5, 2.0)
    step_mc <- step_s
    step_mc$cens.wgt <- rep(1.0, n_s)
    res <- .addTxWgt(make_obj_s("manual.with.censor", tx.wgt.man = wgt), step_mc)
    expect_equal(res$tx.cens.wgt, wgt)
    expect_null(res$cens.wgt)
  })
  
  test_that(".addTxWgt: computed path sets tx.wgt of length n", {
    res <- .addTxWgt(make_obj_s("ipw"), step_s)
    expect_length(res$tx.wgt, n_s)
    expect_true(is.numeric(res$tx.wgt))
  })
})

local({
  n_s <- 20L
  cens <- rep(1.0, n_s)
  delta <- rep(1.0, n_s)
  cc <- list(prob.complete.case = cens)
  
  step_base <- list(n = n_s, delta = delta, A = rep(1.0, n_s))
  
  test_that(".addWgts: 'none' — wgts equals cens.wgt * delta", {
    obj_s <- list(tx.weight = "none", method = "dwols",
                  tx.wgt.man = NA, n.bins = NA_integer_,
                  data = data.frame(X = rnorm(n_s)))
    res <- .addWgts(obj_s, step_base, cc)
    expect_equal(res$wgts, cens * delta)
  })
  
  test_that(".addWgts: 'manual.with.censor' — wgts equals tx.cens.wgt * delta", {
    wgt <- runif(n_s, 0.5, 2.0)
    obj_s <- list(tx.weight = "manual.with.censor", method = "dwols",
                  tx.wgt.man = wgt, n.bins = NA_integer_,
                  data = data.frame(X = rnorm(n_s)))
    step_mc <- step_base
    step_mc$cens.wgt <- cens
    res <- .addWgts(obj_s, step_mc, cc)
    expect_equal(res$wgts, wgt * delta)
  })
})

local({
  obj <- make_obj(method = "dwols")
  cc <- .getCompleteCaseProbability(obj, quiet = TRUE, isSurvival = FALSE)
  # Get a realistic step_obj via .getStage
  step <- .getStage(obj, k = 2L, quiet = TRUE, complete.case.info = cc,
                    isSurvival = FALSE, Y = obj$outcome)
  
  # Build stage-specific obj as .getStage would
  sc <- cc$last.stage >= 2L
  obj_s <- obj
  obj_s$data <- obj$data[sc, , drop = FALSE]
  obj_s$dependent.vars$treat <- obj$dependent.vars$treat[2L]
  obj_s$models <- obj$models[[2L]]
  obj_s$tx.range <- obj$tx.range[[2L]]
  
  test_that(".performMethod DWOLS: returns list with outcome.fit, psi, beta", {
    res <- .performMethod(obj_s, step)
    expect_named(res, c("outcome.fit", "beta", "psi"), ignore.order = TRUE)
  })
  
  test_that(".performMethod DWOLS: psi names contain tx.var", {
    res <- .performMethod(obj_s, step)
    expect_true(all(grepl("A2", names(res$psi))))
  })
  
  test_that(".performMethod DWOLS: beta names do not contain tx.var", {
    res <- .performMethod(obj_s, step)
    expect_false(any(grepl("^A2", names(res$beta))))
  })
  
  test_that(".performMethod DWOLS: psi and beta account for all coefficients", {
    res <- .performMethod(obj_s, step)
    expect_setequal(c(names(res$psi), names(res$beta)),
                    names(coef(res$outcome.fit)))
  })
})

local({
  obj <- make_obj(method = "gest")
  cc <- .getCompleteCaseProbability(obj, quiet = TRUE, isSurvival = FALSE)
  step <- .getStage(obj, k = 2L, quiet = TRUE, complete.case.info = cc,
                    isSurvival = FALSE, Y = obj$outcome)
  
  sc <- cc$last.stage >= 2L
  obj_s <- obj
  obj_s$data <- obj$data[sc, , drop = FALSE]
  obj_s$dependent.vars$treat <- obj$dependent.vars$treat[2L]
  obj_s$models <- obj$models[[2L]]
  obj_s$tx.range <- obj$tx.range[[2L]]
  
  test_that(".performMethod GEST: returns list with outcome.fit, psi, beta", {
    res <- .performMethod(obj_s, step)
    expect_named(res, c("outcome.fit", "beta", "psi"), ignore.order = TRUE)
  })
  
  test_that(".performMethod GEST: outcome.fit has class GEST", {
    res <- .performMethod(obj_s, step)
    expect_true(inherits(res$outcome.fit, "GEST"))
  })
})


local({
  obj <- make_obj()
  cc <- .getCompleteCaseProbability(obj, quiet = TRUE, isSurvival = FALSE)
  
  res_k2 <- .getStage(obj, k = 2L, quiet = TRUE, complete.case.info = cc,
                      isSurvival = FALSE, Y = obj$outcome)
  res_k1 <- .getStage(obj, k = 1L, quiet = TRUE, complete.case.info = cc,
                      isSurvival = FALSE, Y = obj$outcome)
  
  test_that(".getStage: returns list with core structural elements", {
    expect_true(all(c("dp", "tx.var", "models", "n",
                      "A", "A.hat", "cts", "cts.obj",
                      "Y", "delta", "wgts",
                      "outcome.fit", "psi", "beta") %in% names(res_k2)))
  })
  
  test_that(".getStage: dp matches k", {
    expect_equal(res_k2$dp, 2L)
    expect_equal(res_k1$dp, 1L)
  })
  
  test_that(".getStage: tx.var matches dependent.vars$treat[k]", {
    expect_equal(res_k2$tx.var, obj$dependent.vars$treat[2L])
    expect_equal(res_k1$tx.var, obj$dependent.vars$treat[1L])
  })
  
  test_that(".getStage: n equals sum of stage cases", {
    expect_equal(res_k2$n, sum(cc$last.stage >= 2L))
    expect_equal(res_k1$n, sum(cc$last.stage >= 1L))
  })
  
  test_that(".getStage: A is 0/1 integer vector of length n", {
    expect_true(all(res_k2$A %in% c(0L, 1L)))
    expect_length(res_k2$A, res_k2$n)
  })
  
  test_that(".getStage: psi names contain tx.var", {
    expect_true(all(grepl(res_k2$tx.var, names(res_k2$psi))))
  })
  
  test_that(".getStage: beta names do not contain tx.var", {
    expect_false(any(grepl(paste0("^", res_k2$tx.var), names(res_k2$beta))))
  })
  
  test_that(".getStage: delta is all 1s for non-survival", {
    expect_equal(res_k2$delta, rep(1.0, res_k2$n))
  })
  
  test_that(".getStage: wgts is numeric vector of length n", {
    expect_true(is.numeric(res_k2$wgts))
    expect_length(res_k2$wgts, res_k2$n)
  })
})

test_that(".dtrProcedure: stops when required element missing from obj", {
  obj <- make_obj()
  obj$K <- NULL

  expect_error(.dtrProcedure(obj, quiet = TRUE), "`obj`")
  expect_error(.dtrProcedure(make_obj(), quiet = "yes"), "`quiet`")
  expect_error(.dtrProcedure(make_obj(), quiet = TRUE, isSurvival = 0L),
               "`isSurvival`")
  obj <- make_obj(method = "gest")
  obj$tx.type <- "multi"
  expect_error(.dtrProcedure(obj, quiet = TRUE), "multinomial")
})


local({
  obj <- make_obj(method = "dwols")
  res <- suppressMessages(.dtrProcedure(obj, quiet = TRUE))
  
  test_that(".dtrProcedure DWOLS: returns a list", {
    expect_true(is.list(res))
  })
  
  test_that(".dtrProcedure DWOLS: stages element contains K results", {
    expect_length(res$stages, obj$K)
  })
  
  test_that(".dtrProcedure DWOLS: opt.Y is numeric vector of length n", {
    expect_true(is.numeric(res$opt.Y))
    expect_length(res$opt.Y, n)
  })
  
  test_that(".dtrProcedure DWOLS: last.stage is integer vector of length n", {
    expect_true(is.integer(res$last.stage))
    expect_length(res$last.stage, n)
  })
  
  test_that(".dtrProcedure DWOLS: treat.vars matches dependent.vars$treat", {
    expect_equal(res$treat.vars, obj$dependent.vars$treat)
  })
  
  test_that(".dtrProcedure DWOLS: regret is list of length K", {
    expect_length(res$regret, obj$K)
  })
})

local({
  obj <- make_obj(method = "gest")
  res <- suppressMessages(.dtrProcedure(obj, quiet = TRUE))
  
  test_that(".dtrProcedure GEST: opt.Y is numeric vector of length n", {
    expect_true(is.numeric(res$opt.Y))
    expect_length(res$opt.Y, n)
  })
})