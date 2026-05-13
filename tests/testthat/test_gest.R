# Minimal formula and coefficients matching what .gest() would produce
gest_formula <- ~X + A + A:X
gest_coefs <- c("(Intercept)" = 1.0, "X" = 0.5, "A" = 2.0, "A:X" = -0.5)

# Match the structure that .gest() assigns to outcome.fit
gest_obj <- structure(
  list(formula = gest_formula,
       coefficients = gest_coefs),
  class = c("GEST", "list")
)

dat <- data.frame(
  X = c(-1.0, 0.0, 0.5, 1.0, 2.0),
  A = c( 0.0, 1.0, 0.0, 1.0, 0.0)
)

test_that("predict.GEST: matches oracle of explicit model.matrix %*% coefficients", {
  expected <- drop(model.matrix(gest_formula, dat) %*% gest_coefs)
  expect_equal(predict.GEST(gest_obj, newdata = dat), expected)
})

test_that("predict.GEST: returns a vector, not a matrix", {
  result <- predict.GEST(gest_obj, newdata = dat)
  expect_true(is.vector(result))
  expect_false(is.matrix(result))
})

test_that("predict.GEST: result has length equal to nrow(newdata)", {
  result <- predict.GEST(gest_obj, newdata = dat)
  expect_length(result, nrow(dat))
})

test_that("predict.GEST: reflects overwritten treatment column in newdata", {
  dat_opt <- dat
  dat_opt$A <- 1.0  # set all to optimal treatment
  result_obs <- predict.GEST(gest_obj, newdata = dat)
  result_opt <- predict.GEST(gest_obj, newdata = dat_opt)
  # predictions must differ when treatment column differs
  expect_false(isTRUE(all.equal(result_obs, result_opt)))
})

test_that("predict.GEST: works correctly for a single-row newdata", {
  single <- dat[1L, , drop = FALSE]
  expected <- drop(model.matrix(gest_formula, single) %*% gest_coefs)
  result <- predict.GEST(gest_obj, newdata = single)
  expect_equal(result, expected)
  expect_length(result, 1L)
})

n <- 30L
set.seed(42L)

dat <- data.frame(
  X = cos(seq(0.0, pi, length.out = n)),
  A = as.numeric(rbinom(n, 1L, 0.5)),
  Y = rnorm(n)
)

models_bin <- list(blip = ~X, treat = A ~ X, tf = ~X)

cts_bin <- Binary$new(tf.model = models_bin$tf,
                      blip.model = models_bin$blip,
                      tx.var = "A")

fit_tx <- glm(models_bin$treat, dat, family = "binomial")
Ahat <- predict(fit_tx, type = "response")
wgts <- rep(1.0, n)

local({
  res <- .gest(A = dat$A,
               wgts = wgts,
               data = dat,
               Y = dat$Y,
               A.hat = Ahat,
               tx.var = "A",
               cts.obj = cts_bin,
               treat.mod.fitted = fit_tx)
  
  test_that(".gest: returns list with required element names", {
    expect_named(res, c("outcome.fit", "beta", "psi"), ignore.order = TRUE)
  })
  
  test_that(".gest: outcome.fit has class GEST", {
    expect_true(inherits(res$outcome.fit, "GEST"))
  })
  
  test_that(".gest: outcome.fit contains formula, fitted.values, coefficients", {
    expect_named(res$outcome.fit,
                 c("formula", "fitted.values", "coefficients"),
                 ignore.order = TRUE)
  })
  
  test_that(".gest: outcome.fit$formula matches cts.obj$full.model", {
    expect_equal(deparse(res$outcome.fit$formula),
                 deparse(cts_bin$full.model))
  })
  
  test_that(".gest: outcome.fit$fitted.values has length n", {
    expect_length(res$outcome.fit$fitted.values, n)
  })
  
  test_that(".gest: coefficients match oracle of explicit crossprod/solve", {
    Hd <- cts_bin$Hd(dat, dat$A)
    Hw <- cts_bin$Hw(data = dat, A = dat$A, A.hat = Ahat,
                     wgt = wgts, treat.mod.fitted = fit_tx)
    expected <- drop(solve(crossprod(Hw, Hd), crossprod(Hw, dat$Y)))
    names(expected) <- colnames(Hd)
    expect_equal(res$outcome.fit$coefficients, expected)
  })
  
  test_that(".gest: psi contains only blip-related terms (contains tx.var)", {
    expect_true(all(grepl("A", names(res$psi))))
  })
  
  test_that(".gest: beta contains no blip-related terms", {
    expect_false(any(grepl("^A", names(res$beta))))
  })
  
  test_that(".gest: psi and beta together account for all coefficients", {
    expect_setequal(c(names(res$psi), names(res$beta)),
                    names(res$outcome.fit$coefficients))
  })
  
  test_that(".gest: fitted.values match oracle of Hd %*% coefficients", {
    Hd <- cts_bin$Hd(dat, dat$A)
    expected <- drop(Hd %*% res$outcome.fit$coefficients)
    expect_equal(res$outcome.fit$fitted.values, expected)
  })
})

test_that(".gest: errors as expected", {
  expect_error(
    .gest(A = as.character(dat$A), wgts = wgts, data = dat, Y = dat$Y,
          A.hat = Ahat, tx.var = "A", cts.obj = cts_bin,
          treat.mod.fitted = fit_tx),
    "`A`"
  )
  expect_error(
    .gest(A = 1.0, wgts = 1.0, data = dat[1L, ], Y = dat$Y[1L],
          A.hat = Ahat[1L], tx.var = "A", cts.obj = cts_bin,
          treat.mod.fitted = fit_tx),
    "`A`"
  )
  expect_error(
    .gest(A = dat$A, wgts = wgts[1L:5L], data = dat, Y = dat$Y,
          A.hat = Ahat, tx.var = "A", cts.obj = cts_bin,
          treat.mod.fitted = fit_tx),
    "`wgts`"
  )
  expect_error(
    .gest(A = dat$A, wgts = wgts, data = dat,
          Y = as.character(dat$Y),
          A.hat = Ahat, tx.var = "A", cts.obj = cts_bin,
          treat.mod.fitted = fit_tx),
    "`Y`"
  )
  expect_error(
    .gest(A = dat$A, wgts = wgts, data = as.matrix(dat), Y = dat$Y,
          A.hat = Ahat, tx.var = "A", cts.obj = cts_bin,
          treat.mod.fitted = fit_tx),
    "`data`"
  )
  expect_error(
    .gest(A = dat$A, wgts = wgts, data = dat[1L:5L, ], Y = dat$Y,
          A.hat = Ahat, tx.var = "A", cts.obj = cts_bin,
          treat.mod.fitted = fit_tx),
    "`data`"
  )
  expect_error(
    .gest(A = dat$A, wgts = wgts, data = dat, Y = dat$Y,
          A.hat = Ahat[1L:5L], tx.var = "A", cts.obj = cts_bin,
          treat.mod.fitted = fit_tx),
    "`A.hat`"
  )
  expect_error(
    .gest(A = dat$A, wgts = wgts, data = dat, Y = dat$Y,
          A.hat = Ahat, tx.var = c("A", "B"), cts.obj = cts_bin,
          treat.mod.fitted = fit_tx),
    "`tx.var`"
  )
  expect_error(
    .gest(A = dat$A, wgts = wgts, data = dat, Y = dat$Y,
          A.hat = Ahat, tx.var = "A", cts.obj = list(),
          treat.mod.fitted = fit_tx),
    "`cts.obj`"
  )
  expect_error(
    .gest(A = dat$A, wgts = wgts, data = dat, Y = dat$Y,
          A.hat = Ahat, tx.var = "A", cts.obj = cts_bin,
          treat.mod.fitted = "not a model"),
    "`treat.mod.fitted`"
  )
  # Constant X makes design matrix rank-deficient
  dat_sing <- dat
  dat_sing$X <- 0.0
  cts_sing <- Binary$new(tf.model = ~X, blip.model = ~X, tx.var = "A")
  fit_sing <- glm(A ~ X, dat_sing, family = "binomial")
  Ahat_sing  <- predict(fit_sing, type = "response")
  expect_error(
    .gest(A = dat_sing$A, wgts = wgts, data = dat_sing, Y = dat_sing$Y,
          A.hat = Ahat_sing, tx.var = "A", cts.obj = cts_sing,
          treat.mod.fitted = fit_sing),
    "unable to invert matrix"
  )
})