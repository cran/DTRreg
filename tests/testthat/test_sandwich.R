set.seed(42L)
n <- 30L

dat <- data.frame(
  X = cos(seq(0.0, pi, length.out = n)),
  A = as.numeric(rbinom(n, 1L, 0.5)),
  Y = rnorm(n)
)

models <- list(blip = ~X, treat = A ~ X, tf = ~X)

cts_bin <- Binary$new(tf.model = models$tf,
                      blip.model = models$blip,
                      tx.var = "A")

# Fitted treatment model
fit_tx <- glm(models$treat, dat, family = "binomial")
Ahat <- predict(fit_tx, type = "response")
wgts <- rep(1.0, n)

# Fitted outcome models via trusted Layer 4 functions
dwols_fit <- .dwols(Y = dat$Y,
                    A = dat$A,
                    data = dat,
                    wgts = wgts,
                    cts.obj = cts_bin,
                    tx.var = "A")

gest_fit <- .gest(A = dat$A,
                  wgts = wgts,
                  data = dat,
                  Y = dat$Y,
                  A.hat = Ahat,
                  tx.var = "A",
                  cts.obj = cts_bin,
                  treat.mod.fitted = fit_tx)

n_theta <- length(coef(dwols_fit$outcome.fit))


test_that(".sandwich DWOLS: returns a matrix", {
  res <- .sandwich(cts.obj = cts_bin,
                   outcome.fit = dwols_fit$outcome.fit,
                   Y = dat$Y,
                   A = dat$A,
                   A.hat = Ahat,
                   wgts = wgts,
                   tx.var = "A",
                   data = dat,
                   treat.mod.fitted = fit_tx,
                   method = "dwols")
  expect_true(is.matrix(res))
})

test_that(".sandwich DWOLS: dimensions are n_theta x n_theta", {
  res <- .sandwich(cts.obj = cts_bin,
                   outcome.fit = dwols_fit$outcome.fit,
                   Y = dat$Y,
                   A = dat$A,
                   A.hat = Ahat,
                   wgts = wgts,
                   tx.var = "A",
                   data = dat,
                   treat.mod.fitted = fit_tx,
                   method = "dwols")
  expect_equal(dim(res), c(n_theta, n_theta))
})

test_that(".sandwich DWOLS: result is symmetric", {
  res <- .sandwich(cts.obj = cts_bin,
                   outcome.fit = dwols_fit$outcome.fit,
                   Y = dat$Y,
                   A = dat$A,
                   A.hat = Ahat,
                   wgts = wgts,
                   tx.var = "A",
                   data = dat,
                   treat.mod.fitted = fit_tx,
                   method = "dwols")
  expect_equal(res, t(res))
})

test_that(".sandwich GEST: dimensions are n_theta x n_theta", {
  res <- .sandwich(cts.obj = cts_bin,
                   outcome.fit = gest_fit$outcome.fit,
                   Y = dat$Y,
                   A = dat$A,
                   A.hat = Ahat,
                   wgts = wgts,
                   tx.var = "A",
                   data = dat,
                   treat.mod.fitted = fit_tx,
                   method = "gest")
  expect_equal(dim(res), c(n_theta, n_theta))
})

test_that(".sandwich GEST: result is symmetric", {
  res <- .sandwich(cts.obj = cts_bin,
                   outcome.fit = gest_fit$outcome.fit,
                   Y = dat$Y,
                   A = dat$A,
                   A.hat = Ahat,
                   wgts = wgts,
                   tx.var = "A",
                   data = dat,
                   treat.mod.fitted = fit_tx,
                   method = "gest")
  expect_equal(res, t(res))
})


local({
  # DWOLS oracle - write the six steps explicitly
  X_theta <- cts_bin$Hd(dat, dat$A)
  HW <- wgts * X_theta
  dat_pred <- dat
  dat_pred[, "A"] <- dat$A
  dat_pred <- cts_bin$prep(dat_pred, dat$A)
  Yhat <- drop(predict(dwols_fit$outcome.fit, dat_pred))
  U <- drop(dat$Y - Yhat) * HW
  dtheta <- {-1.0 / n} * crossprod(HW, X_theta)
  inv_mat <- solve(dtheta, t(U))
  expected <- {1.0 / n} * var(t(inv_mat))
  
  test_that(".sandwich DWOLS: matches explicit oracle computation", {
    res <- .sandwich(cts.obj = cts_bin,
                     outcome.fit = dwols_fit$outcome.fit,
                     Y = dat$Y,
                     A = dat$A,
                     A.hat = Ahat,
                     wgts = wgts,
                     tx.var = "A",
                     data = dat,
                     treat.mod.fitted = fit_tx,
                     method = "dwols")
    expect_equal(res, expected)
  })
})

local({
  # GEST oracle - Hw replaces wgts * Hd
  X_theta <- cts_bin$Hd(dat, dat$A)
  HW <- cts_bin$Hw(data = dat,
                   A = dat$A,
                   A.hat = Ahat,
                   wgt = wgts,
                   treat.mod.fitted = fit_tx)
  dat_pred <- dat
  dat_pred[, "A"] <- dat$A
  dat_pred <- cts_bin$prep(dat_pred, dat$A)
  Yhat <- drop(predict.GEST(gest_fit$outcome.fit, newdata = dat_pred))
  U <- drop(dat$Y - Yhat) * HW
  dtheta <- {-1.0 / n} * crossprod(HW, X_theta)
  inv_mat <- solve(dtheta, t(U))
  expected <- {1.0 / n} * var(t(inv_mat))
  
  test_that(".sandwich GEST: matches explicit oracle computation", {
    res <- .sandwich(cts.obj = cts_bin,
                     outcome.fit = gest_fit$outcome.fit,
                     Y = dat$Y,
                     A = dat$A,
                     A.hat = Ahat,
                     wgts = wgts,
                     tx.var = "A",
                     data = dat,
                     treat.mod.fitted = fit_tx,
                     method = "gest")
    expect_equal(res, expected)
  })
})


test_that(".sandwich: stops when cts.obj is not R6", {
  expect_error(
    .sandwich(cts.obj = list(),
              outcome.fit = dwols_fit$outcome.fit,
              Y = dat$Y,
              A = dat$A,
              A.hat = Ahat,
              wgts = wgts,
              tx.var = "A",
              data = dat,
              treat.mod.fitted = fit_tx,
              method = "dwols"),
    "`cts.obj`"
  )
  expect_error(
    .sandwich(cts.obj = cts_bin,
              outcome.fit = dwols_fit$outcome.fit,
              Y = as.character(dat$Y),
              A = dat$A,
              A.hat = Ahat,
              wgts = wgts,
              tx.var = "A",
              data = dat,
              treat.mod.fitted = fit_tx,
              method = "dwols"),
    "`Y`"
  )
  expect_error(
    .sandwich(cts.obj = cts_bin,
              outcome.fit = dwols_fit$outcome.fit,
              Y = dat$Y,
              A = dat$A[1L:5L],
              A.hat = Ahat,
              wgts = wgts,
              tx.var = "A",
              data = dat,
              treat.mod.fitted = fit_tx,
              method = "dwols"),
    "`A`"
  )
  expect_error(
    .sandwich(cts.obj = cts_bin,
              outcome.fit = dwols_fit$outcome.fit,
              Y = dat$Y,
              A = dat$A,
              A.hat = Ahat,
              wgts = wgts[1L:5L],
              tx.var = "A",
              data = dat,
              treat.mod.fitted = fit_tx,
              method = "dwols"),
    "`wgts`"
  )
  expect_error(
    .sandwich(cts.obj = cts_bin,
              outcome.fit = dwols_fit$outcome.fit,
              Y = dat$Y,
              A = dat$A,
              A.hat = Ahat,
              wgts = wgts,
              tx.var = c("A", "B"),
              data = dat,
              treat.mod.fitted = fit_tx,
              method = "dwols"),
    "`tx.var`"
  )
  expect_error(
    .sandwich(cts.obj = cts_bin,
              outcome.fit = dwols_fit$outcome.fit,
              Y = dat$Y,
              A = dat$A,
              A.hat = Ahat,
              wgts = wgts,
              tx.var = "A",
              data = dat[1L:5L, ],
              treat.mod.fitted = fit_tx,
              method = "dwols"),
    "`data`"
  )
  expect_error(
    .sandwich(cts.obj = cts_bin,
              outcome.fit = dwols_fit$outcome.fit,
              Y = dat$Y,
              A = dat$A,
              A.hat = Ahat,
              wgts = wgts,
              tx.var = "A",
              data = dat,
              treat.mod.fitted = fit_tx,
              method = 1L),
    "`method`"
  )
  dat_sing <- dat
  dat_sing$X <- 0.0
  cts_sing <- Binary$new(tf.model = ~X, blip.model = ~X, tx.var = "A")
  fit_sing <- .dwols(Y = dat_sing$Y,
                     A = dat_sing$A,
                     data = dat_sing,
                     wgts = wgts,
                     cts.obj = cts_sing,
                     tx.var = "A")
  expect_error(
    .sandwich(cts.obj = cts_sing,
              outcome.fit = fit_sing$outcome.fit,
              Y = dat_sing$Y,
              A = dat_sing$A,
              A.hat = Ahat,
              wgts = wgts,
              tx.var = "A",
              data = dat_sing,
              treat.mod.fitted = fit_tx,
              method = "dwols"),
    "unable to invert matrix"
  )
})