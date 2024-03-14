# testing 12/31/2023 updates

dev <- FALSE

test_that("intercept only is not an issue; censored data", {
  if (!dev) skip("for development only")
  data("twoStageCens")
  
  for (method in c("dwols", "qlearn")) {
    for (weight in c("abs", "ipw", "cipw", "none")) {
      for (var.estim in c("none", "bootstrap", "sandwich")) {
        for (full.cov in c(TRUE, FALSE)) {
          expect_no_error(DWSurv(time = list(~ T1, ~ T2), 
                                 blip.mod = list(~ 1, ~ 1), 
                                 treat.mod = list(A1 ~ 1, A2 ~ 1), 
                                 tf.mod = list(~ 1, ~ 1), 
                                 cens.mod = list(delta ~ 1, delta ~ 1), 
                                 data = twoStageCens, 
                                 method = method,
                                 var.estim = var.estim,
                                 weight = weight))
        }
      }
    }
  }
  
})

test_that("intercept only is not an issue; censored multinomial data", {
  if (!dev) skip("for development only")
  data("twoStageCens")
  twoStageCens$A1 <- withr::with_seed(1234, sample(1:3, nrow(twoStageCens), replace = TRUE))
  twoStageCens$A2[!is.na(twoStageCens$A2)] <- withr::with_seed(2345, 
                                                               sample(1:3, 
                                                                      sum(!is.na(twoStageCens$A2)), 
                                                                      replace = TRUE))
  
  for (method in c("dwols", "qlearn")) {
    for (weight in c("ipw", "cipw", "qpom", "wo", "none")) {
      for (var.estim in c("none", "bootstrap", "sandwich")) {
        for (full.cov in c(TRUE, FALSE)) {
          expect_no_error(DWSurv(time = list(~ T1, ~ T2), 
                                 blip.mod = list(~ 1, ~ 1), 
                                 treat.mod = list(A1 ~ 1, A2 ~ 1), 
                                 tf.mod = list(~ 1, ~ 1), 
                                 cens.mod = list(delta ~ 1, delta ~ 1), 
                                 treat.type = "multi",
                                 data = twoStageCens, 
                                 method = method,
                                 var.estim = var.estim,
                                 weight = weight))
        }
      }
    }
  }
  
})

test_that("intercept only is not an issue; censored continuous data", {
  if (!dev) skip("for development only")
  data("twoStageCens")
  twoStageCens$A1 <- withr::with_seed(1234, stats::runif(nrow(twoStageCens)))
  twoStageCens$A2[!is.na(twoStageCens$A2)] <- 
    withr::with_seed(2345, 
                     stats::runif(sum(!is.na(twoStageCens$A2))))
  
  for (method in c("dwols", "qlearn")) {
    for (weight in c("ipw", "cipw", "qpom", "wo", "none")) {
      for (var.estim in c("none", "bootstrap", "sandwich")) {
        for (full.cov in c(TRUE, FALSE)) {
          expect_no_error(DWSurv(time = list(~ T1, ~ T2), 
                                 blip.mod = list(~ 1 + A1, ~ 1 + A2), 
                                 treat.mod = list(A1 ~ 1, A2 ~ 1), 
                                 tf.mod = list(~ 1, ~ 1), 
                                 cens.mod = list(delta ~ 1, delta ~ 1), 
                                 treat.type = "cont",
                                 data = twoStageCens, 
                                 method = method,
                                 var.estim = var.estim,
                                 weight = weight))
        }
      }
    }
  }
  
})
  
test_that("intercept only is not an issue; survival data", {
  if (!dev) skip("for development only")
  data("twoStageSurv")
  
  for (method in c("dwols", "qlearn")) {
    for (weight in c("abs", "ipw", "cipw", "none")) {
      for (var.estim in c("none", "bootstrap", "sandwich")) {
        for (full.cov in c(TRUE, FALSE)) {
          expect_no_error(DWSurv(time = list(~ T1, ~ T2), 
                                 blip.mod = list(~ 1, ~ 1), 
                                 treat.mod = list(A1 ~ 1, A2 ~ 1), 
                                 tf.mod = list(~ 1, ~ 1), 
                                 cens.mod = list( ~ delta, ~ delta), 
                                 data = twoStageSurv, 
                                 method = method,
                                 var.estim = var.estim,
                                 weight = weight))
        }
      }
    }
  }
  
})

test_that("intercept only is not an issue; survival multinomial data", {
  if (!dev) skip("for development only")
  data("twoStageSurv")
  twoStageSurv$A1 <- withr::with_seed(1234, sample(1:3, nrow(twoStageSurv), replace = TRUE))
  twoStageSurv$A2[!is.na(twoStageSurv$A2)] <- withr::with_seed(2345, 
                                                               sample(1:3, 
                                                                      sum(!is.na(twoStageSurv$A2)), 
                                                                      replace = TRUE))
  
  for (method in c("dwols", "qlearn")) {
    for (weight in c("ipw", "cipw", "qpom", "wo", "none")) {
      for (var.estim in c("none", "bootstrap", "sandwich")) {
        for (full.cov in c(TRUE, FALSE)) {
          expect_no_error(DWSurv(time = list(~ T1, ~ T2), 
                                 blip.mod = list(~ 1, ~ 1), 
                                 treat.mod = list(A1 ~ 1, A2 ~ 1), 
                                 tf.mod = list(~ 1, ~ 1), 
                                 cens.mod = list(~ delta, ~ delta), 
                                 treat.type = "multi",
                                 data = twoStageSurv, 
                                 method = method,
                                 var.estim = var.estim,
                                 weight = weight))
        }
      }
    }
  }
  
})

test_that("intercept only is not an issue; survival continuous data", {
  if (!dev) skip("for development only")
  data("twoStageSurv")
  twoStageSurv$A1 <- withr::with_seed(1234, stats::runif(nrow(twoStageSurv)))
  twoStageSurv$A2[!is.na(twoStageSurv$A2)] <- 
    withr::with_seed(2345, 
                     stats::runif(sum(!is.na(twoStageSurv$A2))))
  
  for (method in c("dwols", "qlearn")) {
    for (weight in c("ipw", "cipw", "qpom", "wo", "none")) {
      for (var.estim in c("none", "bootstrap", "sandwich")) {
        for (full.cov in c(TRUE, FALSE)) {
          expect_no_error(DWSurv(time = list(~ T1, ~ T2), 
                                 blip.mod = list(~ 1 + A1, ~ 1 + A2), 
                                 treat.mod = list(A1 ~ 1, A2 ~ 1), 
                                 tf.mod = list(~ 1, ~ 1), 
                                 cens.mod = list(~ delta, ~ delta), 
                                 treat.type = "cont",
                                 data = twoStageSurv, 
                                 method = method,
                                 var.estim = var.estim,
                                 weight = weight))
        }
      }
    }
  }
  
})

test_that("intercept only is not an issue; cont data", {
  if (!dev) skip("for development only")
  data("twoStageCont")

  for (method in c("gest", "dwols", "qlearn")) {
    for (weight in c("abs", "ipw", "cipw", "none")) {
      for (var.estim in c("none", "bootstrap", "sandwich")) {
        for (full.cov in c(TRUE, FALSE)) {
          for (miss in c("drop", "ipw")) {
            expect_no_error(DTRreg(twoStageCont$Y, 
                                   blip.mod = list(~ 1, ~ 1), 
                                   treat.mod = list(A1 ~ 1, A2 ~ 1), 
                                   tf.mod = list(~ 1, ~ 1), 
                                   weight = weight, 
                                   method = method,
                                   var.estim = var.estim, 
                                   data = twoStageCont,
                                   missing = miss,
                                   missing.mod = list(~ 1, ~ 1)))
          }
        }
      }
    }
  }
  
})

test_that("intercept only is not an issue; cont multinomial data", {
  if (!dev) skip("for development only")
  data("twoStageCont", package = "DTRreg")
  twoStageCont$A1 <- withr::with_seed(1234, sample(1:3, nrow(twoStageCont), replace = TRUE))
  twoStageCont$A2[!is.na(twoStageCont$A2)] <- withr::with_seed(2345, 
                                                               sample(1:3, 
                                                                      sum(!is.na(twoStageCont$A2)), 
                                                                      replace = TRUE))
  
  for (method in c("dwols", "qlearn")) {
    for (weight in c("ipw", "cipw", "qpom", "wo", "none")) {
      for (var.estim in c("none", "bootstrap", "sandwich")) {
        for (full.cov in c(TRUE, FALSE)) {
          for (miss in c("drop", "ipw")) {
            expect_no_error(DTRreg(twoStageCont$Y, 
                                   blip.mod = list(~ 1, ~ 1), 
                                   treat.mod = list(A1 ~ 1, A2 ~ 1), 
                                   treat.type = "multi",
                                   tf.mod = list(~ 1, ~ 1), 
                                   method = method,
                                   weight = weight,
                                   var.estim = var.estim, 
                                   data = twoStageCont,
                                   missing = miss,
                                   missing.mod = list(~ 1, ~ 1)))
          }
        }
      }
    }
  }
  
})

test_that("intercept only is not an issue; cont continuous data", {
  if (!dev) skip("for development only")
  data("twoStageCont", package = "DTRreg")
  twoStageCont$A1 <- withr::with_seed(1234, stats::runif(nrow(twoStageCont)))
  twoStageCont$A2[!is.na(twoStageCont$A2)] <- withr::with_seed(2345, 
                                                               stats::runif(!is.na(twoStageCont$A2)))
  
  for (method in c("gest", "dwols", "qlearn")) {
    for (weight in c("ipw", "cipw", "qpom", "wo", "none")) {
      for (var.estim in c("none", "bootstrap", "sandwich")) {
        for (full.cov in c(TRUE, FALSE)) {
          for (miss in c("drop", "ipw")) {
            expect_no_error(DTRreg(twoStageCont$Y, 
                                   blip.mod = list(~ 1 + A1, ~ 1 + A2), 
                                   treat.mod = list(A1 ~ 1, A2 ~ 1), 
                                   treat.type = "cont",
                                   tf.mod = list(~ 1, ~ 1), 
                                   method = method,
                                   weight = weight,
                                   var.estim = var.estim, 
                                   data = twoStageCont,
                                   missing = miss,
                                   missing.mod = list(~ 1, ~ 1)))
          }
        }
      }
    }
  }
  
})

test_that("testing that qlearn in DWsurv is just setting treatment weights to 1.0", {
  data("twoStageCens")
  
  original <- DWSurv(time = list(~ T1, ~ T2), 
                     blip.mod = list(~ X11, ~ X21), 
                     treat.mod = list(A1 ~ X11, A2 ~ X21), 
                     tf.mod = list(~ X11 + X12, ~ X21), 
                     cens.mod = list(delta ~ X11, delta ~ X21),
                     data = twoStageCens,
                     method = "qlearn")
  
  manual <- DWSurv(time = list(~ T1, ~ T2), 
                   blip.mod = list(~ X11, ~ X21), 
                   treat.mod = list(A1 ~ X11, A2 ~ X21), 
                   tf.mod = list(~ X11 + X12, ~ X21), 
                   cens.mod = list(delta ~ X11, delta ~ X21),
                   data = twoStageCens,
                   weight = "manual",
                   treat.wgt.man = list(rep(1.0, nrow(twoStageCens)),
                                        rep(1.0, nrow(twoStageCens))),
                   method = "dwols")
  
  expect_equal(original$beta, manual$beta)
  expect_equal(original$psi, manual$psi)
  expect_equal(original$analysis$opt.treat, manual$analysis$opt.treat)
  expect_equal(original$analysis$opt.Y, manual$analysis$opt.Y)
  
  # use this to extract the manual weights
  censor_wgts <- original$analysis$cens.wgt
  treat_wgts <- original$analysis$tx.wgt
  
  tx_man <- list()
  tx_man[[1L]] <- censor_wgts[[1L]]
  tx_man[[2L]] <- rep(NA_real_, nrow(twoStageCens))
  tx_man[[2L]][!is.na(twoStageCens$T2)] <- censor_wgts[[2L]]
  
  manual <- DWSurv(time = list(~ T1, ~ T2), 
                   blip.mod = list(~ X11, ~ X21), 
                   treat.mod = list(A1 ~ X11, A2 ~ X21), 
                   tf.mod = list(~ X11 + X12, ~ X21), 
                   cens.mod = list( ~ delta, ~ delta),
                   data = twoStageCens,
                   weight = "manual.with.censor",
                   treat.wgt.man = tx_man,
                   method = "dwols")
  
  expect_equal(original$beta, manual$beta)
  expect_equal(original$psi, manual$psi)
  expect_equal(original$analysis$opt.treat, manual$analysis$opt.treat)
  expect_equal(original$analysis$opt.Y, manual$analysis$opt.Y)
})


test_that("testing that using manual overrides any treatment weight specification", {
  data("twoStageCens")
  
  original <- DWSurv(time = list(~ T1, ~ T2), 
                     blip.mod = list(~ X11, ~ X21), 
                     treat.mod = list(A1 ~ X11, A2 ~ X21), 
                     tf.mod = list(~ X11 + X12, ~ X21), 
                     cens.mod = list(delta ~ X11, delta ~ X21),
                     data = twoStageCens,
                     method = "dwols",
                     weight = "manual",
                     treat.wgt.man = list(rep(1.0, nrow(twoStageCens)),
                                          rep(1.0, nrow(twoStageCens))))
  
  manual <- DWSurv(time = list(~ T1, ~ T2), 
                   blip.mod = list(~ X11, ~ X21), 
                   treat.mod = list(A1 ~ X11, A2 ~ X21), 
                   tf.mod = list(~ X11 + X12, ~ X21), 
                   cens.mod = list(delta ~ X11, delta ~ X21),
                   data = twoStageCens,
                   method = "qlearn")
  
  expect_equal(original$beta, manual$beta)
  expect_equal(original$psi, manual$psi)
  expect_equal(original$analysis$opt.treat, manual$analysis$opt.treat)
  expect_equal(original$analysis$opt.Y, manual$analysis$opt.Y)
  
})

test_that("testing the ability to specify censoring and treatment through manual input", {
  data("twoStageCens")
  
  expect_no_error(DWSurv(time = list(~ T1, ~ T2), 
                         blip.mod = list(~ 1, ~ 1), 
                         treat.mod = list(A1 ~ 1, A2 ~ 1), 
                         tf.mod = list(~ 1, ~ 1), 
                         cens.mod = list(delta ~ 1, delta ~ 1), 
                         var.estim = "sandwich", 
                         data = twoStageCens))
  
  expect_no_error(DWSurv(time = list(~ T1, ~ T2), 
                         blip.mod = list(~ 1, ~ 1), 
                         treat.mod = list(A1 ~ 1, A2 ~ 1), 
                         tf.mod = list(~ 1, ~ 1), 
                         cens.mod = list(delta ~ 1, delta ~ 1), 
                         var.estim = "boot", 
                         data = twoStageCens))
  
  original <- DWSurv(time = list(~ T1, ~ T2), 
                     blip.mod = list(~ X11, ~ 1), 
                     treat.mod = list(A1 ~ X11, A2 ~ 1), 
                     tf.mod = list(~ X11 + X12, ~ 1), 
                     cens.mod = list(delta ~ 1, delta ~ X11), 
                     var.estim = "sandwich", 
                     data = twoStageCens)
  
  # use this to extract the manual weights
  censor_wgts <- original$analysis$cens.wgt
  treat_wgts <- original$analysis$tx.wgt
  
  tx_man <- list()
  tx_man[[1L]] <- censor_wgts[[1L]] * treat_wgts[[1L]]
  tx_man[[2L]] <- treat_wgts[[2L]]
  tx_man[[2L]][!is.na(treat_wgts[[2L]])] <- censor_wgts[[2L]]
  
  manual <- DWSurv(time = list(~ T1, ~ T2), 
                   blip.mod = list(~ X11, ~ 1), 
                   treat.mod = list(A1 ~ X11, A2 ~ 1), 
                   tf.mod = list(~ X11 + X12, ~ 1), 
                   cens.mod = list(~ delta, ~ delta), 
                   var.estim = "sandwich", 
                   weight = "manual.with.censor",
                   treat.wgt.man = tx_man,
                   data = twoStageCens)
  
  expect_equal(original$beta, manual$beta)
  expect_equal(original$psi, manual$psi)
  expect_equal(original$analysis$opt.treat, manual$analysis$opt.treat)
  expect_equal(original$analysis$opt.Y, manual$analysis$opt.Y)
  
})

test_that("testing the ability to specify censoring and treatment through manual input doesn't change DTRreg", {
  data(twoStageCont)
  
  # models to be passed to DTRreg
  # blip model
  blip.mod <- list(~ X1, ~ X2)
  # treatment model (correctly specified)
  treat.mod <- list(A1 ~ X1, A2 ~ X2)
  # treatment-free model (incorrectly specified)
  tf.mod <- list(~ X1, ~ X2)
  
  # perform G-estimation
  expect_no_error(original <- DTRreg(twoStageCont$Y, blip.mod, treat.mod, tf.mod, 
                         data = twoStageCont, method = "dwols"))
  
  treat_wgts <- original$analysis$tx.wgt
  tx_man <- list()
  tx_man[[1L]] <- treat_wgts[[1L]]
  tx_man[[2L]] <- treat_wgts[[2L]]

  expect_no_error(manual <- DTRreg(twoStageCont$Y, blip.mod, treat.mod, tf.mod, 
                                   data = twoStageCont, method = "dwols",
                                   weight = "manual", treat.wgt.man = tx_man))
  
  expect_equal(original$beta, manual$beta)
  expect_equal(original$psi, manual$psi)
  expect_equal(original$analysis$opt.treat, manual$analysis$opt.treat)
  expect_equal(original$analysis$opt.Y, manual$analysis$opt.Y)
  
  expect_error(DTRreg(twoStageCont$Y, blip.mod, treat.mod, tf.mod, 
                      data = twoStageCont, method = "dwols",
                      weight = "manual.with.censor", treat.wgt.man = tx_man))
})

