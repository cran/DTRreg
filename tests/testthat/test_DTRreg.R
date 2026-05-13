# Shared data

data(twoStageCont)
n <- nrow(twoStageCont)
K <- 2L

blip.mod <- list(~X1, ~X2)
treat.mod <- list(A1 ~ X1, A2 ~ X2)
tf.mod <- list(~X1, ~X2)
qblip.mod <- list(~X1+A1, ~X2+A2)


# Output structure - DWOLS, var.estim = "none"


local({
  mod <- suppressMessages(
    DTRreg(twoStageCont$Y, blip.mod, treat.mod, tf.mod,
           data = twoStageCont, method = "dwols", var.estim = "none")
  )
  
  test_that("DTRreg DWOLS: class is 'DTRreg'", {
    expect_true(inherits(mod, "DTRreg"))
  })
  
  test_that("DTRreg DWOLS: K matches number of stages", {
    expect_equal(mod$K, K)
  })
  
  test_that("DTRreg DWOLS: top-level names contain required elements", {
    expect_true(all(c("K", "beta", "psi", "covmat", "nonreg",
                      "setup", "training_data", "analysis", "call") %in% names(mod)))
  })
  
  # - psi / beta -
  
  test_that("DTRreg DWOLS: psi is a list of length K", {
    expect_true(is.list(mod$psi))
    expect_length(mod$psi, K)
  })
  
  test_that("DTRreg DWOLS: psi names contain treatment variable at each stage", {
    expect_true(all(grepl("A1", names(mod$psi[[1L]]))))
    expect_true(all(grepl("A2", names(mod$psi[[2L]]))))
  })
  
  test_that("DTRreg DWOLS: beta is a list of length K", {
    expect_true(is.list(mod$beta))
    expect_length(mod$beta, K)
  })
  
  test_that("DTRreg DWOLS: beta names do not contain treatment variable", {
    expect_false(any(grepl("^A1", names(mod$beta[[1L]]))))
    expect_false(any(grepl("^A2", names(mod$beta[[2L]]))))
  })
  
  test_that("DTRreg DWOLS: covmat is NA when var.estim='none'", {
    expect_true(is.na(mod$covmat))
  })
  
  test_that("DTRreg DWOLS: psi.boot absent when var.estim='none'", {
    expect_null(mod$psi.boot)
  })
  
  # - setup --
  
  test_that("DTRreg DWOLS: setup contains required fields", {
    expect_true(all(c("models", "method", "var.estim", "tx.weight",
                      "tx.type", "type", "tx.vars") %in% names(mod$setup)))
  })
  
  test_that("DTRreg DWOLS: setup$method matches requested method", {
    expect_equal(mod$setup$method, "dwols")
  })
  
  test_that("DTRreg DWOLS: setup$var.estim matches requested var.estim", {
    expect_equal(mod$setup$var.estim, "none")
  })
  
  test_that("DTRreg DWOLS: setup$tx.vars matches treatment variable names", {
    expect_equal(mod$setup$tx.vars, c("A1", "A2"))
  })
  
  test_that("DTRreg DWOLS: setup$models has length K", {
    expect_length(mod$setup$models, K)
  })
  
  # - training_data -
  
  test_that("DTRreg DWOLS: training_data contains data, outcome, A", {
    expect_named(mod$training_data, c("data", "outcome", "A"), ignore.order = TRUE)
  })
  
  test_that("DTRreg DWOLS: training_data$data has n rows", {
    expect_equal(nrow(mod$training_data$data), n)
  })
  
  test_that("DTRreg DWOLS: training_data$outcome matches supplied outcome", {
    expect_equal(mod$training_data$outcome, twoStageCont$Y)
  })
  
  # - analysis -
  
  test_that("DTRreg DWOLS: analysis contains required elements", {
    expect_true(all(c("n", "last.stage", "prob.cc", "cc.mod.fitted",
                      "cc.wgt", "cts", "tx.mod.fitted", "A.hat", "tx.wgt",
                      "outcome.fit", "Y", "regret", "opt.treat",
                      "opt.Y", "fitted.values", "residuals",
                      "blip.data") %in% names(mod$analysis)))
  })
  
  test_that("DTRreg DWOLS: analysis$n is list of length K", {
    expect_length(mod$analysis$n, K)
  })
  
  test_that("DTRreg DWOLS: analysis$last.stage is integer of length n", {
    expect_true(is.integer(mod$analysis$last.stage))
    expect_length(mod$analysis$last.stage, n)
  })
  
  test_that("DTRreg DWOLS: analysis$opt.Y is numeric of length n", {
    expect_true(is.numeric(mod$analysis$opt.Y))
    expect_length(mod$analysis$opt.Y, n)
  })
  
  test_that("DTRreg DWOLS: analysis$cts[[1]] is 'bin' for binary treatment", {
    expect_equal(mod$analysis$cts[[1L]], "bin")
  })
  
  test_that("DTRreg DWOLS: analysis$outcome.fit is list of K fitted models", {
    expect_length(mod$analysis$outcome.fit, K)
    expect_true(inherits(mod$analysis$outcome.fit[[1L]], "lm"))
  })
  
  test_that("DTRreg DWOLS: analysis$opt.treat is list of length K", {
    expect_length(mod$analysis$opt.treat, K)
  })
})



# GEST method


local({
  mod <- suppressMessages(
    DTRreg(twoStageCont$Y, blip.mod, treat.mod, tf.mod,
           data = twoStageCont, method = "gest", var.estim = "none")
  )
  
  test_that("DTRreg GEST: class is 'DTRreg'", {
    expect_true(inherits(mod, "DTRreg"))
  })
  
  test_that("DTRreg GEST: psi names contain treatment variable", {
    expect_true(all(grepl("A1", names(mod$psi[[1L]]))))
    expect_true(all(grepl("A2", names(mod$psi[[2L]]))))
  })
  
  test_that("DTRreg GEST: outcome.fit has class GEST", {
    expect_true(inherits(mod$analysis$outcome.fit[[1L]], "GEST"))
    expect_true(inherits(mod$analysis$outcome.fit[[2L]], "GEST"))
  })
  
  test_that("DTRreg GEST: setup$method is 'gest'", {
    expect_equal(mod$setup$method, "gest")
  })
})



# Sandwich variance


local({
  mod <- suppressMessages(
    DTRreg(twoStageCont$Y, blip.mod, treat.mod, tf.mod,
           data = twoStageCont, method = "dwols", var.estim = "sandwich")
  )
  
  test_that("DTRreg sandwich: covmat is a list of length K", {
    expect_true(is.list(mod$covmat))
    expect_length(mod$covmat, K)
  })
  
  test_that("DTRreg sandwich: each covmat element is a square matrix", {
    for (k in seq_len(K)) {
      expect_true(is.matrix(mod$covmat[[k]]))
      expect_equal(nrow(mod$covmat[[k]]), ncol(mod$covmat[[k]]))
    }
  })
  
  test_that("DTRreg sandwich: covmat row/col names contain psi names", {
    for (k in seq_len(K)) {
      expect_true(all(names(mod$psi[[k]]) %in% rownames(mod$covmat[[k]])))
    }
  })
  
  test_that("DTRreg sandwich: nonreg is numeric of length K for binary treatment", {
    expect_true(is.numeric(mod$nonreg))
    expect_length(mod$nonreg, K)
  })
  
  test_that("DTRreg sandwich: psi.boot absent for sandwich", {
    expect_null(mod$psi.boot)
  })
})



# Bootstrap variance


local({
  withr::with_seed(42L, {
    mod <- suppressMessages(
      DTRreg(twoStageCont$Y, blip.mod, treat.mod, tf.mod,
             data = twoStageCont,
             method = "dwols",
             var.estim = "bootstrap",
             bootstrap.controls = list(B = 10L))
    )
  })
  
  test_that("DTRreg bootstrap: covmat is list of length K", {
    expect_true(is.list(mod$covmat))
    expect_length(mod$covmat, K)
  })
  
  test_that("DTRreg bootstrap: psi.boot is list of length K", {
    expect_true(is.list(mod$psi.boot))
    expect_length(mod$psi.boot, K)
  })
  
  test_that("DTRreg bootstrap: each psi.boot matrix has B rows", {
    for (k in seq_len(K)) {
      expect_equal(nrow(mod$psi.boot[[k]]), 10L)
    }
  })
})



# type = "effect" (blip rather than DTR)


test_that("DTRreg effect type: setup$type is 'effect'", {
  mod <- suppressMessages(
    DTRreg(twoStageCont$Y, blip.mod, treat.mod, tf.mod,
           data = twoStageCont, var.estim = "none", dtr = FALSE)
  )
  expect_equal(mod$setup$type, "effect")
})



# Missing data - IPW


test_that("DTRreg missing='ipw': runs without error", {
  dat_na <- twoStageCont
  dat_na$X2[1L:10L] <- NA_real_
  expect_no_error(suppressMessages(
    DTRreg(dat_na$Y, blip.mod, treat.mod, tf.mod,
           data = dat_na,
           missing = "ipw",
           var.estim = "none")
  ))
})

test_that("DTRreg missing='drop': runs without error", {
  dat_na <- twoStageCont
  dat_na$X2[1L:10L] <- NA_real_
  expect_no_error(suppressMessages(
    DTRreg(dat_na$Y, blip.mod, treat.mod, tf.mod,
           data = dat_na,
           missing = "drop",
           var.estim = "none")
  ))
})



# Weight options


test_that("DTRreg weight='none': runs without error", {
  expect_no_error(suppressMessages(
    DTRreg(twoStageCont$Y, blip.mod, treat.mod, tf.mod,
           data = twoStageCont, weight = "none", var.estim = "none")
  ))
})

test_that("DTRreg weight='cipw': runs without error", {
  expect_no_error(suppressMessages(
    DTRreg(twoStageCont$Y, blip.mod, treat.mod, tf.mod,
           data = twoStageCont, weight = "cipw", var.estim = "none")
  ))
})

test_that("DTRreg weight='overlap': runs without error for binary treatment", {
  expect_no_error(suppressMessages(
    DTRreg(twoStageCont$Y, blip.mod, treat.mod, tf.mod,
           data = twoStageCont, weight = "overlap", var.estim = "none")
  ))
})

test_that("DTRreg weight='abs': runs with message", {
  expect_warning(suppressMessages(
    DTRreg(twoStageCont$Y, blip.mod, treat.mod, tf.mod,
           data = twoStageCont, weight = "abs", var.estim = "none")
  ), "renamed 'overlap'")
})

test_that("DTRreg weight='manual': runs and psi matches ipw result", {
  mod_ipw <- suppressMessages(
    DTRreg(twoStageCont$Y, blip.mod, treat.mod, tf.mod, method = "dwols",
           data = twoStageCont, weight = "ipw", var.estim = "none")
  )
  # Extract computed weights and pass as manual
  tx_man <- list(
    mod_ipw$analysis$tx.wgt[[1L]],
    mod_ipw$analysis$tx.wgt[[2L]]
  )
  mod_man <- suppressMessages(
    DTRreg(twoStageCont$Y, blip.mod, treat.mod, tf.mod,
           data = twoStageCont, weight = "manual", method = "dwols",
           treat.wgt.man = tx_man, var.estim = "none")
  )
  expect_equal(mod_man$psi, mod_ipw$psi)
})



# full.cov = TRUE


test_that("DTRreg full.cov=TRUE: covmat larger than full.cov=FALSE", {
  mod_full <- suppressMessages(
    DTRreg(twoStageCont$Y, blip.mod, treat.mod, tf.mod,
           data = twoStageCont, var.estim = "sandwich", full.cov = TRUE)
  )
  mod_blip <- suppressMessages(
    DTRreg(twoStageCont$Y, blip.mod, treat.mod, tf.mod,
           data = twoStageCont, var.estim = "sandwich", full.cov = FALSE)
  )
  for (k in seq_len(K)) {
    expect_true(nrow(mod_full$covmat[[k]]) >= nrow(mod_blip$covmat[[k]]))
  }
})



# Stopping conditions


test_that("DTRreg: stops when outcome has missing values", {
  Y_na <- twoStageCont$Y
  Y_na[1L] <- NA_real_
  expect_error(
    DTRreg(Y_na, blip.mod, treat.mod, tf.mod, data = twoStageCont),
    "outcome"
  )
})

test_that("DTRreg: stops when data is not a data.frame", {
  expect_error(
    DTRreg(twoStageCont$Y, blip.mod, treat.mod, tf.mod,
           data = as.matrix(twoStageCont)),
    "data.frame"
  )
})

test_that("DTRreg: stops when method is unrecognised", {
  expect_error(
    DTRreg(twoStageCont$Y, blip.mod, treat.mod, tf.mod,
           data = twoStageCont, method = "qbayesian"),
    "should be one of"
  )
})

test_that("DTRreg: deprecated 'abs' weight triggers warning", {
  expect_warning(
    tryCatch(
      suppressMessages(
        DTRreg(twoStageCont$Y, blip.mod, treat.mod, tf.mod,
               data = twoStageCont, weight = "abs", var.estim = "none")
      ),
      error = function(e) NULL
    ),
    "renamed.*overlap|overlap.*renamed"
  )
})

test_that("DTRreg: deprecated 'wo' weight triggers warning", {
  expect_warning(
    tryCatch(
      suppressMessages(
        DTRreg(twoStageCont$Y, blip.mod, treat.mod, tf.mod,
               data = twoStageCont, weight = "wo", var.estim = "none")
      ),
      error = function(e) NULL
    ),
    "renamed.*overlap|overlap.*renamed"
  )
})
test_that("DTRreg: stops when dtr is not logical", {
  expect_error(
    DTRreg(twoStageCont$Y, blip.mod, treat.mod, tf.mod,
           data = twoStageCont, dtr = "yes"),
    "dtr"
  )
})

test_that("DTRreg: stops when blip.mod and treat.mod have different lengths", {
  expect_error(
    DTRreg(twoStageCont$Y, list(~X1), treat.mod, tf.mod,
           data = twoStageCont),
    "same length"
  )
})

test_that("DTRreg: warns when gest used with manual weights", {
  expect_warning(
    suppressMessages(
      DTRreg(twoStageCont$Y, blip.mod, treat.mod, tf.mod,
             data = twoStageCont,
             method = "gest",
             weight = "manual",
             treat.wgt.man = list(rep(1.0, n), rep(1.0, n)),
             var.estim = "none")
    ),
    "not used in g-estimation"
  )
})