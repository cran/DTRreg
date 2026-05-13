dat <- data.frame(
  X = c(1.0, 2.0, 3.0, 4.0, 5.0),
  Z = c(2.0, 4.0, 6.0, 8.0, 10.0),
  Y = c(0.1, 0.2, 0.3, 0.4, 0.5)
)

test_that(".getModelMatrix: returns same matrix as direct model.matrix call", {
  model <- ~X + Z
  expected <- model.matrix(model, dat)
  result <- .getModelMatrix(model, dat)
  expect_equal(result, expected)
})

test_that(".getModelMatrix: column names match formula terms", {
  result <- .getModelMatrix(~X + Z, dat)
  expect_true("X" %in% colnames(result))
  expect_true("Z" %in% colnames(result))
  expect_true("(Intercept)" %in% colnames(result))
})

test_that(".getModelMatrix: returns matrix with correct dimensions", {
  result <- .getModelMatrix(~X + Z, dat)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), nrow(dat))
  expect_equal(ncol(result), 3L) # intercept + X + Z
})

test_that(".getModelMatrix: interaction terms are constructed correctly", {
  model <- ~X + Z + X:Z
  expected <- model.matrix(model, dat)
  result <- .getModelMatrix(model, dat)
  expect_equal(result, expected)
})

test_that(".getModelMatrix: preserves rows with NA rather than dropping them", {
  dat_na <- dat
  dat_na$X[3L] <- NA_real_
  result <- .getModelMatrix(~X + Z, dat_na)
  expect_equal(nrow(result), nrow(dat_na))
  expect_true(any(is.na(result)))
})

test_that(".getModelMatrix: stops with informative message when column is absent", {
  expect_error(.getModelMatrix(~X + W, dat), "unable to extract model frame")
})

local({
  
  set.seed(42L)
  n_u <- 20L
  
  dat_u <- data.frame(
    X = cos(seq(0.0, pi, length.out = n_u)),
    A1 = as.integer(rbinom(n_u, 1L, 0.5)),
    A2 = as.integer(rbinom(n_u, 1L, 0.5)),
    T1 = abs(rnorm(n_u)) + 0.5,
    T2 = abs(rnorm(n_u)) + 0.5,
    Y = rnorm(n_u)
  )
  
  cts1 <- Binary$new(tf.model = ~X, blip.model = ~X, tx.var = "A1")
  cts2 <- Binary$new(tf.model = ~X, blip.model = ~X, tx.var = "A2")
  
  wgts_u <- rep(1.0, n_u)
  fit1 <- .dwols(Y = dat_u$Y, A = dat_u$A1, data = dat_u,
                 wgts = wgts_u, cts.obj = cts1, tx.var = "A1")
  fit2 <- .dwols(Y = dat_u$Y, A = dat_u$A2, data = dat_u,
                 wgts = wgts_u, cts.obj = cts2, tx.var = "A2")
  
  # All participants complete at both stages
  cc_u <- list(
    last.stage = rep(2L, n_u),
    prob.complete.case = matrix(1.0, nrow = n_u, ncol = 2L)
  )
  
  # Minimal steps — opt.treat pre-set so .stageY and .stageRegret can run
  make_steps <- function(opt1 = dat_u$A1, opt2 = dat_u$A2) {
    list(
      list(dp = 1L,
           tx.var = "A1",
           cts = "bin",
           cts.obj = cts1,
           outcome.fit = fit1$outcome.fit,
           opt.treat = opt1,
           A = dat_u$A1),
      list(dp = 2L,
           tx.var = "A2",
           cts = "bin",
           cts.obj = cts2,
           outcome.fit = fit2$outcome.fit,
           opt.treat = opt2,
           A = dat_u$A2)
    )
  }
  
  test_that(".stageTx: returns list of same length as input steps", {
    res <- .stageTx(make_steps(), cc_u, dat_u, quiet = TRUE)
    expect_length(res, 2L)
  })
  
  test_that(".stageTx: sets opt.treat on each stage", {
    res <- .stageTx(make_steps(), cc_u, dat_u, quiet = TRUE)
    expect_false(is.null(res[[1L]]$opt.treat))
    expect_false(is.null(res[[2L]]$opt.treat))
  })
  
  test_that(".stageTx: opt.treat for binary contains only 0/1/NA", {
    res <- .stageTx(make_steps(), cc_u, dat_u, quiet = TRUE)
    expect_true(all(res[[1L]]$opt.treat %in% c(0L, 1L, NA_integer_)))
    expect_true(all(res[[2L]]$opt.treat %in% c(0L, 1L, NA_integer_)))
  })
  
  test_that(".stageTx: oracle — opt.treat at stage_cases matches cts.obj$opt", {
    res <- .stageTx(make_steps(), cc_u, dat_u, quiet = TRUE)
    stage_cases <- cc_u$last.stage >= 1L
    expected <- cts1$opt(fit1$outcome.fit, dat_u[stage_cases, , drop = FALSE],
                         quiet = TRUE) |> drop()
    expect_equal(res[[1L]]$opt.treat[stage_cases], expected)
  })
  
  test_that(".stageTx: NULL steps element is skipped without error", {
    steps_with_null <- make_steps()
    steps_with_null[[1L]] <- NULL
    expect_no_error(.stageTx(steps_with_null, cc_u, dat_u, quiet = TRUE))
  })
  
  test_that(".stageTx: stops when steps is not a list", {
    expect_error(.stageTx("not a list", cc_u, dat_u, quiet = TRUE), "`steps`")
  })
  
  test_that(".stageTx: stops when complete.case.info missing required names", {
    bad_cc <- list(last.stage = cc_u$last.stage)
    expect_error(.stageTx(make_steps(), bad_cc, dat_u, quiet = TRUE),
                 "`complete.case.info`")
  })
  
  test_that(".stageTx: stops when data is not a data.frame", {
    expect_error(.stageTx(make_steps(), cc_u, as.matrix(dat_u), quiet = TRUE),
                 "`data`")
  })
  
  test_that(".stageY: returns numeric vector of length n", {
    res <- .stageY(dat_u$Y, make_steps(), cc_u, dat_u, type = "DTR")
    expect_true(is.numeric(res))
    expect_length(res, n_u)
  })
  
  test_that(".stageY: oracle — when opt equals observed A, regret is 0 and Y unchanged", {
    # opt.treat == A for all participants means regret = 0, Y += 0
    steps_same <- make_steps(opt1 = dat_u$A1, opt2 = dat_u$A2)
    res <- .stageY(dat_u$Y, steps_same, cc_u, dat_u, type = "DTR")
    expect_equal(res, dat_u$Y)
  })
  
  test_that(".stageY: Y changes when opt differs from observed A", {
    steps_diff <- make_steps(opt1 = as.integer(!dat_u$A1),
                             opt2 = as.integer(!dat_u$A2))
    res <- .stageY(dat_u$Y, steps_diff, cc_u, dat_u, type = "DTR")
    expect_false(isTRUE(all.equal(res, dat_u$Y)))
  })
  
  test_that(".stageY: stops when Y is not numeric", {
    expect_error(
      .stageY(as.character(dat_u$Y), make_steps(), cc_u, dat_u, type = "DTR"),
      "`Y`"
    )
  })
  
  test_that(".stageY: stops when steps is not a list", {
    expect_error(
      .stageY(dat_u$Y, "not a list", cc_u, dat_u, type = "DTR"),
      "`steps`"
    )
  })
  
  test_that(".stageY: stops when data nrow != length(Y)", {
    expect_error(
      .stageY(dat_u$Y, make_steps(), cc_u, dat_u[1L:5L, ], type = "DTR"),
      "`data`"
    )
  })
  
  test_that(".stageY: stops when type is not a character", {
    expect_error(
      .stageY(dat_u$Y, make_steps(), cc_u, dat_u, type = 1L),
      "`type`"
    )
  })
  
  test_that(".stageYSurvival: returns numeric vector of length n", {
    times <- c("T1", "T2")
    res <- .stageYSurvival(times, make_steps(), cc_u, dat_u, type = "DTR")
    expect_true(is.numeric(res))
    expect_length(res, n_u)
  })
  
  test_that(".stageYSurvival: oracle — when opt equals observed A, Y = sum of times", {
    # shift = 0 when opt == A, exp(0) = 1, so Y[sc] += time[sc]
    # With all participants complete at both stages:
    # Y = T1 * exp(shift1) + T2 * exp(shift2) = T1 + T2
    times <- c("T1", "T2")
    steps_same <- make_steps(opt1 = dat_u$A1, opt2 = dat_u$A2)
    res <- .stageYSurvival(times, steps_same, cc_u, dat_u, type = "DTR")
    expected <- dat_u$T1 + dat_u$T2
    expect_equal(res, expected)
  })
  
  test_that(".stageYSurvival: stops when times is not a character vector", {
    expect_error(
      .stageYSurvival(1L:2L, make_steps(), cc_u, dat_u, type = "DTR"),
      "`times`"
    )
  })
  
  test_that(".stageYSurvival: stops when steps is not a list", {
    expect_error(
      .stageYSurvival(c("T1", "T2"), "not a list", cc_u, dat_u, type = "DTR"),
      "`steps`"
    )
  })
  
  test_that(".stageYSurvival: stops when type is not a length-1 character", {
    expect_error(
      .stageYSurvival(c("T1", "T2"), make_steps(), cc_u, dat_u, type = 1L),
      "`type`"
    )
  })
  
  test_that(".stageRegret: returns list of length K", {
    res <- .stageRegret(make_steps(), cc_u, dat_u, type = "DTR")
    expect_length(res, 2L)
  })
  
  test_that(".stageRegret: each element is numeric vector", {
    res <- .stageRegret(make_steps(), cc_u, dat_u, type = "DTR")
    expect_true(all(sapply(res, is.numeric)))
  })
  
  test_that(".stageRegret: oracle — when opt equals observed A, all regrets are zero", {
    steps_same <- make_steps(opt1 = dat_u$A1, opt2 = dat_u$A2)
    res <- .stageRegret(steps_same, cc_u, dat_u, type = "DTR")
    expect_equal(unname(res[[1L]]), rep(0.0, n_u))
    expect_equal(unname(res[[2L]]), rep(0.0, n_u))
  })
  
  test_that(".stageRegret: regret is non-zero when opt differs from observed A", {
    steps_diff <- make_steps(opt1 = as.integer(!dat_u$A1),
                             opt2 = as.integer(!dat_u$A2))
    res <- .stageRegret(steps_diff, cc_u, dat_u, type = "DTR")
    expect_false(all(res[[1L]] == 0.0))
  })
  
  test_that(".stageRegret: stops when steps is not a list", {
    expect_error(.stageRegret("not a list", cc_u, dat_u, type = "DTR"), "`steps`")
  })
  
  test_that(".stageRegret: stops when complete.case.info missing required names", {
    bad_cc <- list(last.stage = cc_u$last.stage)
    expect_error(.stageRegret(make_steps(), bad_cc, dat_u, type = "DTR"),
                 "`complete.case.info`")
  })
  
  test_that(".stageRegret: stops when data is not a data.frame", {
    expect_error(
      .stageRegret(make_steps(), cc_u, as.matrix(dat_u), type = "DTR"),
      "`data`"
    )
  })
  
  test_that(".stageRegret: stops when type is not a character", {
    expect_error(.stageRegret(make_steps(), cc_u, dat_u, type = 1L), "`type`")
  })
  
}) # end local
