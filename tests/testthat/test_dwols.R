n <- 30L
set.seed(42L)

dat <- data.frame(
  X = cos(seq(0.0, pi, length.out = n)),
  A = rep(0L:1L, times = n / 2L),
  Y = rnorm(n)
)

wgts <- rep(1.0, n)

models_bin <- list(blip = ~X, treat = A ~ X, tf = ~X)

cts_bin <- Binary$new(tf.model = models_bin$tf,
                      blip.model = models_bin$blip,
                      tx.var = "A")

# Continuous quadratic setup
dat_cont <- data.frame(
  X = cos(seq(0.0, pi, length.out = n)),
  A = seq(0.1, 0.9, length.out = n),
  Y = rnorm(n)
)

models_cont <- list(blip = ~A + X, treat = A ~ X, tf = ~X)

cts_cont <- ContQuadraticBlip$new(tf.model = models_cont$tf,
                                  blip.model = models_cont$blip,
                                  tx.var = "A",
                                  treat.range = c(0.1, 0.9))


local({
  res <- .dwols(Y = dat$Y,
                A = dat$A,
                data = dat,
                wgts = wgts,
                cts.obj = cts_bin,
                tx.var = "A")
  
  test_that(".dwols: returns list with required element names", {
    expect_named(res, c("outcome.fit", "beta", "psi"), ignore.order = TRUE)
  })
  
  test_that(".dwols: outcome.fit is an lm object", {
    expect_true(inherits(res$outcome.fit, "lm"))
  })
  
  test_that(".dwols: psi contains only blip-related terms (contains tx.var)", {
    expect_true(all(grepl("A", names(res$psi))))
  })
  
  test_that(".dwols: beta contains no blip-related terms (no tx.var)", {
    expect_false(any(grepl("^A", names(res$beta))))
  })
  
  test_that(".dwols: psi and beta together account for all coefficients", {
    all_coef_names <- names(coef(res$outcome.fit))
    expect_setequal(c(names(res$psi), names(res$beta)), all_coef_names)
  })
  
  test_that(".dwols: outcome.fit formula LHS is not 'Y' when Y column exists in data", {
    lhs <- as.character(formula(res$outcome.fit))[2L]
    expect_false(lhs == "Y")
  })
  
  test_that(".dwols: outcome.fit formula RHS contains tx.var", {
    rhs_terms <- attr(terms(formula(res$outcome.fit)), "term.labels")
    expect_true(any(grepl("^A", rhs_terms)))
  })
  
  test_that(".dwols: outcome.fit formula RHS contains tf covariate", {
    rhs_terms <- attr(terms(formula(res$outcome.fit)), "term.labels")
    expect_true("X" %in% rhs_terms)
  })
})

test_that(".dwols: handles collision when 'Y' already exists as column in data", {
  dat_collision <- dat
  # "Y" is already in colnames — function must use a different internal name
  expect_no_error(
    .dwols(Y = dat$Y, A = dat$A, data = dat_collision,
           wgts = wgts, cts.obj = cts_bin, tx.var = "A")
  )
})

local({
  res <- .dwols(Y = dat_cont$Y,
                A = dat_cont$A,
                data = dat_cont,
                wgts = wgts,
                cts.obj = cts_cont,
                tx.var = "A")
  
  test_that(".dwols: ContQuadraticBlip — psi contains both linear and quadratic tx terms", {
    expect_true("A" %in% names(res$psi))
    expect_true(length(res$psi) >= 2L)
  })
  
  test_that(".dwols: ContQuadraticBlip — squared tx column appears in outcome.fit terms", {
    fit_terms <- attr(terms(formula(res$outcome.fit)), "term.labels")
    sq_col <- setdiff(names(res$outcome.fit$model),
                      c("X", "A", as.character(formula(res$outcome.fit))[2L]))
    expect_true(length(sq_col) >= 1L)
  })
  
  test_that(".dwols: ContQuadraticBlip — beta contains no tx.var terms", {
    expect_false(any(grepl("^A", names(res$beta))))
  })
})

test_that(".dwols: errors as expected", {
  expect_error(
    .dwols(A = dat$A, data = dat, wgts = wgts, cts.obj = cts_bin, tx.var = "A"),
    "`Y`"
  )
  expect_error(
    .dwols(Y = as.character(dat$Y), A = dat$A, data = dat,
           wgts = wgts, cts.obj = cts_bin, tx.var = "A"),
    "`Y`"
  )
  expect_error(
    .dwols(Y = 1.0, A = dat$A[1L], data = dat[1L, ], wgts = wgts[1L],
           cts.obj = cts_bin, tx.var = "A"),
    "`Y`"
  )
  expect_error(
    .dwols(Y = dat$Y, A = dat$A[1L:5L], data = dat,
           wgts = wgts, cts.obj = cts_bin, tx.var = "A"),
    "`A`"
  )
  expect_error(
    .dwols(Y = dat$Y, A = dat$A, data = as.matrix(dat),
           wgts = wgts, cts.obj = cts_bin, tx.var = "A"),
    "`data`"
  )
  expect_error(
    .dwols(Y = dat$Y, A = dat$A, data = dat[1L:5L, ],
           wgts = wgts, cts.obj = cts_bin, tx.var = "A"),
    "`data`"
  )
  expect_error(
    .dwols(Y = dat$Y, A = dat$A, data = dat,
           wgts = as.character(wgts), cts.obj = cts_bin, tx.var = "A"),
    "`wgts`"
  )
  expect_error(
    .dwols(Y = dat$Y, A = dat$A, data = dat,
           wgts = wgts[1L:5L], cts.obj = cts_bin, tx.var = "A"),
    "`wgts`"
  )
  expect_error(
    .dwols(Y = dat$Y, A = dat$A, data = dat,
           wgts = wgts, cts.obj = list(), tx.var = "A"),
    "`cts.obj`"
  )
  expect_error(
    .dwols(Y = dat$Y, A = dat$A, data = dat,
           wgts = wgts, cts.obj = cts_bin, tx.var = 1L),
    "`tx.var`"
  )
  # Pass a cts.obj whose full.model references a column not in data
  cts_bad <- Binary$new(tf.model = ~Z, blip.model = ~X, tx.var = "A")
  expect_error(
    .dwols(Y = dat$Y, A = dat$A, data = dat,
           wgts = wgts, cts.obj = cts_bad, tx.var = "A"),
    "unable to complete DWOLS"
  )
})