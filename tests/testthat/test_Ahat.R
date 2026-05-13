n <- 50L
set.seed(42L)

dat_bin <- data.frame(
  X  = rnorm(n),
  A  = rbinom(n, 1L, 0.5),
  Y  = rnorm(n)
)

dat_multi <- data.frame(
  X = rnorm(n),
  A = factor(sample(c("low", "med", "high"), n, replace = TRUE),
             levels = c("low", "med", "high")),
  Y = rnorm(n)
)

dat_cont <- data.frame(
  X = cos(seq(0.0, pi, length.out = n)),
  A = seq(0.1, 0.9, length.out = n),
  Y = rnorm(n)
)

# Stage model lists
models_bin <- list(
  blip  = ~X,
  treat = A ~ X,
  tf    = ~X
)

models_multi <- list(
  blip  = ~X,
  treat = A ~ X,
  tf    = ~X
)

# Quadratic blip — tx.var appears in blip model
models_cont_q <- list(
  blip  = ~A + X,
  treat = A ~ X,
  tf    = ~X
)

# Linear blip — tx.var does NOT appear in blip model
models_cont_l <- list(
  blip  = ~X,
  treat = A ~ X,
  tf    = ~X
)


test_that(".Ahat_binary: returns list with required element names", {
  res <- .Ahat_binary(dat_bin$A, dat_bin, "A", models_bin)
  expect_named(res, c("A", "A.hat", "cts", "cts.obj", "tx.mod.fitted"),
               ignore.order = TRUE)
})

test_that(".Ahat_binary: cts is 'bin'", {
  res <- .Ahat_binary(dat_bin$A, dat_bin, "A", models_bin)
  expect_equal(res$cts, "bin")
})

test_that(".Ahat_binary: cts.obj is a Binary R6 object", {
  res <- .Ahat_binary(dat_bin$A, dat_bin, "A", models_bin)
  expect_true(inherits(res$cts.obj, "Binary"))
})

test_that(".Ahat_binary: tx.mod.fitted is a glm object", {
  res <- .Ahat_binary(dat_bin$A, dat_bin, "A", models_bin)
  expect_true(inherits(res$tx.mod.fitted, "glm"))
})

test_that(".Ahat_binary: treatment model formula contains tx.var on LHS", {
  res <- .Ahat_binary(dat_bin$A, dat_bin, "A", models_bin)
  expect_equal(as.character(formula(res$tx.mod.fitted))[2L], "A")
})

test_that(".Ahat_binary: treatment model formula contains covariate on RHS", {
  res <- .Ahat_binary(dat_bin$A, dat_bin, "A", models_bin)
  expect_true("X" %in% names(coef(res$tx.mod.fitted)))
})

test_that(".Ahat_binary: A is recoded to 0/1 integer", {
  res <- .Ahat_binary(dat_bin$A, dat_bin, "A", models_bin)
  expect_true(all(res$A %in% c(0L, 1L)))
  expect_true(is.integer(res$A))
})

test_that(".Ahat_binary: A recoding oracle — non-0/1 binary maps correctly", {
  dat_coded <- dat_bin
  dat_coded$A <- dat_coded$A + 5L  # code as 5/6 instead of 0/1
  res <- .Ahat_binary(dat_coded$A, dat_coded, "A", models_bin)
  expect_true(all(res$A %in% c(0L, 1L)))
})

test_that(".Ahat_binary: A.hat is numeric vector of length n", {
  res <- .Ahat_binary(dat_bin$A, dat_bin, "A", models_bin)
  expect_true(is.numeric(res$A.hat))
  expect_length(res$A.hat, n)
})

test_that(".Ahat_binary: A.hat values are in (0, 1)", {
  res <- .Ahat_binary(dat_bin$A, dat_bin, "A", models_bin)
  expect_true(all(res$A.hat > 0.0 & res$A.hat < 1.0))
})

test_that(".Ahat_binary: stops when A is missing", {
  expect_error(.Ahat_binary(data = dat_bin, tx.var = "A", models = models_bin),
               "`A`")
})

test_that(".Ahat_binary: stops when A has length 1", {
  expect_error(.Ahat_binary(1L, dat_bin[1L, ], "A", models_bin), "`A`")
})

test_that(".Ahat_binary: stops when data is not a data.frame", {
  expect_error(.Ahat_binary(dat_bin$A, as.matrix(dat_bin), "A", models_bin),
               "`data`")
})

test_that(".Ahat_binary: stops when data nrow != length(A)", {
  expect_error(.Ahat_binary(dat_bin$A, dat_bin[1L:10L, ], "A", models_bin),
               "`data`")
})

test_that(".Ahat_binary: stops when tx.var is not a character", {
  expect_error(.Ahat_binary(dat_bin$A, dat_bin, 1L, models_bin), "`tx.var`")
})

test_that(".Ahat_binary: stops when models has no 'treat' element", {
  expect_error(.Ahat_binary(dat_bin$A, dat_bin, "A", list(blip = ~X)),
               "`models`")
})

test_that(".Ahat_binary: stops when treatment model formula references absent column", {
  bad_models <- list(blip = ~X, treat = A ~ Z, tf = ~X)
  expect_error(.Ahat_binary(dat_bin$A, dat_bin, "A", bad_models),
               "unable to fit treatment model")
})


test_that(".Ahat_multinom: returns list with required element names", {
  res <- suppressMessages(
    .Ahat_multinom(dat_multi$A, dat_multi, "A", models_multi)
  )
  expect_named(res, c("A", "A.hat", "cts", "cts.obj", "tx.mod.fitted"),
               ignore.order = TRUE)
})

test_that(".Ahat_multinom: cts is 'multinom'", {
  res <- suppressMessages(
    .Ahat_multinom(dat_multi$A, dat_multi, "A", models_multi)
  )
  expect_equal(res$cts, "multinom")
})

test_that(".Ahat_multinom: cts.obj is a MultiNom R6 object", {
  res <- suppressMessages(
    .Ahat_multinom(dat_multi$A, dat_multi, "A", models_multi)
  )
  expect_true(inherits(res$cts.obj, "MultiNom"))
})

test_that(".Ahat_multinom: treatment model formula contains tx.var on LHS", {
  res <- suppressMessages(
    .Ahat_multinom(dat_multi$A, dat_multi, "A", models_multi)
  )
  expect_equal(as.character(formula(res$tx.mod.fitted))[2L], "A")
})

test_that(".Ahat_multinom: A.hat is a matrix with n rows and n_levels columns", {
  res <- suppressMessages(
    .Ahat_multinom(dat_multi$A, dat_multi, "A", models_multi)
  )
  expect_true(is.matrix(res$A.hat))
  expect_equal(nrow(res$A.hat), n)
  expect_equal(ncol(res$A.hat), 3L)
})

test_that(".Ahat_multinom: A.hat column names match treatment levels", {
  res <- suppressMessages(
    .Ahat_multinom(dat_multi$A, dat_multi, "A", models_multi)
  )
  expect_equal(colnames(res$A.hat), levels(dat_multi$A))
})

test_that(".Ahat_multinom: A.hat rows sum to 1", {
  res <- suppressMessages(
    .Ahat_multinom(dat_multi$A, dat_multi, "A", models_multi)
  )
  expect_equal(unname(rowSums(res$A.hat)), rep(1.0, n), tolerance = 1e-10)
})

test_that(".Ahat_multinom: stops when A has fewer than 3 levels", {
  A_bin <- factor(rbinom(n, 1L, 0.5))
  expect_error(.Ahat_multinom(A_bin, dat_bin, "A", models_multi),
               "not multinomial")
})

test_that(".Ahat_multinom: stops when data nrow != length(A)", {
  expect_error(
    suppressMessages(
      .Ahat_multinom(dat_multi$A, dat_multi[1L:10L, ], "A", models_multi)
    ),
    "`data`"
  )
})

test_that(".Ahat_multinom: stops when tx.var is not a character", {
  expect_error(
    suppressMessages(.Ahat_multinom(dat_multi$A, dat_multi, 1L, models_multi)),
    "`tx.var`"
  )
})

test_that(".Ahat_multinom: stops when models missing required elements", {
  expect_error(
    suppressMessages(
      .Ahat_multinom(dat_multi$A, dat_multi, "A", list(blip = ~X))
    ),
    "`models`"
  )
})

test_that(".Ahat_cont: returns list with required element names", {
  res <- .Ahat_cont(dat_cont$A, dat_cont, "A", models_cont_q,
                    tx.range = c(0.1, 0.9), tx.family = gaussian())
  expect_named(res, c("A", "A.hat", "cts", "cts.obj", "tx.mod.fitted"),
               ignore.order = TRUE)
})

test_that(".Ahat_cont: cts is 'cts.q' when tx.var appears in blip model", {
  res <- .Ahat_cont(dat_cont$A, dat_cont, "A", models_cont_q,
                    tx.range = c(0.1, 0.9), tx.family = gaussian())
  expect_equal(res$cts, "cts.q")
})

test_that(".Ahat_cont: cts.obj is ContQuadraticBlip when tx.var in blip model", {
  res <- .Ahat_cont(dat_cont$A, dat_cont, "A", models_cont_q,
                    tx.range = c(0.1, 0.9), tx.family = gaussian())
  expect_true(inherits(res$cts.obj, "ContQuadraticBlip"))
})

test_that(".Ahat_cont: tx.mod.fitted is a glm object", {
  res <- .Ahat_cont(dat_cont$A, dat_cont, "A", models_cont_q,
                    tx.range = c(0.1, 0.9), tx.family = gaussian())
  expect_true(inherits(res$tx.mod.fitted, "glm"))
})

test_that(".Ahat_cont: treatment model LHS is tx.var", {
  res <- .Ahat_cont(dat_cont$A, dat_cont, "A", models_cont_q,
                    tx.range = c(0.1, 0.9), tx.family = gaussian())
  expect_equal(as.character(formula(res$tx.mod.fitted))[2L], "A")
})

test_that(".Ahat_cont: treatment model RHS contains covariate", {
  res <- .Ahat_cont(dat_cont$A, dat_cont, "A", models_cont_q,
                    tx.range = c(0.1, 0.9), tx.family = gaussian())
  expect_true("X" %in% names(coef(res$tx.mod.fitted)))
})

test_that(".Ahat_cont: A.hat is numeric vector of length n", {
  res <- .Ahat_cont(dat_cont$A, dat_cont, "A", models_cont_q,
                    tx.range = c(0.1, 0.9), tx.family = gaussian())
  expect_true(is.numeric(res$A.hat))
  expect_length(res$A.hat, n)
})

test_that(".Ahat_cont: warns when predicted treatments fall outside tx.range", {
  expect_warning(
    .Ahat_cont(dat_cont$A, dat_cont, "A", models_cont_q,
               tx.range = c(0.5, 0.51), tx.family = gaussian()),
    "outside of allowed"
  )
})

test_that(".Ahat_cont: stops when tx.var not in blip model (linear blip not supported)", {
  expect_error(
    .Ahat_cont(dat_cont$A, dat_cont, "A", models_cont_l,
               tx.range = c(0.1, 0.9), tx.family = gaussian()),
    "not supported"
  )
})

test_that(".Ahat_cont: stops when A is not a vector", {
  expect_error(
    .Ahat_cont(as.list(dat_cont$A), dat_cont, "A", models_cont_q,
               tx.range = c(0.1, 0.9), tx.family = gaussian()),
    "`A`"
  )
})

test_that(".Ahat_cont: stops when tx.range is wrong length", {
  expect_error(
    .Ahat_cont(dat_cont$A, dat_cont, "A", models_cont_q,
               tx.range = c(0.1, 0.5, 0.9), tx.family = gaussian()),
    "tx.range"
  )
})

test_that(".Ahat_cont: stops when tx.family is not gaussian or Gamma", {
  expect_error(
    .Ahat_cont(dat_cont$A, dat_cont, "A", models_cont_q,
               tx.range = c(0.1, 0.9), tx.family = "poisson"),
    "tx.family"
  )
})

test_that(".Ahat_cont: stops when treatment model references absent column", {
  bad_models <- list(blip = ~A + X, treat = A ~ Z, tf = ~X)
  expect_error(
    .Ahat_cont(dat_cont$A, dat_cont, "A", bad_models,
               tx.range = c(0.1, 0.9), tx.family = gaussian()),
    "unable to fit treatment model"
  )
})

pt <- function(stage.models = models_bin,
               data         = dat_bin,
               tx.var       = "A",
               tx.range     = NA_real_,
               tx.type      = "bin",
               tx.family    = gaussian()) {
  .processTreatment(stage.models, data, tx.var, tx.range, tx.type, tx.family)
}

test_that(".processTreatment: dispatches to binary path — cts is 'bin'", {
  res <- pt()
  expect_equal(res$cts, "bin")
})

test_that(".processTreatment: dispatches to multinomial path — cts is 'multinom'", {
  res <- suppressMessages(
    pt(stage.models = models_multi, data = dat_multi, tx.type = "multi")
  )
  expect_equal(res$cts, "multinom")
})

test_that(".processTreatment: dispatches to continuous path — cts is 'cts.q'", {
  res <- pt(stage.models = models_cont_q, data = dat_cont,
            tx.range = c(0.1, 0.9), tx.type = "cont")
  expect_equal(res$cts, "cts.q")
})

test_that(".processTreatment: errors as expected", {
  expect_error(pt(stage.models = list(treat = "A ~ X", blip = ~X, tf = ~X)),
               "stage.models")
  expect_error(pt(data = as.matrix(dat_bin)), "`data`")
  expect_error(pt(tx.var = c("A", "B")), "tx.var")
  expect_error(pt(tx.var = 1L),          "tx.var")
  expect_error(pt(tx.type = "ordered"), "tx.type")
  expect_error(pt(tx.var = "Z"), "treatment variable")
  
  dat_one       <- dat_bin
  dat_one$A     <- 0L
  expect_error(pt(data = dat_one), "binary")

  dat_bin2 <- dat_bin
  dat_bin2$A <- rbinom(n, 1L, 0.5)
  expect_error(
    pt(stage.models = models_multi, data = dat_bin2, tx.type = "multi"),
    "multinomial"
  )

  dat_many   <- dat_bin
  dat_many$A <- factor(sample(letters[1L:12L], n, replace = TRUE))
  models_many <- list(blip = ~X, treat = A ~ X, tf = ~X)
  expect_warning(
    suppressMessages(
      pt(stage.models = models_many, data = dat_many, tx.type = "multi")
    ),
    "large number"
  )
})