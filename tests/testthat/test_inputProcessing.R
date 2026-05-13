# Two-stage model formulae
treat1 <- A1 ~ X1
treat2 <- A2 ~ X2
blip1 <- ~ X1
blip2 <- ~ X2
tf1 <- ~ X1
tf2 <- ~ X2

models_2stage <- list(
  list(blip = blip1, treat = treat1, tf = tf1),
  list(blip = blip2, treat = treat2, tf = tf2)
)

models_1stage <- list(
  list(blip = blip1, treat = treat1, tf = tf1)
)

n <- 50L
dat <- data.frame(
  X1 = rnorm(n),
  A1 = rbinom(n, 1L, 0.5),
  X2 = rnorm(n),
  A2 = rbinom(n, 1L, 0.5),
  Y  = rnorm(n)
)

test_that(".modelsTest: single formulae are coerced to lists", {
  res <- .modelsTest(blip1, treat1, tf1)
  expect_true(is.list(res))
  expect_length(res, 1L)
})

test_that(".modelsTest: output is grouped by stage with correct element names", {
  res <- .modelsTest(list(blip1, blip2), list(treat1, treat2), list(tf1, tf2))
  expect_length(res, 2L)
  expect_named(res[[1L]], c("blip", "treat", "tf"), ignore.order = TRUE)
  expect_named(res[[2L]], c("blip", "treat", "tf"), ignore.order = TRUE)
})

test_that(".modelsTest: formulae are assigned to the correct stage element", {
  res <- .modelsTest(list(blip1, blip2), list(treat1, treat2), list(tf1, tf2))
  expect_identical(res[[1L]]$blip, blip1)
  expect_identical(res[[1L]]$treat, treat1)
  expect_identical(res[[1L]]$tf, tf1)
  expect_identical(res[[2L]]$blip, blip2)
  expect_identical(res[[2L]]$treat, treat2)
  expect_identical(res[[2L]]$tf, tf2)
})

test_that(".modelsTest: errors as expected", {
  expect_error(.modelsTest(blip1, ~X1, tf1), "treat.mod")
  expect_error(.modelsTest(A1 ~ X1, treat1, tf1), "blip.mod")
  expect_error(.modelsTest(blip1, treat1, A1 ~ X1), "tf.mod")
  expect_error(.modelsTest(blip1, list("not a formula"), tf1), "treat.mod")
  expect_error(
    .modelsTest(list(blip1, blip2), list(treat1), list(tf1, tf2)),
    "same length"
  )
})


test_that(".constructMissingness: stage 1 auto-model is intercept-only", {
  res <- suppressMessages(.constructMissingness(models_2stage, NULL))
  expect_equal(deparse(res[[1L]]$cens), "~1")
})

test_that(".constructMissingness: stage 2 auto-model contains stage 1 covariates", {
  res <- suppressMessages(.constructMissingness(models_2stage, NULL))
  vars <- all.vars(res[[2L]]$cens)
  stage1_vars <- lapply(models_2stage[[1L]], all.vars) |> unlist() |> unique()
  expect_true(all(stage1_vars %in% vars))
})

test_that(".constructMissingness: K=1 with both NULL does not error", {
  expect_no_error(suppressMessages(.constructMissingness(models_1stage, NULL)))
})

test_that(".constructMissingness: K=1 auto-model is intercept-only", {
  res <- suppressMessages(.constructMissingness(models_1stage, NULL))
  expect_equal(deparse(res[[1L]]$cens), "~1")
})

test_that(".constructMissingness: adds cens element from user-provided list", {
  mm  <- list(~1, ~X1)
  res <- .constructMissingness(models_2stage, mm)
  expect_identical(res[[1L]]$cens, ~1)
  expect_identical(res[[2L]]$cens, ~X1)
})

test_that(".constructMissingness: single formula is coerced to a list", {
  res <- .constructMissingness(models_1stage, ~1)
  expect_identical(res[[1L]]$cens, ~1)
})

test_that(".constructMissingness: warns when cens already set and missing.mod provided", {
  models_with_cens <- models_2stage
  models_with_cens[[1L]]$cens <- ~1
  models_with_cens[[2L]]$cens <- ~X1
  expect_warning(
    .constructMissingness(models_with_cens, list(~1, ~X1)),
    "will be ignored"
  )
})

test_that(".constructMissingness: errors as expected", {
  expect_error(
    .constructMissingness(models_2stage, list(~1)),
    "appropriate length"
  )
  expect_error(
    .constructMissingness(models_2stage, list(~1, "not a formula")),
    "list of formulae"
  )
  expect_error(
    .constructMissingness(models_2stage, "not a formula or list"),
    "list of models"
  )
})


# Convenience wrapper — fills in the arguments we are not testing
tt <- function(weight = "ipw",
               treat.type = "bin",
               n.bins = NA_integer_,
               treat.wgt.man = NULL,
               treat.range = NULL,
               treat.fam = NULL,
               K = 2L,
               models = models_2stage,
               data = dat) {
  .treatmentTest(weight, treat.type, n.bins, treat.wgt.man,
                 treat.range, treat.fam, K, models, data)
}

test_that(".treatmentTest: returns a list with required elements", {
  res <- tt()
  expect_true(is.list(res))
  expect_true(all(c("n.bins", "tx.wgt.man", "dependent.vars",
                    "tx.range", "tx.family") %in% names(res)))
})

test_that(".treatmentTest: dependent.vars$treat contains treatment variable names", {
  res <- tt()
  expect_equal(res$dependent.vars$treat, c("A1", "A2"))
})

test_that(".treatmentTest: binary treatment sets tx.family to NA", {
  res <- tt(treat.type = "bin")
  expect_true(is.na(res$tx.family))
})

test_that(".treatmentTest: binary treatment sets tx.range to list of NAs", {
  res <- tt(treat.type = "bin")
  expect_true(all(is.na(unlist(res$tx.range))))
})

test_that(".treatmentTest: n.bins stored when weight is qpom and treat.type is cont", {
  res <- tt(weight = "qpom", treat.type = "cont",
            n.bins = 5L, treat.range = c(0, 1))
  expect_equal(res$n.bins, 5L)
})

test_that(".treatmentTest: n.bins stored when weight is overlap and treat.type is cont", {
  res <- tt(weight = "overlap", treat.type = "cont",
            n.bins = 5L, treat.range = c(0, 1))
  expect_equal(res$n.bins, 5L)
})

test_that(".treatmentTest: qpom with multinomial treatment does not require n.bins", {
  expect_no_error(tt(weight = "qpom", treat.type = "multi"))
})

test_that(".treatmentTest: stops when weight is qpom and treat.type is bin", {
  expect_error(tt(weight = "qpom", treat.type = "bin"), "binary treatments")
})

test_that(".treatmentTest: stops when weight is qpom and n.bins is not an integer", {
  expect_error(
    tt(weight = "qpom", treat.type = "cont", n.bins = 2.5, treat.range = c(0, 1)),
    "n.bins"
  )
})

test_that(".treatmentTest: manual weights stored correctly", {
  wgts <- list(rep(1.0, n), rep(1.0, n))
  res  <- tt(weight = "manual", treat.wgt.man = wgts)
  expect_identical(res$tx.wgt.man, wgts)
})

test_that(".treatmentTest: stops when weight is manual but no weights provided", {
  expect_error(tt(weight = "manual", treat.wgt.man = NULL), "must be provided")
})

test_that(".treatmentTest: stops when weights provided but weight is not manual", {
  expect_error(
    tt(weight = "ipw", treat.wgt.man = list(rep(1.0, n), rep(1.0, n))),
    "cannot be provided"
  )
})

test_that(".treatmentTest: stops when treat.wgt.man is wrong length", {
  expect_error(
    tt(weight = "manual", treat.wgt.man = list(rep(1.0, n))),
    "list of K numeric vectors"
  )
})

test_that(".treatmentTest: stops when treat.wgt.man elements are wrong length", {
  expect_error(
    tt(weight = "manual", treat.wgt.man = list(rep(1.0, n - 1L), rep(1.0, n))),
    "all participants"
  )
})

test_that(".treatmentTest: stops when treat.wgt.man contains non-positive values", {
  expect_error(
    tt(weight = "manual", treat.wgt.man = list(rep(-1.0, n), rep(1.0, n))),
    "must be positive"
  )
})

test_that(".treatmentTest: stops when treat.wgt.man is not NULL or a list", {
  expect_error(tt(weight = "manual", treat.wgt.man = 1.0), "NULL or a list")
})

test_that(".treatmentTest: scalar treat.range replicated to list of length K", {
  res <- tt(treat.type = "cont", treat.range = c(0, 1))
  expect_length(res$tx.range, 2L)
  expect_equal(res$tx.range[[1L]], c(0, 1))
  expect_equal(res$tx.range[[2L]], c(0, 1))
})

test_that(".treatmentTest: list treat.range stored as-is", {
  tr  <- list(c(0, 1), c(0.5, 2))
  res <- tt(treat.type = "cont", treat.range = tr)
  expect_identical(res$tx.range, tr)
})

test_that(".treatmentTest: NULL treat.range computed from data for cont treatment", {
  res <- suppressMessages(tt(treat.type = "cont", treat.range = NULL))
  expect_length(res$tx.range, 2L)
  expect_equal(res$tx.range[[1L]], range(dat$A1))
  expect_equal(res$tx.range[[2L]], range(dat$A2))
})

test_that(".treatmentTest: stops when treat.range is invalid for cont treatment", {
  expect_error(
    tt(treat.type = "cont", treat.range = c(0, 1, 2)),
    "treat.range"
  )
})

test_that(".treatmentTest: NULL treat.fam gives gaussian identity for cont", {
  res <- tt(treat.type = "cont", treat.range = c(0, 1), treat.fam = NULL)
  expect_equal(res$tx.family$family, "gaussian")
  expect_equal(res$tx.family$link,   "identity")
})

test_that(".treatmentTest: 'gaussian' character gives gaussian identity", {
  res <- tt(treat.type = "cont", treat.range = c(0, 1), treat.fam = "gaussian")
  expect_equal(res$tx.family$family, "gaussian")
  expect_equal(res$tx.family$link,   "identity")
})

test_that(".treatmentTest: 'Gamma' character gives Gamma log", {
  res <- tt(treat.type = "cont", treat.range = c(0, 1), treat.fam = "Gamma")
  expect_equal(res$tx.family$family, "Gamma")
  expect_equal(res$tx.family$link,   "log")
})

test_that(".treatmentTest: family object is stored directly", {
  fam <- gaussian(link = "identity")
  res <- tt(treat.type = "cont", treat.range = c(0, 1), treat.fam = fam)
  expect_equal(res$tx.family$family, "gaussian")
  expect_equal(res$tx.family$link,   "identity")
})

test_that(".treatmentTest: stops when treat.fam character is unrecognised", {
  expect_error(
    tt(treat.type = "cont", treat.range = c(0, 1), treat.fam = "poisson"),
    "gaussian.*Gamma|Gamma.*gaussian"
  )
})

test_that(".treatmentTest: stops when gaussian family has non-identity link", {
  expect_error(
    tt(treat.type = "cont", treat.range = c(0, 1),
       treat.fam = gaussian(link = "log")),
    "unsupported family"
  )
})

test_that(".treatmentTest: stops when Gamma family has non-log link", {
  expect_error(
    tt(treat.type = "cont", treat.range = c(0, 1),
       treat.fam = Gamma(link = "identity")),
    "unsupported family"
  )
})

test_that(".treatmentTest: stops when treat.fam is not NULL, character, or family", {
  expect_error(
    tt(treat.type = "cont", treat.range = c(0, 1), treat.fam = 42L),
    "treat.fam"
  )
})

dat_status <- data.frame(
  S01  = c(0L, 1L, 0L, 1L, 1L),   # already 0/1 integer
  S12  = c(1L, 2L, 1L, 2L, 1L),   # binary but coded 1/2
  Sna  = c(0L,  NA, 1L, 0L, 1L),  # has NA
  Smul = c(0L, 1L, 2L, 3L, 4L)    # more than 2 values
)

test_that(".statusTest: returns data unchanged when status is already 0/1", {
  res <- .statusTest(rep("S01", 2L), dat_status)
  expect_equal(res$S01, dat_status$S01)
})

test_that(".statusTest: recodes binary non-0/1 status to 0/1", {
  res <- .statusTest(rep("S12", 2L), dat_status)
  expect_true(all(res$S12 %in% c(0L, 1L)))
  expect_equal(as.integer(res$S12), c(0L, 1L, 0L, 1L, 0L))
})

test_that(".statusTest: NA status messages and returns data unchanged", {
  expect_message(
    res <- .statusTest(rep(NA_character_, 2L), dat_status),
    "not censored"
  )
  expect_identical(res, dat_status)
})

test_that(".statusTest: errors as expected", {
  expect_error(.statusTest(c("S01", "S12"), dat_status), "more than 1 status")
  expect_error(.statusTest(rep("Sna", 2L), dat_status), "missing status")
  expect_error(.statusTest(rep("Smul", 2L), dat_status), "binary")
})


n_bc <- 200L

test_that(".bootstrapControlsTest: empty list fills all defaults correctly", {
  res <- .bootstrapControlsTest(list(), n_bc)
  expect_equal(res$B, 100L)
  expect_equal(res$M, n_bc)
  expect_equal(res$type, "standard")
  expect_equal(res$truncate,  0.0)
  expect_false(res$verbose)
  expect_false(res$interrupt)
})

test_that(".bootstrapControlsTest: provided values override defaults", {
  res <- .bootstrapControlsTest(list(B = 500L, type = "empirical"), n_bc)
  expect_equal(res$B, 500L)
  expect_equal(res$type, "empirical")
  expect_equal(res$M, n_bc)   # default still applied
})

test_that(".bootstrapControlsTest: M reset to n when NULL", {
  res <- .bootstrapControlsTest(list(M = NULL), n_bc)
  expect_equal(res$M, n_bc)
})

test_that(".bootstrapControlsTest: M reset to n when zero", {
  res <- .bootstrapControlsTest(list(M = 0L), n_bc)
  expect_equal(res$M, n_bc)
})

test_that(".bootstrapControlsTest: errors as expected", {
  expect_error(.bootstrapControlsTest("not a list", n_bc), "named list")
  expect_error(.bootstrapControlsTest(list(Z = 1L), n_bc), "unrecognized")
  expect_error(.bootstrapControlsTest(list(B = -1L),  n_bc), "B")
  expect_error(.bootstrapControlsTest(list(B = 0L),   n_bc), "B")
  expect_error(.bootstrapControlsTest(list(B = 1.5),  n_bc), "B")
  expect_error(.bootstrapControlsTest(list(M = 1.5),  n_bc), "M")
  expect_error(.bootstrapControlsTest(list(type = "percentile"), n_bc), "type")
  expect_error(.bootstrapControlsTest(list(truncate = -0.1), n_bc), "truncate")
  expect_error(.bootstrapControlsTest(list(truncate =  0.6), n_bc), "truncate")
  expect_error(.bootstrapControlsTest(list(verbose = 1L), n_bc), "verbose")
  expect_error(.bootstrapControlsTest(list(interrupt = 1L), n_bc), "interrupt")
  expect_no_error(.bootstrapControlsTest(list(truncate = 0.0), n_bc))
  expect_no_error(.bootstrapControlsTest(list(truncate = 0.5), n_bc))
})


test_that(".varestimTest: returns NA when var.estim is 'none'", {
  expect_true(is.na(.varestimTest("none", "dwols", list(), n_bc)))
})

test_that(".varestimTest: returns NA when var.estim is 'sandwich'", {
  expect_true(is.na(.varestimTest("sandwich", "dwols", list(), n_bc)))
})

test_that(".varestimTest: returns bootstrap controls list when var.estim is 'bootstrap'", {
  res <- .varestimTest("bootstrap", "dwols", list(B = 200L), n_bc)
  expect_true(is.list(res))
  expect_equal(res$B, 200L)
})

test_that(".varestimTest: propagates stopping conditions from .bootstrapControlsTest", {
  expect_error(
    .varestimTest("bootstrap", "dwols", list(type = "invalid"), n_bc),
    "type"
  )
})
