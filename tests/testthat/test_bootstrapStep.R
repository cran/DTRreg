# Shared minimal data


data(twoStageCont)
n <- nrow(twoStageCont)

withr::with_seed(42L, {
  n_s <- 30L
  dat_s <- data.frame(
    X = cos(seq(0.0, pi, length.out = n_s)),
    A = as.numeric(rbinom(n_s, 1L, 0.5)),
    Y = rnorm(n_s),
    T = abs(rnorm(n_s)) + 0.5
  )
})

# Trusted Layer 0 / Layer 4 objects
cts_b <- Binary$new(tf.model = ~X, blip.model = ~X, tx.var = "A")
fit_dw <- .dwols(Y = dat_s$Y, A = dat_s$A, data = dat_s,
                 wgts = rep(1.0, n_s), cts.obj = cts_b, tx.var = "A")

# Minimal step.obj as .getStage would produce
step_obj_s <- list(
  n = n_s,
  tx.var = "A",
  A = dat_s$A,
  psi = fit_dw$psi,
  cts = "bin",
  cts.obj = cts_b,
  outcome.fit = fit_dw$outcome.fit,
  models = list(blip = ~X, treat = A ~ X, tf = ~X, cens = ~1)
)

# Minimal obj for single-stage
make_ss_obj <- function(isSurvival = FALSE) {
  obj <- list(
    K = 1L,
    outcome = dat_s$Y,
    data = dat_s,
    models = list(list(blip = ~X, treat = A ~ X, tf = ~X, cens = ~1)),
    dependent.vars = list(treat = "A",
                          status = if (isSurvival) "delta" else NA_character_,
                          time = if (isSurvival) "T" else NA_character_),
    tx.range = list(NA_real_),
    method = "dwols",
    var.estim = "none",
    type = "DTR",
    tx.weight = "ipw",
    censoring.modeled = FALSE,
    tx.wgt.man = NA,
    tx.type = "bin",
    n.bins = NA_integer_,
    tx.family = NA,
    boot.controls = list(B = 5L, M = n_s, type = "standard",
                         truncate = 0.0, verbose = FALSE, interrupt = FALSE),
    full.cov = FALSE
  )
  if (isSurvival) {
    obj$manual.censor.weight <- FALSE
    obj$outcome = NULL
  }
  obj
}

# Full obj for .bootstrap() smoke tests — built from trusted .dtrProcedure
make_boot_obj <- function(B = 5L,
                          truncate = 0.0,
                          method = "dwols") {
  obj <- list(
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
    var.estim = "bootstrap",
    type = "DTR",
    tx.weight = "ipw",
    censoring.modeled = FALSE,
    tx.wgt.man = NA,
    tx.type = "bin",
    n.bins = NA_integer_,
    tx.family = NA,
    full.cov = FALSE,
    boot.controls = list(B = B, M = n, type = "standard",
                         truncate = truncate, verbose = FALSE,
                         interrupt = FALSE)
  )
  result <- suppressMessages(.dtrProcedure(obj, quiet = TRUE))
  obj$stages <- result$stages
  obj$last.stage <- result$last.stage
  obj$prob.complete.case <- result$prob.complete.case
  obj$d.hat <- result$d.hat
  obj$cens.mod.fitted <- result$cens.mod.fitted
  obj$opt.Y <- result$opt.Y
  obj$treat.vars <- result$treat.vars
  obj$regret <- result$regret
  obj
}



# .resampleObj


local({
  boot_sample <- sample(seq_len(n_s), n_s, replace = TRUE)
  obj <- make_ss_obj()
  
  test_that(".resampleObj: data is subsampled to boot_sample rows", {
    res <- .resampleObj(obj, boot_sample, isSurvival = FALSE)
    expect_equal(res$data, obj$data[boot_sample, , drop = FALSE])
  })
  
  test_that(".resampleObj: outcome is subsampled for non-survival", {
    res <- .resampleObj(obj, boot_sample, isSurvival = FALSE)
    expect_equal(res$outcome, obj$outcome[boot_sample])
  })
  
  test_that(".resampleObj: outcome not modified for survival", {
    obj_surv <- make_ss_obj(isSurvival = TRUE)
    res <- .resampleObj(obj_surv, boot_sample, isSurvival = TRUE)
    expect_null(res$outcome)
  })
  
  test_that(".resampleObj: var.estim set to 'none'", {
    obj$var.estim <- "bootstrap"
    res <- .resampleObj(obj, boot_sample, isSurvival = FALSE)
    expect_equal(res$var.estim, "none")
  })
  
  test_that(".resampleObj: boot.controls reset to empty list", {
    res <- .resampleObj(obj, boot_sample, isSurvival = FALSE)
    expect_identical(res$boot.controls, list())
  })
  
  test_that(".resampleObj: list tx.wgt.man is subsampled per stage", {
    obj$tx.wgt.man <- list(seq_len(n_s) / n_s, rev(seq_len(n_s) / n_s))
    res <- .resampleObj(obj, boot_sample, isSurvival = FALSE)
    expect_equal(res$tx.wgt.man[[1L]], obj$tx.wgt.man[[1L]][boot_sample])
    expect_equal(res$tx.wgt.man[[2L]], obj$tx.wgt.man[[2L]][boot_sample])
  })
  
  test_that(".resampleObj: non-list tx.wgt.man is left unchanged", {
    res <- .resampleObj(obj, boot_sample, isSurvival = FALSE)
    expect_identical(res$tx.wgt.man, obj$tx.wgt.man)
  })
})



# .trunc


test_that(".trunc: oracle — values clamped to [quantile(x,p), quantile(x,1-p)]", {
  x <- c(-3.0, -1.0, 0.0, 1.0, 3.0, 5.0, 10.0)
  p <- 0.1
  expected <- pmax(pmin(x, quantile(x, 1.0 - p)), quantile(x, p))
  expect_equal(.trunc(x, p), expected)
})

test_that(".trunc: p=0 returns x unchanged", {
  x <- rnorm(20L)
  expect_equal(.trunc(x, 0.0), x)
})

test_that(".trunc: output length equals input length", {
  x <- rnorm(20L)
  expect_length(.trunc(x, 0.1), length(x))
})

test_that(".trunc: max of truncated <= max of original", {
  x <- rnorm(50L)
  expect_true(max(.trunc(x, 0.1)) <= max(x))
})

test_that(".trunc: min of truncated >= min of original", {
  x <- rnorm(50L)
  expect_true(min(.trunc(x, 0.1)) >= min(x))
})



# .stackPsiBoot


local({
  # Hand-crafted: B=3 iterations, K=2 stages, 2 psi per stage
  psi_iter <- function(seed) {
    withr::with_seed(seed,
    list(c(A1 = rnorm(1L), `A1:X` = rnorm(1L)),
         c(A2 = rnorm(1L), `A2:X` = rnorm(1L))))
  }
  psi_boot <- list(psi_iter(1L), psi_iter(2L), psi_iter(3L))
  res <- .stackPsiBoot(psi_boot)
  
  test_that(".stackPsiBoot: returns list of length K", {
    expect_length(res, 2L)
  })
  
  test_that(".stackPsiBoot: each element is a matrix with B rows", {
    expect_equal(nrow(res[[1L]]), 3L)
    expect_equal(nrow(res[[2L]]), 3L)
  })
  
  test_that(".stackPsiBoot: each element has n_psi columns", {
    expect_equal(ncol(res[[1L]]), 2L)
    expect_equal(ncol(res[[2L]]), 2L)
  })
  
  test_that(".stackPsiBoot: oracle — row i of stage k matches psi_boot[[i]][[k]]", {
    for (i in seq_along(psi_boot)) {
      expect_equal(unname(res[[1L]][i, ]), unname(psi_boot[[i]][[1L]]))
      expect_equal(unname(res[[2L]][i, ]), unname(psi_boot[[i]][[2L]]))
    }
  })
})



# .truncPsiBoot


local({
  # Known psi_boot — 10 rows (B), 2 cols (psi) per stage
  withr::with_seed(42L,
  pb <- list(
    matrix(rnorm(20L), nrow = 10L, ncol = 2L),
    matrix(rnorm(20L), nrow = 10L, ncol = 2L)
  ))
  
  test_that(".truncPsiBoot: truncate=0 returns psi_boot unchanged", {
    expect_identical(.truncPsiBoot(0.0, pb), pb)
  })
  
  test_that(".truncPsiBoot: truncate>0 narrows range of each column", {
    res <- .truncPsiBoot(0.1, pb)
    for (k in seq_along(pb)) {
      for (j in seq_len(ncol(pb[[k]]))) {
        expect_true(diff(range(res[[k]][, j])) <= diff(range(pb[[k]][, j])))
      }
    }
  })
  
  test_that(".truncPsiBoot: oracle — matches column-wise .trunc application", {
    p <- 0.1
    res <- .truncPsiBoot(p, pb)
    for (k in seq_along(pb)) {
      expected <- apply(pb[[k]], 2L, .trunc, p)
      expect_equal(res[[k]], expected)
    }
  })
  
  test_that(".truncPsiBoot: dimensions preserved after truncation", {
    res <- .truncPsiBoot(0.1, pb)
    for (k in seq_along(pb)) {
      expect_equal(dim(res[[k]]), dim(pb[[k]]))
    }
  })
})



# .computeCovmat


local({
  withr::with_seed(42L,
  pb <- list(
    matrix(rnorm(20L), nrow = 10L, ncol = 2L),
    matrix(rnorm(20L), nrow = 10L, ncol = 2L)
  ))
  
  test_that(".computeCovmat: returns list of length K", {
    res <- .computeCovmat(pb)
    expect_length(res, 2L)
  })
  
  test_that(".computeCovmat: each element is a square matrix", {
    res <- .computeCovmat(pb)
    for (k in seq_along(pb)) {
      expect_true(is.matrix(res[[k]]))
      expect_equal(nrow(res[[k]]), ncol(res[[k]]))
    }
  })
  
  test_that(".computeCovmat: oracle — matches lapply(psi_boot, var)", {
    expect_equal(.computeCovmat(pb), lapply(pb, var))
  })
})



# .computeHoldY


local({
  dat_pred <- cts_b$prep(dat_s, dat_s$A)
  
  test_that(".computeHoldY: non-survival oracle — equals predict(outcome.fit, data)", {
    expected <- drop(predict(fit_dw$outcome.fit, dat_pred))
    expect_equal(.computeHoldY(step_obj_s, dat_pred, isSurvival = FALSE),
                 expected)
  })
  
  test_that(".computeHoldY: survival oracle — equals exp(predict(outcome.fit, data))", {
    expected <- exp(drop(predict(fit_dw$outcome.fit, dat_pred)))
    expect_equal(.computeHoldY(step_obj_s, dat_pred, isSurvival = TRUE),
                 expected)
  })
  
  test_that(".computeHoldY: returns numeric vector of length n", {
    res <- .computeHoldY(step_obj_s, dat_pred, isSurvival = FALSE)
    expect_true(is.numeric(res))
    expect_length(res, n_s)
  })
  
  test_that(".computeHoldY: survival values are strictly positive", {
    res <- .computeHoldY(step_obj_s, dat_pred, isSurvival = TRUE)
    expect_true(all(res > 0.0))
  })
})



# .buildSingleStageObj


local({
  obj <- make_ss_obj()
  dat_ss <- dat_s
  res_ns <- .buildSingleStageObj(obj, step_obj_s, "A", NA_character_,
                                 NA_character_, dat_ss, NULL, isSurvival = FALSE)
  
  test_that(".buildSingleStageObj: K is set to 1", {
    expect_equal(res_ns$K, 1L)
  })
  
  test_that(".buildSingleStageObj: var.estim is 'none'", {
    expect_equal(res_ns$var.estim, "none")
  })
  
  test_that(".buildSingleStageObj: boot.controls is empty list", {
    expect_identical(res_ns$boot.controls, list())
  })
  
  test_that(".buildSingleStageObj: models is list of one element", {
    expect_length(res_ns$models, 1L)
  })
  
  test_that(".buildSingleStageObj: dependent.vars$treat is tx.var", {
    expect_equal(res_ns$dependent.vars$treat, "A")
  })
  
  test_that(".buildSingleStageObj: data is set to provided data", {
    expect_identical(res_ns$data, dat_ss)
  })
  
  test_that(".buildSingleStageObj: NULL tx.wgt.man stored as NULL list element", {
    expect_null(res_ns$tx.wgt.man[[1L]])
  })
  
  test_that(".buildSingleStageObj: non-NULL tx.wgt.man stored directly", {
    wgt <- runif(n_s)
    res <- .buildSingleStageObj(obj, step_obj_s, "A", NA_character_,
                                NA_character_, dat_ss, wgt, isSurvival = FALSE)
    expect_equal(res$tx.wgt.man, wgt)
  })
  
  test_that(".buildSingleStageObj: survival sets dependent.vars$status and $time", {
    obj_surv <- make_ss_obj(isSurvival = TRUE)
    res_surv <- .buildSingleStageObj(obj_surv, step_obj_s, "A", "delta",
                                     "T", dat_ss, NULL, isSurvival = TRUE)
    expect_equal(res_surv$dependent.vars$status, "delta")
    expect_equal(res_surv$dependent.vars$time, "T")
  })
  
  test_that(".buildSingleStageObj: non-survival does not overwrite status/time", {
    expect_equal(res_ns$dependent.vars$status, NA_character_)
    expect_equal(res_ns$dependent.vars$time, NA_character_)
  })
})



# .sampleResiduals


local({
  residuals <- rnorm(n_s)
  reshist <- hist(residuals, plot = FALSE)
  mean.res <- mean(residuals)
  sd.res <- sd(residuals)
  
  test_that(".sampleResiduals: normal returns numeric vector of length n", {
    res <- .sampleResiduals("normal", reshist, mean.res, sd.res,
                            n = n_s, n_sample = n_s)
    expect_true(is.numeric(res))
    expect_length(res, n_s)
  })
  
  test_that(".sampleResiduals: empirical returns numeric vector of length n", {
    res <- .sampleResiduals("empirical", reshist, mean.res, sd.res,
                            n = n_s, n_sample = n_s)
    expect_true(is.numeric(res))
    expect_length(res, n_s)
  })
  
  test_that(".sampleResiduals: empirical values fall within histogram breaks range", {
    res <- .sampleResiduals("empirical", reshist, mean.res, sd.res,
                            n = n_s, n_sample = n_s)
    expect_true(all(res >= min(reshist$breaks)))
    expect_true(all(res <= max(reshist$breaks)))
  })
  
  test_that(".sampleResiduals: normal produces different results with different seeds", {
    r1 <- withr::with_seed(1L, .sampleResiduals("normal", reshist, mean.res, sd.res, n_s, n_s))
    r2 <- withr::with_seed(2L, .sampleResiduals("normal", reshist, mean.res, sd.res, n_s, n_s))
    expect_false(isTRUE(all.equal(r1, r2)))
  })
})



# .reportProgress
#
# readline() and proc.time() are untestable directly. We test all logic
# branches that don't depend on those — early returns, last initialisation,
# message suppression when eta is small, and the attr(last) update path.
# We mock ptm by setting it far in the past to force large eta values.


local({
  # ptm set 1000 seconds in the past so eta is always large
  ptm_past <- proc.time() - c(0, 0, 1000, 0, 0)
  
  test_that(".reportProgress: returns cont unchanged when verbose=FALSE", {
    res <- .reportProgress(proc.time(), i = 1L, B = 100L, cont = "y",
                           verbose = FALSE, interrupt = FALSE)
    expect_equal(res, "y")
  })
  
  test_that(".reportProgress: returns cont unchanged when i < 10", {
    res <- .reportProgress(proc.time(), i = 5L, B = 100L, cont = "y",
                           verbose = TRUE, interrupt = FALSE)
    expect_equal(res, "y")
  })
  
  test_that(".reportProgress: initialises last attribute when i == 10", {
    res <- .reportProgress(ptm_past, i = 10L, B = 100L, cont = "y",
                           verbose = TRUE, interrupt = FALSE)
    expect_false(is.null(attr(res, "last")))
  })
  
  test_that(".reportProgress: last attribute > 30 when set at i=10", {
    res <- .reportProgress(ptm_past, i = 10L, B = 100L, cont = "y",
                           verbose = TRUE, interrupt = FALSE)
    expect_true(attr(res, "last") > 30.0)
  })
  
  test_that(".reportProgress: issues message when eta > 30 and within last threshold", {
    cont_msg <- "y"
    attr(cont_msg, "last") <- 1e8
    expect_message(
      .reportProgress(ptm_past, i = 15L, B = 100L, cont = cont_msg,
                      verbose = TRUE, interrupt = FALSE),
      "seconds remaining"
    )
  })
  
  test_that(".reportProgress: no message when eta < 30", {
    expect_no_message(
      .reportProgress(proc.time(), i = 15L, B = 100L, cont = "y",
                      verbose = TRUE, interrupt = FALSE)
    )
  })
  
  test_that(".reportProgress: updates last attribute after issuing message", {
    cont_msg <- "y"
    attr(cont_msg, "last") <- 1e8
    res <- suppressMessages(
      .reportProgress(ptm_past, i = 15L, B = 100L, cont = cont_msg,
                      verbose = TRUE, interrupt = FALSE)
    )
    expect_true(attr(res, "last") < 1e8)
  })
})



# .bootstrap — integration smoke test


local({
  withr::with_seed(42L, {
    obj <- make_boot_obj(B = 5L)
    res <- suppressMessages(.bootstrap(obj, isSurvival = FALSE))
  })
  
  test_that(".bootstrap: returns list with covmat and psi.boot", {
    expect_named(res, c("covmat", "psi.boot"), ignore.order = TRUE)
  })
  
  test_that(".bootstrap: covmat is list of length K", {
    expect_true(is.list(res$covmat))
    expect_length(res$covmat, obj$K)
  })
  
  test_that(".bootstrap: each covmat element is a square matrix", {
    for (k in seq_len(obj$K)) {
      expect_true(is.matrix(res$covmat[[k]]))
      expect_equal(nrow(res$covmat[[k]]), ncol(res$covmat[[k]]))
    }
  })
  
  test_that(".bootstrap: covmat dimensions match psi length at each stage", {
    for (k in seq_len(obj$K)) {
      n_psi <- length(obj$stages[[k]]$psi)
      expect_equal(dim(res$covmat[[k]]), c(n_psi, n_psi))
    }
  })
  
  test_that(".bootstrap: psi.boot is list of length K", {
    expect_length(res$psi.boot, obj$K)
  })
  
  test_that(".bootstrap: each psi.boot element is a matrix with B rows", {
    for (k in seq_len(obj$K)) {
      expect_true(is.matrix(res$psi.boot[[k]]))
      expect_equal(nrow(res$psi.boot[[k]]), obj$boot.controls$B)
    }
  })
  
  test_that(".bootstrap: covmat oracle — matches var(psi.boot) per stage", {
    for (k in seq_len(obj$K)) {
      expect_equal(res$covmat[[k]], var(res$psi.boot[[k]]))
    }
  })
})

# --- truncation --------------------------------------------------------------

test_that(".bootstrap: truncation narrows psi.boot range vs no truncation", {
  withr::with_seed(42L,
                   res_plain <- suppressMessages(
                     .bootstrap(make_boot_obj(B = 20L, truncate = 0.0), isSurvival = FALSE)
                   )
  )
  withr::with_seed(42L,
                   res_trunc <- suppressMessages(
                     .bootstrap(make_boot_obj(B = 20L, truncate = 0.1), isSurvival = FALSE)
                   )
  )
  for (k in seq_len(2L)) {
    range_plain <- apply(res_plain$psi.boot[[k]], 2L, function(x) diff(range(x)))
    range_trunc <- apply(res_trunc$psi.boot[[k]], 2L, function(x) diff(range(x)))
    expect_true(all(range_trunc <= range_plain))
  }
})

# --- manual weights ----------------------------------------------------------

test_that(".bootstrap: manual tx.wgt.man is handled without error", {
  withr::with_seed(42L, {
    obj_man <- make_boot_obj(B = 5L)
    obj_man$tx.weight <- "manual"
    obj_man$tx.wgt.man <- list(rep(1.0, n), rep(1.0, n))
    expect_no_error(suppressMessages(.bootstrap(obj_man, isSurvival = FALSE)))
  })
})
