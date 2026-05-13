# Shared setup

data(twoStageCont)
n <- nrow(twoStageCont)

make_obj <- function(method = "dwols",
                     var.estim = "sandwich",
                     full.cov = FALSE) {
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
    type = "DTR",
    tx.weight = "ipw",
    censoring.modeled = FALSE,
    tx.wgt.man = NA,
    tx.type = "bin",
    n.bins = NA_integer_,
    tx.family = NA,
    boot.controls = list(B = 5L, M = n, type = "standard",
                         truncate = 0.0, verbose = FALSE,
                         interrupt = FALSE),
    full.cov = full.cov
  )
}

# Run .dtrProcedure() once - reused across all tests
base_obj <- make_obj()
dtr_result <- suppressMessages(.dtrProcedure(base_obj, quiet = TRUE))

# Attach required fields from dtr_result back to base_obj for .computeVariance()
enrich_obj <- function(obj, result) {
  obj$stages <- result$stages
  obj$last.stage <- result$last.stage
  obj$prob.complete.case <- result$prob.complete.case
  obj$d.hat <- result$d.hat
  obj$cens.mod.fitted <- result$cens.mod.fitted
  obj$opt.Y <- result$opt.Y
  obj$regret <- result$regret
  obj
}

enriched_obj <- enrich_obj(base_obj, dtr_result)



# .finalizeCovmat


local({
  # Hand-crafted: 4x4 covariance matrix, 2 blip vars at positions 3 and 4
  n_theta <- 4L
  covmat <- matrix(seq_len(n_theta^2), nrow = n_theta, ncol = n_theta)
  colnames(covmat) <- rownames(covmat) <- c("(Intercept)", "X", "A1", "A1:X")
  
  blip_vars <- c("A1" = 3L, "A1:X" = 4L)
  
  test_that(".finalizeCovmat: blip row names set correctly", {
    res <- .finalizeCovmat(covmat, blip_vars, full.cov = TRUE)
    expect_equal(rownames(res)[blip_vars], names(blip_vars))
  })
  
  test_that(".finalizeCovmat: blip col names set correctly", {
    res <- .finalizeCovmat(covmat, blip_vars, full.cov = TRUE)
    expect_equal(colnames(res)[blip_vars], names(blip_vars))
  })
  
  test_that(".finalizeCovmat: non-blip names unchanged", {
    res <- .finalizeCovmat(covmat, blip_vars, full.cov = TRUE)
    expect_equal(rownames(res)[1L:2L], c("(Intercept)", "X"))
    expect_equal(colnames(res)[1L:2L], c("(Intercept)", "X"))
  })
  
  test_that(".finalizeCovmat: full.cov=TRUE returns full matrix", {
    res <- .finalizeCovmat(covmat, blip_vars, full.cov = TRUE)
    expect_equal(dim(res), c(n_theta, n_theta))
  })
  
  test_that(".finalizeCovmat: full.cov=FALSE subsets to blip_vars square matrix", {
    res <- .finalizeCovmat(covmat, blip_vars, full.cov = FALSE)
    expect_equal(dim(res), c(length(blip_vars), length(blip_vars)))
  })
  
  test_that(".finalizeCovmat: full.cov=FALSE row and col names are psi names", {
    res <- .finalizeCovmat(covmat, blip_vars, full.cov = FALSE)
    expect_equal(rownames(res), names(blip_vars))
    expect_equal(colnames(res), names(blip_vars))
  })
  
  test_that(".finalizeCovmat: oracle - values at blip positions are correct", {
    res <- .finalizeCovmat(covmat, blip_vars, full.cov = FALSE)
    expected <- covmat[blip_vars, blip_vars, drop = FALSE]
    expect_equal(unname(res), unname(expected))
  })
})



# .getVariance - var.estim = "none"


test_that(".getVariance: var.estim='none' returns step.obj without covmat", {
  obj_none <- enriched_obj
  obj_none$var.estim <- "none"
  cc <- list(
    last.stage = obj_none$last.stage,
    prob.complete.case = obj_none$prob.complete.case,
    d.hat = obj_none$d.hat,
    cens.mod.fitted = obj_none$cens.mod.fitted
  )
  res <- .getVariance(obj_none, obj_none$stages[[2L]], cc, isSurvival = FALSE, k = 2L)
  expect_null(res$covmat)
})



# .getVariance - var.estim = "sandwich"


local({
  cc <- list(
    last.stage = enriched_obj$last.stage,
    prob.complete.case = enriched_obj$prob.complete.case,
    d.hat = enriched_obj$d.hat,
    cens.mod.fitted = enriched_obj$cens.mod.fitted
  )
  
  res_k2 <- .getVariance(enriched_obj, enriched_obj$stages[[2L]],
                         cc, isSurvival = FALSE, k = 2L)
  res_k1 <- .getVariance(enriched_obj, enriched_obj$stages[[1L]],
                         cc, isSurvival = FALSE, k = 1L)
  
  test_that(".getVariance sandwich: covmat is added to step.obj", {
    expect_false(is.null(res_k2$covmat))
    expect_false(is.null(res_k1$covmat))
  })
  
  test_that(".getVariance sandwich: covmat is a square matrix", {
    expect_true(is.matrix(res_k2$covmat))
    expect_equal(nrow(res_k2$covmat), ncol(res_k2$covmat))
    expect_true(is.matrix(res_k1$covmat))
    expect_equal(nrow(res_k1$covmat), ncol(res_k1$covmat))
  })
  
  test_that(".getVariance sandwich: full.cov=FALSE gives psi-only covmat", {
    n_psi_k2 <- length(enriched_obj$stages[[2L]]$psi)
    n_psi_k1 <- length(enriched_obj$stages[[1L]]$psi)
    expect_equal(dim(res_k2$covmat), c(n_psi_k2, n_psi_k2))
    expect_equal(dim(res_k1$covmat), c(n_psi_k1, n_psi_k1))
  })
  
  test_that(".getVariance sandwich: covmat row and col names contain psi names", {
    expect_true(all(names(enriched_obj$stages[[2L]]$psi) %in% rownames(res_k2$covmat)))
    expect_true(all(names(enriched_obj$stages[[2L]]$psi) %in% colnames(res_k2$covmat)))
  })
  
  test_that(".getVariance sandwich: covmat is symmetric", {
    expect_equal(res_k2$covmat, t(res_k2$covmat))
    expect_equal(res_k1$covmat, t(res_k1$covmat))
  })
  
  test_that(".getVariance sandwich: full.cov=TRUE returns larger covmat", {
    obj_full <- enriched_obj
    obj_full$full.cov <- TRUE
    res_full <- .getVariance(obj_full, obj_full$stages[[2L]],
                             cc, isSurvival = FALSE, k = 2L)
    n_theta <- length(coef(enriched_obj$stages[[2L]]$outcome.fit))
    expect_equal(dim(res_full$covmat), c(n_theta, n_theta))
  })
  
  test_that(".getVariance sandwich: full.cov=TRUE covmat larger than full.cov=FALSE", {
    obj_full <- enriched_obj
    obj_full$full.cov <- TRUE
    res_full <- .getVariance(obj_full, obj_full$stages[[2L]],
                             cc, isSurvival = FALSE, k = 2L)
    expect_true(nrow(res_full$covmat) >= nrow(res_k2$covmat))
  })
  
  test_that(".getVariance sandwich: other step.obj elements preserved", {
    expect_equal(res_k2$psi, enriched_obj$stages[[2L]]$psi)
    expect_equal(res_k2$beta, enriched_obj$stages[[2L]]$beta)
  })
})



# .computeVariance - var.estim = "sandwich"


local({
  res <- suppressMessages(
    .computeVariance(enriched_obj, isSurvival = FALSE)
  )
  
  test_that(".computeVariance sandwich: covmat added to each stage", {
    expect_false(is.null(res$stages[[1L]]$covmat))
    expect_false(is.null(res$stages[[2L]]$covmat))
  })
  
  test_that(".computeVariance sandwich: stages length unchanged", {
    expect_length(res$stages, enriched_obj$K)
  })
  
  test_that(".computeVariance sandwich: psi unchanged after variance computation", {
    for (k in seq_len(enriched_obj$K)) {
      expect_equal(res$stages[[k]]$psi, enriched_obj$stages[[k]]$psi)
    }
  })
  
  test_that(".computeVariance sandwich: each stage covmat is square matrix", {
    for (k in seq_len(enriched_obj$K)) {
      expect_true(is.matrix(res$stages[[k]]$covmat))
      expect_equal(nrow(res$stages[[k]]$covmat), ncol(res$stages[[k]]$covmat))
    }
  })
})



# .computeVariance - var.estim = "none"


test_that(".computeVariance none: no covmat added to stages", {
  obj_none <- enrich_obj(make_obj(var.estim = "none"), dtr_result)
  res <- .computeVariance(obj_none, isSurvival = FALSE)
  expect_true(all(sapply(res$stages, function(s) is.null(s$covmat))))
})



# .computeVariance - var.estim = "bootstrap" standard


local({
  withr::with_seed(42L, {
    obj_boot <- enrich_obj(make_obj(var.estim = "bootstrap"), dtr_result)
    res <- suppressMessages(.computeVariance(obj_boot, isSurvival = FALSE))
  })
  
test_that(".computeVariance bootstrap standard: returns covmat and psi.boot", {
  expect_true("covmat" %in% names(res))
  expect_true("psi.boot" %in% names(res))
})

test_that(".computeVariance bootstrap standard: covmat is list of length K", {
  expect_true(is.list(res$covmat))
  expect_length(res$covmat, enriched_obj$K)
})

})