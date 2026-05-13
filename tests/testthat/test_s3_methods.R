# Shared setup - one trusted DTRreg result reused across all tests

data(twoStageCont)
n <- nrow(twoStageCont)

# DWOLS, no variance
mod_none <- suppressMessages(
  DTRreg(twoStageCont$Y,
         blip.mod  = list(~X1, ~X2),
         treat.mod = list(A1 ~ X1, A2 ~ X2),
         tf.mod    = list(~X1, ~X2),
         data      = twoStageCont,
         method    = "dwols",
         var.estim = "none")
)

# DWOLS, sandwich variance
mod_sw <- suppressMessages(
  DTRreg(twoStageCont$Y,
         blip.mod  = list(~X1, ~X2),
         treat.mod = list(A1 ~ X1, A2 ~ X2),
         tf.mod    = list(~X1, ~X2),
         data      = twoStageCont,
         method    = "dwols",
         var.estim = "sandwich")
)

# DWOLS, bootstrap variance
withr::with_seed(42L, {
  mod_boot <- suppressMessages(
    DTRreg(twoStageCont$Y,
           blip.mod       = list(~X1, ~X2),
           treat.mod      = list(A1 ~ X1, A2 ~ X2),
           tf.mod         = list(~X1, ~X2),
           data           = twoStageCont,
           method         = "dwols",
           var.estim      = "bootstrap",
           bootstrap.controls  = list(B = 20L))
  )
})

# GEST, sandwich variance
mod_gest <- suppressMessages(
  DTRreg(twoStageCont$Y,
         blip.mod  = list(~X1, ~X2),
         treat.mod = list(A1 ~ X1, A2 ~ X2),
         tf.mod    = list(~X1, ~X2),
         data      = twoStageCont,
         method    = "gest",
         var.estim = "sandwich")
)


# print.DTRreg

test_that("print.DTRreg: produces output without error (var.estim='none')", {
  expect_output(print.DTRreg(mod_none))
})

test_that("print.DTRreg: produces output without error (sandwich)", {
  expect_output(print.DTRreg(mod_sw))
})

test_that("print.DTRreg: produces output without error (bootstrap)", {
  expect_output(print.DTRreg(mod_boot))
})

test_that("print.DTRreg: output contains 'Stage 1'", {
  expect_output(print.DTRreg(mod_none), "Stage 1")
})

test_that("print.DTRreg: output contains 'Stage 2'", {
  expect_output(print.DTRreg(mod_none), "Stage 2")
})

test_that("print.DTRreg: no-variance output contains psi names", {
  expect_output(print.DTRreg(mod_none), "A1")
  expect_output(print.DTRreg(mod_none), "A2")
})

test_that("print.DTRreg: sandwich output contains Std_Error column", {
  expect_output(print.DTRreg(mod_sw), "Std_Error")
})

test_that("print.DTRreg: var.estim='none' does not print Std_Error", {
  expect_no_match(capture.output(print.DTRreg(mod_none)), "Std_Error")
})

test_that("print.DTRreg: DTR type prints treatment rule", {
  expect_output(print.DTRreg(mod_none), "treat if")
})


# summary.DTRreg

test_that("summary.DTRreg: produces same output as print.DTRreg", {
  expect_equal(capture.output(summary.DTRreg(mod_none)),
               capture.output(print.DTRreg(mod_none)))
})

test_that("summary.DTRreg: produces output without error", {
  expect_output(summary.DTRreg(mod_sw))
})


# coef.DTRreg

test_that("coef.DTRreg: returns a list of length K", {
  res <- coef.DTRreg(mod_none)
  expect_true(is.list(res))
  expect_length(res, mod_none$K)
})

test_that("coef.DTRreg: list names are stage_1, stage_2, ...", {
  res <- coef.DTRreg(mod_none)
  expect_equal(names(res), paste0("stage_", seq_len(mod_none$K)))
})

test_that("coef.DTRreg: each element matches corresponding psi", {
  res <- coef.DTRreg(mod_none)
  for (k in seq_len(mod_none$K)) {
    expect_equal(res[[k]], mod_none$psi[[k]])
  }
})

test_that("coef.DTRreg: psi names preserved within each stage", {
  res <- coef.DTRreg(mod_none)
  expect_true(all(grepl("A1", names(res[["stage_1"]]))))
  expect_true(all(grepl("A2", names(res[["stage_2"]]))))
})

test_that("coef.DTRreg: works for GEST object", {
  res <- coef.DTRreg(mod_gest)
  expect_length(res, mod_gest$K)
  expect_equal(names(res), paste0("stage_", seq_len(mod_gest$K)))
})


# confint.DTRreg

# --- var.estim = "none" early return -----------------------------------------

test_that("confint.DTRreg: returns invisible NULL when var.estim='none'", {
  res <- confint.DTRreg(mod_none)
  expect_null(res)
})

test_that("confint.DTRreg: prints message when var.estim='none'", {
   expect_message(confint.DTRreg(mod_none), "only available")
})

# --- se type (sandwich) ------------------------------------------------------

local({
  res <- confint.DTRreg(mod_sw)
  
  test_that("confint.DTRreg se: returns list of length K", {
    expect_true(is.list(res))
    expect_length(res, mod_sw$K)
  })
  
  test_that("confint.DTRreg se: list names are stage_1, stage_2, ...", {
    expect_equal(names(res), paste0("stage_", seq_len(mod_sw$K)))
  })
  
  test_that("confint.DTRreg se: each element is a matrix with 2 columns", {
    for (k in seq_len(mod_sw$K)) {
      expect_true(is.matrix(res[[k]]))
      expect_equal(ncol(res[[k]]), 2L)
    }
  })
  
  test_that("confint.DTRreg se: row count matches psi length per stage", {
    for (k in seq_len(mod_sw$K)) {
      expect_equal(nrow(res[[k]]), length(mod_sw$psi[[k]]))
    }
  })
  
  test_that("confint.DTRreg se: row names match psi names", {
    for (k in seq_len(mod_sw$K)) {
      expect_equal(rownames(res[[k]]), names(mod_sw$psi[[k]]))
    }
  })
  
  test_that("confint.DTRreg se: lower CI < upper CI for all parameters", {
    for (k in seq_len(mod_sw$K)) {
      expect_true(all(res[[k]][, 1L] < res[[k]][, 2L]))
    }
  })
  
  test_that("confint.DTRreg se: column names reflect default level 0.95", {
    expect_equal(colnames(res[[1L]]), c("2.5 %", "97.5 %"))
  })
  
  test_that("confint.DTRreg se: oracle - CI contains psi estimate", {
    # psi should be within its own CI
    for (k in seq_len(mod_sw$K)) {
      expect_true(all(mod_sw$psi[[k]] >= res[[k]][, 1L] &
                        mod_sw$psi[[k]] <= res[[k]][, 2L]))
    }
  })
  
  test_that("confint.DTRreg se: wider CI with lower level", {
    res_90 <- confint.DTRreg(mod_sw, level = 0.90)
    for (k in seq_len(mod_sw$K)) {
      width_95 <- res[[k]][, 2L]    - res[[k]][, 1L]
      width_90 <- res_90[[k]][, 2L] - res_90[[k]][, 1L]
      expect_true(all(width_95 > width_90))
    }
  })
})

# --- percentile type (bootstrap) ---------------------------------------------

local({
  res <- confint.DTRreg(mod_boot, type = "percentile")
  
  test_that("confint.DTRreg percentile: returns list of length K", {
    expect_length(res, mod_boot$K)
  })
  
  test_that("confint.DTRreg percentile: each element is a matrix with 2 columns", {
    for (k in seq_len(mod_boot$K)) {
      expect_true(is.matrix(res[[k]]))
      expect_equal(ncol(res[[k]]), 2L)
    }
  })
  
  test_that("confint.DTRreg percentile: lower CI < upper CI", {
    for (k in seq_len(mod_boot$K)) {
      expect_true(all(res[[k]][, 1L] < res[[k]][, 2L]))
    }
  })
})

test_that("confint.DTRreg percentile: stops when var.estim is not bootstrap", {
  expect_error(confint.DTRreg(mod_sw, type = "percentile"),
               "Percentile confidence intervals only available with bootstrap")
})


# predict.DTRreg

local({
  # With newdata
  res_new <- predict.DTRreg(mod_none, newdata = twoStageCont)
  # Without newdata (uses training data)
  res_trn <- predict.DTRreg(mod_none)
  
  test_that("predict.DTRreg: returns numeric vector", {
    expect_true(is.numeric(res_new))
    expect_true(is.numeric(res_trn))
  })
  
  test_that("predict.DTRreg: length equals nrow(newdata)", {
    expect_length(res_new, n)
  })
  
  test_that("predict.DTRreg: length equals n when no newdata", {
    expect_length(res_trn, n)
  })
  
  test_that("predict.DTRreg: newdata and training data give same result for same data", {
    expect_equal(res_new, res_trn)
  })
  
  test_that("predict.DTRreg: stops when newdata has missing values", {
    dat_na      <- twoStageCont
    dat_na$X1[1L] <- NA_real_
    expect_error(predict.DTRreg(mod_none, newdata = dat_na),
                 "must be complete")
  })
  
  test_that("predict.DTRreg: treat.range list is accepted without error", {
    expect_no_error(
      predict.DTRreg(mod_none, newdata = twoStageCont,
                     treat.range = NULL)
    )
  })
  
  test_that("predict.DTRreg: stops when treat.range is invalid", {
    expect_error(
      predict.DTRreg(mod_none, newdata = twoStageCont,
                     treat.range = c(0, 1, 2)),
      "treat.range"
    )
  })
})

test_that("predict.DTRreg: works for GEST object", {
  res <- predict.DTRreg(mod_gest, newdata = twoStageCont)
  expect_true(is.numeric(res))
  expect_length(res, n)
})


data(twoStageCont)

mod_dwols <- suppressMessages(
  DTRreg(twoStageCont$Y,
         blip.mod  = list(~X1, ~X2),
         treat.mod = list(A1 ~ X1, A2 ~ X2),
         tf.mod    = list(~X1, ~X2),
         data      = twoStageCont,
         method    = "dwols",
         var.estim = "none")
)

mod_gest <- suppressMessages(
  DTRreg(twoStageCont$Y,
         blip.mod  = list(~X1, ~X2),
         treat.mod = list(A1 ~ X1, A2 ~ X2),
         tf.mod    = list(~X1, ~X2),
         data      = twoStageCont,
         method    = "gest",
         var.estim = "none")
)

# Expected number of residual plots per stage:
# 1 (Fitted) + n_blip_covariates
# blip = ~X1 or ~X2 has 1 covariate → 2 plots per stage → 4 total for 2 stages
n_resid_plots_per_stage <- 2L  # Fitted + X covariate
n_stages <- 2L
n_resid_plots_total <- n_resid_plots_per_stage * n_stages


# autoplot.DTRreg - list structure

local({
  plots <- autoplot.DTRreg(mod_dwols)
  
  test_that("autoplot.DTRreg: returns a list", {
    expect_true(is.list(plots))
  })
  
  test_that("autoplot.DTRreg: contains at least K*2 residual plots", {
    expect_true(length(plots) >= n_resid_plots_total)
  })
  
  test_that("autoplot.DTRreg: names contain stage identifiers", {
    nms <- names(plots)
    expect_true(any(grepl("Stage1", nms)))
    expect_true(any(grepl("Stage2", nms)))
  })
  
  test_that("autoplot.DTRreg: Fitted plots named correctly", {
    expect_true("Stage1_Fitted" %in% names(plots))
    expect_true("Stage2_Fitted" %in% names(plots))
  })
  
  test_that("autoplot.DTRreg: blip covariate plots named correctly", {
    expect_true("Stage1_X1" %in% names(plots))
    expect_true("Stage2_X2" %in% names(plots))
  })
})


# autoplot.DTRreg - individual plot content

local({
  plots <- autoplot.DTRreg(mod_dwols)
  p_fit <- plots[["Stage1_Fitted"]]
  p_cov <- plots[["Stage1_X1"]]
  
  test_that("autoplot.DTRreg: Fitted plot is a ggplot", {
    expect_true(inherits(p_fit, "ggplot"))
  })
  
  test_that("autoplot.DTRreg: covariate plot is a ggplot", {
    expect_true(inherits(p_cov, "ggplot"))
  })
  
  test_that("autoplot.DTRreg: Fitted plot data has correct number of rows", {
    n_stage1 <- sum(mod_dwols$analysis$last.stage >= 1L)
    expect_equal(nrow(p_fit$data), n_stage1)
  })
  
  test_that("autoplot.DTRreg: Fitted plot data has x and y columns", {
    expect_true(all(c("x", "y") %in% names(p_fit$data)))
  })
  
  test_that("autoplot.DTRreg: Fitted plot x-axis label is 'Fitted values'", {
    expect_equal(p_fit$labels$x, "Fitted values")
  })
  
  test_that("autoplot.DTRreg: Fitted plot y-axis label is 'Residuals'", {
    expect_equal(p_fit$labels$y, "Residuals")
  })
  
  test_that("autoplot.DTRreg: Fitted plot title contains stage number", {
    expect_match(p_fit$labels$title, "Stage 1")
  })
  
  test_that("autoplot.DTRreg: covariate plot x-axis label is covariate name", {
    expect_equal(p_cov$labels$x, "X1")
  })
  
  test_that("autoplot.DTRreg: plots contain point layer", {
    layer_geoms <- sapply(p_fit$layers, function(l) class(l$geom)[1L])
    expect_true("GeomPoint" %in% layer_geoms)
  })
  
  test_that("autoplot.DTRreg: plots contain hline layer", {
    layer_geoms <- sapply(p_fit$layers, function(l) class(l$geom)[1L])
    expect_true("GeomHline" %in% layer_geoms)
  })
  
  test_that("autoplot.DTRreg: plots contain smooth layer", {
    layer_geoms <- sapply(p_fit$layers, function(l) class(l$geom)[1L])
    expect_true("GeomSmooth" %in% layer_geoms)
  })
  
  test_that("autoplot.DTRreg: x values in Fitted plot match fitted.values", {
    fitK <- mod_dwols$analysis$outcome.fit[[1L]]$fitted.values
    expect_equal(p_fit$data$x, unname(fitK))
  })
  
  test_that("autoplot.DTRreg: y values in Fitted plot are residuals oracle", {
    obsK  <- mod_dwols$analysis$Y[[1L]]
    fitK  <- mod_dwols$analysis$outcome.fit[[1L]]$fitted.values
    expected_resids <- obsK - fitK
    expect_equal(p_fit$data$y, unname(expected_resids))
  })
})


# autoplot.DTRreg - stage 2 correctness

local({
  plots   <- autoplot.DTRreg(mod_dwols)
  p_fit2  <- plots[["Stage2_Fitted"]]
  n_stage2 <- sum(mod_dwols$analysis$last.stage >= 2L)
  
  test_that("autoplot.DTRreg: Stage 2 Fitted plot has correct number of rows", {
    expect_equal(nrow(p_fit2$data), n_stage2)
  })
  
  test_that("autoplot.DTRreg: Stage 2 title contains 'Stage 2'", {
    expect_match(p_fit2$labels$title, "Stage 2")
  })
  
  test_that("autoplot.DTRreg: Stage 2 residual oracle matches Y - fitted", {
    obsK     <- mod_dwols$analysis$Y[[2L]]
    fitK     <- mod_dwols$analysis$outcome.fit[[2L]]$fitted.values
    expected <- obsK - fitK
    expect_equal(p_fit2$data$y, unname(expected))
  })
})


# autoplot.DTRreg - intercept-only blip produces no covariate plots

test_that("autoplot.DTRreg: intercept-only blip produces no covariate plot for that stage", {
  mod_int <- suppressMessages(
    DTRreg(twoStageCont$Y,
           blip.mod = list(~1, ~X2),
           treat.mod = list(A1 ~ X1, A2 ~ X2),
           tf.mod = list(~X1, ~X2),
           data = twoStageCont,
           method = "dwols",
           var.estim = "none")
  )
  plots <- autoplot.DTRreg(mod_int)
  # Stage 1 has intercept-only blip - no X covariate plot
  expect_false(any(grepl("Stage1_X", names(plots))))
  # Stage 1 still has Fitted plot
  expect_true("Stage1_Fitted" %in% names(plots))
  # Stage 2 still has X2 covariate plot
  expect_true("Stage2_X2" %in% names(plots))
})


# autoplot.DTRreg - GEST object

test_that("autoplot.DTRreg: works for GEST object", {
  plots <- autoplot.DTRreg(mod_gest)
  expect_true(is.list(plots))
  expect_true("Stage1_Fitted" %in% names(plots))
  expect_true(inherits(plots[["Stage1_Fitted"]], "ggplot"))
})


# plot.DTRreg

test_that("plot.DTRreg: returns invisibly", {
  pdf(NULL)
  on.exit(dev.off())
  res <- withVisible(plot.DTRreg(mod_dwols))
  expect_false(res$visible)
})

test_that("plot.DTRreg: invisible return value is a list", {
  pdf(NULL)
  on.exit(dev.off())
  res <- plot.DTRreg(mod_dwols)
  expect_true(is.list(res))
})

test_that("plot.DTRreg: return value matches autoplot.DTRreg output", {
  pdf(NULL)
  on.exit(dev.off())
  res_plot <- plot.DTRreg(mod_dwols)
  res_autoplot <- autoplot.DTRreg(mod_dwols)
  expect_equal(names(res_plot), names(res_autoplot))
})

test_that("plot.DTRreg: each returned element is a ggplot", {
  pdf(NULL)
  on.exit(dev.off())
  res <- plot.DTRreg(mod_dwols)
  ggplot_elements <- sapply(res, function(p) inherits(p, "ggplot"))
  # At minimum the residual plots should all be ggplot objects
  resid_plots <- res[grepl("Fitted|_X", names(res))]
  expect_true(all(sapply(resid_plots, function(p) inherits(p, "ggplot"))))
})

test_that("plot.DTRreg: works for GEST object without error", {
  pdf(NULL)
  on.exit(dev.off())
  expect_no_error(plot.DTRreg(mod_gest))
})


# autoplot.DTRreg - submodel plots (treatment, complete cases)

local({
  plots <- autoplot.DTRreg(mod_dwols)
  nms   <- names(plots)
  
  test_that("autoplot.DTRreg: treatment model plots present for each stage", {
    expect_true("Stage1_Treatment" %in% nms)
    expect_true("Stage2_Treatment" %in% nms)
  })
  
  test_that("autoplot.DTRreg: no CompleteCases plots when censoring not modeled", {
    expect_false(any(grepl("CompleteCases", names(plots))))
  })
  
  test_that("autoplot.DTRreg: treatment plots are ggplot objects", {
    tx_plots <- plots[grepl("_Treatment$", nms)]
    expect_true(all(sapply(tx_plots, inherits, what = "ggplot")))
  })
  
  n_treatment_plots_total <- n_stages   # 1 per stage
  n_total_plots <- n_resid_plots_total + n_treatment_plots_total  # 6

  test_that("autoplot.DTRreg: total plot count equals residual plus treatment plots", {
    expect_equal(length(plots), n_total_plots)
  })
})

test_that("autoplot.DTRreg: no CompleteCases plots for DWSurv object", {
  data(twoStageCens)
  mod_surv <- suppressMessages(
    DWSurv(time = list(~T1, ~T2),
           blip.mod = list(~X11, ~X21),
           treat.mod = list(A1 ~ X11, A2 ~ X21),
           tf.mod = list(~X11, ~X21),
           cens.mod = list(delta~X11, delta~X21),
           var.estim = "none",
           data = twoStageCens)
  )
  plots <- autoplot.DTRreg(mod_surv)
  nms   <- names(plots)
  expect_false(any(grepl("CompleteCases", nms)))
  expect_true(any(grepl("Censoring", nms)))
})
