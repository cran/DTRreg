# Shared minimal data
n_obs <- 20L

dat_bin <- data.frame(
  X = seq(-1.0, 1.0, length.out = n_obs),
  A = rep(0L:1L, times = n_obs / 2L),
  Y = seq_len(n_obs) / n_obs
)

dat_multi <- data.frame(
  X = seq(-1.0, 1.0, length.out = n_obs),
  A = factor(rep(c("low", "med", "high"), length.out = n_obs),
             levels = c("low", "med", "high")),
  Y = seq_len(n_obs) / n_obs
)

dat_cont <- data.frame(
  X = cos(seq(0.0, pi, length.out = n_obs)),  # not linear in A
  A = seq(0.1, 0.9, length.out = n_obs),
  Y = seq_len(n_obs) / n_obs
)

test_that("Binary$initialize: blip.model has no intercept and contains tx.var", {
  obj <- Binary$new(tf.model = ~X, blip.model = ~X, tx.var = "A")
  blip_terms <- attr(terms(obj$blip.model), "term.labels")
  expect_equal(attr(terms(obj$blip.model), "intercept"), 0L)
  expect_true("A" %in% blip_terms)
})

test_that("Binary$initialize: blip.model with covariate adds interaction term A:X", {
  obj <- Binary$new(tf.model = ~X, blip.model = ~X, tx.var = "A")
  blip_terms <- attr(terms(obj$blip.model), "term.labels")
  expect_true("A:X" %in% blip_terms)
})

test_that("Binary$initialize: intercept-only blip produces only tx.var in blip.model", {
  obj <- Binary$new(tf.model = ~X, blip.model = ~1, tx.var = "A")
  blip_terms <- attr(terms(obj$blip.model), "term.labels")
  expect_identical(blip_terms, "A")
})

test_that("Binary$initialize: full.model contains tf and blip terms", {
  obj <- Binary$new(tf.model = ~X, blip.model = ~X, tx.var = "A")
  full_terms <- attr(terms(obj$full.model), "term.labels")
  expect_true("X" %in% full_terms)
  expect_true("A" %in% full_terms)
  expect_true("A:X" %in% full_terms || "X:A" %in% full_terms)
})

test_that("Binary$blip_params: returns named indices for coefficients containing tx.var", {
  obj <- Binary$new(tf.model = ~X, blip.model = ~X, tx.var = "A")
  coefs <- c("(Intercept)" = 1.0, "X" = 2.0, "A" = 3.0, "A:X" = 4.0)
  res <- obj$blip_params(coefs)
  expect_named(res, c("A", "A:X"), ignore.order = FALSE)
  expect_equal(unname(res), c(3L, 4L))
})

test_that("Binary$blip_params: treatment-free coefficients are excluded", {
  obj <- Binary$new(tf.model = ~X, blip.model = ~X, tx.var = "A")
  coefs <- c("(Intercept)" = 1.0, "X" = 2.0, "A" = 3.0, "A:X" = 4.0)
  res <- obj$blip_params(coefs)
  expect_false("(Intercept)" %in% names(res))
  expect_false("X" %in% names(res))
})

test_that("Binary$ipw: matches 1 / (A*Ahat + (1-A)*(1-Ahat))", {
  obj <- Binary$new(tf.model = ~X, blip.model = ~X, tx.var = "A")
  A <- c(1L, 0L, 1L, 0L, 1L)
  Ahat <- c(0.8, 0.3, 0.6, 0.4, 0.9)
  expected <- 1.0 / (A * Ahat + (1.0 - A) * (1.0 - Ahat))
  expect_equal(obj$ipw(A, Ahat), expected)
})

test_that("Binary$.pom: stops - not appropriate for binary treatments", {
  obj <- Binary$new(tf.model = ~X, blip.model = ~X, tx.var = "A")
  expect_error(obj$.pom(), "not appropriate for binary")
})

test_that("Binary$prep: returns data unchanged", {
  obj <- Binary$new(tf.model = ~X, blip.model = ~X, tx.var = "A")
  expect_identical(obj$prep(dat_bin, dat_bin$A), dat_bin)
})

test_that("Binary$Hd: returns matrix with correct dimensions and expected column names", {
  obj <- Binary$new(tf.model = ~X, blip.model = ~X, tx.var = "A")
  mm <- obj$Hd(dat_bin, dat_bin$A)
  expect_true(is.matrix(mm))
  expect_equal(nrow(mm), n_obs)
  expect_true(all(c("X", "A", "A:X") %in% colnames(mm)) ||
              all(c("X", "A", "X:A") %in% colnames(mm)))
})

test_that("Binary$Hpsi: returns no-intercept matrix containing tx.var term", {
  obj <- Binary$new(tf.model = ~X, blip.model = ~X, tx.var = "A")
  mm <- obj$Hpsi(dat_bin, dat_bin$A)
  expect_true(is.matrix(mm))
  expect_false("(Intercept)" %in% colnames(mm))
  expect_true("A" %in% colnames(mm))
})

test_that("Binary$Hw: returns matrix with same dimensions as Hd", {
  obj <- Binary$new(tf.model = ~X, blip.model = ~X, tx.var = "A")
  wgt <- rep(1.0, n_obs)
  Ahat <- rep(0.5, n_obs)
  expect_equal(dim(obj$Hw(dat_bin, dat_bin$A, Ahat, wgt)),
               dim(obj$Hd(dat_bin, dat_bin$A)))
})

local({
  obj <- Binary$new(tf.model = ~X, blip.model = ~X, tx.var = "A")
  dat <- dat_bin
  dat$Y_int <- dat$Y
  fit <- lm(update(obj$full.model, Y_int ~ .), dat)
  
  test_that("Binary$regret: returns numeric vector of length n", {
    res <- obj$regret(fit, dat, opt = dat$A, A = dat$A)
    expect_true(is.numeric(res))
    expect_length(res, n_obs)
  })
  
  test_that("Binary$regret: is zero when opt equals observed treatment", {
    res <- obj$regret(fit, dat, opt = dat$A, A = dat$A)
    expect_equal(unname(res), rep(0.0, n_obs))
  })
  
  test_that("Binary$opt: returns integer 0/1 vector of length n", {
    res <- obj$opt(fit, dat)
    expect_true(is.integer(res))
    expect_true(all(res %in% c(0L, 1L)))
    expect_length(res, n_obs)
  })
  
  test_that("Binary$opt: matches oracle - treat if predicted(A=1) - predicted(A=0) > 1e-8", {
    dat0 <- dat; dat0$A <- 0L
    dat1 <- dat; dat1$A <- 1L
    expected <- as.integer((predict(fit, dat1) - predict(fit, dat0)) > 1e-8)
    expect_equal(obj$opt(fit, dat), expected)
  })
  
  test_that("Binary$shiftY: 'DTR' passes opt through to regret", {
    opt_vec <- rep(0L, n_obs)
    expect_equal(obj$shiftY("DTR", fit, dat, opt = opt_vec, A = dat$A),
                 obj$regret(fit, dat, opt = opt_vec, A = dat$A))
  })
  
  test_that("Binary$shiftY: non-DTR type uses 0L as opt (effect parameterisation)", {
    expect_equal(obj$shiftY("effect", fit, dat, opt = rep(1L, n_obs), A = dat$A),
                 obj$regret(fit, dat, opt = 0L, A = dat$A))
  })
})


# MultiNom
tx_levels <- c("low", "med", "high")

test_that("MultiNom$initialize: blip.model has no intercept and contains tx.var terms", {
  obj<- MultiNom$new(tf.model = ~X, blip.model = ~X, tx.var = "A", tx.levels = tx_levels)
  blip_terms <- attr(terms(obj$blip.model), "term.labels")
  expect_equal(attr(terms(obj$blip.model), "intercept"), 0L)
  expect_true(any(grepl("^A", blip_terms)))
})

test_that("MultiNom$initialize: full.model contains tf covariate and blip terms", {
  obj <- MultiNom$new(tf.model = ~X, blip.model = ~X, tx.var = "A", tx.levels = tx_levels)
  full_terms <- attr(terms(obj$full.model), "term.labels")
  expect_true("X" %in% full_terms)
  expect_true(any(grepl("^A", full_terms)))
})

test_that("MultiNom$blip_params: returns indices for non-reference level tx terms", {
  obj <- MultiNom$new(tf.model = ~X, blip.model = ~X, tx.var = "A", tx.levels = tx_levels)
  coefs <- c("(Intercept)" = 1.0, "X" = 2.0,
             "Amed" = 3.0, "Ahigh" = 4.0, "Amed:X" = 5.0, "Ahigh:X" = 6.0)
  res <- obj$blip_params(coefs)
  expect_equal(sort(unname(res)), c(3L, 4L, 5L, 6L))
})

test_that("MultiNom$blip_params: reference level and tf coefficients are excluded", {
  obj <- MultiNom$new(tf.model = ~X, blip.model = ~X, tx.var = "A", tx.levels = tx_levels)
  coefs <- c("(Intercept)" = 1.0, "X" = 2.0,
             "Amed" = 3.0, "Ahigh" = 4.0, "Amed:X" = 5.0, "Ahigh:X" = 6.0)
  res <- obj$blip_params(coefs)
  expect_false("(Intercept)" %in% names(res))
  expect_false("X"           %in% names(res))
  expect_false(any(grepl("Alow", names(res))))
})

test_that("MultiNom$ipw: matches 1 / prob(observed treatment level) for each observation", {
  obj <- MultiNom$new(tf.model = ~X, blip.model = ~1, tx.var = "A", tx.levels = tx_levels)
  A <- factor(c("low", "med", "high", "med", "low"), levels = tx_levels)
  Ahat <- matrix(c(0.5, 0.3, 0.2,
                   0.2, 0.6, 0.2,
                   0.1, 0.2, 0.7,
                   0.3, 0.5, 0.2,
                   0.6, 0.2, 0.2), nrow = 5L, byrow = TRUE)
  expected <- sapply(seq_along(A), function(i) {
    1.0 / Ahat[i, match(as.character(A[i]), tx_levels)]
  })
  expect_equal(obj$ipw(A, Ahat), expected)
})

test_that("MultiNom$Hpsi: stops - g-estimation not supported for multinomial", {
  obj <- MultiNom$new(tf.model = ~X, blip.model = ~1, tx.var = "A",
                      tx.levels = tx_levels)
  expect_error(obj$Hpsi(dat_multi, dat_multi$A), "not yet supported")
})

test_that("MultiNom$Hw: stops - g-estimation not supported for multinomial", {
  obj  <- MultiNom$new(tf.model = ~X, blip.model = ~1, tx.var = "A",
                       tx.levels = tx_levels)
  Ahat <- matrix(1.0 / 3.0, nrow = n_obs, ncol = 3L)
  expect_error(obj$Hw(dat_multi, dat_multi$A, Ahat, rep(1.0, n_obs)),
               "not yet supported")
})

test_that("MultiNom$prep: returns data unchanged", {
  obj <- MultiNom$new(tf.model = ~X, blip.model = ~1, tx.var = "A",
                      tx.levels = tx_levels)
  expect_identical(obj$prep(dat_multi, dat_multi$A), dat_multi)
})

test_that("MultiNom$Hd: returns matrix with n rows", {
  obj <- MultiNom$new(tf.model = ~X, blip.model = ~1, tx.var = "A",
                      tx.levels = tx_levels)
  mm <- obj$Hd(dat_multi, dat_multi$A)
  expect_true(is.matrix(mm))
  expect_equal(nrow(mm), n_obs)
})

local({
  obj <- MultiNom$new(tf.model = ~X, blip.model = ~1, tx.var = "A",
                            tx.levels = tx_levels)
  dat <- dat_multi
  dat$Y_int <- dat$Y
  fit <- lm(update(obj$full.model, Y_int ~ .), dat)
  
  test_that("MultiNom$regret: returns numeric vector of length n", {
    res <- obj$regret(fit, dat, opt = dat$A, A = dat$A)
    expect_true(is.numeric(res))
    expect_length(res, n_obs)
  })
  
  test_that("MultiNom$regret: is zero when opt equals observed treatment", {
    res <- obj$regret(fit, dat, opt = dat$A, A = dat$A)
    expect_equal(unname(res), rep(0.0, n_obs))
  })
  
  test_that("MultiNom$opt: returns values drawn from tx.levels", {
    res <- obj$opt(fit, dat)
    expect_true(all(res %in% tx_levels))
    expect_length(res, n_obs)
  })
  
  test_that("MultiNom$opt: matches oracle - level with highest predicted outcome", {
    pred_mat <- sapply(tx_levels, function(lv) {
      d <- dat
      d$A <- factor(lv, levels = tx_levels)
      predict(fit, d)
    })
    expected <- factor(tx_levels[apply(pred_mat, 1L, which.max)], levels = tx_levels)
    expect_equal(obj$opt(fit, dat), expected)
  })
  
  test_that("MultiNom$shiftY: 'DTR' passes opt through to regret", {
    opt_vec <- factor(rep("low", n_obs), levels = tx_levels)
    expect_equal(obj$shiftY("DTR",fit, dat, opt = opt_vec, A = dat$A),
                 obj$regret(fit, dat, opt = opt_vec, A = dat$A))
  })
  
  test_that("MultiNom$shiftY: non-DTR type uses first tx.level as reference", {
    opt_vec <- factor(rep("med", n_obs), levels = tx_levels)
    expect_equal(obj$shiftY("effect", fit, dat, opt = opt_vec, A = dat$A),
                 obj$regret(fit, dat, opt = tx_levels[1L], A = dat$A))
  })
})

# ContLinearBlip

test_that("ContLinearBlip$initialize: stops immediately - not supported", {
  expect_error(
    ContLinearBlip$new(tf.model = ~X, blip.model = ~X, tx.var = "A",
                       treat.range = c(0.0, 1.0)),
    "not supported"
  )
})


# ContQuadraticBlip

test_that("ContQuadraticBlip$initialize: stops on intercept-only blip model", {
  expect_error(
    ContQuadraticBlip$new(tf.model = ~X, blip.model = ~1, tx.var = "A",
                          treat.range = c(0.0, 1.0)),
    "intercept only"
  )
})

test_that("ContQuadraticBlip$initialize: blip.model has no intercept and contains tx.var", {
  obj <- ContQuadraticBlip$new(tf.model = ~X, blip.model = ~A, tx.var = "A",
                               treat.range = c(0.0, 1.0))
  blip_terms <- attr(terms(obj$blip.model), "term.labels")
  expect_equal(attr(terms(obj$blip.model), "intercept"), 0L)
  expect_true("A" %in% blip_terms)
})

test_that("ContQuadraticBlip$initialize: full.model contains tf covariate and blip terms", {
  obj <- ContQuadraticBlip$new(tf.model = ~X, blip.model = ~A, tx.var = "A",
                               treat.range = c(0.0, 1.0))
  full_terms <- attr(terms(obj$full.model), "term.labels")
  expect_true("X" %in% full_terms)
  expect_true("A" %in% full_terms)
})

test_that("ContQuadraticBlip$blip_params: identifies both linear and quadratic tx terms", {
  obj <- ContQuadraticBlip$new(tf.model = ~X, blip.model = ~A, tx.var = "A",
                               treat.range = c(0.0, 1.0))
  dat <- dat_cont
  prepped <- obj$prep(dat, dat$A)
  prepped$Y_int <- prepped$Y
  fit <- lm(update(obj$full.model, Y_int ~ .), prepped)
  res <- obj$blip_params(coef(fit))
  # Treatment-free terms must be absent; tx terms must be present
  expect_false("(Intercept)" %in% names(res))
  expect_false("X" %in% names(res))
  expect_true(length(res) >= 1L)
  expect_true("A" %in% names(res))
})

test_that("ContQuadraticBlip$ipw gaussian: matches 1 / dnorm(A, fitted, sqrt(dispersion))", {
  obj <- ContQuadraticBlip$new(tf.model = ~X, blip.model = ~A, tx.var = "A",
                               treat.range = c(0.0, 1.0))
  fit <- glm(A ~ X, dat_cont, family = gaussian())
  disp <- summary(fit)$dispersion
  expected <- 1.0 / dnorm(dat_cont$A, mean = fit$fitted.values, sd = sqrt(disp))
  expect_equal(obj$ipw(dat_cont$A, fit), expected)
})

test_that("ContQuadraticBlip$ipw Gamma: matches 1 / dgamma(A, shape, scale)", {
  dat_pos <- data.frame(X = seq(-1.0, 1.0, length.out = n_obs),
                        A = seq(0.1, 1.0, length.out = n_obs))
  obj <- ContQuadraticBlip$new(tf.model = ~X, blip.model = ~A, tx.var = "A",
                               treat.range = c(0.1, 1.0))
  fit <- glm(A ~ X, dat_pos, family = Gamma(link = "log"))
  disp <- summary(fit)$dispersion
  expected <- 1.0 / dgamma(dat_pos$A, shape = 1.0 / disp,
                           scale = fit$fitted.values * disp)
  expect_equal(obj$ipw(dat_pos$A, fit), expected)
})

test_that("ContQuadraticBlip$ipw: stops on unsupported glm family", {
  obj <- ContQuadraticBlip$new(tf.model = ~X, blip.model = ~A, tx.var = "A",
                               treat.range = c(0.0, 1.0))
  fit <- glm(A ~ X, dat_cont, family = gaussian())
  fit$family$family <- "poisson"
  expect_error(obj$ipw(dat_cont$A, fit), "unsupported family")
})

test_that("ContQuadraticBlip$.pom: stops when tx.mod.fitted is NULL", {
  obj <- ContQuadraticBlip$new(tf.model = ~X, blip.model = ~A, tx.var = "A",
                               treat.range = c(0.0, 1.0))
  expect_error(obj$.pom(dat_cont$A, NULL, dat_cont, 4L))
})

test_that("ContQuadraticBlip$.pom: stops when tx.mod.fitted is not a glm", {
  obj <- ContQuadraticBlip$new(tf.model = ~X, blip.model = ~A, tx.var = "A",
                               treat.range = c(0.0, 1.0))
  expect_error(obj$.pom(dat_cont$A, list(not = "a glm"), dat_cont, 4L))
})

test_that("ContQuadraticBlip$prep: adds exactly one new column with values A^2", {
  obj <- ContQuadraticBlip$new(tf.model = ~X, blip.model = ~A, tx.var = "A",
                               treat.range = c(0.0, 1.0))
  A <- dat_cont$A
  result <- obj$prep(dat_cont, A)
  new_col <- setdiff(names(result), names(dat_cont))
  expect_length(new_col, 1L)
  expect_equal(result[[new_col]], A^2)
})

test_that("ContQuadraticBlip$prep: original columns are unchanged", {
  obj <- ContQuadraticBlip$new(tf.model = ~X, blip.model = ~A, tx.var = "A",
                               treat.range = c(0.0, 1.0))
  result <- obj$prep(dat_cont, dat_cont$A)
  expect_equal(result[, names(dat_cont)], dat_cont)
})

test_that("ContQuadraticBlip$Hd: returns matrix with n rows and tx.var column", {
  obj <- ContQuadraticBlip$new(tf.model = ~X, blip.model = ~A, tx.var = "A",
                               treat.range = c(0.0, 1.0))
  mm <- obj$Hd(dat_cont, dat_cont$A)
  expect_true(is.matrix(mm))
  expect_equal(nrow(mm), n_obs)
  expect_true("A" %in% colnames(mm))
})

test_that("ContQuadraticBlip$Hw: returns matrix with same dimensions as Hd", {
  obj <- ContQuadraticBlip$new(tf.model = ~X, blip.model = ~A, tx.var = "A",
                                treat.range = c(0.0, 1.0))
  fit <- glm(A ~ X, dat_cont, family = gaussian())
  Ahat <- fit$fitted.values
  wgt <- rep(1.0, n_obs)
  expect_equal(dim(obj$Hw(dat_cont, dat_cont$A, Ahat, wgt, treat.mod.fitted = fit)),
               dim(obj$Hd(dat_cont, dat_cont$A)))
})

local({
  obj <- ContQuadraticBlip$new(tf.model = ~X, blip.model = ~A, tx.var = "A",
                               treat.range = c(0.0, 1.0))
  dat <- dat_cont
  prepped <- obj$prep(dat, dat$A)
  prepped$Y_int <- prepped$Y
  fit <- lm(update(obj$full.model, Y_int ~ .), prepped)
  
  test_that("ContQuadraticBlip$regret: returns numeric vector of length n", {
    res <- obj$regret(fit, dat, opt = dat$A, A = dat$A)
    expect_true(is.numeric(res))
    expect_length(res, n_obs)
  })
  
  test_that("ContQuadraticBlip$regret: is zero when opt equals observed treatment", {
    res <- obj$regret(fit, dat, opt = dat$A, A = dat$A)
    expect_equal(unname(res), rep(0.0, n_obs))
  })
  
  test_that("ContQuadraticBlip$opt: result is within treat.range for all observations", {
    res <- obj$opt(fit, dat)
    expect_true(is.numeric(res))
    expect_length(res, n_obs)
    expect_true(all(res >= 0.0 & res <= 1.0))
  })
  
  test_that("ContQuadraticBlip$opt: no warning issued when quiet = TRUE", {
    obj_inf <- ContQuadraticBlip$new(tf.model = ~X, blip.model = ~A, tx.var = "A",
                                     treat.range = c(-Inf, Inf))
    expect_no_warning(obj_inf$opt(fit, dat, quiet = TRUE))
  })
  
  test_that("ContQuadraticBlip$shiftY: 'DTR' passes opt through to regret", {
    opt_vec <- rep(0.5, n_obs)
    expect_equal(obj$shiftY("DTR", fit, dat, opt = opt_vec, A = dat$A),
                 obj$regret(fit, dat, opt = opt_vec, A = dat$A))
  })
  
  test_that("ContQuadraticBlip$shiftY: non-DTR type uses 0.0 as reference treatment", {
    expect_equal(obj$shiftY("effect", fit, dat, opt = rep(0.5, n_obs), A = dat$A),
                 obj$regret(fit, dat, opt = 0.0, A = dat$A))
  })
})

local({
  # Y = A^2 is perfectly U-shaped: quadratic coefficient >> 0, guaranteeing
  # y_quad > 0 and the warning branch fires when treat.range is infinite
  A_seq <- seq(0.1, 0.9, length.out = n_obs)
  dat_cup <- data.frame(X = rep(0.0, n_obs), A = A_seq, Y = A_seq^2)
  obj_inf <- ContQuadraticBlip$new(tf.model = ~X, blip.model = ~A, tx.var = "A",
                                   treat.range = c(-Inf, Inf))
  prepped <- obj_inf$prep(dat_cup, dat_cup$A)
  prepped$Y_int <- prepped$Y
  fit_cup <- lm(update(obj_inf$full.model, Y_int ~ .), prepped)
  
  test_that("ContQuadraticBlip$opt: warns about potential minimum when quiet = FALSE", {
    expect_warning(obj_inf$opt(fit_cup, dat_cup, quiet = FALSE),
                   "may be a minimum")
  })
})