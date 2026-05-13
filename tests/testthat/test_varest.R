set.seed(42L)
n_v <- 30L

dat_v <- data.frame(
  X = cos(seq(0.0, pi, length.out = n_v))
)

# Two-stage models
models_v <- list(
  list(blip = ~X),
  list(blip = ~X)
)

# Known psi and covmat for two stages
psi_v <- list(
  c("A1" = 1.0, "A1:X" = 0.5),
  c("A2" = 2.0, "A2:X" = -1.0)
)

# Tight covariance - narrow CIs, few participants straddle zero
covmat_tight <- list(
  diag(c(0.001, 0.001)),
  diag(c(0.001, 0.001))
)

# Wide covariance - CIs so wide everything straddles zero
covmat_wide <- list(
  diag(c(1e6, 1e6)),
  diag(c(1e6, 1e6))
)

make_obj_v <- function(covmat = covmat_tight,
                       psi = psi_v,
                       tx.type = "bin",
                       K = 2L) {
  list(
    K = K,
    data = dat_v,
    tx.type = tx.type,
    models = models_v[seq_len(K)],
    covmat = covmat[seq_len(K)],
    psi = psi[seq_len(K)]
  )
}

test_that(".varest: returns NA when covmat is NULL", {
  obj <- make_obj_v()
  obj$covmat <- NULL
  expect_true(is.na(.varest(obj)))
})

test_that(".varest: returns NA when tx.type is not 'bin'", {
  expect_true(is.na(.varest(make_obj_v(tx.type = "cont"))))
  expect_true(is.na(.varest(make_obj_v(tx.type = "multi"))))
})

test_that(".varest: returns numeric vector of length K", {
  res <- .varest(make_obj_v())
  expect_true(is.numeric(res))
  expect_length(res, 2L)
})

test_that(".varest: returns numeric vector of length 1 for K=1", {
  res <- .varest(make_obj_v(
    covmat = covmat_tight[1L],
    psi = psi_v[1L],
    K = 1L
  ))
  expect_true(is.numeric(res))
  expect_length(res, 1L)
})

test_that(".varest: all values are in [0, 1]", {
  res <- .varest(make_obj_v())
  expect_true(all(res >= 0.0 & res <= 1.0))
})

test_that(".varest: tight CIs produce near-zero non-regularity", {
  # With very small variance, CIs are narrow and unlikely to straddle zero
  # when psi is clearly non-zero - nonreg should be 0 or very small
  res <- .varest(make_obj_v(covmat = covmat_tight))
  expect_true(all(res < 0.1))
})

test_that(".varest: wide CIs produce non-regularity of 1", {
  # With enormous variance, CIs always span zero - every participant
  # is non-regular
  res <- .varest(make_obj_v(covmat = covmat_wide))
  expect_true(all(res == 1.0))
})

test_that(".varest: psi of all zeros with any variance produces nonreg of 1", {
  # If psi = 0, CI always straddles 0 regardless of variance
  psi_zero <- list(c("A1" = 0.0, "A1:X" = 0.0),
                   c("A2" = 0.0, "A2:X" = 0.0))
  res <- .varest(make_obj_v(covmat = covmat_tight, psi = psi_zero))
  expect_true(all(res == 1.0))
})

local({
  # Oracle: explicit .each_stage computation for stage 1
  .each_stage_oracle <- function(models, data, covmat, psi) {
    H_psi <- model.matrix(models$blip, data)
    H_psi <- H_psi[complete.cases(H_psi), , drop = FALSE]
    tmp <- sqrt(diag(covmat))
    psi_l <- psi - 1.96 * tmp
    psi_u <- psi + 1.96 * tmp
    H_psi_p <- pmax(H_psi, 0.0)
    H_psi_n <- pmin(H_psi, 0.0)
    blip_l <- H_psi_p %*% psi_l + H_psi_n %*% psi_u
    blip_u <- H_psi_p %*% psi_u + H_psi_n %*% psi_l
    sum(blip_l < 0 & blip_u > 0) / nrow(H_psi)
  }
  
  # Intermediate covariance - some but not all participants non-regular
  covmat_mid <- list(
    diag(c(0.5, 0.5)),
    diag(c(0.5, 0.5))
  )
  
  obj <- make_obj_v(covmat = covmat_mid)
  res <- .varest(obj)
  
  test_that(".varest: stage 1 matches oracle computation", {
    expected <- .each_stage_oracle(models_v[[1L]], dat_v,
                                   covmat_mid[[1L]], psi_v[[1L]])
    expect_equal(res[1L], expected)
  })
  
  test_that(".varest: stage 2 matches oracle computation", {
    expected <- .each_stage_oracle(models_v[[2L]], dat_v,
                                   covmat_mid[[2L]], psi_v[[2L]])
    expect_equal(res[2L], expected)
  })
  
  test_that(".varest: both stages oracle simultaneously", {
    expected <- c(
      .each_stage_oracle(models_v[[1L]], dat_v, covmat_mid[[1L]], psi_v[[1L]]),
      .each_stage_oracle(models_v[[2L]], dat_v, covmat_mid[[2L]], psi_v[[2L]])
    )
    expect_equal(res, expected)
  })
})

test_that(".varest: NA rows in data are excluded from proportion calculation", {
  dat_na <- dat_v
  dat_na$X[5L] <- NA_real_
  
  obj_na <- make_obj_v(K = 1L, covmat = covmat_wide[1L], psi = psi_v[1L])
  obj_na$data <- dat_na
  obj_clean <- make_obj_v(K = 1L, covmat = covmat_wide[1L], psi = psi_v[1L])
  
  res_na <- .varest(obj_na)
  res_clean <- .varest(obj_clean)
  
  # Wide CIs: both should be 1.0 - NA exclusion doesn't change the proportion
  # when all non-NA participants are non-regular
  expect_equal(res_na, 1.0)
  expect_equal(res_clean, 1.0)
})