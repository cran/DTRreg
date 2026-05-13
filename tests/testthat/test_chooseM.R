data(twoStageCont)
n <- nrow(twoStageCont)

blip.mod <- list(~X1, ~X2)
treat.mod <- list(A1 ~ X1, A2 ~ X2)
tf.mod <- list(~X1, ~X2)


# Output structure


local({
  skip_on_cran()
  
  withr::with_seed(42L, {
    res <- suppressMessages(
      chooseM(twoStageCont$Y, blip.mod, treat.mod, tf.mod,
              data = twoStageCont,
              method = "dwols",
              B1 = 3L,
              B2 = 3L)
    )
  })
  
  test_that("chooseM: returns a list", {
    expect_true(is.list(res))
  })
  
  test_that("chooseM: returns list with single element 'm'", {
    expect_named(res, "m")
  })
  
  test_that("chooseM: m is numeric", {
    expect_true(is.numeric(res$m))
  })
  
  test_that("chooseM: m is between 1 and n", {
    expect_true(res$m >= 1.0)
    expect_true(res$m <= n)
  })
  
  test_that("chooseM: m is a whole number", {
    expect_equal(res$m, floor(res$m))
  })
})



# m formula bounds


test_that("chooseM: m formula gives n when p=0 and alpha=0.5", {
  skip_on_cran()
  # When p=0 (maximally non-regular): m = n^((1 + 0.5) / (1 + 0.5)) = n^1 = n
  # When p=1 (regular): m = n^((1 + 0) / (1 + 0.5)) = n^(2/3)
  # So m should always be <= n
  withr::with_seed(42L, {
    res <- suppressMessages(
      chooseM(twoStageCont$Y, blip.mod, treat.mod, tf.mod,
              data = twoStageCont,
              method = "dwols",
              B1 = 3L,
              B2 = 3L)
    )
  })
  expect_true(res$m <= n)
})



# Messages produced during run


test_that("chooseM: produces alpha/coverage messages during run", {
  skip_on_cran()
  withr::with_seed(42L,
                   expect_message(
                     chooseM(twoStageCont$Y, blip.mod, treat.mod, tf.mod,
                             data = twoStageCont,
                             method = "dwols",
                             B1 = 3L,
                             B2 = 3L),
                     "alpha"
                   )
  )
})



# Stopping conditions


test_that("chooseM: stops when outcome is missing", {
  skip_on_cran()
  expect_error(
    chooseM(blip.mod = blip.mod, treat.mod = treat.mod, tf.mod = tf.mod,
            data = twoStageCont, B1 = 3L, B2 = 3L),
    "must be provided"
  )
})

test_that("chooseM: stops when blip.mod is missing", {
  skip_on_cran()
  expect_error(
    chooseM(twoStageCont$Y, treat.mod = treat.mod, tf.mod = tf.mod,
            data = twoStageCont, B1 = 3L, B2 = 3L),
    "must be provided"
  )
})

test_that("chooseM: stops when K != 2", {
  skip_on_cran()
  # Single stage - K=1
  expect_error(
    suppressMessages(
      chooseM(twoStageCont$Y,
              blip.mod = list(~X1),
              treat.mod = list(A1 ~ X1),
              tf.mod = list(~X1),
              data = twoStageCont,
              B1 = 3L,
              B2 = 3L)
    ),
    "K = 2"
  )
})



# Reproducibility


test_that("chooseM: same seed produces same m", {
  skip_on_cran()
  withr::with_seed(42L,
                   res1 <- suppressMessages(
                     chooseM(twoStageCont$Y, blip.mod, treat.mod, tf.mod,
                             data = twoStageCont,
                             method = "dwols",
                             B1 = 3L,
                             B2 = 3L)
                   )
  )
  withr::with_seed(42L,
                   res2 <- suppressMessages(
                     chooseM(twoStageCont$Y, blip.mod, treat.mod, tf.mod,
                             data = twoStageCont,
                             method = "dwols",
                             B1 = 3L,
                             B2 = 3L)
                   )
  )
  expect_equal(res1$m, res2$m)
})