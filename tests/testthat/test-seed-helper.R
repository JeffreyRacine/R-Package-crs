test_that("seed helper restores pre-existing global RNG state", {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    original_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }

  on.exit({
    if (had_seed) {
      assign(".Random.seed", original_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  set.seed(123)
  baseline_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)

  seed_state <- crs:::.crs_capture_seed()
  set.seed(456)
  crs:::.crs_restore_seed(seed_state)

  expect_equal(get(".Random.seed", envir = .GlobalEnv, inherits = FALSE), baseline_seed)
})

test_that("seed helper is a no-op when no pre-existing seed is present", {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    original_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }

  on.exit({
    if (had_seed) {
      assign(".Random.seed", original_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    rm(".Random.seed", envir = .GlobalEnv)
  }

  seed_state <- crs:::.crs_capture_seed()
  set.seed(789)
  generated_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)

  crs:::.crs_restore_seed(seed_state)

  expect_equal(get(".Random.seed", envir = .GlobalEnv, inherits = FALSE), generated_seed)
})
