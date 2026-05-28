native_nomad_solve <- function(spec) {
  .Call("crs_nomad_native_test_solve", spec, PACKAGE = "crs")
}

native_nomad_spec <- function(...,
                              mode = "c",
                              scenario = "quadratic",
                              x0 = 0,
                              lower = -1,
                              upper = 1,
                              input_type = 0L,
                              output_type = 0L,
                              max_eval = 20L,
                              random_seed = 42L,
                              quiet = TRUE,
                              options = character(),
                              start_count = 0L,
                              read_nomad_opt_file = FALSE) {
  options <- {
    out <- as.character(options)
    names(out) <- names(options)
    out
  }
  list(
    mode = mode,
    scenario = scenario,
    x0 = as.double(x0),
    lower = as.double(lower),
    upper = as.double(upper),
    input_type = as.integer(input_type),
    output_type = as.integer(output_type),
    max_eval = as.integer(max_eval),
    random_seed = as.integer(random_seed),
    quiet = isTRUE(quiet),
    options = options,
    start_count = as.integer(start_count),
    read_nomad_opt_file = isTRUE(read_nomad_opt_file),
    ...
  )
}

test_that("native R callback helper distinguishes success, errors, and interrupts", {
  ok <- crs:::.crs_nomad_native_eval_r(function(x) sum(x), c(1, 2))
  expect_identical(ok[[1]], 0L)
  expect_equal(ok[[2]], 3)

  err <- crs:::.crs_nomad_native_eval_r(function(x) stop("boom"), c(1, 2))
  expect_identical(err[[1]], 2L)
  expect_match(err[[3]], "boom")

  intr <- crs:::.crs_nomad_native_eval_r(
    function(x) {
      stop(structure(
        list(message = "simulated interrupt"),
        class = c("interrupt", "condition")
      ))
    },
    c(1, 2)
  )
  expect_identical(intr[[1]], 1L)
  expect_match(intr[[3]], "simulated interrupt")
})

test_that("crs_nomad_solve succeeds through the direct C callback bridge", {
  res <- native_nomad_solve(native_nomad_spec(
    x0 = c(0, 0),
    lower = c(-1, -1),
    upper = c(1, 1),
    input_type = c(0L, 0L),
    output_type = 0L,
    options = c(MAX_BB_EVAL = "15")
  ))

  expect_equal(res$status, 0L)
  expect_equal(res$result_status, 0L)
  expect_true(is.finite(res$objective))
  expect_length(res$solution, 2L)
  expect_true(res$callback_evaluations >= 1L)
  expect_true(res$total_evaluations >= res$callback_evaluations)
})

test_that("crs_nomad_solve succeeds through the direct R callback bridge", {
  eval_f <- function(x) sum((x - 0.25)^2)

  res <- native_nomad_solve(native_nomad_spec(
    mode = "r",
    eval_f = eval_f,
    x0 = c(0, 0),
    lower = c(-1, -1),
    upper = c(1, 1),
    input_type = c(0L, 0L),
    output_type = 0L,
    options = c(MAX_BB_EVAL = "15", NB_THREADS_PARALLEL_EVAL = "1")
  ))

  expect_equal(res$status, 0L)
  expect_equal(res$result_status, 0L)
  expect_true(is.finite(res$objective))
  expect_length(res$solution, 2L)
  expect_true(res$callback_evaluations >= 1L)
})

test_that("crs_nomad_solve has an explicit native EVAL_USE_CACHE contract", {
  cache_off_message <- "requires EVAL_USE_CACHE = TRUE"
  eval_f <- function(x) sum((x - 0.25)^2)

  c_on <- native_nomad_solve(native_nomad_spec(
    x0 = c(0, 0),
    lower = c(-1, -1),
    upper = c(1, 1),
    input_type = c(0L, 0L),
    output_type = 0L,
    options = c(MAX_BB_EVAL = "15", EVAL_USE_CACHE = "true")
  ))
  expect_equal(c_on$status, 0L)
  expect_equal(c_on$result_status, 0L)
  expect_true(c_on$callback_evaluations >= 1L)

  c_on_upper <- native_nomad_solve(native_nomad_spec(
    x0 = c(0, 0),
    lower = c(-1, -1),
    upper = c(1, 1),
    input_type = c(0L, 0L),
    output_type = 0L,
    options = c(MAX_BB_EVAL = "15", EVAL_USE_CACHE = "TRUE")
  ))
  expect_equal(c_on_upper$status, 0L)
  expect_equal(c_on_upper$result_status, 0L)

  c_off <- native_nomad_solve(native_nomad_spec(
    x0 = c(0, 0),
    lower = c(-1, -1),
    upper = c(1, 1),
    input_type = c(0L, 0L),
    output_type = 0L,
    options = c(MAX_BB_EVAL = "15", EVAL_USE_CACHE = "false")
  ))
  expect_equal(c_off$status, 1L)
  expect_equal(c_off$result_status, 1L)
  expect_match(c_off$message, cache_off_message)
  expect_equal(c_off$callback_evaluations, 0L)

  r_on <- native_nomad_solve(native_nomad_spec(
    mode = "r",
    eval_f = eval_f,
    x0 = c(0, 0),
    lower = c(-1, -1),
    upper = c(1, 1),
    input_type = c(0L, 0L),
    output_type = 0L,
    options = c(MAX_BB_EVAL = "15", NB_THREADS_PARALLEL_EVAL = "1",
                EVAL_USE_CACHE = "true")
  ))
  expect_equal(r_on$status, 0L)
  expect_equal(r_on$result_status, 0L)
  expect_true(r_on$callback_evaluations >= 1L)

  r_off <- native_nomad_solve(native_nomad_spec(
    mode = "r",
    eval_f = eval_f,
    x0 = c(0, 0),
    lower = c(-1, -1),
    upper = c(1, 1),
    input_type = c(0L, 0L),
    output_type = 0L,
    options = c(MAX_BB_EVAL = "15", NB_THREADS_PARALLEL_EVAL = "1",
                EVAL_USE_CACHE = "FALSE")
  ))
  expect_equal(r_off$status, 1L)
  expect_equal(r_off$result_status, 1L)
  expect_match(r_off$message, cache_off_message)
  expect_equal(r_off$callback_evaluations, 0L)

  oldwd <- setwd(tempdir())
  on.exit(setwd(oldwd), add = TRUE)
  writeLines("EVAL_USE_CACHE false", "nomad.opt")
  file_off <- native_nomad_solve(native_nomad_spec(
    x0 = c(0, 0),
    lower = c(-1, -1),
    upper = c(1, 1),
    input_type = c(0L, 0L),
    output_type = 0L,
    options = c(MAX_BB_EVAL = "15"),
    read_nomad_opt_file = TRUE
  ))
  unlink("nomad.opt")
  expect_equal(file_off$status, 1L)
  expect_equal(file_off$result_status, 1L)
  expect_match(file_off$message, cache_off_message)
  expect_equal(file_off$callback_evaluations, 0L)
})

test_that("R callback errors do not poison the next native solve", {
  bad_eval <- function(x) stop("intentional native R callback failure")
  good_eval <- function(x) sum((x - 0.25)^2)

  bad <- native_nomad_solve(native_nomad_spec(
    mode = "r",
    eval_f = bad_eval,
    x0 = 0,
    lower = -1,
    upper = 1,
    input_type = 0L,
    output_type = 0L,
    max_eval = 1L,
    options = c(MAX_BB_EVAL = "1")
  ))

  expect_equal(bad$status, 2L)
  expect_equal(bad$result_status, 2L)
  expect_match(bad$message, "callback", ignore.case = TRUE)

  good <- native_nomad_solve(native_nomad_spec(
    mode = "r",
    eval_f = good_eval,
    x0 = 0,
    lower = -1,
    upper = 1,
    input_type = 0L,
    output_type = 0L,
    options = c(MAX_BB_EVAL = "10")
  ))

  expect_equal(good$status, 0L)
  expect_equal(good$result_status, 0L)
  expect_true(is.finite(good$objective))
})

test_that("R callbacks reject unsafe NOMAD parallel-evaluation options", {
  eval_f <- function(x) sum((x - 0.25)^2)

  rejected <- native_nomad_solve(native_nomad_spec(
    mode = "r",
    eval_f = eval_f,
    x0 = 0,
    lower = -1,
    upper = 1,
    input_type = 0L,
    output_type = 0L,
    options = c(NB_THREADS_PARALLEL_EVAL = "2")
  ))

  expect_equal(rejected$status, 1L)
  expect_equal(rejected$result_status, 1L)
  expect_match(rejected$message, "NB_THREADS_PARALLEL_EVAL > 1")

  accepted <- native_nomad_solve(native_nomad_spec(
    mode = "r",
    eval_f = eval_f,
    x0 = 0,
    lower = -1,
    upper = 1,
    input_type = 0L,
    output_type = 0L,
    options = c(NB_THREADS_PARALLEL_EVAL = "1", MAX_BB_EVAL = "10")
  ))

  expect_equal(accepted$status, 0L)
  expect_equal(accepted$result_status, 0L)
})

test_that("nomad.opt input is trimmed and screened for R callback threading", {
  eval_f <- function(x) sum((x - 0.25)^2)
  old <- getwd()
  dir <- tempfile("crs-native-nomad-opt-")
  dir.create(dir)
  on.exit(setwd(old), add = TRUE)
  setwd(dir)

  writeLines("   MAX_BB_EVAL 10   # inline comment", "nomad.opt")
  ok <- native_nomad_solve(native_nomad_spec(
    mode = "r",
    eval_f = eval_f,
    x0 = 0,
    lower = -1,
    upper = 1,
    input_type = 0L,
    output_type = 0L,
    read_nomad_opt_file = TRUE
  ))
  expect_equal(ok$status, 0L)

  writeLines("   NB_THREADS_PARALLEL_EVAL 2", "nomad.opt")
  rejected <- native_nomad_solve(native_nomad_spec(
    mode = "r",
    eval_f = eval_f,
    x0 = 0,
    lower = -1,
    upper = 1,
    input_type = 0L,
    output_type = 0L,
    read_nomad_opt_file = TRUE
  ))
  expect_equal(rejected$status, 1L)
  expect_match(rejected$message, "NB_THREADS_PARALLEL_EVAL > 1")
})

test_that("callback output validation matches snomadr Inf/NaN contracts", {
  nan_obj <- native_nomad_solve(native_nomad_spec(
    scenario = "nan_obj",
    output_type = 0L,
    max_eval = 1L,
    options = c(MAX_BB_EVAL = "1")
  ))
  expect_equal(nan_obj$status, 2L)
  expect_match(nan_obj$message, "callback", ignore.case = TRUE)

  inf_obj <- native_nomad_solve(native_nomad_spec(
    scenario = "inf_obj",
    output_type = 0L,
    max_eval = 1L,
    options = c(MAX_BB_EVAL = "1")
  ))
  expect_equal(inf_obj$status, 0L)
  expect_true(is.finite(inf_obj$objective))
  expect_equal(inf_obj$objective, .Machine$double.xmax)

  negative_inf_obj <- native_nomad_solve(native_nomad_spec(
    scenario = "negative_inf_obj",
    output_type = 0L,
    max_eval = 1L,
    options = c(MAX_BB_EVAL = "1")
  ))
  expect_equal(negative_inf_obj$status, 0L)
  expect_true(is.finite(negative_inf_obj$objective))
  expect_equal(negative_inf_obj$objective, -.Machine$double.xmax)

  r_inf_obj <- native_nomad_solve(native_nomad_spec(
    mode = "r",
    eval_f = function(x) Inf,
    x0 = 0,
    lower = -1,
    upper = 1,
    input_type = 0L,
    output_type = 0L,
    max_eval = 1L,
    options = c(MAX_BB_EVAL = "1")
  ))
  expect_equal(r_inf_obj$status, 0L)
  expect_true(is.finite(r_inf_obj$objective))
  expect_equal(r_inf_obj$objective, .Machine$double.xmax)

  inf_pb <- native_nomad_solve(native_nomad_spec(
    scenario = "inf_constraint",
    output_type = c(0L, 1L),
    options = c(MAX_BB_EVAL = "5")
  ))
  expect_false(grepl("native callback reported evaluation failure",
                     inf_pb$message,
                     fixed = TRUE))

  inf_eb <- native_nomad_solve(native_nomad_spec(
    scenario = "inf_constraint",
    output_type = c(0L, 2L),
    options = c(MAX_BB_EVAL = "5")
  ))
  expect_false(grepl("native callback reported evaluation failure",
                     inf_eb$message,
                     fixed = TRUE))
})

test_that("sequential native solves do not reuse stale objective values", {
  first <- native_nomad_solve(native_nomad_spec(
    scenario = "quadratic",
    x0 = 0,
    lower = -1,
    upper = 1,
    input_type = 0L,
    output_type = 0L,
    options = c(MAX_BB_EVAL = "10"),
    target = 0.25
  ))
  second <- native_nomad_solve(native_nomad_spec(
    scenario = "shifted",
    x0 = 0,
    lower = -1,
    upper = 1,
    input_type = 0L,
    output_type = 0L,
    options = c(MAX_BB_EVAL = "10"),
    target = 0.25
  ))

  expect_equal(first$status, 0L)
  expect_equal(second$status, 0L)
  expect_true(first$callback_evaluations >= 1L)
  expect_true(second$callback_evaluations >= 1L)
  expect_true(second$objective > first$objective + 1)
})

test_that("direct NOMAD start capture preserves and respects RNG state", {
  set.seed(314159)
  before <- get(".Random.seed", envir = .GlobalEnv)
  starts <- crs:::.crs_nomad_capture_start_matrix(
    x0 = c(0.5, 1),
    nstart = 3L,
    bbin = c(0L, 0L),
    bbout = 0L,
    lb = c(0.1, 0.1),
    ub = c(2, 2),
    random.seed = 42L,
    opts = list()
  )
  after <- get(".Random.seed", envir = .GlobalEnv)

  expect_identical(after, before)
  expect_equal(dim(starts), c(3L, 2L))

  first <- crs:::.crs_nomad_capture_start_matrix(
    x0 = c(0.5, 1),
    nstart = 3L,
    bbin = c(0L, 0L),
    bbout = 0L,
    lb = c(0.1, 0.1),
    ub = c(2, 2),
    random.seed = 99L,
    opts = list()
  )
  second <- crs:::.crs_nomad_capture_start_matrix(
    x0 = c(0.5, 1),
    nstart = 3L,
    bbin = c(0L, 0L),
    bbout = 0L,
    lb = c(0.1, 0.1),
    ub = c(2, 2),
    random.seed = 99L,
    opts = list()
  )

  expect_equal(first, second, tolerance = 0)
})
