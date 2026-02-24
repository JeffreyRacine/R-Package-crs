test_that("snomadr reports missing required eval.f arguments from ...", {
  eval.f <- function(x, shift) sum(x) + shift

  expect_error(
    crs::snomadr(
      n = 2,
      eval.f = eval.f,
      x0 = c(0, 0),
      lb = c(-1, -1),
      ub = c(1, 1),
      display.nomad.progress = FALSE
    ),
    "requires argument 'shift'"
  )
})

test_that("snomadr reports extraneous ... arguments", {
  eval.f <- function(x, shift) sum(x) + shift

  expect_error(
    crs::snomadr(
      n = 2,
      eval.f = eval.f,
      x0 = c(0, 0),
      lb = c(-1, -1),
      ub = c(1, 1),
      shift = 0,
      extra_arg = 123,
      display.nomad.progress = FALSE
    ),
    "'extra_arg' passed to \\(\\.\\.\\.\\) in 'snomadr'"
  )
})

test_that("snomadr rejects ... arguments when eval.f only takes x", {
  eval.f <- function(x) sum(x)

  expect_error(
    crs::snomadr(
      n = 2,
      eval.f = eval.f,
      x0 = c(0, 0),
      lb = c(-1, -1),
      ub = c(1, 1),
      extra_arg = 123,
      display.nomad.progress = FALSE
    ),
    "'extra_arg' passed to \\(\\.\\.\\.\\) in 'snomadr'"
  )
})

test_that("snomadr does not require defaulted eval.f arguments", {
  eval.f <- function(x, shift = 0) sum(x) + shift

  expect_error(
    crs::snomadr(
      n = 2,
      eval.f = eval.f,
      x0 = c(0, 0),
      lb = c(-1, -1),
      ub = c(1, 1),
      extra_arg = 123,
      display.nomad.progress = FALSE
    ),
    "'extra_arg' passed to \\(\\.\\.\\.\\) in 'snomadr'"
  )
})

test_that("snomadr accepts extra arguments when eval.f has ellipsis", {
  eval.f <- function(x, shift, ...) sum(x) + shift

  expect_error(
    crs::snomadr(
      n = 2,
      eval.f = eval.f,
      x0 = c(0, 0),
      lb = c(-1, -1),
      ub = c(1, 1),
      extra_arg = 123,
      display.nomad.progress = FALSE
    ),
    "requires argument 'shift'"
  )
})
