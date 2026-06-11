test_that("CRS plot controls validate their inputs", {
  boot <- getFromNamespace("np_boot_control", "crs")(
    nonfixed = "frozen",
    wild = "mammen",
    blocklen = 3
  )
  grid <- getFromNamespace("np_grid_control", "crs")(
    xtrim = c(0.1, 0.9),
    xq = 0.25
  )
  render <- getFromNamespace("np_render_control", "crs")(
    style = "bar",
    bar = "I",
    bar_num = 10
  )
  crs.boot <- getFromNamespace(".crs_boot_control", "crs")()

  expect_s3_class(boot, "np_boot_control")
  expect_s3_class(grid, "np_grid_control")
  expect_s3_class(render, "np_render_control")
  expect_s3_class(crs.boot, "crs_boot_control")
  expect_s3_class(crs.boot, "np_boot_control")
  expect_error(getFromNamespace("np_boot_control", "crs")(blocklen = 0),
               "positive numeric scalar")
  expect_error(getFromNamespace("np_grid_control", "crs")(xtrim = c(0.9, 0.1)),
               "xtrim")
  expect_error(getFromNamespace("np_render_control", "crs")(bar_num = 0),
               "positive numeric scalar")
})

test_that("CRS plot public dot validation accepts canonical and prefixed args", {
  validate <- getFromNamespace(".crs_plot_validate_public_dots", "crs")

  expect_true(validate(list(output = "data",
                            errors = "bootstrap",
                            plot.col = "red",
                            persp.theta = 30,
                            rgl.persp3d.alpha = 0.5)))
  expect_error(validate(list(neva = 10)), "unused plot argument")
})

test_that("CRS plot public dot normalization maps np-style aliases", {
  normalize <- getFromNamespace(".crs_plot_normalize_public_dots", "crs")
  boot <- getFromNamespace("np_boot_control", "crs")(blocklen = 4)
  grid <- getFromNamespace("np_grid_control", "crs")(xtrim = c(0.05, 0.95),
                                                       xq = 0.4)
  render <- getFromNamespace("np_render_control", "crs")(style = "bar",
                                                           bar = "I",
                                                           bar_num = 8)

  dots <- normalize(list(errors = "bootstrap",
                         band = "pointwise",
                         alpha = 0.1,
                         B = 7,
                         output = "data",
                         data_overlay = TRUE,
                         data_rug = FALSE,
                         layout = "current",
                         gradient = TRUE,
                         neval = 12,
                         renderer = "base",
                         boot_control = boot,
                         grid_control = grid,
                         render_control = render),
                    context = "plot.crs")

  expect_identical(dots$plot.errors.method, "bootstrap")
  expect_identical(dots$plot.errors.type, "pointwise")
  expect_equal(dots$plot.errors.alpha, 0.1)
  expect_identical(dots$plot.errors.boot.num, 7L)
  expect_identical(dots$plot.behavior, "data")
  expect_true(dots$plot.data.overlay)
  expect_false(dots$plot.rug)
  expect_false(dots$plot.par.mfrow)
  expect_true(dots$gradients)
  expect_identical(dots$num.eval, 12L)
  expect_identical(dots$renderer, "base")
  expect_identical(dots$plot.errors.boot.nonfixed, "exact")
  expect_identical(dots$plot.errors.boot.wild, "rademacher")
  expect_equal(dots$plot.errors.boot.blocklen, 4)
  expect_equal(dots$xtrim, c(0.05, 0.95))
  expect_equal(dots$xq, 0.4)
  expect_identical(dots$plot.errors.style, "bar")
  expect_identical(dots$plot.errors.bar, "I")
  expect_equal(dots$plot.errors.bar.num, 8)
})

test_that("CRS plot normalization rejects ambiguous or inconsistent controls", {
  normalize <- getFromNamespace(".crs_plot_normalize_public_dots", "crs")

  expect_error(normalize(list(output = "data", plot.behavior = "plot")),
               "conflicting plot arguments")
  expect_error(normalize(list(B = 9)),
               "bootstrap controls require errors")
  expect_error(normalize(list(alpha = 1.5)), "alpha")
  expect_error(normalize(list(data_overlay = NA)), "TRUE or FALSE")
  expect_error(normalize(list(intervals = TRUE)), "did you mean errors")
})

test_that("CRS plot colors and user arg helpers are deterministic", {
  color <- getFromNamespace(".crs_plot_color", "crs")
  surface_colors <- getFromNamespace(".crs_plot_surface_colors", "crs")
  user_args <- getFromNamespace(".crs_plot_user_args", "crs")
  merge_args <- getFromNamespace(".crs_plot_merge_user_args", "crs")

  expect_match(color("fit"), "^#")
  expect_match(color("fit", alpha = 0.5), "^#")
  expect_length(surface_colors(matrix(1:4, 2, 2), num.colors = 5), 4)
  expect_equal(user_args(list(plot = list(col = "red"),
                              plot.lwd = 2,
                              lines.col = "blue"),
                         "plot"),
               list(col = "red", lwd = 2))
  expect_equal(merge_args(list(col = "black", lwd = 1),
                          list(col = "red")),
               list(col = "red", lwd = 1))
})
