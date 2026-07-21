## Public-interface support shared by crsiv() and crsivderiv().

.crs_iv_pipe_parts <- function(expr) {
  if (is.call(expr) && identical(expr[[1L]], as.name("|"))) {
    return(c(.crs_iv_pipe_parts(expr[[2L]]), list(expr[[3L]])))
  }
  list(expr)
}

.crs_iv_expr_contains_symbol <- function(expr, symbol) {
  if (is.symbol(expr)) {
    return(identical(as.character(expr), symbol))
  }
  if (!is.call(expr)) {
    return(FALSE)
  }
  any(vapply(as.list(expr)[-1L], .crs_iv_expr_contains_symbol,
             logical(1L), symbol = symbol))
}

.crs_iv_add_exprs <- function(exprs) {
  if (length(exprs) == 1L) {
    return(exprs[[1L]])
  }
  Reduce(function(lhs, rhs) call("+", lhs, rhs), exprs)
}

.crs_iv_make_formula <- function(lhs = NULL, rhs, env) {
  ans <- if (is.null(lhs)) {
    as.call(list(as.name("~"), rhs))
  } else {
    as.call(list(as.name("~"), lhs, rhs))
  }
  class(ans) <- "formula"
  environment(ans) <- env
  ans
}

.crs_iv_term_sources <- function(labels) {
  lapply(labels, function(label) {
    tryCatch(all.vars(str2lang(label)), error = function(e) character())
  })
}

.crs_iv_parse_formula <- function(formula, where) {
  if (!inherits(formula, "formula") || length(formula) != 3L) {
    stop(sprintf("%s requires a two-sided formula", where), call. = FALSE)
  }

  env <- environment(formula)
  if (is.null(env)) {
    env <- parent.frame()
  }
  parts <- .crs_iv_pipe_parts(formula[[3L]])
  if (!length(parts) %in% c(2L, 3L)) {
    stop(sprintf(
      "%s formula must be 'y ~ z | w' or 'y ~ z | w | x'",
      where
    ), call. = FALSE)
  }

  role.names <- c("z", "w", "x")[seq_along(parts)]
  roles <- setNames(vector("list", length(parts)), role.names)
  for (i in seq_along(parts)) {
    expr <- parts[[i]]
    if (.crs_iv_expr_contains_symbol(expr, ".")) {
      stop(sprintf(
        "%s does not support '.' in IV formula partitions; name each role explicitly",
        where
      ), call. = FALSE)
    }
    role.formula <- .crs_iv_make_formula(rhs = expr, env = env)
    role.terms <- terms(role.formula)
    labels <- attr(role.terms, "term.labels")
    if (length(attr(role.terms, "offset"))) {
      stop(sprintf("%s does not support offset() in IV formula partitions",
                   where), call. = FALSE)
    }
    if (!length(labels)) {
      stop(sprintf("%s formula partition '%s' is empty", where,
                   role.names[[i]]), call. = FALSE)
    }
    if (any(attr(role.terms, "order") > 1L)) {
      stop(sprintf(
        "%s does not support interaction operators in IV formula partitions; use an explicit transformation such as I(a * b)",
        where
      ), call. = FALSE)
    }
    roles[[i]] <- list(
      expr = expr,
      formula = role.formula,
      terms = role.terms,
      labels = labels,
      source_variables = .crs_iv_term_sources(labels)
    )
  }

  response.vars <- all.vars(formula[[2L]])
  if (length(response.vars) != 1L) {
    stop(sprintf("%s requires one response variable", where), call. = FALSE)
  }
  combined <- .crs_iv_make_formula(
    lhs = formula[[2L]],
    rhs = .crs_iv_add_exprs(lapply(roles, `[[`, "expr")),
    env = env
  )

  list(
    formula = formula,
    combined = combined,
    roles = roles,
    response = response.vars[[1L]],
    environment = env
  )
}

.crs_iv_validate_data_columns <- function(data, parsed, where) {
  if (is.null(data)) {
    return(invisible(NULL))
  }
  required <- unique(all.vars(parsed[["combined"]]))
  available <- if (is.environment(data)) {
    ls(data, all.names = TRUE)
  } else {
    names(data)
  }
  if (is.null(available)) {
    stop(sprintf("%s 'data' must be a data frame, list, or environment",
                 where), call. = FALSE)
  }
  missing.columns <- setdiff(required, available)
  if (length(missing.columns)) {
    stop(sprintf(
      "%s 'data' is missing required variable%s: %s",
      where,
      if (length(missing.columns) == 1L) "" else "s",
      paste(shQuote(missing.columns), collapse = ", ")
    ), call. = FALSE)
  }
  invisible(NULL)
}

.crs_iv_dot_expressions <- function(method.call, where) {
  dots <- method.call[["..."]]
  if (is.null(dots)) {
    return(list())
  }
  dots <- as.list(dots)
  dot.names <- names(dots)
  if (is.null(dot.names) || any(!nzchar(dot.names))) {
    stop(sprintf("%s requires all arguments in '...' to be named", where),
         call. = FALSE)
  }
  duplicated.names <- unique(dot.names[duplicated(dot.names)])
  if (length(duplicated.names)) {
    stop(sprintf("%s received duplicated argument%s in '...': %s",
                 where,
                 if (length(duplicated.names) == 1L) "" else "s",
                 paste(shQuote(duplicated.names), collapse = ", ")),
         call. = FALSE)
  }
  dots
}

.crs_iv_reject_formula_evaluation <- function(dot.expressions, where) {
  evaluation.args <- intersect(
    names(dot.expressions), c("newdata", "zeval", "weval", "xeval")
  )
  if (length(evaluation.args)) {
    stop(sprintf(
      "%s formula interface currently supports training-row evaluation only; use the native vector interface for %s",
      where, paste(shQuote(evaluation.args), collapse = ", ")
    ), call. = FALSE)
  }
  invisible(NULL)
}

.crs_iv_training_frame <- function(parsed, method.call, data, where,
                                   eval.env, dot.expressions) {
  .crs_iv_validate_data_columns(data, parsed, where)

  mf.call <- list(quote(stats::model.frame), formula = parsed[["combined"]])
  if (!is.null(method.call[["data"]])) {
    mf.call[["data"]] <- method.call[["data"]]
  }
  if (!is.null(method.call[["subset"]])) {
    mf.call[["subset"]] <- method.call[["subset"]]
  }
  if (!is.null(method.call[["na.action"]])) {
    mf.call[["na.action"]] <- method.call[["na.action"]]
  }
  for (arg in c("weights", "starting.values")) {
    if (!is.null(dot.expressions[[arg]])) {
      mf.call[[arg]] <- dot.expressions[[arg]]
    }
  }
  mf.call[["drop.unused.levels"]] <- TRUE
  frame <- eval(as.call(mf.call), envir = eval.env)

  omit <- attr(frame, "na.action")
  n.input <- NROW(frame) + length(omit)
  omitted <- as.integer(omit)
  retained <- if (length(omitted)) {
    setdiff(seq_len(n.input), omitted)
  } else {
    seq_len(n.input)
  }

  list(
    frame = frame,
    weights = model.weights(frame),
    starting.values = model.extract(frame, "starting.values"),
    rows = list(
      n_input = n.input,
      retained = retained,
      omitted = omitted,
      na_action = omit
    )
  )
}

.crs_iv_role_frame <- function(frame, role, role.name, where) {
  labels <- role[["labels"]]
  missing.labels <- setdiff(labels, names(frame))
  if (length(missing.labels)) {
    stop(sprintf(
      "%s could not construct IV formula term%s: %s",
      where,
      if (length(missing.labels) == 1L) "" else "s",
      paste(shQuote(missing.labels), collapse = ", ")
    ), call. = FALSE)
  }

  ans <- frame[, labels, drop = FALSE]
  for (j in seq_along(ans)) {
    if (is.matrix(ans[[j]])) {
      stop(sprintf(
        "%s formula term %s produces a matrix; supply its columns as explicit variables",
        where, shQuote(labels[[j]])
      ), call. = FALSE)
    }
    if (inherits(ans[[j]], "AsIs")) {
      cls <- setdiff(class(ans[[j]]), "AsIs")
      if (length(cls)) {
        class(ans[[j]]) <- cls
      } else {
        class(ans[[j]]) <- NULL
      }
    }
  }
  internal.names <- sprintf("crsiv_%s_%03d", role.name, seq_along(ans))
  names(ans) <- internal.names

  list(
    data = ans,
    metadata = list(
      labels = labels,
      internal_names = internal.names,
      source_variables = role[["source_variables"]],
      classes = lapply(ans, class),
      xlevels = lapply(ans, function(x) if (is.factor(x)) levels(x) else NULL)
    )
  )
}

.crs_iv_compile_training <- function(parsed, training, where) {
  frame <- training[["frame"]]
  response <- model.response(frame)
  if (is.matrix(response) || NCOL(response) != 1L) {
    stop(sprintf("%s response must be univariate", where), call. = FALSE)
  }
  response <- as.vector(response)

  role.names <- names(parsed[["roles"]])
  roles <- setNames(vector("list", length(role.names)), role.names)
  metadata <- setNames(vector("list", length(role.names)), role.names)
  for (role.name in role.names) {
    compiled <- .crs_iv_role_frame(
      frame = frame,
      role = parsed[["roles"]][[role.name]],
      role.name = role.name,
      where = where
    )
    roles[[role.name]] <- compiled[["data"]]
    metadata[[role.name]] <- compiled[["metadata"]]
  }

  list(response = response, roles = roles, metadata = metadata)
}

.crs_iv_eval_dots <- function(dot.expressions, eval.env,
                              exclude = c("weights", "starting.values")) {
  keep <- setdiff(names(dot.expressions), exclude)
  out <- lapply(dot.expressions[keep], eval, envir = eval.env)
  names(out) <- keep
  out
}

.crs_iv_compact_formula <- function(formula) {
  ans <- formula
  environment(ans) <- emptyenv()
  ans
}

.crs_iv_compact_public_call <- function(call) {
  ans <- call
  observation.args <- c(
    "y", "z", "w", "x", "data", "subset", "weights",
    "starting.values", "zeval", "weval", "xeval", "newdata"
  )
  for (i in seq_along(ans)) {
    call.names <- names(ans)
    arg.name <- if (is.null(call.names)) "" else call.names[[i]]
    if (inherits(ans[[i]], "formula")) {
      ans[[i]] <- .crs_iv_compact_formula(ans[[i]])
    } else if (i > 1L && arg.name %in% observation.args &&
               !is.symbol(ans[[i]]) && !is.call(ans[[i]])) {
      ans[[i]] <- as.name(sprintf("<supplied-%s>", arg.name))
    }
  }
  ans
}

.crs_iv_selected_state <- function(object) {
  selected <- object[["num.iterations"]]
  evaluated <- length(object[["norm.stop"]])
  stopping <- if (!is.null(selected) && length(selected) == 1L &&
                  !is.na(selected) && selected >= 1L &&
                  evaluated >= selected) {
    object[["norm.stop"]][[selected]]
  } else {
    NULL
  }
  list(
    iteration = selected,
    evaluated_state_count = evaluated,
    stopping_value = stopping,
    convergence = object[["convergence"]]
  )
}

.crs_iv_native_role_metadata <- function(data) {
  if (is.null(data)) {
    return(NULL)
  }
  data <- data.frame(data, check.names = FALSE)
  role.names <- names(data)
  list(
    labels = role.names,
    internal_names = role.names,
    source_variables = lapply(role.names, function(x) x),
    classes = lapply(data, class),
    xlevels = lapply(data, function(x) if (is.factor(x)) levels(x) else NULL)
  )
}

.crs_iv_attach_native_metadata <- function(object, estimator, public.call,
                                           y, z, w, x, zeval,
                                           is.training.grid) {
  y <- as.vector(y)
  z <- data.frame(z, check.names = FALSE)
  w <- data.frame(w, check.names = FALSE)
  if (!is.null(x)) {
    x <- data.frame(x, check.names = FALSE)
  }
  zeval <- data.frame(zeval, check.names = FALSE)

  object[["iv"]] <- list(
    schema_version = 1L,
    estimator = estimator,
    public_call = public.call,
    public_formula = NULL,
    roles = list(
      z = .crs_iv_native_role_metadata(z),
      w = .crs_iv_native_role_metadata(w),
      x = .crs_iv_native_role_metadata(x)
    ),
    rows = list(
      n_input = NROW(y),
      retained = seq_len(NROW(y)),
      omitted = integer(),
      na_action = NULL
    ),
    training = list(response_original = y),
    evaluation = list(z = zeval, is_training_grid = isTRUE(is.training.grid)),
    selected_state = .crs_iv_selected_state(object)
  )
  object
}

.crs_iv_attach_formula_metadata <- function(object, public.call,
                                             public.formula, compiled,
                                             training) {
  iv <- object[["iv"]]
  iv[["public_call"]] <- public.call
  iv[["public_formula"]] <- .crs_iv_compact_formula(public.formula)
  roles <- compiled[["metadata"]]
  if (!"x" %in% names(roles)) {
    roles <- c(roles, list(x = NULL))
  }
  iv[["roles"]] <- roles
  iv[["rows"]] <- training[["rows"]]
  iv[["training"]][["response_original"]] <- compiled[["response"]]
  iv[["evaluation"]][["z"]] <- compiled[["roles"]][["z"]]
  iv[["evaluation"]][["is_training_grid"]] <- TRUE
  object[["iv"]] <- iv
  object
}

.crs_iv_formula_fit <- function(formula, data, subset, na.action,
                                method.call, public.call, where, estimator,
                                eval.env) {
  parsed <- .crs_iv_parse_formula(formula, where = where)
  dot.expressions <- .crs_iv_dot_expressions(method.call, where = where)
  .crs_iv_reject_formula_evaluation(dot.expressions, where = where)
  training <- .crs_iv_training_frame(
    parsed = parsed,
    method.call = method.call,
    data = data,
    where = where,
    eval.env = eval.env,
    dot.expressions = dot.expressions
  )
  compiled <- .crs_iv_compile_training(parsed, training, where = where)
  dots <- .crs_iv_eval_dots(dot.expressions, eval.env = eval.env)

  args <- list(
    y = compiled[["response"]],
    z = compiled[["roles"]][["z"]],
    w = compiled[["roles"]][["w"]]
  )
  if (!is.null(compiled[["roles"]][["x"]])) {
    args[["x"]] <- compiled[["roles"]][["x"]]
  }
  if (!is.null(training[["starting.values"]])) {
    args[["starting.values"]] <- as.vector(training[["starting.values"]])
  }
  if (!is.null(training[["weights"]])) {
    dots[["weights"]] <- as.vector(training[["weights"]])
  }

  fit.fun <- if (identical(estimator, "crsiv")) {
    crsiv.default
  } else {
    crsivderiv.default
  }
  fit <- do.call(fit.fun, c(args, dots))
  .crs_iv_attach_formula_metadata(
    object = fit,
    public.call = public.call,
    public.formula = formula,
    compiled = compiled,
    training = training
  )
}

crsiv.formula <- function(y, data = NULL, subset, na.action, ...) {
  mc <- match.call(expand.dots = FALSE)
  public.call <- .crs_iv_compact_public_call(match.call())
  public.call[[1L]] <- quote(crsiv)
  .crs_iv_formula_fit(
    formula = y,
    data = data,
    subset = if (missing(subset)) NULL else subset,
    na.action = if (missing(na.action)) NULL else na.action,
    method.call = mc,
    public.call = public.call,
    where = "crsiv()",
    estimator = "crsiv",
    eval.env = parent.frame()
  )
}

crsivderiv.formula <- function(y, data = NULL, subset, na.action, ...) {
  mc <- match.call(expand.dots = FALSE)
  public.call <- .crs_iv_compact_public_call(match.call())
  public.call[[1L]] <- quote(crsivderiv)
  .crs_iv_formula_fit(
    formula = y,
    data = data,
    subset = if (missing(subset)) NULL else subset,
    na.action = if (missing(na.action)) NULL else na.action,
    method.call = mc,
    public.call = public.call,
    where = "crsivderiv()",
    estimator = "crsivderiv",
    eval.env = parent.frame()
  )
}

.crs_iv_metadata <- function(object) {
  object[["iv"]]
}

.crs_iv_public_call <- function(object) {
  iv <- .crs_iv_metadata(object)
  if (!is.null(iv) && !is.null(iv[["public_call"]])) {
    return(iv[["public_call"]])
  }
  object[["call"]]
}

.crs_iv_is_formula_object <- function(object) {
  iv <- .crs_iv_metadata(object)
  !is.null(iv) && !is.null(iv[["public_formula"]])
}

.crs_iv_is_training_grid <- function(object) {
  iv <- .crs_iv_metadata(object)
  if (is.null(iv)) {
    return(NROW(object[["phi"]]) == NROW(object[["y"]]))
  }
  isTRUE(iv[["evaluation"]][["is_training_grid"]])
}

.crs_iv_evaluation_z <- function(object) {
  iv <- .crs_iv_metadata(object)
  if (!is.null(iv) && !is.null(iv[["evaluation"]][["z"]])) {
    return(data.frame(iv[["evaluation"]][["z"]], check.names = FALSE))
  }
  data.frame(object[["xz"]][, 1L, drop = FALSE], check.names = FALSE)
}

.crs_iv_plot_xname <- function(object) {
  iv <- .crs_iv_metadata(object)
  if (!is.null(iv) && !is.null(iv[["roles"]][["z"]][["labels"]])) {
    return(iv[["roles"]][["z"]][["labels"]][[1L]])
  }
  object[["xnames"]][[1L]]
}

.crs_iv_plot_training_data <- function(object) {
  iv <- .crs_iv_metadata(object)
  if (is.null(iv)) {
    return(list(z = object[["xz"]][, 1L], y = object[["y"]]))
  }
  if (!.crs_iv_is_training_grid(object)) {
    stop("training-data overlays are available only when the structural function was evaluated on the training rows",
         call. = FALSE)
  }
  z <- .crs_iv_evaluation_z(object)[[1L]]
  y <- iv[["training"]][["response_original"]]
  if (length(z) != length(y)) {
    stop("training-data overlay rows are not aligned with the selected IV state",
         call. = FALSE)
  }
  list(z = z, y = y)
}

.crs_iv_restore_omitted <- function(object, value, residual = FALSE) {
  iv <- .crs_iv_metadata(object)
  omit <- if (is.null(iv)) NULL else iv[["rows"]][["na_action"]]
  if (is.null(omit)) {
    return(value)
  }
  if (isTRUE(residual)) {
    naresid(omit, value)
  } else {
    napredict(omit, value)
  }
}

fitted.crsiv <- function(object, ...) {
  .crs_iv_restore_omitted(object, object[["phi"]])
}

fitted.crsivderiv <- function(object, ...) {
  .crs_iv_restore_omitted(object, object[["phi"]])
}

.crs_iv_residuals <- function(object) {
  iv <- .crs_iv_metadata(object)
  if (is.null(iv)) {
    return(object[["residuals"]])
  }
  response <- iv[["training"]][["response_original"]]
  if (!.crs_iv_is_training_grid(object) ||
      NROW(response) != length(object[["phi"]])) {
    stop("residuals are available only when the structural function was evaluated on the training rows",
         call. = FALSE)
  }
  ans <- as.vector(response) - as.vector(object[["phi"]])
  .crs_iv_restore_omitted(object, ans, residual = TRUE)
}

residuals.crsiv <- function(object, ...) {
  .crs_iv_residuals(object)
}

residuals.crsivderiv <- function(object, ...) {
  .crs_iv_residuals(object)
}

.crs_iv_validate_deriv <- function(deriv) {
  if (!is.numeric(deriv) || length(deriv) != 1L || is.na(deriv) ||
      !is.finite(deriv) || deriv < 0 ||
      abs(deriv - round(deriv)) > sqrt(.Machine$double.eps)) {
    stop("deriv must be a non-negative integer scalar", call. = FALSE)
  }
  as.integer(round(deriv))
}

.crs_iv_reject_formula_prediction <- function(object, newdata) {
  if (.crs_iv_is_formula_object(object) && !is.null(newdata)) {
    stop("predict() with newdata is not yet supported for formula-created crsiv objects; use the native vector interface for post-fit CRS projection",
         call. = FALSE)
  }
  invisible(NULL)
}

predict.crsiv <- function(object, newdata = NULL, deriv = 0, ...) {
  deriv.supplied <- !missing(deriv)
  if (is.null(newdata)) {
    if (!deriv.supplied || .crs_iv_validate_deriv(deriv) == 0L) {
      return(fitted(object))
    }
    return(predict.crs(object, newdata = NULL, deriv = deriv, ...))
  }
  .crs_iv_reject_formula_prediction(object, newdata)
  if (deriv.supplied) {
    predict.crs(object, newdata = newdata, deriv = deriv, ...)
  } else {
    predict.crs(object, newdata = newdata, ...)
  }
}

predict.crsivderiv <- function(object, newdata = NULL, deriv = 0, ...) {
  deriv.supplied <- !missing(deriv)
  if (is.null(newdata)) {
    if (!deriv.supplied) {
      return(fitted(object))
    }
    deriv <- .crs_iv_validate_deriv(deriv)
    if (deriv == 0L) {
      return(fitted(object))
    }
    if (deriv == 1L) {
      return(.crs_iv_restore_omitted(object, object[["phi.prime"]]))
    }
    return(predict.crs(object, newdata = NULL, deriv = deriv, ...))
  }
  .crs_iv_reject_formula_prediction(object, newdata)
  if (deriv.supplied) {
    predict.crs(object, newdata = newdata, deriv = deriv, ...)
  } else {
    predict.crs(object, newdata = newdata, ...)
  }
}

.crs_iv_display_names <- function(object, internal.names) {
  iv <- .crs_iv_metadata(object)
  if (is.null(iv) || is.null(internal.names)) {
    return(internal.names)
  }
  roles <- iv[["roles"]]
  internal <- unlist(lapply(roles, `[[`, "internal_names"), use.names = FALSE)
  labels <- unlist(lapply(roles, `[[`, "labels"), use.names = FALSE)
  index <- match(internal.names, internal)
  out <- internal.names
  out[!is.na(index)] <- labels[index[!is.na(index)]]
  out
}

.crs_iv_summary_payload <- function(object, derivative = FALSE) {
  iv <- .crs_iv_metadata(object)
  roles <- if (is.null(iv)) NULL else iv[["roles"]]
  role.counts <- if (is.null(roles)) {
    list(z = NA_integer_, w = NA_integer_, x = NA_integer_)
  } else {
    setNames(lapply(c("z", "w", "x"), function(role.name) {
      role <- roles[[role.name]]
      if (is.null(role)) 0L else length(role[["labels"]])
    }), c("z", "w", "x"))
  }
  selected.state <- if (is.null(iv)) {
    .crs_iv_selected_state(object)
  } else {
    iv[["selected_state"]]
  }

  list(
    call = .crs_iv_public_call(object),
    formula = if (is.null(iv)) NULL else iv[["public_formula"]],
    title = if (derivative) {
      if (isTRUE(object[["kernel"]])) {
        if (is.null(object[["tau"]]))
          "Nonparametric Instrumental Spline Derivative Estimation (Kernel Weighting)"
        else
          "Nonparametric Instrumental Spline Quantile Derivative Estimation (Kernel Weighting)"
      } else if (is.null(object[["tau"]])) {
        "Nonparametric Instrumental Spline Derivative Estimation"
      } else {
        "Nonparametric Instrumental Spline Quantile Derivative Estimation"
      }
    } else if (isTRUE(object[["kernel"]])) {
      if (is.null(object[["tau"]]))
        "Nonparametric Instrumental Spline Regression (Kernel Weighting)"
      else
        "Nonparametric Instrumental Spline Quantile Regression (Kernel Weighting)"
    } else if (is.null(object[["tau"]])) {
      "Nonparametric Instrumental Spline Regression"
    } else {
      "Nonparametric Instrumental Spline Quantile Regression"
    },
    derivative = derivative,
    role_counts = role.counts,
    role_labels = if (is.null(roles)) NULL else
      setNames(lapply(c("z", "w", "x"), function(role.name) {
        roles[[role.name]][["labels"]]
      }), c("z", "w", "x")),
    tau = object[["tau"]],
    degree = object[["degree"]],
    segments = object[["segments"]],
    lambda = object[["lambda"]],
    include = object[["include"]],
    xnames = .crs_iv_display_names(object, object[["xnames"]]),
    znames = .crs_iv_display_names(object, object[["znames"]]),
    degree_w = object[["degree.w"]],
    segments_w = object[["segments.w"]],
    lambda_w = object[["lambda.w"]],
    include_w = object[["include.w"]],
    xnames_w = .crs_iv_display_names(object, object[["xnames.w"]]),
    znames_w = .crs_iv_display_names(object, object[["znames.w"]]),
    num_x = object[["num.x"]],
    num_z = object[["num.z"]],
    num_x_w = object[["num.x.w"]],
    num_z_w = object[["num.z.w"]],
    complexity = object[["complexity"]],
    knots = object[["knots"]],
    basis = object[["basis"]],
    ntrain = if (is.null(iv)) object[["nobs"]] else length(iv[["training"]][["response_original"]]),
    n_input = if (is.null(iv)) object[["nobs"]] else iv[["rows"]][["n_input"]],
    n_omitted = if (is.null(iv)) 0L else length(iv[["rows"]][["omitted"]]),
    alpha = object[["alpha"]],
    selected = selected.state[["iteration"]],
    evaluated = selected.state[["evaluated_state_count"]],
    stopping = selected.state[["stopping_value"]],
    convergence = selected.state[["convergence"]],
    nmulti = object[["nmulti"]],
    nomad_summary = object[["nomad.summary"]],
    elapsed = .crs_elapsed_seconds(object[["ptm"]]),
    capabilities = list(
      fitted = !is.null(object[["phi"]]),
      residuals = .crs_iv_is_training_grid(object),
      selected_derivative = !is.null(object[["phi.prime"]]),
      formula_newdata = FALSE
    )
  )
}

.crs_iv_print_summary <- function(x) {
  cat("Call:\n")
  print(x[["call"]])
  cat("\n", x[["title"]], "\n", sep = "")

  counts <- x[["role_counts"]]
  if (!is.null(counts) && !is.na(counts[["z"]])) {
    cat("\nEndogenous terms: ", format(counts[["z"]]), sep = "")
    cat("\nInstrument terms: ", format(counts[["w"]]), sep = "")
    if (counts[["x"]] > 0L) {
      cat("\nExogenous terms: ", format(counts[["x"]]), sep = "")
    }
  }
  if (!is.null(x[["tau"]])) {
    cat("\nQuantile estimated: tau = ", format(x[["tau"]]), sep = "")
  }

  if (!is.null(x[["degree"]])) {
    for (j in seq_along(x[["degree"]])) {
      label <- if (length(x[["xnames"]]) >= j) x[["xnames"]][[j]] else j
      cat("\nSpline degree/number of segments for ", label, ": ",
          format(x[["degree"]][[j]]), "/",
          format(x[["segments"]][[j]]), sep = "")
    }
  }
  if (!is.null(x[["include"]])) {
    for (j in seq_along(x[["include"]])) {
      label <- if (length(x[["znames"]]) >= j) x[["znames"]][[j]] else j
      cat("\nInclusion indicator for ", label, ": ",
          format(x[["include"]][[j]]), sep = "")
    }
  }
  if (!is.null(x[["lambda"]])) {
    for (j in seq_along(x[["lambda"]])) {
      label <- if (length(x[["znames"]]) >= j) x[["znames"]][[j]] else j
      cat("\nBandwidth for ", label, ": ",
          format(x[["lambda"]][[j]]), sep = "")
    }
  }
  if (!is.null(x[["degree_w"]])) {
    for (j in seq_along(x[["degree_w"]])) {
      label <- if (length(x[["xnames_w"]]) >= j) x[["xnames_w"]][[j]] else j
      cat("\nInstrument spline degree/number of segments for ", label, ": ",
          format(x[["degree_w"]][[j]]), "/",
          format(x[["segments_w"]][[j]]), sep = "")
    }
  }
  if (!is.null(x[["include_w"]])) {
    for (j in seq_along(x[["include_w"]])) {
      label <- if (length(x[["znames_w"]]) >= j) x[["znames_w"]][[j]] else j
      cat("\nInstrument inclusion indicator for ", label, ": ",
          format(x[["include_w"]][[j]]), sep = "")
    }
  }
  if (!is.null(x[["lambda_w"]])) {
    for (j in seq_along(x[["lambda_w"]])) {
      label <- if (length(x[["znames_w"]]) >= j) x[["znames_w"]][[j]] else j
      cat("\nInstrument bandwidth for ", label, ": ",
          format(x[["lambda_w"]][[j]]), sep = "")
    }
  }

  cat("\nModel complexity proxy: ", format(x[["complexity"]]), sep = "")
  cat("\nKnot type: ", format(x[["knots"]]), sep = "")
  if (!is.null(x[["basis"]])) {
    cat("\nBasis type: ", format(x[["basis"]]), sep = "")
  }
  cat("\nTraining observations: ", format(x[["ntrain"]]), sep = "")
  if (x[["n_omitted"]] > 0L) {
    cat("\nObservations omitted by the formula frame: ",
        format(x[["n_omitted"]]), sep = "")
  }

  if (!is.null(x[["alpha"]])) {
    cat("\n\nRegularization method: Tikhonov")
    cat("\nTikhonov parameter (alpha): ",
        format(x[["alpha"]], digits = 8), sep = "")
  } else {
    cat("\n\nRegularization method: Landweber-Fridman")
    cat("\nSelected iteration: ", format(x[["selected"]]), sep = "")
    cat("\nStates evaluated: ", format(x[["evaluated"]]), sep = "")
    cat("\nStopping rule value: ",
        format(x[["stopping"]], digits = 8), sep = "")
    if (!is.null(x[["convergence"]])) {
      cat("\nConvergence status: ", x[["convergence"]], sep = "")
    }
  }
  cat("\nNumber of multistarts: ", format(x[["nmulti"]]), sep = "")
  .crs_nomad_summary_print(list(nomad.summary = x[["nomad_summary"]]))
  if (is.finite(x[["elapsed"]])) {
    cat("\nEstimation time: ",
        formatC(x[["elapsed"]], digits = 1, format = "f"),
        " seconds", sep = "")
  }
  cat("\n\n")
  invisible(x)
}

print.summary.crsiv <- function(x, ...) {
  .crs_iv_print_summary(x)
}

print.summary.crsivderiv <- function(x, ...) {
  .crs_iv_print_summary(x)
}
