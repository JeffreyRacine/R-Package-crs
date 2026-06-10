.crs_nomad_restart_is_better <- function(candidate,
                                         incumbent,
                                         tol = 0) {
  candidate <- as.numeric(candidate)[1L]
  incumbent <- as.numeric(incumbent)[1L]

  if (!is.finite(candidate)) {
    return(FALSE)
  }
  if (!is.finite(incumbent)) {
    return(TRUE)
  }

  candidate < incumbent - abs(tol)
}

.crs_nomad_restart_cache_key <- function(x) {
  x <- as.numeric(x)
  paste(c(length(x), as.character(writeBin(x, raw(), size = 8L))), collapse = ":")
}

.crs_nomad_restart_cached_eval <- function(eval.f,
                                           params) {
  cache.env <- new.env(parent = emptyenv())
  eval.env <- new.env(parent = environment(eval.f))
  eval.env$.crs_nomad_cached_eval_f <- eval.f
  eval.env$.crs_nomad_cache <- cache.env
  eval.env$.crs_nomad_cache_hits <- 0L
  eval.env$.crs_nomad_cache_misses <- 0L
  eval.env$.crs_nomad_restart_cache_key <- .crs_nomad_restart_cache_key

  cached.eval <- function(x, params) {
    current.env <- parent.env(environment())
    key.fun <- get(".crs_nomad_restart_cache_key", envir = current.env, inherits = FALSE)
    cache.env <- get(".crs_nomad_cache", envir = current.env, inherits = FALSE)
    key <- key.fun(x)
    if (exists(key, envir = cache.env, inherits = FALSE)) {
      assign(
        ".crs_nomad_cache_hits",
        get(".crs_nomad_cache_hits", envir = current.env, inherits = FALSE) + 1L,
        envir = current.env
      )
      return(get(key, envir = cache.env, inherits = FALSE))
    }

    eval.fun <- get(".crs_nomad_cached_eval_f", envir = current.env, inherits = FALSE)
    value <- eval.fun(x, params)
    assign(key, value, envir = cache.env)
    assign(
      ".crs_nomad_cache_misses",
      get(".crs_nomad_cache_misses", envir = current.env, inherits = FALSE) + 1L,
      envir = current.env
    )
    value
  }
  environment(cached.eval) <- eval.env

  list(eval.f = cached.eval, environment = eval.env)
}

.crs_nomad_restart_sweep <- function(eval.f,
                                     n,
                                     starts,
                                     bbin,
                                     bbout,
                                     lb,
                                     ub,
                                     random.seed,
                                     opts,
                                     display.nomad.progress,
                                     params,
                                     use.cache = !isTRUE(display.nomad.progress)) {
  starts <- as.matrix(starts)
  nstart <- NROW(starts)
  n <- as.integer(n)[1L]

  if (!is.finite(n) || n <= 0L) {
    stop("NOMAD restart sweep received invalid problem dimension", call. = FALSE)
  }
  if (NCOL(starts) != n) {
    stop("NOMAD restart sweep received a start matrix with the wrong number of columns",
         call. = FALSE)
  }
  if (nstart < 1L) {
    stop("NOMAD restart sweep requires at least one start", call. = FALSE)
  }

  restart.results <- vector("list", nstart)
  restart.solutions <- matrix(NA_real_, nrow = nstart, ncol = n)
  best.index <- NA_integer_
  best.objective <- Inf
  best.solution <- NULL
  best.result <- NULL

  total.bbe <- 0
  total.cache.hits <- 0
  total.callback.evaluations <- 0
  total.evaluations <- 0
  cached.wrapper <- NULL
  callback.cache.hits.before <- 0L
  callback.cache.misses.before <- 0L
  callback.cache.entries.before <- 0L

  if (isTRUE(use.cache)) {
    cached.wrapper <- .crs_nomad_restart_cached_eval(eval.f, params)
    eval.f <- cached.wrapper$eval.f
  }

  for (i in seq_len(nstart)) {
    start.i <- as.numeric(starts[i, ])
    if (!is.null(cached.wrapper)) {
      callback.cache.hits.before <- get(
        ".crs_nomad_cache_hits",
        envir = cached.wrapper$environment,
        inherits = FALSE
      )
      callback.cache.misses.before <- get(
        ".crs_nomad_cache_misses",
        envir = cached.wrapper$environment,
        inherits = FALSE
      )
      callback.cache.entries.before <- length(ls(
        envir = get(".crs_nomad_cache", envir = cached.wrapper$environment, inherits = FALSE),
        all.names = TRUE
      ))
    }

    elapsed <- system.time({
      snomadr.args <- list(
        eval.f = eval.f,
        n = n,
        x0 = start.i,
        bbin = bbin,
        bbout = bbout,
        lb = lb,
        ub = ub,
        nmulti = 0L,
        random.seed = random.seed,
        opts = opts,
        display.nomad.progress = display.nomad.progress,
        params = params
      )
      if (!is.null(cached.wrapper)) {
        snomadr.args$snomadr.environment <- cached.wrapper$environment
      }
      result.i <- do.call(snomadr, snomadr.args)
    })[["elapsed"]]

    callback.cache.hits.i <- NA_integer_
    callback.cache.misses.i <- NA_integer_
    callback.cache.entries.i <- NA_integer_
    if (!is.null(cached.wrapper)) {
      callback.cache.hits.now <- get(
        ".crs_nomad_cache_hits",
        envir = cached.wrapper$environment,
        inherits = FALSE
      )
      callback.cache.misses.now <- get(
        ".crs_nomad_cache_misses",
        envir = cached.wrapper$environment,
        inherits = FALSE
      )
      callback.cache.entries.now <- length(ls(
        envir = get(".crs_nomad_cache", envir = cached.wrapper$environment, inherits = FALSE),
        all.names = TRUE
      ))
      callback.cache.hits.i <- callback.cache.hits.now - callback.cache.hits.before
      callback.cache.misses.i <- callback.cache.misses.now - callback.cache.misses.before
      callback.cache.entries.i <- callback.cache.entries.now - callback.cache.entries.before
    }

    objective.i <- as.numeric(result.i$objective)[1L]
    solution.i <- as.numeric(result.i$solution)
    if (length(solution.i) == n) {
      restart.solutions[i, ] <- solution.i
    }

    bbe.i <- as.numeric(result.i$bbe)[1L]
    cache.hits.i <- as.numeric(result.i$cache.hits)[1L]
    callback.evaluations.i <- as.numeric(result.i$callback.evaluations)[1L]
    total.evaluations.i <- as.numeric(result.i$total.evaluations)[1L]

    restart.results[[i]] <- list(
      restart = as.integer(i),
      objective = objective.i,
      status = as.integer(result.i$status)[1L],
      message = as.character(result.i$message)[1L],
      bbe = bbe.i,
      cache.hits = cache.hits.i,
      cache.size = as.numeric(result.i$cache.size)[1L],
      callback.evaluations = callback.evaluations.i,
      total.evaluations = total.evaluations.i,
      iterations = as.numeric(result.i$iterations)[1L],
      callback.cache.hits = callback.cache.hits.i,
      callback.cache.misses = callback.cache.misses.i,
      callback.cache.entries.added = callback.cache.entries.i,
      elapsed = as.numeric(elapsed)
    )

    if (is.finite(bbe.i)) {
      total.bbe <- total.bbe + bbe.i
    }
    if (is.finite(cache.hits.i)) {
      total.cache.hits <- total.cache.hits + cache.hits.i
    }
    if (is.finite(callback.evaluations.i)) {
      total.callback.evaluations <- total.callback.evaluations + callback.evaluations.i
    }
    if (is.finite(total.evaluations.i)) {
      total.evaluations <- total.evaluations + total.evaluations.i
    }

    if (.crs_nomad_restart_is_better(objective.i, best.objective)) {
      best.index <- as.integer(i)
      best.objective <- objective.i
      best.solution <- solution.i
      best.result <- result.i
    }
  }

  if (is.null(best.result) || !is.finite(best.objective)) {
    stop("NOMAD restart sweep did not produce a finite objective from any restart",
         call. = FALSE)
  }

  ledger <- do.call(
    rbind.data.frame,
    c(restart.results, list(stringsAsFactors = FALSE, make.row.names = FALSE))
  )

  best.result$objective <- best.objective
  best.result$solution <- best.solution
  best.result$bbe <- as.integer(total.bbe)
  best.result$cache.hits <- as.integer(total.cache.hits)
  best.result$callback.evaluations <- as.integer(total.callback.evaluations)
  best.result$total.evaluations <- as.integer(total.evaluations)
  best.result$restart.results <- ledger
  best.result$restart.fval <- ledger$objective
  best.result$best.restart <- best.index
  best.result$restart.starts <- starts
  best.result$restart.solutions <- restart.solutions
  best.result$restart.contract <- "independent-single-start-best"
  best.result$restart.callback.cache.enabled <- isTRUE(use.cache)
  if (!is.null(cached.wrapper)) {
    best.result$restart.callback.cache.hits <- get(
      ".crs_nomad_cache_hits",
      envir = cached.wrapper$environment,
      inherits = FALSE
    )
    best.result$restart.callback.cache.misses <- get(
      ".crs_nomad_cache_misses",
      envir = cached.wrapper$environment,
      inherits = FALSE
    )
    best.result$restart.callback.cache.size <- length(ls(
      envir = get(".crs_nomad_cache", envir = cached.wrapper$environment, inherits = FALSE),
      all.names = TRUE
    ))
  } else {
    best.result$restart.callback.cache.hits <- NA_integer_
    best.result$restart.callback.cache.misses <- NA_integer_
    best.result$restart.callback.cache.size <- NA_integer_
  }

  best.result
}
