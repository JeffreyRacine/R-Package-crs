.crs_nomad_geometry_user_has <- function(opts, key) {
  opt.names <- names(opts)
  if (is.null(opt.names))
    return(FALSE)
  key %in% opt.names
}

.crs_nomad_geometry_values <- function(roles,
                                       integer.value,
                                       real.value) {
  roles <- as.character(roles)
  if (!length(roles))
    return(character())
  if (anyNA(roles))
    stop("NOMAD geometry roles must be complete", call. = FALSE)
  ifelse(roles %in% c("integer", "degree_integer", "segments_integer",
                      "include_integer", "lambda_integer_scaled"),
         integer.value,
         real.value)
}

.crs_nomad_apply_source_geometry <- function(opts,
                                             roles,
                                             initial.mesh.size.integer = "1",
                                             initial.mesh.size.real = "r0.1",
                                             min.mesh.size.integer = "1",
                                             min.mesh.size.real = paste("r", sqrt(.Machine$double.eps), sep = ""),
                                             min.frame.size.integer = NULL,
                                             min.frame.size.real = NULL,
                                             initial.frame.size.integer = NULL,
                                             initial.frame.size.real = NULL) {
  roles <- as.character(roles)
  if (!length(roles))
    return(opts)

  generated <- list(
    INITIAL_MESH_SIZE = .crs_nomad_geometry_values(
      roles, initial.mesh.size.integer, initial.mesh.size.real
    ),
    MIN_MESH_SIZE = .crs_nomad_geometry_values(
      roles, min.mesh.size.integer, min.mesh.size.real
    )
  )

  if (!is.null(min.frame.size.integer) || !is.null(min.frame.size.real)) {
    generated$MIN_FRAME_SIZE <- .crs_nomad_geometry_values(
      roles,
      if (is.null(min.frame.size.integer)) min.mesh.size.integer else min.frame.size.integer,
      if (is.null(min.frame.size.real)) min.mesh.size.real else min.frame.size.real
    )
  }

  if (!is.null(initial.frame.size.integer) || !is.null(initial.frame.size.real)) {
    generated$INITIAL_FRAME_SIZE <- .crs_nomad_geometry_values(
      roles,
      if (is.null(initial.frame.size.integer)) initial.mesh.size.integer else initial.frame.size.integer,
      if (is.null(initial.frame.size.real)) initial.mesh.size.real else initial.frame.size.real
    )
  }

  for (key in names(generated)) {
    if (!.crs_nomad_geometry_user_has(opts, key))
      opts[[key]] <- generated[[key]]
  }

  opts
}

