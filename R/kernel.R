## Faster kernel: same signature & behavior as original
kernel <- function(Z,
                   z,
                   lambda,
                   is.ordered.z = NULL) {

  if (is.null(is.ordered.z) || missing(Z) || missing(z) || missing(lambda))
    stop(" must provide is.ordered.z, Z, z, and lambda")

  # Ensure vector behavior is consistent
  Z <- as.vector(Z)

  if (!is.ordered.z) {
    res <- rep.int(lambda, length(Z))
    res[Z == z] <- 1
    res
  } else {
    res <- lambda^abs(Z - z)
    res[Z == z] <- 1
    res
  }
}

## Faster product kernel: same signature & behavior as original
prod.kernel <- function(Z,
                        z,
                        lambda,
                        is.ordered.z = NULL,
                        ...,
                        na.rm) {

  if (!is.matrix(Z)) Z <- as.matrix(Z)

  if (is.null(is.ordered.z) || missing(Z) || missing(z) || missing(lambda))
    stop(" must provide is.ordered.z, Z, z, and lambda")
  if (length(is.ordered.z) != NCOL(Z))
    stop(" is.ordered.z and Z incompatible")

  num.z <- NCOL(Z)

  ## Follow original API: z and lambda can be column vectors (p x 1)
  ## NROW() is used in the original; keep the same checks
  if (num.z != NROW(z) || num.z != NROW(lambda))
    stop(paste(" incompatible dimensions for Z, z, and lambda (",
               num.z, ",", NROW(z), ",", NROW(lambda), ")", sep = ""))

  ## Inline the kernel logic per column (avoids function-call overhead)
  n <- NROW(Z)
  prodker <- rep.int(1, n)

  ## Access z/lambda as vectors to mirror original indexing semantics
  z_vec <- as.vector(z)
  lambda_vec <- as.vector(lambda)

  for (j in seq_len(num.z)) {
    if (!is.ordered.z[j]) {
      tmp <- rep.int(lambda_vec[j], n)
      eq <- Z[, j] == z_vec[j]
      tmp[eq] <- 1
    } else {
      tmp <- lambda_vec[j]^abs(Z[, j] - z_vec[j])
      eq <- Z[, j] == z_vec[j]
      tmp[eq] <- 1
    }
    prodker <- prodker * tmp
  }

  prodker
}

## Fast-path for matrix Z (avoids repeated is.matrix/as.matrix checks).
## Behavior matches prod.kernel for valid inputs.
prod.kernel.matrix <- function(Z,
                               z,
                               lambda,
                               is.ordered.z = NULL,
                               ...,
                               na.rm) {
  if (is.null(is.ordered.z) || missing(Z) || missing(z) || missing(lambda))
    stop(" must provide is.ordered.z, Z, z, and lambda")
  if (length(is.ordered.z) != NCOL(Z))
    stop(" is.ordered.z and Z incompatible")
  num.z <- NCOL(Z)
  if (num.z != NROW(z) || num.z != NROW(lambda))
    stop(paste(" incompatible dimensions for Z, z, and lambda (",
               num.z, ",", NROW(z), ",", NROW(lambda), ")", sep = ""))

  n <- NROW(Z)
  prodker <- rep.int(1, n)
  z_vec <- as.vector(z)
  lambda_vec <- as.vector(lambda)

  for (j in seq_len(num.z)) {
    if (!is.ordered.z[j]) {
      tmp <- rep.int(lambda_vec[j], n)
      eq <- Z[, j] == z_vec[j]
      tmp[eq] <- 1
    } else {
      tmp <- lambda_vec[j]^abs(Z[, j] - z_vec[j])
      eq <- Z[, j] == z_vec[j]
      tmp[eq] <- 1
    }
    prodker <- prodker * tmp
  }

  prodker
}
