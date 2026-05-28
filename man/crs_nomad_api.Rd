\name{crs_nomad_api}
\alias{crs_nomad_api}
\title{Native NOMAD C API for package authors}
\description{
  \pkg{crs} exposes a native C interface to the embedded NOMAD solver for
  package authors who need to call NOMAD from compiled code without routing
  each black-box evaluation through the \code{snomadr()} R interface.
}
\details{
  The installed header \file{crs_nomad_native.h} defines the final API surface:
  \code{crs_nomad_problem}, \code{crs_nomad_result},
  \code{crs_nomad_option}, \code{crs_nomad_eval_fn}, and
  \code{crs_nomad_solve}.  The callable is registered with
  \code{R_RegisterCCallable()} under the package name \code{"crs"} and symbol
  name \code{"crs_nomad_solve"}.

  The final API supports two callback modes.  \code{CRS_NOMAD_CALLBACK_C}
  calls a package-supplied C callback and is the preferred mode when the
  objective is available in compiled code.  \code{CRS_NOMAD_CALLBACK_R}
  is an explicit bridge for routes whose objective still lives in R; callers
  pass a protected R function and environment through \code{crs_nomad_r_callback}.

  This is the only native NOMAD C callable exported by \pkg{crs}; legacy
  v1/v2 development callables are intentionally not registered.
}
\section{Callback Contract}{
  C callbacks return \code{0} on successful evaluation, fill exactly
  \code{problem->m} black-box outputs, and return nonzero on evaluation
  failure.  They must not retain pointers supplied by \pkg{crs}.

  R callback bridge users must ensure the R function and environment stored in
  \code{crs_nomad_r_callback} remain protected for the duration of
  \code{crs_nomad_solve}.  R errors are trapped and reported as callback
  failures.  R callbacks are main-thread callbacks; explicit
  \code{NB_THREADS_PARALLEL_EVAL > 1} is rejected for R callbacks, including
  when supplied through \file{nomad.opt}.  Use the default or set it explicitly
  to \code{1}.
}
\section{Numeric And Option Semantics}{
  \code{random_seed > 0} sets NOMAD's \code{SEED}.  When \pkg{crs} generates
  multi-start points, the same seed also controls that start generation;
  \code{random_seed = 0} uses time-varying generated starts, matching the
  historical \code{snomadr()} posture.

  Bounds may use \code{-Inf} and \code{Inf} for unbounded coordinates.
  \code{NaN} is invalid in starting points and bounds.  The embedded NOMAD4 C
  interface represents categorical inputs through its integer input type, so
  callers remain responsible for category encoding and decoding.

  \code{NaN} black-box outputs are callback failures for all output types.
  Infinite objective outputs are callback failures.  Infinite \code{PB} or
  \code{EB} constraint outputs are passed through to NOMAD as infeasibility
  signals.

  The native API currently requires NOMAD's evaluation cache
  (\code{EVAL_USE_CACHE = TRUE}, the NOMAD default), because result recovery is
  based on NOMAD's cache-backed incumbent state.  Supplying
  \code{EVAL_USE_CACHE = FALSE} to \code{crs_nomad_solve} is rejected early
  with a clear error.  The public R-level \code{snomadr()} interface retains
  its separate cache option behavior.
}
\section{Ownership}{
  The problem structure and all pointer fields are borrowed and read-only.
  The caller owns \code{result->solution} and \code{result->outputs}; \pkg{crs}
  does not return memory that callers must free.
}
\seealso{
  \code{\link{snomadr}}
}
\keyword{interface}
