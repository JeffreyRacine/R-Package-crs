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

  The callable itself must be invoked from the main R interpreter thread of the
  current process.  It uses R's protection stack, unwind protection, and user
  interrupt machinery, and must not be called from native worker threads, signal
  handlers, or contexts that are not executing through R or a registered native
  call on that thread.
}
\section{Quick Start For Package Authors}{
  A package that calls the native API typically needs:

  \enumerate{
    \item \code{LinkingTo: crs} in \file{DESCRIPTION}, so
      \file{crs_nomad_native.h} is available at compile time.
    \item A runtime dependency path that loads the \pkg{crs} namespace before
      calling \code{R_GetCCallable("crs", "crs_nomad_solve")}.  If the route is
      mandatory, use \code{Imports: crs}; if it is optional, guard the route
      with an explicit \code{requireNamespace("crs")} check before entering the
      compiled solver path.
    \item A callback that fills \code{m} black-box outputs and returns \code{0}
      on success.
    \item A zero-initialized \code{crs_nomad_problem} and
      \code{crs_nomad_result}.  Always set \code{api_version} and
      \code{struct_size}; this lets future \pkg{crs} versions reject
      incompatible structures cleanly.
    \item Caller-owned buffers for \code{result.solution} and
      \code{result.outputs}.  \pkg{crs} does not allocate result memory for the
      caller.
  }

  In ordinary package code, prefer explicit options over \file{nomad.opt}
  files, and leave \code{read_nomad_opt_file = 0} unless reading a working
  directory file is a deliberate part of your user contract.
}
\section{Minimal C Callback Skeleton}{
  The following abridged skeleton shows the direct C-callback route used by
  packages whose objective is available in compiled code.

  \preformatted{
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <string.h>
#include <crs_nomad_native.h>

typedef struct {
  double target;
} my_eval_data;

static int my_eval(int n, const double *x, int m, double *bb, void *user_data)
{
  my_eval_data *data = (my_eval_data *) user_data;
  double value = 0.0;
  int i;

  if (m < 1)
    return 1;

  for (i = 0; i < n; ++i) {
    const double d = x[i] - data->target;
    value += d * d;
  }

  bb[0] = value;      /* one CRS_NOMAD_OUTPUT_OBJ output */
  return 0;
}

SEXP my_nomad_solve(SEXP x0_s)
{
  SEXP x0;
  crs_nomad_solve_fn solve;
  crs_nomad_problem problem;
  crs_nomad_result result;
  crs_nomad_option options[2];
  int input_type[2] = {CRS_NOMAD_INPUT_REAL, CRS_NOMAD_INPUT_REAL};
  int output_type[1] = {CRS_NOMAD_OUTPUT_OBJ};
  double lower[2] = {-1.0, -1.0};
  double upper[2] = { 1.0,  1.0};
  double solution[2] = {R_NaN, R_NaN};
  double outputs[1] = {R_NaN};
  my_eval_data data;
  int status;

  PROTECT(x0 = coerceVector(x0_s, REALSXP));
  if (XLENGTH(x0) != 2) {
    UNPROTECT(1);
    error("x0 must have length 2");
  }

  /*
   * Ensure the crs namespace is loaded before R_GetCCallable().
   * Many packages do this in an R wrapper via requireNamespace("crs").
   */
  solve = (crs_nomad_solve_fn)
    R_GetCCallable("crs", "crs_nomad_solve");
  if (solve == NULL) {
    UNPROTECT(1);
    error("failed to resolve crs_nomad_solve");
  }

  options[0].name = "MAX_BB_EVAL";
  options[0].value = "100";
  options[1].name = "DISPLAY_DEGREE";
  options[1].value = "0";

  memset(&problem, 0, sizeof(problem));
  problem.api_version = CRS_NOMAD_API_VERSION;
  problem.struct_size = sizeof(problem);
  problem.callback_mode = CRS_NOMAD_CALLBACK_C;
  problem.n = 2;
  problem.m = 1;
  problem.x0 = REAL(x0);
  problem.bb_input_type = input_type;
  problem.bb_output_type = output_type;
  problem.lower = lower;
  problem.upper = upper;
  problem.random_seed = 42;
  problem.quiet = 1;
  problem.option_count = 2;
  problem.options = options;
  problem.start_count = 1;
  problem.read_nomad_opt_file = 0;

  memset(&result, 0, sizeof(result));
  result.api_version = CRS_NOMAD_API_VERSION;
  result.struct_size = sizeof(result);
  result.solution = solution;
  result.solution_len = 2;
  result.outputs = outputs;
  result.outputs_len = 1;

  data.target = 0.25;
  status = solve(&problem, my_eval, &data, &result);

  if (status != CRS_NOMAD_OK || result.status != CRS_NOMAD_OK) {
    UNPROTECT(1);
    error("NOMAD failed: %s", result.message);
  }

  UNPROTECT(1);
  return ScalarReal(result.objective);
}
  }

  Real package code should validate all R inputs before constructing the
  borrowed pointer fields, should protect any R objects whose memory is used
  during the call, and should return the solution, objective, evaluation
  counts, and diagnostic message in an object shape appropriate for the
  package.
}
\section{R Callback Bridge Skeleton}{
  The R-callback bridge is useful when the search geometry is managed in
  compiled code but the black-box objective still lives in R.  This is the
  pattern used by the \pkg{np} and \pkg{npRmpi} NOMAD routes while their
  objective functions remain in R.

  \preformatted{
SEXP my_nomad_solve_r_callback(SEXP eval_f, SEXP eval_env, SEXP x0_s)
{
  crs_nomad_solve_fn solve;
  crs_nomad_problem problem;
  crs_nomad_result result;
  crs_nomad_r_callback callback;
  crs_nomad_option options[1];
  int input_type[2] = {CRS_NOMAD_INPUT_REAL, CRS_NOMAD_INPUT_REAL};
  int output_type[1] = {CRS_NOMAD_OUTPUT_OBJ};
  double lower[2] = {-1.0, -1.0};
  double upper[2] = { 1.0,  1.0};
  double solution[2] = {R_NaN, R_NaN};
  double outputs[1] = {R_NaN};
  int status;

  PROTECT(eval_f);
  PROTECT(eval_env);
  PROTECT(x0_s = coerceVector(x0_s, REALSXP));
  if (XLENGTH(x0_s) != 2) {
    UNPROTECT(3);
    error("x0 must have length 2");
  }

  solve = (crs_nomad_solve_fn)
    R_GetCCallable("crs", "crs_nomad_solve");
  if (solve == NULL) {
    UNPROTECT(3);
    error("failed to resolve crs_nomad_solve");
  }

  options[0].name = "NB_THREADS_PARALLEL_EVAL";
  options[0].value = "1";

  memset(&callback, 0, sizeof(callback));
  callback.eval_f = (void *) eval_f;
  callback.environment = (void *) eval_env;

  memset(&problem, 0, sizeof(problem));
  problem.api_version = CRS_NOMAD_API_VERSION;
  problem.struct_size = sizeof(problem);
  problem.callback_mode = CRS_NOMAD_CALLBACK_R;
  problem.n = 2;
  problem.m = 1;
  problem.x0 = REAL(x0_s);
  problem.bb_input_type = input_type;
  problem.bb_output_type = output_type;
  problem.lower = lower;
  problem.upper = upper;
  problem.random_seed = 42;
  problem.quiet = 1;
  problem.option_count = 1;
  problem.options = options;
  problem.start_count = 1;

  memset(&result, 0, sizeof(result));
  result.api_version = CRS_NOMAD_API_VERSION;
  result.struct_size = sizeof(result);
  result.solution = solution;
  result.solution_len = 2;
  result.outputs = outputs;
  result.outputs_len = 1;

  status = solve(&problem, NULL, &callback, &result);

  if (status != CRS_NOMAD_OK || result.status != CRS_NOMAD_OK) {
    UNPROTECT(3);
    error("NOMAD failed: %s", result.message);
  }

  UNPROTECT(3);
  return ScalarReal(result.objective);
}
  }

  For R callbacks, \code{eval_f} must accept a numeric vector and return a
  numeric vector of length at least \code{problem.m}.  Because this bridge calls
  back into R, do not request parallel black-box evaluation:
  \code{NB_THREADS_PARALLEL_EVAL > 1} is rejected.
}
\section{Callback Contract}{
  C callbacks return \code{0} on successful evaluation, fill exactly
  \code{problem->m} black-box outputs, and return nonzero on evaluation
  failure.  They must not retain pointers supplied by \pkg{crs}.

  R callback bridge users must ensure the R function and environment stored in
  \code{crs_nomad_r_callback} remain protected for the duration of
  \code{crs_nomad_solve}.  R errors are trapped and reported as callback
  failures, with the first trapped R diagnostic retained in the fixed-capacity
  result message.  R callbacks are main-thread callbacks; explicit
  \code{NB_THREADS_PARALLEL_EVAL > 1} is rejected for R callbacks, including
  when supplied through \file{nomad.opt}.  Use the default or set it explicitly
  to \code{1}.
}
\section{Interrupt Semantics}{
  Serial consumers normally leave interrupt ownership at
  \code{CRS_NOMAD_INTERRUPT_OWNERSHIP_DEFAULT}.  The solve-scoped handler
  records only signal-safe state; NOMAD consumes it from ordinary solver
  control flow, performs native cleanup, restores the previous process handler,
  and returns \code{CRS_NOMAD_INTERRUPTED}.  If a callback failure occurs
  before the same solve is interrupted, the interrupt status takes precedence.

  Lockstep consumers such as MPI-backed packages may use
  \code{CRS_NOMAD_INTERRUPT_OWNERSHIP_EXTERNAL}.  They must own protected
  polling and defer any R unwind until a rank-common synchronization boundary.
  This option does not add synchronization, messages, or payloads inside
  \pkg{crs}.
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
  Infinite objective outputs are coerced to \code{+/-DBL_MAX}, matching
  \code{snomadr()} compatibility.  Infinite \code{PB} or \code{EB} constraint
  outputs are passed through to NOMAD as infeasibility signals.

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
\section{Common Result Fields}{
  The most useful result fields for package authors are:

  \itemize{
    \item \code{status}: \pkg{crs} API status, where \code{CRS_NOMAD_OK}
      indicates that the native wrapper completed.
    \item \code{nomad_run_flag}: NOMAD's native run flag.
    \item \code{objective}: best objective reported by NOMAD.
    \item \code{solution}: caller-owned solution buffer filled by \pkg{crs}.
    \item \code{outputs}: caller-owned output buffer for the best solution.
    \item \code{callback_evaluations}: number of true package callback
      evaluations.
    \item \code{total_evaluations} and \code{cache_hits}: NOMAD point lookups
      and native cache hits.  These are lookup diagnostics; cache hits are not
      expensive objective recomputations.
    \item \code{message}: diagnostic text when the wrapper reports a failure.
  }
}
\seealso{
  \code{\link{snomadr}}
}
\keyword{interface}
