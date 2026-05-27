#ifndef CRS_NOMAD_NATIVE_H
#define CRS_NOMAD_NATIVE_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

#define CRS_NOMAD_API_VERSION 1

typedef enum {
  CRS_NOMAD_OK = 0,
  CRS_NOMAD_INVALID_INPUT = 1,
  CRS_NOMAD_CALLBACK_FAILURE = 2,
  CRS_NOMAD_INTERNAL_ERROR = 3
} crs_nomad_status;

typedef enum {
  CRS_NOMAD_CALLBACK_C = 0,
  CRS_NOMAD_CALLBACK_R = 1
} crs_nomad_callback_mode;

typedef enum {
  CRS_NOMAD_INPUT_REAL = 0,
  CRS_NOMAD_INPUT_INTEGER = 1,
  CRS_NOMAD_INPUT_CATEGORICAL = 2,
  CRS_NOMAD_INPUT_BINARY = 3
} crs_nomad_input_type;

typedef enum {
  CRS_NOMAD_OUTPUT_OBJ = 0,
  CRS_NOMAD_OUTPUT_PB = 1,
  CRS_NOMAD_OUTPUT_EB = 2
} crs_nomad_output_type;

typedef int (*crs_nomad_eval_fn)(
  int n,
  const double *x,
  int m,
  double *bb_outputs,
  void *user_data
);

typedef struct {
  const char *name;
  const char *value;
} crs_nomad_option;

typedef struct {
  void *eval_f;
  void *environment;
} crs_nomad_r_callback;

typedef struct {
  int api_version;
  size_t struct_size;
  int callback_mode;
  int n;
  int m;
  const double *x0;
  const int *bb_input_type;
  const int *bb_output_type;
  const double *lower;
  const double *upper;
  /*
   * 0 leaves MAX_BB_EVAL/MAX_EVAL unset; positive applies a convenience
   * budget through MAX_BB_EVAL and MAX_EVAL.
   */
  int max_eval;
  unsigned int random_seed;
  int quiet;
  int option_count;
  const crs_nomad_option *options;
  /*
   * Optional row-major start matrix with start_count rows and n columns.
   * If starts is NULL and start_count > 1, crs generates snomadr-compatible
   * random starts from x0, bounds, input types, and random_seed. In that
   * generated-start mode, random_seed applies to start generation; pass an
   * explicit SEED option if NOMAD's internal search seed must also be set.
   * If start_count <= 1 and starts is NULL, x0 is used exactly as the
   * single-start point.
   */
  int start_count;
  const double *starts;
  /*
   * When nonzero, read and apply NOMAD parameters from a path named
   * nomad.opt in the current working directory before applying explicit
   * options. The default native package-consumer contract is 0.
   */
  int read_nomad_opt_file;
  const void *reserved_ptr1;
  const void *reserved_ptr2;
  int reserved_int1;
  int reserved_int2;
} crs_nomad_problem;

typedef struct {
  int api_version;
  size_t struct_size;
  int status;
  int nomad_run_flag;
  int blackbox_evaluations;
  int callback_evaluations;
  int cache_hits;
  int cache_size;
  int total_evaluations;
  int iterations;
  int solution_count;
  int feasible_solution;
  double objective;
  double *solution;
  int solution_len;
  double *outputs;
  int outputs_len;
  char message[512];
  void *reserved_ptr1;
  void *reserved_ptr2;
  int reserved_int1;
  int reserved_int2;
} crs_nomad_result;

typedef int (*crs_nomad_solve_fn)(
  const crs_nomad_problem *problem,
  crs_nomad_eval_fn eval,
  void *user_data,
  crs_nomad_result *result
);

/*
 * Package-author native NOMAD API.
 *
 * Ownership:
 * - problem and its pointer fields are borrowed and read-only.
 * - result->solution is caller-owned storage of length result->solution_len.
 * - result->outputs is caller-owned storage of length result->outputs_len.
 * - crs never returns memory that callers must free.
 * - callback pointers are borrowed for the duration of this call only.
 * - problem->options is a borrowed name/value array of NOMAD option values
 *   expressed as strings. Options are applied after core bounds, optional
 *   convenience budget, seed, and quiet defaults, so explicit caller options
 *   take precedence. If max_eval is 0, crs does not set MAX_BB_EVAL or
 *   MAX_EVAL unless provided in options. A MAX_BB_EVAL option without
 *   MAX_EVAL follows snomadr() compatibility by also setting MAX_EVAL.
 * - random_seed > 0 sets NOMAD's SEED. For generated multi-start points,
 *   random_seed also controls crs-side start generation; random_seed = 0 uses
 *   time-varying generated starts, matching the historical snomadr() posture.
 * - lower/upper may contain -Inf/Inf to express unbounded coordinates. NaN is
 *   invalid in x0, bounds, and generated/explicit starts.
 * - CRS_NOMAD_INPUT_CATEGORICAL is currently represented through NOMAD's
 *   integer input type in the embedded NOMAD4 C interface; callers remain
 *   responsible for category-to-numeric encoding and decoding.
 * - Black-box outputs of NaN are callback failures for all output types.
 *   Infinite OBJ values are callback failures. Infinite PB/EB constraint
 *   values are passed through to NOMAD as infeasibility signals.
 *
 * Callback contract for CRS_NOMAD_CALLBACK_C:
 * - return 0 on successful evaluation, nonzero on evaluation failure.
 * - fill exactly m black-box outputs on success.
 * - do not retain x or bb_outputs after returning.
 * - do not call back into R, longjmp, or throw across this C ABI.
 *
 * Callback contract for CRS_NOMAD_CALLBACK_R:
 * - eval is ignored and may be NULL.
 * - user_data must point to a crs_nomad_r_callback whose eval_f and
 *   environment fields are R SEXP pointers protected by the caller for the
 *   duration of this call.
 * - eval_f must be an R function accepting a numeric vector and returning a
 *   numeric vector of length at least m.
 * - R errors are trapped and returned as callback failures.
 * - R callbacks run on R's main thread. Explicit
 *   NB_THREADS_PARALLEL_EVAL > 1 is rejected for R callbacks, including when
 *   supplied through nomad.opt. Use the default or set it explicitly to 1.
 */
int crs_nomad_solve(
  const crs_nomad_problem *problem,
  crs_nomad_eval_fn eval,
  void *user_data,
  crs_nomad_result *result
);

#ifdef __cplusplus
}
#endif

#endif
