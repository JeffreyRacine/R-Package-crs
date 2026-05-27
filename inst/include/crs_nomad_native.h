#ifndef CRS_NOMAD_NATIVE_H
#define CRS_NOMAD_NATIVE_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

#define CRS_NOMAD_NATIVE_API_VERSION 1
#define CRS_NOMAD_NATIVE_API_VERSION_V2 2

typedef enum {
  CRS_NOMAD_OK = 0,
  CRS_NOMAD_INVALID_INPUT = 1,
  CRS_NOMAD_CALLBACK_FAILURE = 2,
  CRS_NOMAD_INTERNAL_ERROR = 3
} crs_nomad_native_status_v1;

typedef int (*crs_nomad_native_eval_fn)(
  int n,
  const double *x,
  int m,
  double *bb_outputs,
  void *user_data
);

typedef struct {
  const char *name;
  const char *value;
} crs_nomad_native_option_v1;

typedef struct {
  int api_version;
  size_t struct_size;
  int n;
  int m;
  const double *x0;
  const int *bb_input_type;
  const double *lower;
  const double *upper;
  /* 0 leaves MAX_BB_EVAL/MAX_EVAL unset; positive preserves v1 convenience budget. */
  int max_eval;
  unsigned int random_seed;
  int quiet;
  int option_count;
  const crs_nomad_native_option_v1 *options;
} crs_nomad_native_problem_v1;

typedef struct {
  int api_version;
  size_t struct_size;
  int n;
  int m;
  const double *x0;
  const int *bb_input_type;
  const double *lower;
  const double *upper;
  /* 0 leaves MAX_BB_EVAL/MAX_EVAL unset; positive preserves v1 convenience budget. */
  int max_eval;
  unsigned int random_seed;
  int quiet;
  int option_count;
  const crs_nomad_native_option_v1 *options;
  /*
   * Optional row-major start matrix with start_count rows and n columns.
   * If starts is NULL and start_count > 1, crs generates snomadr-compatible
   * random starts from x0, bounds, input types, and random_seed. In that
   * generated-start mode, random_seed applies to start generation; pass an
   * explicit SEED option if NOMAD's internal search seed must also be set.
   * If start_count <= 1 and starts is NULL, x0 is used exactly as in v1.
   */
  int start_count;
  const double *starts;
} crs_nomad_native_problem_v2;

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
  char message[512];
} crs_nomad_native_result_v1;

typedef int (*crs_nomad_native_solve_fn_v1)(
  const crs_nomad_native_problem_v1 *problem,
  crs_nomad_native_eval_fn eval,
  void *user_data,
  crs_nomad_native_result_v1 *result
);

typedef crs_nomad_native_result_v1 crs_nomad_native_result_v2;

typedef int (*crs_nomad_native_solve_fn_v2)(
  const crs_nomad_native_problem_v2 *problem,
  crs_nomad_native_eval_fn eval,
  void *user_data,
  crs_nomad_native_result_v2 *result
);

/*
 * Experimental package-author API.
 *
 * Ownership:
 * - problem and its pointer fields are borrowed and read-only.
 * - result->solution is caller-owned storage of length result->solution_len.
 * - crs never returns memory that callers must free.
 * - callback pointers are borrowed for the duration of this call only.
 * - problem->options is a borrowed name/value array of NOMAD option values
 *   expressed as strings. Options are applied after core bounds, optional
 *   convenience budget, seed, and quiet defaults, so explicit caller options
 *   take precedence. If max_eval is 0, crs does not set MAX_BB_EVAL or
 *   MAX_EVAL unless provided in options. A MAX_BB_EVAL option without
 *   MAX_EVAL follows snomadr() compatibility by also setting MAX_EVAL.
 *
 * Callback contract:
 * - return 0 on successful evaluation, nonzero on evaluation failure.
 * - fill exactly m black-box outputs on success.
 * - do not retain x or bb_outputs after returning.
 * - do not call back into R, longjmp, or throw across this C ABI.
 */
int crs_nomad_native_solve_v1(
  const crs_nomad_native_problem_v1 *problem,
  crs_nomad_native_eval_fn eval,
  void *user_data,
  crs_nomad_native_result_v1 *result
);

/*
 * Native v2 extends v1 by accepting either an explicit row-major start matrix
 * or snomadr-compatible generated starts for inner multistart without
 * changing v1.
 * Result output is unchanged and uses the v1/v2-compatible result layout.
 */
int crs_nomad_native_solve_v2(
  const crs_nomad_native_problem_v2 *problem,
  crs_nomad_native_eval_fn eval,
  void *user_data,
  crs_nomad_native_result_v2 *result
);

#ifdef __cplusplus
}
#endif

#endif
