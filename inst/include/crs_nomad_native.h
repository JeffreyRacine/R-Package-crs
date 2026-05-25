#ifndef CRS_NOMAD_NATIVE_H
#define CRS_NOMAD_NATIVE_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

#define CRS_NOMAD_NATIVE_API_VERSION 1

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
  int max_eval;
  unsigned int random_seed;
  int quiet;
  int option_count;
  const crs_nomad_native_option_v1 *options;
} crs_nomad_native_problem_v1;

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

/*
 * Experimental package-author API.
 *
 * Ownership:
 * - problem and its pointer fields are borrowed and read-only.
 * - result->solution is caller-owned storage of length result->solution_len.
 * - crs never returns memory that callers must free.
 * - callback pointers are borrowed for the duration of this call only.
 * - problem->options is a borrowed name/value array of NOMAD option values
 *   expressed as strings. Options are applied after core bounds, budget, seed,
 *   and quiet defaults, so explicit caller options take precedence.
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

#ifdef __cplusplus
}
#endif

#endif
