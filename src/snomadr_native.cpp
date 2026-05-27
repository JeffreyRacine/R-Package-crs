#include "../inst/include/crs_nomad_native.h"

#include "NomadStdCInterface.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <vector>

namespace {

std::atomic_flag native_solver_busy = ATOMIC_FLAG_INIT;

struct NativeEvalContext {
  crs_nomad_native_eval_fn eval;
  void *user_data;
  int callback_evaluations;
  int callback_failures;
};

void set_message(crs_nomad_native_result_v1 *result, const char *message) {
  if (result == nullptr) {
    return;
  }
  std::strncpy(result->message, message, sizeof(result->message) - 1);
  result->message[sizeof(result->message) - 1] = '\0';
}

void initialize_result(crs_nomad_native_result_v1 *result) {
  if (result == nullptr) {
    return;
  }
  result->status = CRS_NOMAD_INVALID_INPUT;
  result->nomad_run_flag = 0;
  result->blackbox_evaluations = 0;
  result->callback_evaluations = 0;
  result->cache_hits = 0;
  result->cache_size = 0;
  result->total_evaluations = 0;
  result->iterations = 0;
  result->solution_count = 0;
  result->feasible_solution = 0;
  result->objective = std::numeric_limits<double>::quiet_NaN();
  set_message(result, "");
}

int fail(crs_nomad_native_result_v1 *result,
         crs_nomad_native_status_v1 status,
         const char *message) {
  if (result != nullptr) {
    result->status = status;
    set_message(result, message);
  }
  return status;
}

crs_nomad_native_problem_v1 problem_v1_view(const crs_nomad_native_problem_v2 *problem) {
  crs_nomad_native_problem_v1 view;
  view.api_version = CRS_NOMAD_NATIVE_API_VERSION;
  view.struct_size = sizeof(crs_nomad_native_problem_v1);
  view.n = problem->n;
  view.m = problem->m;
  view.x0 = problem->x0;
  view.bb_input_type = problem->bb_input_type;
  view.lower = problem->lower;
  view.upper = problem->upper;
  view.max_eval = problem->max_eval;
  view.random_seed = problem->random_seed;
  view.quiet = problem->quiet;
  view.option_count = problem->option_count;
  view.options = problem->options;
  return view;
}

int size_to_int_saturated(std::size_t value) {
  const std::size_t max_int =
    static_cast<std::size_t>(std::numeric_limits<int>::max());
  if (value > max_int) {
    return std::numeric_limits<int>::max();
  }
  return static_cast<int>(value);
}

const char *run_flag_message(int run_flag) {
  switch (run_flag) {
    case 1:
      return "NOMAD converged to a feasible point or reached objective target";
    case 0:
      return "NOMAD returned a feasible point after spending its evaluation budget";
    case -1:
      return "NOMAD mesh converged with no feasible point";
    case -2:
      return "NOMAD evaluation budget was spent with no feasible point";
    case -3:
      return "NOMAD initial point failed to evaluate";
    case -4:
      return "NOMAD time limit reached";
    case -5:
      return "NOMAD stopped by callback/user stop";
    case -6:
      return "NOMAD stopped on feasible point";
    case -7:
      return "NOMAD rejected the supplied parameters";
    case -8:
      return "NOMAD reported an evaluation failure";
    default:
      return "NOMAD returned an unrecognized run flag";
  }
}

void suppress_native_nomad_display(NomadProblem pb) {
  addNomadValParam(pb, "DISPLAY_DEGREE", 0);
  addNomadBoolParam(pb, "DISPLAY_ALL_EVAL", false);
  addNomadBoolParam(pb, "DISPLAY_UNSUCCESSFUL", false);
  addNomadBoolParam(pb, "DISPLAY_INFEASIBLE", false);
}

bool map_input_type(int code, char *out) {
  switch (code) {
    case 0:
      *out = 'R';
      return true;
    case 1:
    case 2:
      *out = 'I';
      return true;
    case 3:
      *out = 'B';
      return true;
    default:
      *out = '\0';
      return false;
  }
}

std::string input_type_string(const crs_nomad_native_problem_v1 *problem) {
  std::ostringstream oss;
  oss << "( ";
  for (int i = 0; i < problem->n; ++i) {
    char code = '\0';
    map_input_type(problem->bb_input_type[i], &code);
    oss << code;
    if (i + 1 < problem->n) {
      oss << ' ';
    }
  }
  oss << " )";
  return oss.str();
}

bool is_bad_number(double x) {
  return std::isnan(x);
}

double coerce_x0_value(double x, int input_type, double lower, double upper) {
  if (input_type == 3) {
    x = (x >= 0.5) ? 1.0 : 0.0;
  } else if (input_type == 1 || input_type == 2) {
    x = std::round(x);
  }
  if (!std::isnan(lower)) {
    x = std::max(x, lower);
  }
  if (!std::isnan(upper)) {
    x = std::min(x, upper);
  }
  return x;
}

double finite_or_default(double value, double fallback) {
  return std::isfinite(value) ? value : fallback;
}

void fill_random_start(std::vector<double> &starts,
                       int point_index,
                       const crs_nomad_native_problem_v2 *problem,
                       std::mt19937_64 &rng) {
  const int n = problem->n;
  for (int i = 0; i < n; ++i) {
    double lo = finite_or_default(problem->lower[i], -1.0);
    double hi = finite_or_default(problem->upper[i], 1.0);
    if (hi < lo) {
      std::swap(lo, hi);
    }
    std::uniform_real_distribution<double> unif(lo, hi);
    const double raw = unif(rng);
    starts[static_cast<std::size_t>(point_index) * static_cast<std::size_t>(n) +
           static_cast<std::size_t>(i)] = coerce_x0_value(
      raw,
      problem->bb_input_type[i],
      problem->lower[i],
      problem->upper[i]
    );
  }
}

int validate_problem(const crs_nomad_native_problem_v1 *problem,
                     crs_nomad_native_eval_fn eval,
                     crs_nomad_native_result_v1 *result) {
  if (problem == nullptr) {
    return fail(result, CRS_NOMAD_INVALID_INPUT, "problem pointer is null");
  }
  if (eval == nullptr) {
    return fail(result, CRS_NOMAD_INVALID_INPUT, "native callback pointer is null");
  }
  if (problem->api_version != CRS_NOMAD_NATIVE_API_VERSION) {
    return fail(result, CRS_NOMAD_INVALID_INPUT, "unsupported problem api_version");
  }
  if (problem->struct_size < sizeof(crs_nomad_native_problem_v1)) {
    return fail(result, CRS_NOMAD_INVALID_INPUT, "problem struct_size is too small");
  }
  if (problem->n <= 0) {
    return fail(result, CRS_NOMAD_INVALID_INPUT, "problem n must be positive");
  }
  if (problem->m != 1) {
    return fail(result, CRS_NOMAD_INVALID_INPUT, "problem m must be 1 in native v1");
  }
  if (problem->x0 == nullptr) {
    return fail(result, CRS_NOMAD_INVALID_INPUT, "problem x0 pointer is null");
  }
  if (problem->bb_input_type == nullptr) {
    return fail(result, CRS_NOMAD_INVALID_INPUT, "problem bb_input_type pointer is null");
  }
  if (problem->lower == nullptr) {
    return fail(result, CRS_NOMAD_INVALID_INPUT, "problem lower pointer is null");
  }
  if (problem->upper == nullptr) {
    return fail(result, CRS_NOMAD_INVALID_INPUT, "problem upper pointer is null");
  }
  if (problem->max_eval < 0) {
    return fail(result, CRS_NOMAD_INVALID_INPUT, "problem max_eval must be nonnegative");
  }
  if (result->solution == nullptr) {
    return fail(result, CRS_NOMAD_INVALID_INPUT, "result solution pointer is null");
  }
  if (result->solution_len < problem->n) {
    return fail(result, CRS_NOMAD_INVALID_INPUT, "result solution_len is smaller than problem n");
  }
  for (int i = 0; i < problem->n; ++i) {
    char input_type = '\0';
    if (!map_input_type(problem->bb_input_type[i], &input_type)) {
      return fail(result, CRS_NOMAD_INVALID_INPUT, "unsupported bb_input_type value");
    }
    if (is_bad_number(problem->x0[i])) {
      return fail(result, CRS_NOMAD_INVALID_INPUT, "problem x0 contains NaN");
    }
    if (is_bad_number(problem->lower[i])) {
      return fail(result, CRS_NOMAD_INVALID_INPUT, "problem lower contains NaN");
    }
    if (is_bad_number(problem->upper[i])) {
      return fail(result, CRS_NOMAD_INVALID_INPUT, "problem upper contains NaN");
    }
    if (problem->lower[i] > problem->upper[i]) {
      return fail(result, CRS_NOMAD_INVALID_INPUT, "problem lower bound exceeds upper bound");
    }
  }
  if (problem->random_seed >
      static_cast<unsigned int>(std::numeric_limits<int>::max())) {
    return fail(result, CRS_NOMAD_INVALID_INPUT, "problem random_seed exceeds NOMAD integer range");
  }
  if (problem->option_count < 0) {
    return fail(result, CRS_NOMAD_INVALID_INPUT, "problem option_count must be nonnegative");
  }
  if (problem->option_count > 0 && problem->options == nullptr) {
    return fail(result, CRS_NOMAD_INVALID_INPUT, "problem options pointer is null");
  }
  for (int i = 0; i < problem->option_count; ++i) {
    const crs_nomad_native_option_v1 *option = problem->options + i;
    if (option->name == nullptr || option->name[0] == '\0') {
      return fail(result, CRS_NOMAD_INVALID_INPUT, "problem option name is null or empty");
    }
    if (option->value == nullptr) {
      return fail(result, CRS_NOMAD_INVALID_INPUT, "problem option value is null");
    }
  }
  return CRS_NOMAD_OK;
}

int validate_problem(const crs_nomad_native_problem_v2 *problem,
                     crs_nomad_native_eval_fn eval,
                     crs_nomad_native_result_v2 *result) {
  if (problem == nullptr) {
    return fail(result, CRS_NOMAD_INVALID_INPUT, "problem pointer is null");
  }
  if (problem->api_version != CRS_NOMAD_NATIVE_API_VERSION_V2) {
    return fail(result, CRS_NOMAD_INVALID_INPUT, "unsupported problem api_version");
  }
  if (problem->struct_size < sizeof(crs_nomad_native_problem_v2)) {
    return fail(result, CRS_NOMAD_INVALID_INPUT, "problem struct_size is too small");
  }
  const crs_nomad_native_problem_v1 view = problem_v1_view(problem);
  int status = validate_problem(&view, eval, result);
  if (status != CRS_NOMAD_OK) {
    return status;
  }
  if (problem->start_count < 0) {
    return fail(result, CRS_NOMAD_INVALID_INPUT, "problem start_count must be nonnegative");
  }
  if (problem->start_count > 0) {
    const std::size_t nstart = static_cast<std::size_t>(problem->start_count);
    const std::size_t n = static_cast<std::size_t>(problem->n);
    if (nstart > std::numeric_limits<std::size_t>::max() / n) {
      return fail(result, CRS_NOMAD_INVALID_INPUT, "problem start matrix is too large");
    }
    if (problem->starts != nullptr) {
      for (std::size_t j = 0; j < nstart; ++j) {
        for (std::size_t i = 0; i < n; ++i) {
          if (is_bad_number(problem->starts[j * n + i])) {
            return fail(result, CRS_NOMAD_INVALID_INPUT, "problem starts contains NaN");
          }
        }
      }
    }
  }
  return CRS_NOMAD_OK;
}

bool native_eval_single(int nb_inputs,
                        double *x,
                        int nb_outputs,
                        double *bb_outputs,
                        bool *count_eval,
                        NomadUserDataPtr user_data) {
  NativeEvalContext *context = static_cast<NativeEvalContext *>(user_data);
  if (count_eval != nullptr) {
    *count_eval = false;
  }
  if (context == nullptr || context->eval == nullptr ||
      x == nullptr || bb_outputs == nullptr ||
      nb_inputs <= 0 || nb_outputs <= 0) {
    return false;
  }

  ++context->callback_evaluations;
  int status = 1;
  try {
    status = context->eval(
      nb_inputs,
      static_cast<const double *>(x),
      nb_outputs,
      bb_outputs,
      context->user_data
    );
  } catch (...) {
    ++context->callback_failures;
    return false;
  }

  if (status != 0) {
    ++context->callback_failures;
    return false;
  }
  for (int i = 0; i < nb_outputs; ++i) {
    if (!std::isfinite(bb_outputs[i])) {
      ++context->callback_failures;
      return false;
    }
  }
  if (count_eval != nullptr) {
    *count_eval = true;
  }
  return true;
}

int apply_native_options(const crs_nomad_native_problem_v1 *problem,
                         NomadProblem pb,
                         crs_nomad_native_result_v1 *result) {
  bool seen_max_eval = false;
  bool seen_max_bb_eval = false;
  int max_bb_eval_value = 0;

  for (int i = 0; i < problem->option_count; ++i) {
    const crs_nomad_native_option_v1 *option = problem->options + i;
    if (std::strcmp(option->name, "MAX_EVAL") == 0) {
      seen_max_eval = true;
    }
    if (std::strcmp(option->name, "MAX_BB_EVAL") == 0) {
      char *endptr = nullptr;
      const long parsed = std::strtol(option->value, &endptr, 10);
      if (endptr != option->value && *endptr == '\0' &&
          parsed > 0 && parsed <= std::numeric_limits<int>::max()) {
        seen_max_bb_eval = true;
        max_bb_eval_value = static_cast<int>(parsed);
      }
    }
    const std::string line = std::string(option->name) + " " + option->value;
    if (!addNomadParam(pb, line.c_str())) {
      const std::string message = std::string("failed to set NOMAD option ") + option->name;
      return fail(result, CRS_NOMAD_INVALID_INPUT, message.c_str());
    }
  }

  if (!seen_max_eval && seen_max_bb_eval && max_bb_eval_value > 0) {
    if (!addNomadValParam(pb, "MAX_EVAL", max_bb_eval_value)) {
      return fail(result, CRS_NOMAD_INVALID_INPUT, "failed to synthesize NOMAD MAX_EVAL from MAX_BB_EVAL");
    }
  }

  return CRS_NOMAD_OK;
}

int apply_problem_parameters(const crs_nomad_native_problem_v1 *problem,
                             NomadProblem pb,
                             crs_nomad_native_result_v1 *result) {
  if (!addNomadValParam(pb, "DIMENSION", problem->n)) {
    return fail(result, CRS_NOMAD_INTERNAL_ERROR, "failed to set NOMAD DIMENSION");
  }
  const std::string bbin = input_type_string(problem);
  if (!addNomadStringParam(pb, "BB_INPUT_TYPE", bbin.c_str())) {
    return fail(result, CRS_NOMAD_INVALID_INPUT, "failed to set NOMAD BB_INPUT_TYPE");
  }
  if (!addNomadStringParam(pb, "BB_OUTPUT_TYPE", "OBJ")) {
    return fail(result, CRS_NOMAD_INTERNAL_ERROR, "failed to set NOMAD BB_OUTPUT_TYPE");
  }
  if (!addNomadArrayOfDoubleParam(pb, "LOWER_BOUND", problem->lower)) {
    return fail(result, CRS_NOMAD_INVALID_INPUT, "failed to set NOMAD LOWER_BOUND");
  }
  if (!addNomadArrayOfDoubleParam(pb, "UPPER_BOUND", problem->upper)) {
    return fail(result, CRS_NOMAD_INVALID_INPUT, "failed to set NOMAD UPPER_BOUND");
  }
  if (problem->max_eval > 0) {
    if (!addNomadValParam(pb, "MAX_BB_EVAL", problem->max_eval)) {
      return fail(result, CRS_NOMAD_INVALID_INPUT, "failed to set NOMAD MAX_BB_EVAL");
    }
    if (!addNomadValParam(pb, "MAX_EVAL", problem->max_eval)) {
      return fail(result, CRS_NOMAD_INVALID_INPUT, "failed to set NOMAD MAX_EVAL");
    }
  }
  if (problem->random_seed > 0) {
    if (!addNomadValParam(pb, "SEED", static_cast<int>(problem->random_seed))) {
      return fail(result, CRS_NOMAD_INVALID_INPUT, "failed to set NOMAD SEED");
    }
  }
  if (problem->quiet) {
    suppress_native_nomad_display(pb);
  }
  return apply_native_options(problem, pb, result);
}

std::vector<double> build_starting_point(const crs_nomad_native_problem_v1 *problem) {
  std::vector<double> x0(static_cast<std::size_t>(problem->n), 0.0);
  for (int i = 0; i < problem->n; ++i) {
    x0[static_cast<std::size_t>(i)] = coerce_x0_value(
      problem->x0[i],
      problem->bb_input_type[i],
      problem->lower[i],
      problem->upper[i]
    );
  }
  return x0;
}

std::vector<double> build_starting_points(const crs_nomad_native_problem_v2 *problem,
                                          int *start_count) {
  if (problem->start_count <= 0) {
    if (start_count != nullptr) {
      *start_count = 1;
    }
    const crs_nomad_native_problem_v1 view = problem_v1_view(problem);
    return build_starting_point(&view);
  }

  if (start_count != nullptr) {
    *start_count = std::max(1, problem->start_count);
  }
  const std::size_t nstart = static_cast<std::size_t>(std::max(1, problem->start_count));
  const std::size_t n = static_cast<std::size_t>(problem->n);
  std::vector<double> starts(nstart * n, 0.0);

  if (problem->starts == nullptr) {
    std::mt19937_64 rng(
      problem->random_seed == 0 ?
        static_cast<unsigned int>(std::time(nullptr)) :
        problem->random_seed
    );
    for (int j = 0; j < std::max(1, problem->start_count); ++j) {
      fill_random_start(starts, j, problem, rng);
    }
    for (std::size_t i = 0; i < n; ++i) {
      starts[i] = coerce_x0_value(
        problem->x0[i],
        problem->bb_input_type[i],
        problem->lower[i],
        problem->upper[i]
      );
    }
    return starts;
  }

  for (std::size_t j = 0; j < nstart; ++j) {
    for (std::size_t i = 0; i < n; ++i) {
      starts[j * n + i] = coerce_x0_value(
        problem->starts[j * n + i],
        problem->bb_input_type[i],
        problem->lower[i],
        problem->upper[i]
      );
    }
  }
  return starts;
}

int solve_native_problem_with_starts(const crs_nomad_native_problem_v1 *problem,
                                     const double *starts,
                                     int start_count,
                                     crs_nomad_native_eval_fn eval,
                                     void *user_data,
                                     crs_nomad_native_result_v1 *result) {
  NomadProblem pb = createNomadProblem(native_eval_single, nullptr, problem->n, problem->m);
  if (pb == nullptr) {
    return fail(result, CRS_NOMAD_INTERNAL_ERROR, "unable to initialize NOMAD problem");
  }

  int status = apply_problem_parameters(problem, pb, result);
  if (status != CRS_NOMAD_OK) {
    freeNomadProblem(pb);
    return status;
  }

  NomadResult nomad_result = createNomadResult();
  if (nomad_result == nullptr) {
    freeNomadProblem(pb);
    return fail(result, CRS_NOMAD_INTERNAL_ERROR, "unable to allocate NOMAD result");
  }

  NativeEvalContext context;
  context.eval = eval;
  context.user_data = user_data;
  context.callback_evaluations = 0;
  context.callback_failures = 0;

  const int run_flag = solveNomadProblem(nomad_result, pb, start_count, starts, &context);
  result->nomad_run_flag = run_flag;
  result->blackbox_evaluations = context.callback_evaluations;
  result->callback_evaluations = context.callback_evaluations;
  result->cache_hits = cacheHitsNomadResult(nomad_result);
  result->cache_size = cacheSizeNomadResult(nomad_result);
  result->total_evaluations = totalEvaluationsNomadResult(nomad_result);
  if (result->total_evaluations <= 0) {
    result->total_evaluations = size_to_int_saturated(
      static_cast<std::size_t>(std::max(0, context.callback_evaluations)) +
      static_cast<std::size_t>(std::max(0, result->cache_hits))
    );
  }

  int nb_solutions = nbSolutionsNomadResult(nomad_result);
  if (nb_solutions < 0) {
    nb_solutions = 0;
  }
  result->solution_count = nb_solutions;
  result->feasible_solution = feasibleSolutionsFoundNomadResult(nomad_result) ? 1 : 0;
  if (nb_solutions > 0) {
    std::vector<double> best_x(static_cast<std::size_t>(problem->n), 0.0);
    std::vector<double> best_out(static_cast<std::size_t>(problem->m), 0.0);
    if (!loadInputSolutionsNomadResult(best_x.data(), 1, nomad_result) ||
        !loadOutputSolutionsNomadResult(best_out.data(), 1, nomad_result)) {
      freeNomadResult(nomad_result);
      freeNomadProblem(pb);
      return fail(result, CRS_NOMAD_INTERNAL_ERROR, "failed to load NOMAD solution");
    }
    for (int i = 0; i < problem->n; ++i) {
      result->solution[i] = best_x[static_cast<std::size_t>(i)];
    }
    result->objective = best_out[0];
  }

  if (context.callback_failures > 0) {
    freeNomadResult(nomad_result);
    freeNomadProblem(pb);
    return fail(result, CRS_NOMAD_CALLBACK_FAILURE, "native callback reported evaluation failure");
  }

  if (nb_solutions > 0 && (run_flag == 1 || run_flag == 0 || run_flag == -6)) {
    result->status = CRS_NOMAD_OK;
    set_message(result, run_flag_message(run_flag));
    freeNomadResult(nomad_result);
    freeNomadProblem(pb);
    return CRS_NOMAD_OK;
  }

  const crs_nomad_native_status_v1 failure_status =
    (run_flag == -3 || run_flag == -5 || run_flag == -8) ?
      CRS_NOMAD_CALLBACK_FAILURE : CRS_NOMAD_INTERNAL_ERROR;
  status = fail(result, failure_status, run_flag_message(run_flag));
  freeNomadResult(nomad_result);
  freeNomadProblem(pb);
  return status;
}

int solve_native_problem(const crs_nomad_native_problem_v1 *problem,
                         crs_nomad_native_eval_fn eval,
                         void *user_data,
                         crs_nomad_native_result_v1 *result) {
  const std::vector<double> x0 = build_starting_point(problem);
  return solve_native_problem_with_starts(problem, x0.data(), 1, eval, user_data, result);
}

int solve_native_problem(const crs_nomad_native_problem_v2 *problem,
                         crs_nomad_native_eval_fn eval,
                         void *user_data,
                         crs_nomad_native_result_v2 *result) {
  int start_count = 1;
  const std::vector<double> starts = build_starting_points(problem, &start_count);
  crs_nomad_native_problem_v1 view = problem_v1_view(problem);
  if (problem->starts == nullptr && problem->start_count > 1) {
    view.random_seed = 0;
  }
  return solve_native_problem_with_starts(&view, starts.data(), start_count, eval, user_data, result);
}

} // namespace

extern "C" int crs_nomad_native_solve_v1(
  const crs_nomad_native_problem_v1 *problem,
  crs_nomad_native_eval_fn eval,
  void *user_data,
  crs_nomad_native_result_v1 *result
) {
  (void) user_data;

  if (result == nullptr) {
    return CRS_NOMAD_INVALID_INPUT;
  }
  if (result->api_version != CRS_NOMAD_NATIVE_API_VERSION) {
    return fail(result, CRS_NOMAD_INVALID_INPUT, "unsupported result api_version");
  }
  if (result->struct_size < sizeof(crs_nomad_native_result_v1)) {
    return fail(result, CRS_NOMAD_INVALID_INPUT, "result struct_size is too small");
  }

  initialize_result(result);

  if (native_solver_busy.test_and_set(std::memory_order_acquire)) {
    return fail(result, CRS_NOMAD_INVALID_INPUT, "native NOMAD solver is already active");
  }

  int status = CRS_NOMAD_INTERNAL_ERROR;
  try {
    status = validate_problem(problem, eval, result);
    if (status == CRS_NOMAD_OK) {
      status = solve_native_problem(problem, eval, user_data, result);
    }
  } catch (...) {
    status = fail(result, CRS_NOMAD_INTERNAL_ERROR, "unexpected exception in native NOMAD solve");
  }

  native_solver_busy.clear(std::memory_order_release);
  return status;
}

extern "C" int crs_nomad_native_solve_v2(
  const crs_nomad_native_problem_v2 *problem,
  crs_nomad_native_eval_fn eval,
  void *user_data,
  crs_nomad_native_result_v2 *result
) {
  (void) user_data;

  if (result == nullptr) {
    return CRS_NOMAD_INVALID_INPUT;
  }
  if (result->api_version != CRS_NOMAD_NATIVE_API_VERSION &&
      result->api_version != CRS_NOMAD_NATIVE_API_VERSION_V2) {
    return fail(result, CRS_NOMAD_INVALID_INPUT, "unsupported result api_version");
  }
  if (result->struct_size < sizeof(crs_nomad_native_result_v2)) {
    return fail(result, CRS_NOMAD_INVALID_INPUT, "result struct_size is too small");
  }

  initialize_result(result);

  if (native_solver_busy.test_and_set(std::memory_order_acquire)) {
    return fail(result, CRS_NOMAD_INVALID_INPUT, "native NOMAD solver is already active");
  }

  int status = CRS_NOMAD_INTERNAL_ERROR;
  try {
    status = validate_problem(problem, eval, result);
    if (status == CRS_NOMAD_OK) {
      status = solve_native_problem(problem, eval, user_data, result);
    }
  } catch (...) {
    status = fail(result, CRS_NOMAD_INTERNAL_ERROR, "unexpected exception in native NOMAD solve");
  }

  native_solver_busy.clear(std::memory_order_release);
  return status;
}
