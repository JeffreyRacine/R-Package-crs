#include "../inst/include/crs_nomad_native.h"

#include "NomadStdCInterface.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Error.h>

#ifdef length
#undef length
#endif

#ifdef error
#undef error
#endif

#ifdef match
#undef match
#endif

namespace {

std::atomic_flag native_solver_busy = ATOMIC_FLAG_INIT;

struct NativeEvalContext {
  int callback_mode;
  crs_nomad_eval_fn eval;
  void *user_data;
  int callback_evaluations;
  int callback_failures;
};

void set_message(crs_nomad_result *result, const char *message) {
  if (result == nullptr) {
    return;
  }
  std::strncpy(result->message, message, sizeof(result->message) - 1);
  result->message[sizeof(result->message) - 1] = '\0';
}

void initialize_result(crs_nomad_result *result) {
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
  if (result->solution != nullptr) {
    for (int i = 0; i < result->solution_len; ++i) {
      result->solution[i] = std::numeric_limits<double>::quiet_NaN();
    }
  }
  if (result->outputs != nullptr) {
    for (int i = 0; i < result->outputs_len; ++i) {
      result->outputs[i] = std::numeric_limits<double>::quiet_NaN();
    }
  }
  set_message(result, "");
}

int fail(crs_nomad_result *result, crs_nomad_status status, const char *message) {
  if (result != nullptr) {
    result->status = status;
    set_message(result, message);
  }
  return status;
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
    case CRS_NOMAD_INPUT_REAL:
      *out = 'R';
      return true;
    case CRS_NOMAD_INPUT_INTEGER:
    case CRS_NOMAD_INPUT_CATEGORICAL:
      *out = 'I';
      return true;
    case CRS_NOMAD_INPUT_BINARY:
      *out = 'B';
      return true;
    default:
      *out = '\0';
      return false;
  }
}

bool map_output_type(int code, const char **out) {
  switch (code) {
    case CRS_NOMAD_OUTPUT_OBJ:
      *out = "OBJ";
      return true;
    case CRS_NOMAD_OUTPUT_PB:
      *out = "PB";
      return true;
    case CRS_NOMAD_OUTPUT_EB:
      *out = "EB";
      return true;
    default:
      *out = "BBO_UNDEFINED";
      return false;
  }
}

std::string input_type_string(const crs_nomad_problem *problem) {
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

std::string output_type_string(const crs_nomad_problem *problem) {
  std::ostringstream oss;
  for (int i = 0; i < problem->m; ++i) {
    const char *code = nullptr;
    map_output_type(problem->bb_output_type[i], &code);
    if (i > 0) {
      oss << ' ';
    }
    oss << code;
  }
  return oss.str();
}

bool is_bad_number(double x) {
  return std::isnan(x);
}

double coerce_x0_value(double x, int input_type, double lower, double upper) {
  if (input_type == CRS_NOMAD_INPUT_BINARY) {
    x = (x >= 0.5) ? 1.0 : 0.0;
  } else if (input_type == CRS_NOMAD_INPUT_INTEGER ||
             input_type == CRS_NOMAD_INPUT_CATEGORICAL) {
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
                       const crs_nomad_problem *problem,
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

bool native_eval_c(int nb_inputs,
                   double *x,
                   int nb_outputs,
                   double *bb_outputs,
                   NativeEvalContext *context) {
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
  return true;
}

bool native_eval_r(int nb_inputs,
                   double *x,
                   int nb_outputs,
                   double *bb_outputs,
                   NativeEvalContext *context) {
  crs_nomad_r_callback *callback =
    static_cast<crs_nomad_r_callback *>(context->user_data);
  SEXP eval_f = (callback == nullptr) ?
    R_NilValue : reinterpret_cast<SEXP>(callback->eval_f);
  if (callback == nullptr || eval_f == nullptr || eval_f == R_NilValue ||
      !Rf_isFunction(eval_f)) {
    ++context->callback_failures;
    return false;
  }

  SEXP rho = reinterpret_cast<SEXP>(callback->environment);
  if (rho == R_NilValue || rho == nullptr) {
    rho = R_GlobalEnv;
  }

  int nprotect = 0;
  SEXP rargs = PROTECT(Rf_allocVector(REALSXP, nb_inputs));
  ++nprotect;
  for (int i = 0; i < nb_inputs; ++i) {
    REAL(rargs)[i] = x[i];
  }

  SEXP call = PROTECT(Rf_lang2(eval_f, rargs));
  ++nprotect;

  int error = 0;
  SEXP result = R_tryEval(call, rho, &error);
  if (error || result == R_NilValue) {
    UNPROTECT(nprotect);
    ++context->callback_failures;
    return false;
  }
  result = PROTECT(result);
  ++nprotect;

  if (!Rf_isNumeric(result)) {
    UNPROTECT(nprotect);
    ++context->callback_failures;
    return false;
  }

  SEXP rnum = result;
  if (TYPEOF(rnum) != REALSXP) {
    rnum = PROTECT(Rf_coerceVector(result, REALSXP));
    ++nprotect;
  }

  if (Rf_length(rnum) < nb_outputs) {
    UNPROTECT(nprotect);
    ++context->callback_failures;
    return false;
  }

  for (int i = 0; i < nb_outputs; ++i) {
    bb_outputs[i] = REAL(rnum)[i];
  }

  UNPROTECT(nprotect);
  return true;
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
  if (context == nullptr || x == nullptr || bb_outputs == nullptr ||
      nb_inputs <= 0 || nb_outputs <= 0) {
    return false;
  }

  ++context->callback_evaluations;
  bool ok = false;
  if (context->callback_mode == CRS_NOMAD_CALLBACK_R) {
    ok = native_eval_r(nb_inputs, x, nb_outputs, bb_outputs, context);
  } else {
    if (context->eval == nullptr) {
      ++context->callback_failures;
      return false;
    }
    ok = native_eval_c(nb_inputs, x, nb_outputs, bb_outputs, context);
  }
  if (!ok) {
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

int validate_problem(const crs_nomad_problem *problem,
                     crs_nomad_eval_fn eval,
                     crs_nomad_result *result) {
  if (problem == nullptr) {
    return fail(result, CRS_NOMAD_INVALID_INPUT, "problem pointer is null");
  }
  if (problem->api_version != CRS_NOMAD_API_VERSION) {
    return fail(result, CRS_NOMAD_INVALID_INPUT, "unsupported problem api_version");
  }
  if (problem->struct_size < sizeof(crs_nomad_problem)) {
    return fail(result, CRS_NOMAD_INVALID_INPUT, "problem struct_size is too small");
  }
  if (problem->callback_mode != CRS_NOMAD_CALLBACK_C &&
      problem->callback_mode != CRS_NOMAD_CALLBACK_R) {
    return fail(result, CRS_NOMAD_INVALID_INPUT, "unsupported callback_mode");
  }
  if (problem->callback_mode == CRS_NOMAD_CALLBACK_C && eval == nullptr) {
    return fail(result, CRS_NOMAD_INVALID_INPUT, "native callback pointer is null");
  }
  if (problem->n <= 0) {
    return fail(result, CRS_NOMAD_INVALID_INPUT, "problem n must be positive");
  }
  if (problem->m <= 0) {
    return fail(result, CRS_NOMAD_INVALID_INPUT, "problem m must be positive");
  }
  if (problem->x0 == nullptr) {
    return fail(result, CRS_NOMAD_INVALID_INPUT, "problem x0 pointer is null");
  }
  if (problem->bb_input_type == nullptr) {
    return fail(result, CRS_NOMAD_INVALID_INPUT, "problem bb_input_type pointer is null");
  }
  if (problem->bb_output_type == nullptr) {
    return fail(result, CRS_NOMAD_INVALID_INPUT, "problem bb_output_type pointer is null");
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
  if (result->outputs == nullptr) {
    return fail(result, CRS_NOMAD_INVALID_INPUT, "result outputs pointer is null");
  }
  if (result->outputs_len < problem->m) {
    return fail(result, CRS_NOMAD_INVALID_INPUT, "result outputs_len is smaller than problem m");
  }

  bool has_objective = false;
  for (int i = 0; i < problem->m; ++i) {
    const char *code = nullptr;
    if (!map_output_type(problem->bb_output_type[i], &code)) {
      return fail(result, CRS_NOMAD_INVALID_INPUT, "unsupported bb_output_type value");
    }
    if (problem->bb_output_type[i] == CRS_NOMAD_OUTPUT_OBJ) {
      has_objective = true;
    }
  }
  if (!has_objective) {
    return fail(result, CRS_NOMAD_INVALID_INPUT, "problem bb_output_type must contain an OBJ output");
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
    const crs_nomad_option *option = problem->options + i;
    if (option->name == nullptr || option->name[0] == '\0') {
      return fail(result, CRS_NOMAD_INVALID_INPUT, "problem option name is null or empty");
    }
    if (option->value == nullptr) {
      return fail(result, CRS_NOMAD_INVALID_INPUT, "problem option value is null");
    }
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

int apply_native_options(const crs_nomad_problem *problem,
                         NomadProblem pb,
                         crs_nomad_result *result) {
  bool seen_max_eval = false;
  bool seen_max_bb_eval = false;
  int max_bb_eval_value = 0;

  for (int i = 0; i < problem->option_count; ++i) {
    const crs_nomad_option *option = problem->options + i;
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

int apply_nomad_opt_file(NomadProblem pb, crs_nomad_result *result) {
  std::ifstream fin("nomad.opt");
  if (!fin.good()) {
    return CRS_NOMAD_OK;
  }
  std::string line;
  while (std::getline(fin, line)) {
    std::size_t first = line.find_first_not_of(" \t\r\n");
    if (first == std::string::npos || line[first] == '#') {
      continue;
    }
    if (!addNomadParam(pb, line.c_str())) {
      return fail(result, CRS_NOMAD_INVALID_INPUT, "failed to set NOMAD option from nomad.opt");
    }
  }
  return CRS_NOMAD_OK;
}

int apply_problem_parameters(const crs_nomad_problem *problem,
                             NomadProblem pb,
                             crs_nomad_result *result) {
  if (!addNomadValParam(pb, "DIMENSION", problem->n)) {
    return fail(result, CRS_NOMAD_INTERNAL_ERROR, "failed to set NOMAD DIMENSION");
  }
  const std::string bbin = input_type_string(problem);
  if (!addNomadStringParam(pb, "BB_INPUT_TYPE", bbin.c_str())) {
    return fail(result, CRS_NOMAD_INVALID_INPUT, "failed to set NOMAD BB_INPUT_TYPE");
  }
  const std::string bbout = output_type_string(problem);
  if (!addNomadStringParam(pb, "BB_OUTPUT_TYPE", bbout.c_str())) {
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
  if (problem->read_nomad_opt_file) {
    const int file_status = apply_nomad_opt_file(pb, result);
    if (file_status != CRS_NOMAD_OK) {
      return file_status;
    }
  }
  return apply_native_options(problem, pb, result);
}

std::vector<double> build_starting_point(const crs_nomad_problem *problem) {
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

std::vector<double> build_starting_points(const crs_nomad_problem *problem,
                                          int *start_count) {
  if (problem->start_count <= 0) {
    if (start_count != nullptr) {
      *start_count = 1;
    }
    return build_starting_point(problem);
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

int solve_final_problem_with_starts(const crs_nomad_problem *problem,
                                    const double *starts,
                                    int start_count,
                                    crs_nomad_eval_fn eval,
                                    void *user_data,
                                    crs_nomad_result *result) {
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
  context.callback_mode = problem->callback_mode;
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
    for (int i = 0; i < problem->m; ++i) {
      result->outputs[i] = best_out[static_cast<std::size_t>(i)];
    }
    int objective_index = 0;
    for (int i = 0; i < problem->m; ++i) {
      if (problem->bb_output_type[i] == CRS_NOMAD_OUTPUT_OBJ) {
        objective_index = i;
        break;
      }
    }
    result->objective = best_out[static_cast<std::size_t>(objective_index)];
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

  const crs_nomad_status failure_status =
    (run_flag == -3 || run_flag == -5 || run_flag == -8) ?
      CRS_NOMAD_CALLBACK_FAILURE : CRS_NOMAD_INTERNAL_ERROR;
  status = fail(result, failure_status, run_flag_message(run_flag));
  freeNomadResult(nomad_result);
  freeNomadProblem(pb);
  return status;
}

int solve_final_problem(const crs_nomad_problem *problem,
                        crs_nomad_eval_fn eval,
                        void *user_data,
                        crs_nomad_result *result) {
  int start_count = 1;
  const std::vector<double> starts = build_starting_points(problem, &start_count);
  crs_nomad_problem run_problem = *problem;
  if (problem->starts == nullptr && problem->start_count > 1) {
    run_problem.random_seed = 0;
  }
  return solve_final_problem_with_starts(&run_problem, starts.data(), start_count, eval, user_data, result);
}

} // namespace

extern "C" int crs_nomad_solve(
  const crs_nomad_problem *problem,
  crs_nomad_eval_fn eval,
  void *user_data,
  crs_nomad_result *result
) {
  if (result == nullptr) {
    return CRS_NOMAD_INVALID_INPUT;
  }
  if (result->api_version != CRS_NOMAD_API_VERSION) {
    return fail(result, CRS_NOMAD_INVALID_INPUT, "unsupported result api_version");
  }
  if (result->struct_size < sizeof(crs_nomad_result)) {
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
      status = solve_final_problem(problem, eval, user_data, result);
    }
  } catch (...) {
    status = fail(result, CRS_NOMAD_INTERNAL_ERROR, "unexpected exception in native NOMAD solve");
  }

  native_solver_busy.clear(std::memory_order_release);
  return status;
}
