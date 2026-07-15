#include "../inst/include/crs_nomad_native.h"

#include "NomadStdCInterface.h"
#include "Algos/Step.hpp"

#include <algorithm>
#include <atomic>
#include <cerrno>
#include <cctype>
#include <chrono>
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
#include <R_ext/Utils.h>

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

struct NativeObserverRuntime;
NativeObserverRuntime *active_observer_runtime = nullptr;

struct NativeEvalContext {
  int callback_mode;
  crs_nomad_eval_fn eval;
  void *user_data;
  const int *bb_output_type;
  int callback_evaluations;
  int callback_failures;
  bool user_interrupted;
  NativeObserverRuntime *observer_runtime;
};

struct NativeObserverRuntime {
  crs_nomad_observer *observer;
  NativeEvalContext *eval_context;
  const double *current_x;
  const double *current_outputs;
  int current_n;
  int current_m;
  int current_evaluation;
  bool outputs_available;
  bool active;
  std::chrono::steady_clock::time_point started;
  std::chrono::steady_clock::time_point last_emit;

  NativeObserverRuntime()
    : observer(nullptr),
      eval_context(nullptr),
      current_x(nullptr),
      current_outputs(nullptr),
      current_n(0),
      current_m(0),
      current_evaluation(0),
      outputs_available(false),
      active(false),
      started(),
      last_emit() {}
};

struct NativeSolveHandles {
  NomadProblem problem;
  NomadResult result;

  NativeSolveHandles() : problem(nullptr), result(nullptr) {}
};

void release_nomad_handles(NativeSolveHandles *handles) {
  if (handles == nullptr) {
    return;
  }
  if (handles->result != nullptr) {
    freeNomadResult(handles->result);
    handles->result = nullptr;
  }
  if (handles->problem != nullptr) {
    freeNomadProblem(handles->problem);
    handles->problem = nullptr;
  }
}

class NomadHandleGuard {
 public:
  explicit NomadHandleGuard(NativeSolveHandles *handles) : handles_(handles) {}
  ~NomadHandleGuard() {
    release_nomad_handles(handles_);
  }

  NomadHandleGuard(const NomadHandleGuard&) = delete;
  NomadHandleGuard& operator=(const NomadHandleGuard&) = delete;

 private:
  NativeSolveHandles *handles_;
};

struct CrsNomadSolveState {
  const crs_nomad_problem *problem;
  crs_nomad_eval_fn eval;
  void *user_data;
  crs_nomad_observer *observer;
  crs_nomad_result *result;
  NativeObserverRuntime observer_runtime;
  NativeSolveHandles handles;
  int status;
  bool busy_acquired;
  bool preserved_eval;
  bool preserved_env;
  SEXP eval_f;
  SEXP environment;

  CrsNomadSolveState(const crs_nomad_problem *problem_,
                     crs_nomad_eval_fn eval_,
                     void *user_data_,
                     crs_nomad_observer *observer_,
                     crs_nomad_result *result_)
    : problem(problem_),
      eval(eval_),
      user_data(user_data_),
      observer(observer_),
      result(result_),
      status(CRS_NOMAD_INTERNAL_ERROR),
      busy_acquired(false),
      preserved_eval(false),
      preserved_env(false),
      eval_f(R_NilValue),
      environment(R_NilValue) {}
};

void set_message(crs_nomad_result *result, const char *message) {
  if (result == nullptr) {
    return;
  }
  std::strncpy(result->message, message, sizeof(result->message) - 1);
  result->message[sizeof(result->message) - 1] = '\0';
}

int fail(crs_nomad_result *result, crs_nomad_status status, const char *message);

void set_observer_message(crs_nomad_observer *observer, const char *message) {
  if (observer == nullptr) {
    return;
  }
  std::strncpy(observer->message, message, sizeof(observer->message) - 1);
  observer->message[sizeof(observer->message) - 1] = '\0';
}

void initialize_observer(crs_nomad_observer *observer) {
  if (observer == nullptr) {
    return;
  }
  observer->state = CRS_NOMAD_OBSERVER_ACTIVE;
  observer->outcome = CRS_NOMAD_OBSERVER_OUTCOME_OK;
  observer->events_emitted = 0;
  observer->poll_calls = 0;
  set_observer_message(observer, "");
}

int validate_observer(crs_nomad_observer *observer,
                      crs_nomad_result *result) {
  if (observer == nullptr) {
    return CRS_NOMAD_OK;
  }
  if (observer->api_version != CRS_NOMAD_OBSERVER_API_VERSION) {
    return fail(result,
                CRS_NOMAD_INVALID_INPUT,
                "unsupported observer api_version");
  }
  if (observer->struct_size < sizeof(crs_nomad_observer)) {
    return fail(result,
                CRS_NOMAD_INVALID_INPUT,
                "observer struct_size is too small");
  }
  if (observer->observe == nullptr) {
    return fail(result,
                CRS_NOMAD_INVALID_INPUT,
                "observer callback is null");
  }
  if (!std::isfinite(observer->interval_sec) || observer->interval_sec < 0.0) {
    return fail(result,
                CRS_NOMAD_INVALID_INPUT,
                "observer interval_sec must be finite and nonnegative");
  }
  return CRS_NOMAD_OK;
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

std::string trim_copy(const std::string &value) {
  const std::size_t first = value.find_first_not_of(" \t\r\n");
  if (first == std::string::npos) {
    return std::string();
  }
  const std::size_t last = value.find_last_not_of(" \t\r\n");
  return value.substr(first, last - first + 1);
}

std::string strip_inline_comment(const std::string &value) {
  return trim_copy(value.substr(0, value.find('#')));
}

std::string ascii_upper_copy(const char *value) {
  std::string out = trim_copy(value == nullptr ? std::string() : std::string(value));
  for (std::size_t i = 0; i < out.size(); ++i) {
    out[i] = static_cast<char>(std::toupper(static_cast<unsigned char>(out[i])));
  }
  return out;
}

bool parse_positive_int_like(const char *value, int *out) {
  if (value == nullptr) {
    return false;
  }
  const std::string clean = trim_copy(value);
  if (clean.empty()) {
    return false;
  }
  errno = 0;
  char *endptr = nullptr;
  const double parsed = std::strtod(clean.c_str(), &endptr);
  if (endptr == clean.c_str() || *endptr != '\0' || errno == ERANGE ||
      !std::isfinite(parsed) || parsed <= 0.0 ||
      parsed > static_cast<double>(std::numeric_limits<int>::max())) {
    return false;
  }
  const double rounded = std::round(parsed);
  if (std::fabs(parsed - rounded) > 1e-9) {
    return false;
  }
  if (out != nullptr) {
    *out = static_cast<int>(rounded);
  }
  return true;
}

bool parse_bool_like(const char *value, bool *out) {
  if (value == nullptr) {
    return false;
  }
  const std::string clean = ascii_upper_copy(value);
  if (clean == "TRUE" || clean == "T" || clean == "YES" ||
      clean == "Y" || clean == "1") {
    if (out != nullptr) {
      *out = true;
    }
    return true;
  }
  if (clean == "FALSE" || clean == "F" || clean == "NO" ||
      clean == "N" || clean == "0") {
    if (out != nullptr) {
      *out = false;
    }
    return true;
  }
  return false;
}

std::string first_token_upper(const std::string &line) {
  const std::string clean = trim_copy(line);
  const std::size_t split = clean.find_first_of(" \t\r\n");
  return ascii_upper_copy((split == std::string::npos) ?
                          clean.c_str() :
                          clean.substr(0, split).c_str());
}

std::string second_token(const std::string &line) {
  const std::string clean = trim_copy(line);
  const std::size_t split = clean.find_first_of(" \t\r\n");
  if (split == std::string::npos) {
    return std::string();
  }
  return trim_copy(clean.substr(split + 1));
}

bool is_array_option_key(const std::string &key) {
  return key == "INITIAL_MESH_SIZE" || key == "MIN_MESH_SIZE" ||
         key == "INITIAL_FRAME_SIZE" || key == "MIN_FRAME_SIZE";
}

bool parse_mesh_frame_scalar(const std::string &token,
                             int index,
                             const crs_nomad_problem *problem,
                             double *value) {
  const std::string clean = trim_copy(token);
  if (clean.empty()) {
    return false;
  }
  const bool relative = clean[0] == 'r' || clean[0] == 'R';
  const char *start = relative ? clean.c_str() + 1 : clean.c_str();
  if (relative && *start == '\0') {
    return false;
  }
  errno = 0;
  char *endptr = nullptr;
  const double parsed = std::strtod(start, &endptr);
  if (endptr == start || *endptr != '\0' || errno == ERANGE ||
      !std::isfinite(parsed) || parsed <= 0.0) {
    return false;
  }

  double out = parsed;
  if (relative) {
    double span = 1.0;
    if (problem != nullptr && index >= 0 && index < problem->n &&
        std::isfinite(problem->lower[index]) &&
        std::isfinite(problem->upper[index])) {
      const double raw_span = problem->upper[index] - problem->lower[index];
      if (std::isfinite(raw_span) && raw_span > 0.0) {
        span = raw_span;
      }
    }
    out *= span;
  }

  if (value != nullptr) {
    *value = out;
  }
  return true;
}

bool parse_mesh_frame_values(const char *raw_value,
                             const crs_nomad_problem *problem,
                             std::vector<double> *values) {
  if (raw_value == nullptr || problem == nullptr || values == nullptr ||
      problem->n <= 0) {
    return false;
  }

  std::string clean = trim_copy(raw_value);
  if (clean.empty()) {
    return false;
  }
  if (!clean.empty() && clean[0] == '*') {
    clean = trim_copy(clean.substr(1));
  }
  if (clean.size() >= 2 && clean.front() == '(' && clean.back() == ')') {
    clean = trim_copy(clean.substr(1, clean.size() - 2));
  }
  if (clean.empty()) {
    return false;
  }

  std::istringstream iss(clean);
  std::vector<std::string> tokens;
  std::string token;
  while (iss >> token) {
    tokens.push_back(token);
  }
  if (tokens.empty()) {
    return false;
  }
  if (!(tokens.size() == 1 ||
        tokens.size() == static_cast<std::size_t>(problem->n))) {
    return false;
  }

  values->assign(static_cast<std::size_t>(problem->n), 0.0);
  if (tokens.size() == 1) {
    double scalar = 0.0;
    if (!parse_mesh_frame_scalar(tokens[0], 0, problem, &scalar)) {
      return false;
    }
    std::fill(values->begin(), values->end(), scalar);
  } else {
    for (int i = 0; i < problem->n; ++i) {
      double parsed = 0.0;
      if (!parse_mesh_frame_scalar(tokens[static_cast<std::size_t>(i)], i, problem, &parsed)) {
        return false;
      }
      (*values)[static_cast<std::size_t>(i)] = parsed;
    }
  }

  for (int i = 0; i < problem->n; ++i) {
    if (problem->bb_input_type[i] != CRS_NOMAD_INPUT_REAL &&
        (*values)[static_cast<std::size_t>(i)] < 1.0) {
      (*values)[static_cast<std::size_t>(i)] = 1.0;
    }
  }
  return true;
}

bool apply_mesh_frame_array_option(NomadProblem pb,
                                   const crs_nomad_problem *problem,
                                   const std::string &key,
                                   const char *raw_value) {
  std::vector<double> values;
  if (!parse_mesh_frame_values(raw_value, problem, &values)) {
    return false;
  }
  return addNomadArrayOfDoubleParam(pb, key.c_str(), values.data());
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
    return false;
  }
  if (status != 0) {
    return false;
  }
  return true;
}

struct RCallbackEvalState {
  int nb_inputs;
  double *x;
  int nb_outputs;
  double *bb_outputs;
  NativeEvalContext *context;
  bool ok;
  bool user_interrupted;
};

void native_eval_r_toplevel(void *data) {
  RCallbackEvalState *state = static_cast<RCallbackEvalState *>(data);
  state->ok = false;
  const int nb_inputs = state->nb_inputs;
  double *x = state->x;
  const int nb_outputs = state->nb_outputs;
  double *bb_outputs = state->bb_outputs;
  NativeEvalContext *context = state->context;

  R_CheckUserInterrupt();

  crs_nomad_r_callback *callback =
    static_cast<crs_nomad_r_callback *>(context->user_data);
  SEXP eval_f = (callback == nullptr) ?
    R_NilValue : reinterpret_cast<SEXP>(callback->eval_f);
  if (callback == nullptr || eval_f == nullptr || eval_f == R_NilValue ||
      !Rf_isFunction(eval_f)) {
    return;
  }

  int nprotect = 0;
  SEXP rargs = PROTECT(Rf_allocVector(REALSXP, nb_inputs));
  ++nprotect;
  for (int i = 0; i < nb_inputs; ++i) {
    REAL(rargs)[i] = x[i];
  }

  SEXP ns_name = PROTECT(Rf_mkString("crs"));
  ++nprotect;
  SEXP ns = PROTECT(R_FindNamespace(ns_name));
  ++nprotect;
  SEXP helper = PROTECT(Rf_findFun(Rf_install(".crs_nomad_native_eval_r"), ns));
  ++nprotect;
  if (helper == R_UnboundValue || helper == R_NilValue || !Rf_isFunction(helper)) {
    UNPROTECT(nprotect);
    return;
  }

  SEXP call = PROTECT(Rf_lang3(helper, eval_f, rargs));
  ++nprotect;

  int error = 0;
  SEXP wrapped = R_tryEvalSilent(call, ns, &error);
  if (error || wrapped == R_NilValue) {
    if (error) {
      R_CheckUserInterrupt();
    }
    UNPROTECT(nprotect);
    return;
  }
  wrapped = PROTECT(wrapped);
  ++nprotect;
  if (!Rf_isNewList(wrapped) || Rf_length(wrapped) < 2) {
    UNPROTECT(nprotect);
    return;
  }

  SEXP status_sexp = VECTOR_ELT(wrapped, 0);
  int wrapped_status = NA_INTEGER;
  if (Rf_isInteger(status_sexp) && Rf_length(status_sexp) >= 1) {
    wrapped_status = INTEGER(status_sexp)[0];
  } else if (Rf_isReal(status_sexp) && Rf_length(status_sexp) >= 1) {
    wrapped_status = static_cast<int>(REAL(status_sexp)[0]);
  }
  if (wrapped_status == 1) {
    state->user_interrupted = true;
    UNPROTECT(nprotect);
    return;
  }
  if (wrapped_status != 0) {
    UNPROTECT(nprotect);
    return;
  }

  SEXP result = VECTOR_ELT(wrapped, 1);

  if (!Rf_isNumeric(result)) {
    UNPROTECT(nprotect);
    return;
  }

  SEXP rnum = result;
  if (TYPEOF(rnum) != REALSXP) {
    rnum = PROTECT(Rf_coerceVector(result, REALSXP));
    ++nprotect;
  }

  if (Rf_length(rnum) < nb_outputs) {
    UNPROTECT(nprotect);
    return;
  }

  for (int i = 0; i < nb_outputs; ++i) {
    bb_outputs[i] = REAL(rnum)[i];
  }

  UNPROTECT(nprotect);
  state->ok = true;
}

bool native_eval_r(int nb_inputs,
                   double *x,
                   int nb_outputs,
                   double *bb_outputs,
                   NativeEvalContext *context) {
  RCallbackEvalState state;
  state.nb_inputs = nb_inputs;
  state.x = x;
  state.nb_outputs = nb_outputs;
  state.bb_outputs = bb_outputs;
  state.context = context;
  state.ok = false;
  state.user_interrupted = false;

  const bool toplevel_ok = R_ToplevelExec(native_eval_r_toplevel, &state);
  if ((!toplevel_ok || state.user_interrupted) && context != nullptr) {
    context->user_interrupted = true;
  }
  if (!toplevel_ok || !state.ok) {
    return false;
  }
  return true;
}

struct NativeObserverCallState {
  NativeObserverRuntime *runtime;
  int phase;
  int status;
  char message[256];
};

void native_observer_toplevel(void *data) {
  NativeObserverCallState *state = static_cast<NativeObserverCallState *>(data);
  state->status = CRS_NOMAD_OBSERVER_OUTCOME_ERROR;
  state->message[0] = '\0';
  if (state->runtime == nullptr || state->runtime->observer == nullptr) {
    std::strncpy(state->message,
                 "observer runtime is unavailable",
                 sizeof(state->message) - 1);
    state->message[sizeof(state->message) - 1] = '\0';
    return;
  }

  R_CheckUserInterrupt();

  crs_nomad_observer *observer = state->runtime->observer;
  if (observer->observe == nullptr) {
    std::strncpy(state->message,
                 "observer callback is unavailable",
                 sizeof(state->message) - 1);
    state->message[sizeof(state->message) - 1] = '\0';
    return;
  }

  crs_nomad_observer_event event;
  event.phase = state->phase;
  event.evaluation = state->runtime->current_evaluation;
  event.n = state->runtime->current_n;
  event.x = state->runtime->current_x;
  event.m = state->runtime->outputs_available ? state->runtime->current_m : 0;
  event.outputs = state->runtime->outputs_available ?
    state->runtime->current_outputs : nullptr;
  event.elapsed_sec = std::chrono::duration<double>(
    std::chrono::steady_clock::now() - state->runtime->started
  ).count();
  state->status = observer->observe(&event,
                                    observer->user_data,
                                    state->message,
                                    sizeof(state->message));
}

bool native_observer_due(const NativeObserverRuntime *runtime) {
  if (runtime == nullptr || runtime->observer == nullptr || !runtime->active ||
      runtime->observer->state != CRS_NOMAD_OBSERVER_ACTIVE ||
      runtime->current_x == nullptr || runtime->current_n <= 0 ||
      runtime->current_evaluation <= 0) {
    return false;
  }
  if (runtime->observer->interval_sec <= 0.0) {
    return true;
  }
  const double elapsed = std::chrono::duration<double>(
    std::chrono::steady_clock::now() - runtime->last_emit
  ).count();
  return elapsed >= runtime->observer->interval_sec;
}

bool native_observer_emit(NativeObserverRuntime *runtime, int phase) {
  if (!native_observer_due(runtime)) {
    return false;
  }

  NativeObserverCallState state;
  state.runtime = runtime;
  state.phase = phase;
  state.status = CRS_NOMAD_OBSERVER_OUTCOME_ERROR;
  state.message[0] = '\0';

  const bool toplevel_ok = R_ToplevelExec(native_observer_toplevel, &state);
  crs_nomad_observer *observer = runtime->observer;
  if (!toplevel_ok || state.status == CRS_NOMAD_OBSERVER_OUTCOME_INTERRUPT) {
    observer->state = CRS_NOMAD_OBSERVER_INTERRUPT_PENDING;
    observer->outcome = CRS_NOMAD_OBSERVER_OUTCOME_INTERRUPT;
    set_observer_message(observer,
                         state.message[0] == '\0' ?
                           "observer interrupted by user" : state.message);
    if (runtime->eval_context != nullptr) {
      runtime->eval_context->user_interrupted = true;
    }
    return false;
  }
  if (state.status != CRS_NOMAD_OBSERVER_OUTCOME_OK) {
    observer->state = CRS_NOMAD_OBSERVER_DISABLED_ERROR;
    observer->outcome = CRS_NOMAD_OBSERVER_OUTCOME_ERROR;
    set_observer_message(observer,
                         state.message[0] == '\0' ?
                           "observer callback failed" : state.message);
    return false;
  }

  ++observer->events_emitted;
  runtime->last_emit = std::chrono::steady_clock::now();
  return true;
}

bool mark_callback_failure(int nb_outputs,
                           double *bb_outputs,
                           bool *count_eval,
                           NativeEvalContext *context) {
  if (context != nullptr) {
    ++context->callback_failures;
  }
  if (bb_outputs != nullptr) {
    for (int i = 0; i < nb_outputs; ++i) {
      bb_outputs[i] = 1.0e6;
    }
  }
  if (count_eval != nullptr) {
    *count_eval = true;
  }
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

  if (context->callback_mode == CRS_NOMAD_CALLBACK_R &&
      NOMAD::Step::getUserInterrupt()) {
    context->user_interrupted = true;
    return false;
  }

  ++context->callback_evaluations;
  if (context->observer_runtime != nullptr && context->observer_runtime->active) {
    context->observer_runtime->eval_context = context;
    context->observer_runtime->current_x = x;
    context->observer_runtime->current_outputs = nullptr;
    context->observer_runtime->current_n = nb_inputs;
    context->observer_runtime->current_m = nb_outputs;
    context->observer_runtime->current_evaluation = context->callback_evaluations;
    context->observer_runtime->outputs_available = false;
    native_observer_emit(context->observer_runtime,
                         CRS_NOMAD_OBSERVER_EVALUATION_BEGIN);
    if (context->user_interrupted) {
      context->observer_runtime->current_x = nullptr;
      context->observer_runtime->eval_context = nullptr;
      return false;
    }
  }
  bool ok = false;
  if (context->callback_mode == CRS_NOMAD_CALLBACK_R) {
    ok = native_eval_r(nb_inputs, x, nb_outputs, bb_outputs, context);
  } else {
    if (context->eval == nullptr) {
      return mark_callback_failure(nb_outputs, bb_outputs, count_eval, context);
    }
    ok = native_eval_c(nb_inputs, x, nb_outputs, bb_outputs, context);
  }
  if (context->user_interrupted) {
    if (count_eval != nullptr) {
      *count_eval = false;
    }
    if (context->observer_runtime != nullptr) {
      context->observer_runtime->current_x = nullptr;
      context->observer_runtime->current_outputs = nullptr;
      context->observer_runtime->outputs_available = false;
      context->observer_runtime->eval_context = nullptr;
    }
    return false;
  }
  if (!ok) {
    if (context->observer_runtime != nullptr) {
      context->observer_runtime->current_x = nullptr;
      context->observer_runtime->current_outputs = nullptr;
      context->observer_runtime->outputs_available = false;
      context->observer_runtime->eval_context = nullptr;
    }
    return mark_callback_failure(nb_outputs, bb_outputs, count_eval, context);
  }

  for (int i = 0; i < nb_outputs; ++i) {
    if (std::isnan(bb_outputs[i])) {
      if (context->observer_runtime != nullptr) {
        context->observer_runtime->current_x = nullptr;
        context->observer_runtime->current_outputs = nullptr;
        context->observer_runtime->outputs_available = false;
        context->observer_runtime->eval_context = nullptr;
      }
      return mark_callback_failure(nb_outputs, bb_outputs, count_eval, context);
    }
    if (!std::isfinite(bb_outputs[i]) &&
        (context->bb_output_type == nullptr ||
         context->bb_output_type[i] == CRS_NOMAD_OUTPUT_OBJ)) {
      bb_outputs[i] = std::copysign(std::numeric_limits<double>::max(),
                                    bb_outputs[i]);
    }
  }
  if (context->observer_runtime != nullptr && context->observer_runtime->active) {
    context->observer_runtime->current_outputs = bb_outputs;
    context->observer_runtime->outputs_available = true;
    native_observer_emit(context->observer_runtime,
                         CRS_NOMAD_OBSERVER_EVALUATION_END);
    if (context->user_interrupted) {
      context->observer_runtime->current_x = nullptr;
      context->observer_runtime->current_outputs = nullptr;
      context->observer_runtime->outputs_available = false;
      context->observer_runtime->eval_context = nullptr;
      if (count_eval != nullptr) {
        *count_eval = false;
      }
      return false;
    }
    context->observer_runtime->current_x = nullptr;
    context->observer_runtime->current_outputs = nullptr;
    context->observer_runtime->outputs_available = false;
    context->observer_runtime->eval_context = nullptr;
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
    const std::string key = ascii_upper_copy(option->name);
    if (key == "MAX_EVAL") {
      seen_max_eval = true;
    }
    if (key == "MAX_BB_EVAL") {
      int parsed = 0;
      if (parse_positive_int_like(option->value, &parsed)) {
        seen_max_bb_eval = true;
        max_bb_eval_value = parsed;
      }
    }
    if (key == "EVAL_USE_CACHE") {
      bool parsed = true;
      if (parse_bool_like(option->value, &parsed) && !parsed) {
        return fail(result,
                    CRS_NOMAD_INVALID_INPUT,
                    "crs_nomad_solve currently requires EVAL_USE_CACHE = TRUE; cache-off native solves are not supported");
      }
    }
    if (key == "NB_THREADS_PARALLEL_EVAL" &&
        problem->callback_mode == CRS_NOMAD_CALLBACK_R) {
      int parsed = 0;
      if (parse_positive_int_like(option->value, &parsed) && parsed > 1) {
        return fail(result,
                    CRS_NOMAD_INVALID_INPUT,
                    "NB_THREADS_PARALLEL_EVAL > 1 is not supported for native R callbacks");
      }
    }
    if (is_array_option_key(key)) {
      if (!apply_mesh_frame_array_option(pb, problem, key, option->value)) {
        const std::string message =
          std::string("failed to set NOMAD array option ") + option->name;
        return fail(result, CRS_NOMAD_INVALID_INPUT, message.c_str());
      }
      continue;
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

int apply_nomad_opt_file(const crs_nomad_problem *problem,
                         NomadProblem pb,
                         crs_nomad_result *result) {
  std::ifstream fin("nomad.opt");
  if (!fin.good()) {
    return CRS_NOMAD_OK;
  }
  std::string line;
  while (std::getline(fin, line)) {
    const std::string clean = strip_inline_comment(line);
    if (clean.empty()) {
      continue;
    }
    if (problem != nullptr && problem->callback_mode == CRS_NOMAD_CALLBACK_R &&
        first_token_upper(clean) == "NB_THREADS_PARALLEL_EVAL") {
      int parsed = 0;
      const std::string value = second_token(clean);
      if (parse_positive_int_like(value.c_str(), &parsed) && parsed > 1) {
        return fail(result,
                    CRS_NOMAD_INVALID_INPUT,
                    "NB_THREADS_PARALLEL_EVAL > 1 is not supported for native R callbacks");
      }
    }
    if (problem != nullptr && first_token_upper(clean) == "EVAL_USE_CACHE") {
      bool parsed = true;
      const std::string value = second_token(clean);
      if (parse_bool_like(value.c_str(), &parsed) && !parsed) {
        return fail(result,
                    CRS_NOMAD_INVALID_INPUT,
                    "crs_nomad_solve currently requires EVAL_USE_CACHE = TRUE; cache-off native solves are not supported");
      }
    }
    if (!addNomadParam(pb, clean.c_str())) {
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
    const int file_status = apply_nomad_opt_file(problem, pb, result);
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
                                    crs_nomad_result *result,
                                    NativeSolveHandles *handles,
                                    NativeObserverRuntime *observer_runtime) {
  NativeSolveHandles local_handles;
  if (handles == nullptr) {
    handles = &local_handles;
  }
  NomadHandleGuard handle_guard(handles);

  handles->problem = createNomadProblem(native_eval_single, nullptr, problem->n, problem->m);
  NomadProblem pb = handles->problem;
  if (pb == nullptr) {
    return fail(result, CRS_NOMAD_INTERNAL_ERROR, "unable to initialize NOMAD problem");
  }

  int status = apply_problem_parameters(problem, pb, result);
  if (status != CRS_NOMAD_OK) {
    return status;
  }

  handles->result = createNomadResult();
  NomadResult nomad_result = handles->result;
  if (nomad_result == nullptr) {
    return fail(result, CRS_NOMAD_INTERNAL_ERROR, "unable to allocate NOMAD result");
  }

  NativeEvalContext context;
  context.callback_mode = problem->callback_mode;
  context.eval = eval;
  context.user_data = user_data;
  context.bb_output_type = problem->bb_output_type;
  context.callback_evaluations = 0;
  context.callback_failures = 0;
  context.user_interrupted = false;
  context.observer_runtime = observer_runtime;

  NOMAD::Step::resetUserInterrupt();
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
  if (context.user_interrupted) {
    NOMAD::Step::resetUserInterrupt();
    return fail(result, CRS_NOMAD_CALLBACK_FAILURE, "native R callback interrupted by user");
  }
  if (nb_solutions > 0) {
    std::vector<double> best_x(static_cast<std::size_t>(problem->n), 0.0);
    std::vector<double> best_out(static_cast<std::size_t>(problem->m), 0.0);
    if (!loadInputSolutionsNomadResult(best_x.data(), 1, nomad_result) ||
        !loadOutputSolutionsNomadResult(best_out.data(), 1, nomad_result)) {
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
    return fail(result, CRS_NOMAD_CALLBACK_FAILURE, "native callback reported evaluation failure");
  }

  if (nb_solutions > 0 && (run_flag == 1 || run_flag == 0 || run_flag == -6)) {
    result->status = CRS_NOMAD_OK;
    set_message(result, run_flag_message(run_flag));
    return CRS_NOMAD_OK;
  }

  const crs_nomad_status failure_status =
    (run_flag == -3 || run_flag == -5 || run_flag == -8) ?
      CRS_NOMAD_CALLBACK_FAILURE : CRS_NOMAD_INTERNAL_ERROR;
  status = fail(result, failure_status, run_flag_message(run_flag));
  return status;
}

int solve_final_problem(const crs_nomad_problem *problem,
                        crs_nomad_eval_fn eval,
                        void *user_data,
                        crs_nomad_result *result,
                        NativeSolveHandles *handles,
                        NativeObserverRuntime *observer_runtime) {
  int start_count = 1;
  const std::vector<double> starts = build_starting_points(problem, &start_count);
  crs_nomad_problem run_problem = *problem;
  if (problem->starts == nullptr && problem->start_count > 1) {
    run_problem.random_seed = 0;
  }
  return solve_final_problem_with_starts(&run_problem,
                                         starts.data(),
                                         start_count,
                                         eval,
                                         user_data,
                                         result,
                                         handles,
                                         observer_runtime);
}

void preserve_r_callback_objects(CrsNomadSolveState *state) {
  if (state == nullptr) {
    return;
  }

  if (state->problem != nullptr &&
      state->problem->callback_mode == CRS_NOMAD_CALLBACK_R &&
      state->user_data != nullptr) {
    crs_nomad_r_callback *callback =
      static_cast<crs_nomad_r_callback *>(state->user_data);
    state->eval_f = reinterpret_cast<SEXP>(callback->eval_f);
    state->environment = reinterpret_cast<SEXP>(callback->environment);
    if (state->eval_f != nullptr && state->eval_f != R_NilValue) {
      R_PreserveObject(state->eval_f);
      state->preserved_eval = true;
    }
    if (state->environment != nullptr && state->environment != R_NilValue) {
      R_PreserveObject(state->environment);
      state->preserved_env = true;
    }
  }

}

void release_r_callback_objects(CrsNomadSolveState *state) {
  if (state == nullptr) {
    return;
  }
  if (state->preserved_env) {
    R_ReleaseObject(state->environment);
    state->preserved_env = false;
  }
  if (state->preserved_eval) {
    R_ReleaseObject(state->eval_f);
    state->preserved_eval = false;
  }
}

SEXP crs_nomad_solve_unwind_body(void *data) {
  CrsNomadSolveState *state = static_cast<CrsNomadSolveState *>(data);
  state->status = CRS_NOMAD_INTERNAL_ERROR;
  try {
    state->status = validate_problem(state->problem, state->eval, state->result);
    if (state->status == CRS_NOMAD_OK) {
      state->status = solve_final_problem(state->problem,
                                          state->eval,
                                          state->user_data,
                                          state->result,
                                          &state->handles,
                                          state->observer == nullptr ?
                                            nullptr : &state->observer_runtime);
    }
  } catch (...) {
    state->status = fail(state->result,
                         CRS_NOMAD_INTERNAL_ERROR,
                         "unexpected exception in native NOMAD solve");
  }
  return R_NilValue;
}

void crs_nomad_solve_unwind_cleanup(void *data, Rboolean jump) {
  (void) jump;
  CrsNomadSolveState *state = static_cast<CrsNomadSolveState *>(data);
  if (active_observer_runtime == &state->observer_runtime) {
    active_observer_runtime = nullptr;
  }
  state->observer_runtime.active = false;
  state->observer_runtime.eval_context = nullptr;
  state->observer_runtime.current_x = nullptr;
  state->observer_runtime.current_outputs = nullptr;
  state->observer_runtime.outputs_available = false;
  if (state->observer != nullptr) {
    state->observer->state = CRS_NOMAD_OBSERVER_CLOSED;
  }
  release_nomad_handles(&state->handles);
  release_r_callback_objects(state);
  if (state->busy_acquired) {
    native_solver_busy.clear(std::memory_order_release);
    state->busy_acquired = false;
  }
}

SEXP list_element(SEXP list, const char *name) {
  if (list == R_NilValue || !Rf_isNewList(list)) {
    return R_NilValue;
  }
  SEXP names = Rf_getAttrib(list, R_NamesSymbol);
  if (names == R_NilValue) {
    return R_NilValue;
  }
  if (XLENGTH(names) != XLENGTH(list)) {
    return R_NilValue;
  }
  const R_xlen_t n = XLENGTH(list);
  for (R_xlen_t i = 0; i < n; ++i) {
    if (std::strcmp(CHAR(STRING_ELT(names, i)), name) == 0) {
      return VECTOR_ELT(list, i);
    }
  }
  return R_NilValue;
}

int scalar_int_or(SEXP value, int fallback) {
  if (value == R_NilValue || XLENGTH(value) < 1) {
    return fallback;
  }
  return Rf_asInteger(value);
}

bool scalar_logical_or(SEXP value, bool fallback) {
  if (value == R_NilValue || XLENGTH(value) < 1) {
    return fallback;
  }
  const int v = Rf_asLogical(value);
  return v == NA_LOGICAL ? fallback : (v != 0);
}

std::string scalar_string_or(SEXP value, const char *fallback) {
  if (value == R_NilValue || XLENGTH(value) < 1) {
    return std::string(fallback);
  }
  SEXP chr = Rf_asChar(value);
  if (chr == NA_STRING) {
    return std::string(fallback);
  }
  return std::string(CHAR(chr));
}

struct NativeTestEvalData {
  std::string scenario;
  double target;
  int evaluations;
};

struct NativeTestObserverData {
  SEXP observe_f;
  SEXP environment;
};

int native_test_observe(const crs_nomad_observer_event *event,
                        void *user_data,
                        char *message,
                        size_t message_size) {
  NativeTestObserverData *data =
    static_cast<NativeTestObserverData *>(user_data);
  if (event == nullptr || data == nullptr ||
      data->observe_f == R_NilValue || !Rf_isFunction(data->observe_f)) {
    if (message != nullptr && message_size > 0) {
      std::snprintf(message, message_size, "native test observer is unavailable");
    }
    return CRS_NOMAD_OBSERVER_OUTCOME_ERROR;
  }

  const char *phase = "activity";
  if (event->phase == CRS_NOMAD_OBSERVER_EVALUATION_BEGIN) {
    phase = "evaluation_begin";
  } else if (event->phase == CRS_NOMAD_OBSERVER_EVALUATION_END) {
    phase = "evaluation_end";
  }

  int nprotect = 0;
  SEXP event_r = PROTECT(Rf_allocVector(VECSXP, 5)); ++nprotect;
  SEXP names = PROTECT(Rf_allocVector(STRSXP, 5)); ++nprotect;
  SEXP x = PROTECT(Rf_allocVector(REALSXP, event->n)); ++nprotect;
  for (int i = 0; i < event->n; ++i) {
    REAL(x)[i] = event->x[i];
  }
  SEXP outputs = R_NilValue;
  if (event->outputs != nullptr && event->m > 0) {
    outputs = PROTECT(Rf_allocVector(REALSXP, event->m)); ++nprotect;
    for (int i = 0; i < event->m; ++i) {
      REAL(outputs)[i] = event->outputs[i];
    }
  }
  SET_VECTOR_ELT(event_r, 0, Rf_mkString(phase));
  SET_VECTOR_ELT(event_r, 1, Rf_ScalarInteger(event->evaluation));
  SET_VECTOR_ELT(event_r, 2, x);
  SET_VECTOR_ELT(event_r, 3, outputs);
  SET_VECTOR_ELT(event_r, 4, Rf_ScalarReal(event->elapsed_sec));
  const char *event_names[] = {"phase", "evaluation", "x", "outputs", "elapsed"};
  for (int i = 0; i < 5; ++i) {
    SET_STRING_ELT(names, i, Rf_mkChar(event_names[i]));
  }
  Rf_setAttrib(event_r, R_NamesSymbol, names);

  SEXP call = PROTECT(Rf_lang2(data->observe_f, event_r)); ++nprotect;
  int error = 0;
  SEXP value = R_tryEvalSilent(call, data->environment, &error);
  if (error) {
    if (message != nullptr && message_size > 0) {
      std::snprintf(message,
                    message_size,
                    "native test observer evaluation failed");
    }
    UNPROTECT(nprotect);
    return CRS_NOMAD_OBSERVER_OUTCOME_ERROR;
  }
  value = PROTECT(value); ++nprotect;

  int status = CRS_NOMAD_OBSERVER_OUTCOME_OK;
  SEXP status_value = value;
  SEXP message_value = R_NilValue;
  if (Rf_isNewList(value) && XLENGTH(value) >= 1) {
    status_value = VECTOR_ELT(value, 0);
    if (XLENGTH(value) >= 2) {
      message_value = VECTOR_ELT(value, 1);
    }
  }
  if ((TYPEOF(status_value) == INTSXP || TYPEOF(status_value) == REALSXP) &&
      XLENGTH(status_value) >= 1) {
    const int candidate = Rf_asInteger(status_value);
    if (candidate >= CRS_NOMAD_OBSERVER_OUTCOME_OK &&
        candidate <= CRS_NOMAD_OBSERVER_OUTCOME_INTERRUPT) {
      status = candidate;
    }
  }
  if (message != nullptr && message_size > 0 &&
      message_value != R_NilValue && XLENGTH(message_value) >= 1) {
    SEXP message_char = Rf_asChar(message_value);
    if (message_char != NA_STRING) {
      std::snprintf(message, message_size, "%s", CHAR(message_char));
    }
  }
  UNPROTECT(nprotect);
  return status;
}

int native_test_eval(int n,
                     const double *x,
                     int m,
                     double *bb_outputs,
                     void *user_data) {
  NativeTestEvalData *data = static_cast<NativeTestEvalData *>(user_data);
  if (data == nullptr || x == nullptr || bb_outputs == nullptr || n <= 0 || m <= 0) {
    return 1;
  }
  ++data->evaluations;
  if (data->scenario == "poll") {
    crs_nomad_observer_poll();
  } else if (data->scenario == "many_polls") {
    for (int i = 0; i < 8; ++i) {
      crs_nomad_observer_poll();
    }
  }
  double objective = 0.0;
  for (int i = 0; i < n; ++i) {
    const double diff = x[i] - data->target;
    objective += diff * diff;
  }
  for (int i = 0; i < m; ++i) {
    bb_outputs[i] = 0.0;
  }
  bb_outputs[0] = objective;
  if (data->scenario == "nan_obj") {
    bb_outputs[0] = std::numeric_limits<double>::quiet_NaN();
  } else if (data->scenario == "inf_obj") {
    bb_outputs[0] = std::numeric_limits<double>::infinity();
  } else if (data->scenario == "negative_inf_obj") {
    bb_outputs[0] = -std::numeric_limits<double>::infinity();
  } else if (data->scenario == "inf_constraint" && m > 1) {
    bb_outputs[1] = std::numeric_limits<double>::infinity();
  } else if (data->scenario == "negative_inf_constraint" && m > 1) {
    bb_outputs[1] = -std::numeric_limits<double>::infinity();
  } else if (data->scenario == "shifted") {
    bb_outputs[0] = objective + 10.0;
  }
  return 0;
}

} // namespace

int crs_nomad_solve_impl(
  const crs_nomad_problem *problem,
  crs_nomad_eval_fn eval,
  void *user_data,
  crs_nomad_observer *observer,
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
  const int observer_status = validate_observer(observer, result);
  if (observer_status != CRS_NOMAD_OK) {
    return observer_status;
  }
  initialize_observer(observer);

  CrsNomadSolveState state(problem, eval, user_data, observer, result);
  preserve_r_callback_objects(&state);

  if (native_solver_busy.test_and_set(std::memory_order_acquire)) {
    release_r_callback_objects(&state);
    if (observer != nullptr) {
      observer->state = CRS_NOMAD_OBSERVER_CLOSED;
    }
    return fail(result, CRS_NOMAD_INVALID_INPUT, "native NOMAD solver is already active");
  }
  state.busy_acquired = true;

  if (observer != nullptr) {
    state.observer_runtime.observer = observer;
    state.observer_runtime.active = true;
    state.observer_runtime.started = std::chrono::steady_clock::now();
    state.observer_runtime.last_emit = state.observer_runtime.started;
    active_observer_runtime = &state.observer_runtime;
  }

  SEXP unwind_token = PROTECT(R_MakeUnwindCont());
  R_UnwindProtect(crs_nomad_solve_unwind_body,
                  &state,
                  crs_nomad_solve_unwind_cleanup,
                  &state,
                  unwind_token);
  SETCAR(unwind_token, R_NilValue);
  UNPROTECT(1);
  return state.status;
}

extern "C" int crs_nomad_solve(
  const crs_nomad_problem *problem,
  crs_nomad_eval_fn eval,
  void *user_data,
  crs_nomad_result *result
) {
  return crs_nomad_solve_impl(problem, eval, user_data, nullptr, result);
}

extern "C" int crs_nomad_solve_observed(
  const crs_nomad_problem *problem,
  crs_nomad_eval_fn eval,
  void *user_data,
  crs_nomad_observer *observer,
  crs_nomad_result *result
) {
  return crs_nomad_solve_impl(problem, eval, user_data, observer, result);
}

extern "C" int crs_nomad_observer_poll(void) {
  NativeObserverRuntime *runtime = active_observer_runtime;
  if (runtime == nullptr || !runtime->active || runtime->observer == nullptr) {
    return 0;
  }
  ++runtime->observer->poll_calls;
  if (runtime->current_x == nullptr || runtime->current_n <= 0 ||
      runtime->current_evaluation <= 0) {
    return 0;
  }
  return native_observer_emit(runtime, CRS_NOMAD_OBSERVER_ACTIVITY) ? 1 : 0;
}

extern "C" SEXP crs_nomad_native_test_solve(SEXP spec) {
  if (!Rf_isNewList(spec)) {
    Rf_error("crs_nomad_native_test_solve requires a list spec");
  }

  const std::string mode = scalar_string_or(list_element(spec, "mode"), "c");
  const std::string scenario = scalar_string_or(list_element(spec, "scenario"), "quadratic");

  SEXP x0_s = PROTECT(Rf_coerceVector(list_element(spec, "x0"), REALSXP));
  SEXP lower_s = PROTECT(Rf_coerceVector(list_element(spec, "lower"), REALSXP));
  SEXP upper_s = PROTECT(Rf_coerceVector(list_element(spec, "upper"), REALSXP));
  SEXP input_type_s = PROTECT(Rf_coerceVector(list_element(spec, "input_type"), INTSXP));
  SEXP output_type_s = PROTECT(Rf_coerceVector(list_element(spec, "output_type"), INTSXP));
  int nprotect = 5;

  const R_xlen_t n_x = XLENGTH(x0_s);
  const R_xlen_t m_x = XLENGTH(output_type_s);
  if (n_x <= 0 || m_x <= 0 ||
      XLENGTH(lower_s) != n_x ||
      XLENGTH(upper_s) != n_x ||
      XLENGTH(input_type_s) != n_x ||
      n_x > std::numeric_limits<int>::max() ||
      m_x > std::numeric_limits<int>::max()) {
    UNPROTECT(nprotect);
    Rf_error("crs_nomad_native_test_solve received inconsistent dimensions");
  }
  const int n = static_cast<int>(n_x);
  const int m = static_cast<int>(m_x);

  std::vector<std::string> option_names;
  std::vector<std::string> option_values;
  std::vector<crs_nomad_option> options;
  SEXP options_s = list_element(spec, "options");
  if (options_s != R_NilValue) {
    if (TYPEOF(options_s) != STRSXP) {
      UNPROTECT(nprotect);
      Rf_error("crs_nomad_native_test_solve options must be a named character vector");
    }
    const R_xlen_t n_options = XLENGTH(options_s);
    if (n_options == 0) {
      options_s = R_NilValue;
    }
  }
  if (options_s != R_NilValue) {
    SEXP names_s = Rf_getAttrib(options_s, R_NamesSymbol);
    if (names_s == R_NilValue || XLENGTH(names_s) != XLENGTH(options_s)) {
      UNPROTECT(nprotect);
      Rf_error("crs_nomad_native_test_solve options must be named");
    }
    const R_xlen_t n_options = XLENGTH(options_s);
    option_names.reserve(static_cast<std::size_t>(n_options));
    option_values.reserve(static_cast<std::size_t>(n_options));
    options.reserve(static_cast<std::size_t>(n_options));
    for (R_xlen_t i = 0; i < n_options; ++i) {
      option_names.push_back(CHAR(STRING_ELT(names_s, i)));
      option_values.push_back(CHAR(STRING_ELT(options_s, i)));
    }
    for (std::size_t i = 0; i < option_names.size(); ++i) {
      crs_nomad_option option;
      option.name = option_names[i].c_str();
      option.value = option_values[i].c_str();
      options.push_back(option);
    }
  }

  std::vector<double> solution(static_cast<std::size_t>(n), R_NaN);
  std::vector<double> outputs(static_cast<std::size_t>(m), R_NaN);

  crs_nomad_problem problem;
  std::memset(&problem, 0, sizeof(problem));
  problem.api_version = CRS_NOMAD_API_VERSION;
  problem.struct_size = sizeof(problem);
  problem.callback_mode = (mode == "r") ? CRS_NOMAD_CALLBACK_R : CRS_NOMAD_CALLBACK_C;
  problem.n = n;
  problem.m = m;
  problem.x0 = REAL(x0_s);
  problem.bb_input_type = INTEGER(input_type_s);
  problem.bb_output_type = INTEGER(output_type_s);
  problem.lower = REAL(lower_s);
  problem.upper = REAL(upper_s);
  problem.max_eval = scalar_int_or(list_element(spec, "max_eval"), 20);
  problem.random_seed = static_cast<unsigned int>(
    std::max(0, scalar_int_or(list_element(spec, "random_seed"), 42))
  );
  problem.quiet = scalar_logical_or(list_element(spec, "quiet"), true) ? 1 : 0;
  problem.option_count = static_cast<int>(options.size());
  problem.options = options.empty() ? nullptr : options.data();
  problem.start_count = scalar_int_or(list_element(spec, "start_count"), 0);
  problem.starts = nullptr;
  problem.read_nomad_opt_file = scalar_logical_or(list_element(spec, "read_nomad_opt_file"), false) ? 1 : 0;

  crs_nomad_result result;
  std::memset(&result, 0, sizeof(result));
  result.api_version = CRS_NOMAD_API_VERSION;
  result.struct_size = sizeof(result);
  result.solution = solution.data();
  result.solution_len = n;
  result.outputs = outputs.data();
  result.outputs_len = m;

  NativeTestEvalData test_data;
  test_data.scenario = scenario;
  SEXP target_s = list_element(spec, "target");
  test_data.target = (target_s == R_NilValue || XLENGTH(target_s) < 1) ?
    0.25 : Rf_asReal(target_s);
  if (!std::isfinite(test_data.target)) {
    test_data.target = 0.25;
  }
  test_data.evaluations = 0;

  crs_nomad_r_callback r_callback;
  std::memset(&r_callback, 0, sizeof(r_callback));
  crs_nomad_eval_fn eval = native_test_eval;
  void *user_data = &test_data;

  if (mode == "r") {
    SEXP eval_f = list_element(spec, "eval_f");
    if (!Rf_isFunction(eval_f)) {
      UNPROTECT(nprotect);
      Rf_error("crs_nomad_native_test_solve R mode requires eval_f");
    }
    SEXP eval_env = list_element(spec, "eval_env");
    if (eval_env == R_NilValue || TYPEOF(eval_env) != ENVSXP) {
      eval_env = R_GlobalEnv;
    }
    r_callback.eval_f = reinterpret_cast<void *>(eval_f);
    r_callback.environment = reinterpret_cast<void *>(eval_env);
    eval = nullptr;
    user_data = &r_callback;
  }

  crs_nomad_observer observer;
  std::memset(&observer, 0, sizeof(observer));
  NativeTestObserverData observer_data;
  observer_data.observe_f = R_NilValue;
  observer_data.environment = R_GlobalEnv;
  bool observed = false;
  SEXP observe_f = list_element(spec, "observe_f");
  if (Rf_isFunction(observe_f)) {
    SEXP observe_env = list_element(spec, "observe_env");
    if (observe_env == R_NilValue || TYPEOF(observe_env) != ENVSXP) {
      observe_env = R_GlobalEnv;
    }
    SEXP interval_s = list_element(spec, "observe_interval");
    double interval = (interval_s == R_NilValue || XLENGTH(interval_s) < 1) ?
      0.0 : Rf_asReal(interval_s);
    observer.api_version = CRS_NOMAD_OBSERVER_API_VERSION;
    observer.struct_size = sizeof(observer);
    observer_data.observe_f = observe_f;
    observer_data.environment = observe_env;
    observer.observe = native_test_observe;
    observer.user_data = &observer_data;
    observer.interval_sec = interval;
    observed = true;
  }

  const int status = observed ?
    crs_nomad_solve_observed(&problem, eval, user_data, &observer, &result) :
    crs_nomad_solve(&problem, eval, user_data, &result);

  SEXP out = PROTECT(Rf_allocVector(VECSXP, 19)); ++nprotect;
  SEXP names = PROTECT(Rf_allocVector(STRSXP, 19)); ++nprotect;
  SEXP sol = PROTECT(Rf_allocVector(REALSXP, n)); ++nprotect;
  SEXP bbout = PROTECT(Rf_allocVector(REALSXP, m)); ++nprotect;
  for (int i = 0; i < n; ++i) {
    REAL(sol)[i] = solution[static_cast<std::size_t>(i)];
  }
  for (int i = 0; i < m; ++i) {
    REAL(bbout)[i] = outputs[static_cast<std::size_t>(i)];
  }

  SET_VECTOR_ELT(out, 0, Rf_ScalarInteger(status));
  SET_VECTOR_ELT(out, 1, Rf_ScalarInteger(result.status));
  SET_VECTOR_ELT(out, 2, Rf_mkString(result.message));
  SET_VECTOR_ELT(out, 3, Rf_ScalarInteger(result.nomad_run_flag));
  SET_VECTOR_ELT(out, 4, Rf_ScalarReal(result.objective));
  SET_VECTOR_ELT(out, 5, sol);
  SET_VECTOR_ELT(out, 6, bbout);
  SET_VECTOR_ELT(out, 7, Rf_ScalarInteger(result.blackbox_evaluations));
  SET_VECTOR_ELT(out, 8, Rf_ScalarInteger(result.callback_evaluations));
  SET_VECTOR_ELT(out, 9, Rf_ScalarInteger(result.cache_hits));
  SET_VECTOR_ELT(out, 10, Rf_ScalarInteger(result.cache_size));
  SET_VECTOR_ELT(out, 11, Rf_ScalarInteger(result.total_evaluations));
  SET_VECTOR_ELT(out, 12, Rf_ScalarInteger(result.solution_count));
  SET_VECTOR_ELT(out, 13, Rf_ScalarInteger(test_data.evaluations));
  SET_VECTOR_ELT(out, 14, Rf_ScalarInteger(observed ? observer.state : CRS_NOMAD_OBSERVER_OFF));
  SET_VECTOR_ELT(out, 15, Rf_ScalarInteger(observed ? observer.outcome : CRS_NOMAD_OBSERVER_OUTCOME_OK));
  SET_VECTOR_ELT(out, 16, Rf_ScalarInteger(observed ? observer.events_emitted : 0));
  SET_VECTOR_ELT(out, 17, Rf_ScalarInteger(observed ? observer.poll_calls : 0));
  SET_VECTOR_ELT(out, 18, Rf_mkString(observed ? observer.message : ""));

  const char *out_names[] = {
    "status",
    "result_status",
    "message",
    "nomad_run_flag",
    "objective",
    "solution",
    "outputs",
    "blackbox_evaluations",
    "callback_evaluations",
    "cache_hits",
    "cache_size",
    "total_evaluations",
    "solution_count",
    "test_evaluations",
    "observer_state",
    "observer_outcome",
    "observer_events",
    "observer_polls",
    "observer_message"
  };
  for (int i = 0; i < 19; ++i) {
    SET_STRING_ELT(names, i, Rf_mkChar(out_names[i]));
  }
  Rf_setAttrib(out, R_NamesSymbol, names);

  UNPROTECT(nprotect);
  return out;
}
