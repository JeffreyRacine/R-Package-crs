#include "snomadr.h"

#include "NomadStdCInterface.h"
#include "Algos/EvcInterface.hpp"
#include "nomad_version.hpp"

#include <R_ext/Boolean.h>

#include <algorithm>
#include <atomic>
#include <cctype>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <fstream>
#include <limits>
#include <ostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

static SEXP thefun = R_NilValue;
static SEXP theenv = R_NilValue;
static Routbuf routbuf;
static std::ostream rout(&routbuf);
static std::atomic<int> g_bbe_counter{0};

static SEXP getListElement(SEXP list, const char* str) {
  SEXP names = Rf_getAttrib(list, R_NamesSymbol);
  const int n = Rf_length(list);
  for (int i = 0; i < n; ++i) {
    if (std::strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
      return VECTOR_ELT(list, i);
    }
  }
  return R_NilValue;
}

static std::string trim(const std::string& s) {
  std::size_t b = 0;
  while (b < s.size() && std::isspace(static_cast<unsigned char>(s[b]))) {
    ++b;
  }
  std::size_t e = s.size();
  while (e > b && std::isspace(static_cast<unsigned char>(s[e - 1]))) {
    --e;
  }
  return s.substr(b, e - b);
}

static bool r_eval_single(int nb_inputs,
                          double* x,
                          int nb_outputs,
                          double* bb_outputs,
                          bool* count_eval,
                          NomadUserDataPtr) {
  R_CheckUserInterrupt();

  int nprotect = 0;
  SEXP rargs = PROTECT(Rf_allocVector(REALSXP, nb_inputs));
  ++nprotect;
  for (int i = 0; i < nb_inputs; ++i) {
    REAL(rargs)[i] = x[i];
  }

  SEXP call = PROTECT(Rf_lang2(thefun, rargs));
  ++nprotect;

  int error = 0;
  SEXP result = R_tryEval(call, theenv, &error);
  if (error || result == R_NilValue) {
    UNPROTECT(nprotect);
    return false;
  }

  if (!Rf_isNumeric(result)) {
    UNPROTECT(nprotect);
    return false;
  }

  SEXP rnum = result;
  if (TYPEOF(rnum) != REALSXP) {
    rnum = PROTECT(Rf_coerceVector(result, REALSXP));
    ++nprotect;
  }

  if (Rf_length(rnum) < nb_outputs) {
    UNPROTECT(nprotect);
    return false;
  }

  for (int i = 0; i < nb_outputs; ++i) {
    bb_outputs[i] = REAL(rnum)[i];
  }
  *count_eval = true;
  g_bbe_counter.fetch_add(1, std::memory_order_relaxed);

  UNPROTECT(nprotect);
  return true;
}

static std::string map_bbin_to_string(SEXP sbbin) {
  std::ostringstream oss;
  const int n = Rf_length(sbbin);
  oss << "(";
  bool warned_categorical = false;
  for (int i = 0; i < n; ++i) {
    const int type = INTEGER(sbbin)[i];
    oss << ' ';
    switch (type) {
      case 0:
        oss << "R";
        break;
      case 1:
        oss << "I";
        break;
      case 3:
        oss << "B";
        break;
      case 2:
      default:
        // NOMAD4 C API does not expose categorical input type; degrade to integer.
        if (type == 2 && !warned_categorical) {
          Rprintf("Warning: categorical input type is mapped to integer in embedded NOMAD4.\n");
          warned_categorical = true;
        }
        oss << "I";
        break;
    }
  }
  oss << " )";
  return oss.str();
}

static std::string map_bbout_to_string(SEXP sbbout) {
  std::ostringstream oss;
  const int n = Rf_length(sbbout);
  for (int i = 0; i < n; ++i) {
    if (i > 0) {
      oss << ' ';
    }
    switch (INTEGER(sbbout)[i]) {
      case 0:
        oss << "OBJ";
        break;
      case 1:
        oss << "PB";
        break;
      case 2:
        oss << "EB";
        break;
      default:
        oss << "BBO_UNDEFINED";
        break;
    }
  }
  return oss.str();
}

static bool apply_named_option_line(NomadProblem pb, const std::string& key, const std::string& value) {
  std::string line = key + " " + value;
  return addNomadParam(pb, line.c_str());
}

static bool is_array_option_key(const std::string& key) {
  return key == "INITIAL_MESH_SIZE" || key == "MIN_MESH_SIZE" || key == "MIN_POLL_SIZE" ||
         key == "INITIAL_POLL_SIZE" || key == "INITIAL_FRAME_SIZE" || key == "MIN_FRAME_SIZE";
}

static bool is_ignored_legacy_option_key(const std::string& key) {
  (void)key;
  return false;
}

static std::string canonical_option_key(const std::string& key) {
  if (key == "MIN_POLL_SIZE") {
    return "MIN_FRAME_SIZE";
  }
  if (key == "INITIAL_POLL_SIZE") {
    return "INITIAL_FRAME_SIZE";
  }
  return key;
}

static bool is_mesh_or_frame_option_key(const std::string& key) {
  return key == "INITIAL_MESH_SIZE" || key == "MIN_MESH_SIZE" || key == "INITIAL_FRAME_SIZE" ||
         key == "MIN_FRAME_SIZE";
}

static std::string normalize_array_option_value(const std::string& raw) {
  const std::string v = trim(raw);
  if (v.empty()) {
    return v;
  }
  if (v[0] == '(' || v[0] == '*') {
    return v;
  }
  if (v.find_first_of(" \t") == std::string::npos) {
    return std::string("* ") + v;
  }
  return v;
}

static std::string array_values_to_literal(const std::vector<double>& values) {
  std::ostringstream oss;
  oss.precision(17);
  oss << "(";
  for (std::size_t i = 0; i < values.size(); ++i) {
    oss << " " << values[i];
  }
  oss << " )";
  return oss.str();
}

static bool parse_scalar_option_value(const std::string& token,
                                      int index,
                                      const std::vector<double>& lb,
                                      const std::vector<double>& ub,
                                      double* value) {
  const std::string t = trim(token);
  if (t.empty()) {
    return false;
  }
  if (t[0] == 'r' || t[0] == 'R') {
    if (t.size() == 1) {
      return false;
    }
    char* endptr = nullptr;
    const double rel = std::strtod(t.c_str() + 1, &endptr);
    if (endptr == t.c_str() + 1 || *endptr != '\0') {
      return false;
    }
    double span = 1.0;
    if (index >= 0 && static_cast<std::size_t>(index) < lb.size() &&
        static_cast<std::size_t>(index) < ub.size() && std::isfinite(lb[static_cast<std::size_t>(index)]) &&
        std::isfinite(ub[static_cast<std::size_t>(index)])) {
      span = ub[static_cast<std::size_t>(index)] - lb[static_cast<std::size_t>(index)];
      if (!std::isfinite(span) || span <= 0.0) {
        span = 1.0;
      }
    }
    *value = rel * span;
    return true;
  }
  char* endptr = nullptr;
  const double parsed = std::strtod(t.c_str(), &endptr);
  if (endptr == t.c_str() || *endptr != '\0') {
    return false;
  }
  *value = parsed;
  return true;
}

static void enforce_integer_mesh_frame_floor(std::vector<double>& values,
                                             const std::string& canonical_key,
                                             const std::vector<int>& bbin) {
  if (!is_mesh_or_frame_option_key(canonical_key)) {
    return;
  }
  const std::size_t n = std::min(values.size(), bbin.size());
  for (std::size_t i = 0; i < n; ++i) {
    if (bbin[i] == 1 || bbin[i] == 2 || bbin[i] == 3) {
      if (values[i] < 1.0) {
        values[i] = 1.0;
      }
    }
  }
}

static bool apply_array_option_values(NomadProblem pb,
                                      const std::string& key,
                                      const std::vector<double>& raw_values,
                                      int n_inputs,
                                      const std::vector<int>& bbin) {
  if (n_inputs <= 0 || raw_values.empty()) {
    return false;
  }
  std::vector<double> values(static_cast<std::size_t>(n_inputs), raw_values[0]);
  if (raw_values.size() == static_cast<std::size_t>(n_inputs)) {
    values = raw_values;
  }
  const std::string canonical_key = canonical_option_key(key);
  enforce_integer_mesh_frame_floor(values, canonical_key, bbin);
  if (addNomadArrayOfDoubleParam(pb, canonical_key.c_str(), values.data())) {
    return true;
  }
  return apply_named_option_line(pb, canonical_key, array_values_to_literal(values));
}

static bool apply_array_option_from_string(NomadProblem pb,
                                           const std::string& key,
                                           SEXP val_sexp,
                                           int n_inputs,
                                           const std::vector<int>& bbin,
                                           const std::vector<double>& lb,
                                           const std::vector<double>& ub) {
  if (TYPEOF(val_sexp) != STRSXP || Rf_length(val_sexp) <= 0) {
    return false;
  }
  const int nval = Rf_length(val_sexp);
  if (nval == 1) {
    const std::string token = trim(CHAR(STRING_ELT(val_sexp, 0)));
    if (token.empty()) {
      return false;
    }
    if (token[0] == '(' || token[0] == '*') {
      return apply_named_option_line(pb, canonical_option_key(key), normalize_array_option_value(token));
    }
    double scalar = NA_REAL;
    if (!parse_scalar_option_value(token, 0, lb, ub, &scalar)) {
      return apply_named_option_line(pb, canonical_option_key(key), normalize_array_option_value(token));
    }
    return apply_array_option_values(pb, key, std::vector<double>{scalar}, n_inputs, bbin);
  }

  std::vector<double> raw_values;
  raw_values.reserve(static_cast<std::size_t>(nval));
  for (int j = 0; j < nval; ++j) {
    double parsed = NA_REAL;
    if (!parse_scalar_option_value(CHAR(STRING_ELT(val_sexp, j)), j, lb, ub, &parsed)) {
      return false;
    }
    raw_values.push_back(parsed);
  }
  return apply_array_option_values(pb, key, raw_values, n_inputs, bbin);
}

static double sanitize_epsilon(double v) {
  const double eps = std::numeric_limits<double>::epsilon();
  if (v <= eps) {
    return std::nextafter(eps, std::numeric_limits<double>::infinity());
  }
  return v;
}

static void apply_options(NomadProblem pb,
                          SEXP opts,
                          int n_inputs,
                          const std::vector<int>& bbin,
                          const std::vector<double>& lb,
                          const std::vector<double>& ub) {
  if (Rf_isNull(opts)) {
    return;
  }

  SEXP opts_integer = getListElement(opts, "integer");
  SEXP opts_numeric = getListElement(opts, "numeric");
  SEXP opts_string = getListElement(opts, "string");
  bool seen_max_eval = false;
  bool seen_max_bb_eval = false;
  int max_bb_eval_value = 0;

  if (!Rf_isNull(opts_integer)) {
    SEXP names = Rf_getAttrib(opts_integer, R_NamesSymbol);
    const int n = Rf_length(opts_integer);
    for (int i = 0; i < n; ++i) {
      const std::string key = CHAR(STRING_ELT(names, i));
      SEXP val_sexp = VECTOR_ELT(opts_integer, i);
      const int val = INTEGER(val_sexp)[0];
      if (key == "MAX_EVAL") {
        seen_max_eval = true;
      }
      if (key == "MAX_BB_EVAL") {
        seen_max_bb_eval = true;
        max_bb_eval_value = val;
      }
      if (is_ignored_legacy_option_key(key)) {
        continue;
      }
      if (is_array_option_key(key)) {
        std::vector<double> raw_values;
        const int nval = Rf_length(val_sexp);
        for (int j = 0; j < nval; ++j) {
          raw_values.push_back(static_cast<double>(INTEGER(val_sexp)[j]));
        }
        apply_array_option_values(pb, key, raw_values, n_inputs, bbin);
        continue;
      }
      if (!addNomadValParam(pb, key.c_str(), val)) {
        std::ostringstream oss;
        oss << val;
        apply_named_option_line(pb, key, oss.str());
      }
    }
  }

  if (!Rf_isNull(opts_numeric)) {
    SEXP names = Rf_getAttrib(opts_numeric, R_NamesSymbol);
    const int n = Rf_length(opts_numeric);
    for (int i = 0; i < n; ++i) {
      const std::string key = CHAR(STRING_ELT(names, i));
      SEXP val_sexp = VECTOR_ELT(opts_numeric, i);
      const double val = REAL(val_sexp)[0];
      if (key == "MAX_EVAL") {
        seen_max_eval = true;
      }
      if (key == "MAX_BB_EVAL") {
        seen_max_bb_eval = true;
        max_bb_eval_value = static_cast<int>(std::lround(val));
      }
      if (is_ignored_legacy_option_key(key)) {
        continue;
      }
      if (is_array_option_key(key)) {
        std::vector<double> raw_values;
        const int nval = Rf_length(val_sexp);
        for (int j = 0; j < nval; ++j) {
          raw_values.push_back(REAL(val_sexp)[j]);
        }
        apply_array_option_values(pb, key, raw_values, n_inputs, bbin);
        continue;
      }
      const double safe_val = (key == "EPSILON") ? sanitize_epsilon(val) : val;
      if (!addNomadDoubleParam(pb, key.c_str(), safe_val)) {
        std::ostringstream oss;
        oss.precision(17);
        oss << safe_val;
        apply_named_option_line(pb, key, oss.str());
      }
    }
  }

  if (!Rf_isNull(opts_string)) {
    SEXP names = Rf_getAttrib(opts_string, R_NamesSymbol);
    const int n = Rf_length(opts_string);
    for (int i = 0; i < n; ++i) {
      const std::string key = CHAR(STRING_ELT(names, i));
      SEXP val_sexp = VECTOR_ELT(opts_string, i);
      const char* val = CHAR(STRING_ELT(val_sexp, 0));
      if (key == "MAX_EVAL") {
        seen_max_eval = true;
      }
      if (key == "MAX_BB_EVAL") {
        char* endptr = nullptr;
        const long v = std::strtol(val, &endptr, 10);
        if (endptr != val && *endptr == '\0' && v > 0 && v <= std::numeric_limits<int>::max()) {
          seen_max_bb_eval = true;
          max_bb_eval_value = static_cast<int>(v);
        }
      }
      if (is_ignored_legacy_option_key(key)) {
        continue;
      }
      if (is_array_option_key(key)) {
        if (!apply_array_option_from_string(pb, key, val_sexp, n_inputs, bbin, lb, ub)) {
          apply_named_option_line(pb, canonical_option_key(key), normalize_array_option_value(val));
        }
        continue;
      }
      if (key == "EPSILON") {
        char* endptr = nullptr;
        const double parsed = std::strtod(val, &endptr);
        if (endptr != val && *endptr == '\0') {
          const double safe_val = sanitize_epsilon(parsed);
          if (!addNomadDoubleParam(pb, key.c_str(), safe_val)) {
            std::ostringstream oss;
            oss.precision(17);
            oss << safe_val;
            apply_named_option_line(pb, key, oss.str());
          }
          continue;
        }
      }
      if (!addNomadStringParam(pb, key.c_str(), val)) {
        apply_named_option_line(pb, key, val);
      }
    }
  }

  if (!seen_max_eval && seen_max_bb_eval && max_bb_eval_value > 0) {
    if (!addNomadValParam(pb, "MAX_EVAL", max_bb_eval_value)) {
      std::ostringstream oss;
      oss << max_bb_eval_value;
      apply_named_option_line(pb, "MAX_EVAL", oss.str());
    }
  }
}

static void apply_nomad_opt_file(NomadProblem pb) {
  std::ifstream fin("nomad.opt");
  if (!fin.good()) {
    return;
  }
  std::string line;
  while (std::getline(fin, line)) {
    const std::string clean = trim(line);
    if (clean.empty()) {
      continue;
    }
    if (clean[0] == '#') {
      continue;
    }
    addNomadParam(pb, clean.c_str());
  }
}

static double finite_or_default(double x, double fallback) {
  return std::isfinite(x) ? x : fallback;
}

static double clamp_to_bounds(double x, double lb, double ub) {
  if (std::isfinite(lb) && x < lb) {
    x = lb;
  }
  if (std::isfinite(ub) && x > ub) {
    x = ub;
  }
  return x;
}

static double coerce_by_input_type(double x, int type, double lb, double ub) {
  if (type == 3) {
    x = (x >= 0.5) ? 1.0 : 0.0;
  } else if (type == 1 || type == 2) {
    x = std::round(x);
  }
  return clamp_to_bounds(x, lb, ub);
}

static void fill_random_start(std::vector<double>& x0s,
                              int point_index,
                              int n,
                              const std::vector<int>& bbin,
                              const std::vector<double>& lb,
                              const std::vector<double>& ub,
                              std::mt19937_64& rng) {
  for (int i = 0; i < n; ++i) {
    double lo = finite_or_default(lb[i], -1.0);
    double hi = finite_or_default(ub[i], 1.0);
    if (hi < lo) {
      std::swap(lo, hi);
    }
    std::uniform_real_distribution<double> unif(lo, hi);
    double v = unif(rng);
    x0s[point_index * n + i] = coerce_by_input_type(v, bbin[i], lb[i], ub[i]);
  }
}

static void build_starting_points(SEXP sx0,
                                  int n,
                                  int nstart,
                                  const std::vector<int>& bbin,
                                  const std::vector<double>& lb,
                                  const std::vector<double>& ub,
                                  unsigned int seed,
                                  std::vector<double>& x0s) {
  std::mt19937_64 rng(seed == 0 ? static_cast<unsigned int>(std::time(nullptr)) : seed);
  x0s.assign(static_cast<std::size_t>(nstart) * static_cast<std::size_t>(n), 0.0);

  for (int j = 0; j < nstart; ++j) {
    fill_random_start(x0s, j, n, bbin, lb, ub, rng);
  }

  if (Rf_isNull(sx0)) {
    return;
  }

  const int len = Rf_length(sx0);
  if (len == n) {
    for (int i = 0; i < n; ++i) {
      x0s[i] = coerce_by_input_type(REAL(sx0)[i], bbin[i], lb[i], ub[i]);
    }
    return;
  }

  if (len >= n * nstart) {
    for (int j = 0; j < nstart; ++j) {
      for (int i = 0; i < n; ++i) {
        // Keep legacy R column-major convention used by old implementation.
        const double val = REAL(sx0)[j + i * nstart];
        x0s[j * n + i] = coerce_by_input_type(val, bbin[i], lb[i], ub[i]);
      }
    }
  }
}

static int map_run_flag_to_status(int run_flag) {
  switch (run_flag) {
    case 1:
      return 8;   // mesh/target reached on feasible point
    case 0:
      return 13;  // budget/iter reached with feasible point
    case -1:
      return 8;   // mesh convergence (infeasible)
    case -2:
      return 13;  // budget/iter reached (infeasible)
    case -3:
      return 6;   // x0 fail
    case -4:
      return 12;  // max time
    case -5:
      return 3;   // user stop / CTRL-C
    case -6:
      return 18;  // stop on feasible
    case -7:
    case -8:
    default:
      return 1;   // generic error
  }
}

static const char* run_flag_message(int run_flag) {
  switch (run_flag) {
    case 1:
      return "NOMAD4 run flag 1: target reached or mesh converged on feasible point";
    case 0:
      return "NOMAD4 run flag 0: feasible point found; evaluation/iteration budget reached";
    case -1:
      return "NOMAD4 run flag -1: mesh converged without feasible point";
    case -2:
      return "NOMAD4 run flag -2: no feasible point; evaluation/iteration budget reached";
    case -3:
      return "NOMAD4 run flag -3: starting point failed to evaluate";
    case -4:
      return "NOMAD4 run flag -4: time limit reached";
    case -5:
      return "NOMAD4 run flag -5: user interruption";
    case -6:
      return "NOMAD4 run flag -6: stop on feasible point";
    case -7:
      return "NOMAD4 run flag -7: invalid parameters";
    case -8:
      return "NOMAD4 run flag -8: evaluation/runtime failure";
    default:
      return "NOMAD4 run flag: unknown";
  }
}

static SEXP build_solution(int status,
                           const char* message,
                           int bbe,
                           int iterations,
                           double objective,
                           const std::vector<double>& solution) {
  const int n = static_cast<int>(solution.size());
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 6));
  SEXP names = PROTECT(Rf_allocVector(STRSXP, 6));

  SET_STRING_ELT(names, 0, Rf_mkChar("status"));
  SET_STRING_ELT(names, 1, Rf_mkChar("message"));
  SET_STRING_ELT(names, 2, Rf_mkChar("bbe"));
  SET_STRING_ELT(names, 3, Rf_mkChar("objective"));
  SET_STRING_ELT(names, 4, Rf_mkChar("solution"));
  SET_STRING_ELT(names, 5, Rf_mkChar("iterations"));
  Rf_setAttrib(out, R_NamesSymbol, names);

  SEXP r_status = PROTECT(Rf_allocVector(INTSXP, 1));
  INTEGER(r_status)[0] = status;
  SET_VECTOR_ELT(out, 0, r_status);

  SEXP r_message = PROTECT(Rf_allocVector(STRSXP, 1));
  SET_STRING_ELT(r_message, 0, Rf_mkChar(message));
  SET_VECTOR_ELT(out, 1, r_message);

  SEXP r_bbe = PROTECT(Rf_allocVector(INTSXP, 1));
  INTEGER(r_bbe)[0] = bbe;
  SET_VECTOR_ELT(out, 2, r_bbe);

  SEXP r_obj = PROTECT(Rf_allocVector(REALSXP, 1));
  REAL(r_obj)[0] = objective;
  SET_VECTOR_ELT(out, 3, r_obj);

  SEXP r_sol = PROTECT(Rf_allocVector(REALSXP, n));
  for (int i = 0; i < n; ++i) {
    REAL(r_sol)[i] = solution[static_cast<std::size_t>(i)];
  }
  SET_VECTOR_ELT(out, 4, r_sol);

  SEXP r_iter = PROTECT(Rf_allocVector(INTSXP, 1));
  INTEGER(r_iter)[0] = iterations;
  SET_VECTOR_ELT(out, 5, r_iter);

  UNPROTECT(8);
  return out;
}

static SEXP solve_nomad4(SEXP args, bool use_multi) {
  theenv = getListElement(args, "snomadr.environment");
  thefun = getListElement(args, "eval.f");

  if (Rf_isNull(thefun) || !Rf_isFunction(thefun)) {
    Rf_error("Invalid eval.f supplied to snomadr.");
  }
  if (Rf_isNull(theenv)) {
    theenv = R_GlobalEnv;
  }

  SEXP sN = getListElement(args, "n");
  const int N = INTEGER(sN)[0];
  if (N <= 0) {
    Rf_error("Invalid dimension n supplied to snomadr.");
  }

  SEXP sbbin = getListElement(args, "bbin");
  SEXP sbbout = getListElement(args, "bbout");
  SEXP slb = getListElement(args, "lower.bounds");
  SEXP sub = getListElement(args, "upper.bounds");
  SEXP sx0 = getListElement(args, "x0");
  SEXP sopts = getListElement(args, "options");
  SEXP snmulti = getListElement(args, "nmulti");
  SEXP sseed = getListElement(args, "random.seed");
  SEXP sprint = getListElement(args, "print.output");

  int nstart = 1;
  if (use_multi && !Rf_isNull(snmulti)) {
    nstart = std::max(1, INTEGER(snmulti)[0]);
  }

  unsigned int seed = 0;
  if (!Rf_isNull(sseed) && Rf_length(sseed) > 0) {
    seed = static_cast<unsigned int>(INTEGER(sseed)[0]);
  }
  if (seed == 0) {
    seed = static_cast<unsigned int>(std::time(nullptr));
  }

  int print_output = 1;
  if (!Rf_isNull(sprint) && Rf_length(sprint) > 0) {
    print_output = (Rf_asLogical(sprint) == TRUE) ? 1 : 0;
  }

  std::vector<int> bbin(static_cast<std::size_t>(N), 0);
  std::vector<double> lb(static_cast<std::size_t>(N), -std::numeric_limits<double>::infinity());
  std::vector<double> ub(static_cast<std::size_t>(N), std::numeric_limits<double>::infinity());
  for (int i = 0; i < N; ++i) {
    if (!Rf_isNull(sbbin) && Rf_length(sbbin) > i) {
      bbin[static_cast<std::size_t>(i)] = INTEGER(sbbin)[i];
    }
    if (!Rf_isNull(slb) && Rf_length(slb) > i) {
      lb[static_cast<std::size_t>(i)] = REAL(slb)[i];
    }
    if (!Rf_isNull(sub) && Rf_length(sub) > i) {
      ub[static_cast<std::size_t>(i)] = REAL(sub)[i];
    }
  }

  const int M = Rf_length(sbbout);
  if (M <= 0) {
    Rf_error("Invalid bbout supplied to snomadr.");
  }

  NomadProblem pb = createNomadProblem(r_eval_single, nullptr, N, M);
  if (pb == nullptr) {
    Rf_error("Unable to initialize NOMAD4 problem.");
  }
  if (!addNomadValParam(pb, "DIMENSION", N)) {
    freeNomadProblem(pb);
    Rf_error("Failed to set DIMENSION for NOMAD4.");
  }

  if (!addNomadStringParam(pb, "BB_INPUT_TYPE", map_bbin_to_string(sbbin).c_str())) {
    freeNomadProblem(pb);
    Rf_error("Failed to set BB_INPUT_TYPE for NOMAD4.");
  }
  if (!addNomadStringParam(pb, "BB_OUTPUT_TYPE", map_bbout_to_string(sbbout).c_str())) {
    freeNomadProblem(pb);
    Rf_error("Failed to set BB_OUTPUT_TYPE for NOMAD4.");
  }
  if (!addNomadArrayOfDoubleParam(pb, "LOWER_BOUND", lb.data())) {
    freeNomadProblem(pb);
    Rf_error("Failed to set LOWER_BOUND for NOMAD4.");
  }
  if (!addNomadArrayOfDoubleParam(pb, "UPPER_BOUND", ub.data())) {
    freeNomadProblem(pb);
    Rf_error("Failed to set UPPER_BOUND for NOMAD4.");
  }

  // Keep old behavior: hide progress when print.output is FALSE.
  if (!print_output) {
    addNomadValParam(pb, "DISPLAY_DEGREE", 0);
    addNomadBoolParam(pb, "DISPLAY_ALL_EVAL", false);
    addNomadBoolParam(pb, "DISPLAY_UNSUCCESSFUL", false);
  }

  apply_options(pb, sopts, N, bbin, lb, ub);
  apply_nomad_opt_file(pb);

  std::vector<double> x0s;
  build_starting_points(sx0, N, nstart, bbin, lb, ub, seed, x0s);

  NomadResult result = createNomadResult();
  if (result == nullptr) {
    freeNomadProblem(pb);
    Rf_error("Unable to allocate NOMAD4 result container.");
  }

  g_bbe_counter.store(0, std::memory_order_relaxed);
  const int run_flag = solveNomadProblem(result, pb, nstart, x0s.data(), nullptr);

  int nb_solutions = nbSolutionsNomadResult(result);
  if (nb_solutions < 0) {
    nb_solutions = 0;
  }

  std::vector<double> best_x(static_cast<std::size_t>(N), NA_REAL);
  std::vector<double> best_out(static_cast<std::size_t>(M), NA_REAL);
  if (nb_solutions > 0) {
    loadInputSolutionsNomadResult(best_x.data(), 1, result);
    loadOutputSolutionsNomadResult(best_out.data(), 1, result);
  }

  int objective_index = 0;
  for (int i = 0; i < M; ++i) {
    if (INTEGER(sbbout)[i] == 0) {
      objective_index = i;
      break;
    }
  }
  const double objective = (nb_solutions > 0) ? best_out[static_cast<std::size_t>(objective_index)] : NA_REAL;

  int bbe = g_bbe_counter.load(std::memory_order_relaxed);
  try {
    auto evc = NOMAD::EvcInterface::getEvaluatorControl();
    if (evc != nullptr) {
      const std::size_t nbe = evc->getNbEval();
      const int bbe_evc = (nbe > static_cast<std::size_t>(std::numeric_limits<int>::max()))
                              ? std::numeric_limits<int>::max()
                              : static_cast<int>(nbe);
      if (bbe_evc > 0) {
        bbe = bbe_evc;
      }
    }
  } catch (...) {
    // Keep a safe scalar to preserve legacy downstream comparisons.
    bbe = 0;
  }

  const int iterations = 0;
  const int status = map_run_flag_to_status(run_flag);
  const char* message = run_flag_message(run_flag);

  SEXP ret = build_solution(status, message, bbe, iterations, objective, best_x);

  freeNomadResult(result);
  freeNomadProblem(pb);

  return ret;
}

extern "C" {

SEXP snomadRInfo(SEXP args) {
  SEXP solution = PROTECT(args);

  SEXP sinfo = getListElement(args, "info");
  SEXP sversion = getListElement(args, "version");
  SEXP shelp = getListElement(args, "help");

  if (!Rf_isNull(sinfo)) {
    const char* s = CHAR(STRING_ELT(sinfo, 0));
    if (s[0] == '-' && (s[1] == 'i' || s[1] == 'I')) {
      rout << "NOMAD (embedded in crs) -- C interface backend\n";
      rout << "See https://www.gerad.ca/en/software/nomad\n";
    }
  }

  if (!Rf_isNull(sversion)) {
    const char* s = CHAR(STRING_ELT(sversion, 0));
    if (s[0] == '-' && (s[1] == 'v' || s[1] == 'V')) {
      rout << "NOMAD version " << NOMAD_VERSION_NUMBER << '\n';
    }
  }

  if (!Rf_isNull(shelp)) {
    const char* s = CHAR(STRING_ELT(shelp, 0));
    if (s[0] == '-' && (s[1] == 'h' || s[1] == 'H')) {
      rout << "NOMAD parameter help is available in the user guide:\n";
      rout << "https://nomad-4-user-guide.readthedocs.io/en/latest/\n";
    }
  }

  UNPROTECT(1);
  return solution;
}

SEXP snomadRSolve(SEXP args) {
  return solve_nomad4(args, false);
}

SEXP smultinomadRSolve(SEXP args) {
  return solve_nomad4(args, true);
}

} // extern "C"
