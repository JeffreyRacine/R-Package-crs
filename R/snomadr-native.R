.crs_nomad_native_eval_r <- function(eval.f, x) {
  tryCatch(
    list(0L, eval.f(x), ""),
    interrupt = function(e) {
      list(1L, NULL, conditionMessage(e))
    },
    error = function(e) {
      list(2L, NULL, conditionMessage(e))
    }
  )
}
