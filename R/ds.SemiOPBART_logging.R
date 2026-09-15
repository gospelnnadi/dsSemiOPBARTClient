# dsSemiOPBART_logging.R
# ---------------------------------------------------------------------------
# OUTPUT-STREAM LOGGING
#
# Several server-side functions in this package have (or have had) stray
# print()/cat() debug calls left in from development -- a couple were
# cleaned up directly (localFitDS.R, semiOPBARTLocalSwapDS() in
# semiopbartDS.R), but DSLite/DSI's own aggregate-call progress output
# also prints, and more debug prints could get added back by anyone
# working on this code later. Rather than relying on every one of those
# being caught by review, semiOPBART_withLogging() wraps a whole
# orchestrator run (or any expression) and tees EVERYTHING printed during
# it -- print(), cat(), message(), warning() -- to a log file, in addition
# to (or instead of) the console.
# ---------------------------------------------------------------------------

#' Run `expr` with all console output teed to `log_file`.
#'
#' HOW: print()/cat() go through R's stdout connection, which sink(...,
#' split = TRUE) natively tees to both the file and the console in one
#' call -- no extra work needed for those. message()/warning() do NOT go
#' through stdout at all (they go to a separate message connection/
#' stderr), and sink(type = "message") cannot be split, so those are
#' instead intercepted directly via withCallingHandlers(), re-emitted
#' through cat() (which -- because the stdout sink above is already
#' active for the whole duration -- automatically gets teed the same way),
#' and then muffled so R's default handler doesn't ALSO print them a
#' second time.
#'
#' Errors are NOT intercepted (only message/warning conditions are) --
#' they propagate normally. on.exit() guarantees the sink is removed and
#' the log file closed on both normal return AND an error inside `expr`,
#' so a failed run never leaves the console silently redirected to a
#' file. One known gap: a FATAL top-level error's own final printed
#' message goes through R's own error-printing path after the sink is
#' already torn down by on.exit(), so it will show on the console but
#' NOT land in the log file -- wrap the call site in tryCatch() yourself
#' if you need that captured too.
#'
#' @param expr  an expression (use `{ ... }` for multiple statements),
#'   evaluated via substitute()/eval() so it runs in the caller's
#'   environment, not this function's
#' @param log_file  path to write to
#' @param append  append to an existing file instead of overwriting (default FALSE)
#' @param also_console  tee to the console too (default TRUE); FALSE sends
#'   output to the file ONLY, nothing printed live
#' @return whatever `expr` returns, invisibly is NOT forced -- same
#'   visibility as a normal top-level call
#' @export
#'
#' @examples
#' \dontrun{
#' result <- semiOPBART_withLogging({
#'   run_semiOPBART_IP_EP_orchestrator(conns = conns, ...)
#' }, log_file = "output/logs/run_2026-08-14.log")
#' }
semiOPBART_withLogging <- function(expr, log_file, append = FALSE, also_console = TRUE) {
  expr <- substitute(expr)
  parent <- parent.frame()

  log_dir <- dirname(log_file)
  if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)

  con <- file(log_file, open = if (append) "a" else "w")
  on.exit(close(con), add = TRUE)

  cat(paste0("===== semiOPBART log started ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " =====\n"),
      file = con)

  sink(con, split = also_console)
  # after = FALSE inserts THIS on.exit() ahead of the close(con) one above,
  # so it runs first on exit -- sink must be removed before the connection
  # it points at is closed, in both the normal-return and error-unwind case.
  on.exit(sink(), add = TRUE, after = FALSE)

  withCallingHandlers(
    eval(expr, envir = parent),
    message = function(m) {
      cat(sub("\n$", "", conditionMessage(m)), "\n", sep = "")
      invokeRestart("muffleMessage")
    },
    warning = function(w) {
      cat("Warning: ", conditionMessage(w), "\n", sep = "")
      invokeRestart("muffleWarning")
    }
  )
}
