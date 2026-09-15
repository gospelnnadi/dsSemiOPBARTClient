# dsSemiOPBARTExplain.R
# ---------------------------------------------------------------------------
# MODEL EXPLAINABILITY MODULE -- Architectures D, E, F
#
# What this is: two kinds of explainability, both built ONLY from
# summaries the pipeline already produces or can produce disclosure-safely
# -- nothing here ever touches row-level data.
#
#   1. LINEAR COEFFICIENTS (theta) and THRESHOLDS (us) -- the semi-
#      parametric part of the model. Already fully available client-side
#      the moment run_semiOPBART_orchestrator() returns:
#        - Architecture D:  fit_D$theta / fit_D$us are full posterior draw
#          matrices (every kept sweep) -- real credible intervals.
#        - Architecture E:  fit_E$theta_draws/us_draws if present, else
#          fit_E$theta_mean/us_mean (point only, no interval).
#        - Architecture F:  fit_F$theta/us are POINT estimates from
#          inverse-variance meta-analysis; fit_F$theta_cov (now returned
#          by ds.semiOPBARTCombineF() -- see ds_localCombine.R) gives a
#          normal-approximation interval. No equivalent covariance exists
#          for us in F, so its threshold plot has no interval.
#      No new server-side calls needed for any of this.
#
#   2. VARIABLE IMPORTANCE for the BART part (f(x)) -- how often each X
#      feature is used as a splitting variable, via the forest's
#      $get_counts() method. This already existed and was already being
#      computed, just never returned to the client:
#        - Architecture F:  smopbart() (the local, non-federated fit
#          Architecture F calls per site) already computes var_counts
#          every post-burn-in sweep (semiOPBART.R). It just wasn't in
#          semiOPBARTLocalFitFDS()'s return list -- now it is
#          (var_counts_mean, localFitDS.R), and ds.semiOPBARTCombineF()
#          pools it the same n-weighted way as theta/us.
#        - Architectures D/E: no local smopbart() call exists to expose
#          this from, so a new (small, disclosure-safe) aggregate function
#          was added -- semiOPBARTLocalVarCountsDS() (semiopbartDS.R),
#          called via the new ds.semiOPBARTVarImportance() (ds_semiopbart.R).
#          It's a SNAPSHOT of the final trained forest's split counts, not
#          a full per-sweep posterior like F gets, specifically to avoid
#          touching the already-tested semiOPBARTLocalMCMC(Improved)() hot
#          loop for this -- see that function's docstring for the trade-off.
#
# What this is NOT: no partial dependence / model-surface plots. The
# orchestrator's own docstring records that a grid-combined f(x) capability
# (Architecture D's old "Option C") was deliberately REMOVED because it
# required per-site model-surface data to leave a site. Reintroducing
# anything that evaluates f(x) at grid points and returns it centrally
# would be the same category of capability -- that's a data-governance
# decision for whoever owns this pipeline, not something to quietly add
# back in via an "explainability" label. Variable importance (a count, not
# a curve) does not carry that risk.
# ---------------------------------------------------------------------------

# ============================================================================
# INTERNAL HELPERS
# ============================================================================

.semiOPBART_summarizeFromDraws <- function(draws) {
  data.frame(
    mean  = colMeans(draws),
    sd    = apply(draws, 2, stats::sd),
    ci_lo = apply(draws, 2, stats::quantile, probs = 0.025, names = FALSE),
    ci_hi = apply(draws, 2, stats::quantile, probs = 0.975, names = FALSE)
  )
}

.semiOPBART_summarizeFromPoint <- function(point, cov_mat = NULL) {
  point <- as.numeric(point)
  sd <- if (!is.null(cov_mat)) sqrt(diag(as.matrix(cov_mat))) else rep(NA_real_, length(point))
  data.frame(
    mean  = point,
    sd    = sd,
    ci_lo = point - 1.96 * sd,   # normal approximation -- only meaningful when sd is non-NA
    ci_hi = point + 1.96 * sd
  )
}

# Label theta's p coefficients from w_features when the lengths line up
# (p == length(w_features): no intercept in the design; p == length(w_features)
# + 1: intercept first, R's usual model.matrix() column order) -- generic
# W1..Wp names otherwise, rather than guessing wrong or erroring.
.semiOPBART_thetaTermNames <- function(p, w_features) {
  if (!is.null(w_features) && length(w_features) == p) return(w_features)
  if (!is.null(w_features) && length(w_features) == p - 1) return(c("(Intercept)", w_features))
  paste0("W", seq_len(p))
}

.semiOPBART_thresholdTermNames <- function(p) {
  nm <- paste0("u", seq_len(p))
  if (p >= 1) nm[1] <- "u1 (fixed = 0)"   # us[1] is always fixed at 0 by construction
  nm
}


# ============================================================================
# 1. LINEAR COEFFICIENTS (theta) AND THRESHOLDS (us)
# ============================================================================

#' Posterior summary of the linear (W) coefficients for one architecture's
#' fit -- see this file's header for exactly what's available per
#' architecture.
#'
#' @param fit  fit_D, fit_E, or fit_F from a run_semiOPBART_orchestrator() result
#' @param w_features  the SAME w_features passed to that orchestrator call,
#'   for labeling -- falls back to generic names if lengths don't match
#'   (see .semiOPBART_thetaTermNames())
#' @param architecture  "D", "E", or "F"
#' @return data.frame(architecture, term, mean, sd, ci_lo, ci_hi), or NULL
#'   if fit is NULL
#' @export
semiOPBART_explainCoefficients <- function(fit, w_features = NULL,
                                            architecture = c("D", "E", "F")) {
  architecture <- match.arg(architecture)
  if (is.null(fit)) return(NULL)

  summary_df <- switch(architecture,
    D = {
      if (is.null(fit$theta) || !is.matrix(fit$theta))
        stop("semiOPBART_explainCoefficients(): Architecture D fit has no theta draw matrix.")
      .semiOPBART_summarizeFromDraws(fit$theta)
    },
    E = {
      if (!is.null(fit$theta_draws)) {
        .semiOPBART_summarizeFromDraws(fit$theta_draws)
      } else if (!is.null(fit$theta_mean)) {
        message("semiOPBART_explainCoefficients(): Architecture E fit has no ",
                "theta_draws -- plotting theta_mean with no uncertainty band.")
        .semiOPBART_summarizeFromPoint(fit$theta_mean)
      } else {
        stop("semiOPBART_explainCoefficients(): Architecture E fit has neither ",
             "theta_draws nor theta_mean.")
      }
    },
    F = {
      if (is.null(fit$theta))
        stop("semiOPBART_explainCoefficients(): Architecture F fit has no theta.")
      .semiOPBART_summarizeFromPoint(fit$theta, fit$theta_cov)
    }
  )

  summary_df$term <- .semiOPBART_thetaTermNames(nrow(summary_df), w_features)
  summary_df$architecture <- architecture
  summary_df[, c("architecture", "term", "mean", "sd", "ci_lo", "ci_hi")]
}

#' Same idea as semiOPBART_explainCoefficients(), for the ordinal
#' thresholds (us) instead of the linear coefficients.
#' @export
semiOPBART_explainThresholds <- function(fit, architecture = c("D", "E", "F")) {
  architecture <- match.arg(architecture)
  if (is.null(fit)) return(NULL)

  summary_df <- switch(architecture,
    D = {
      if (is.null(fit$us) || !is.matrix(fit$us))
        stop("semiOPBART_explainThresholds(): Architecture D fit has no us draw matrix.")
      .semiOPBART_summarizeFromDraws(fit$us)
    },
    E = {
      if (!is.null(fit$us_draws)) {
        .semiOPBART_summarizeFromDraws(fit$us_draws)
      } else if (!is.null(fit$us_mean)) {
        .semiOPBART_summarizeFromPoint(fit$us_mean)
      } else {
        stop("semiOPBART_explainThresholds(): Architecture E fit has neither ",
             "us_draws nor us_mean.")
      }
    },
    F = {
      if (is.null(fit$us))
        stop("semiOPBART_explainThresholds(): Architecture F fit has no us.")
      .semiOPBART_summarizeFromPoint(fit$us)   # no us covariance anywhere -> sd/CI are NA
    }
  )

  summary_df$term <- .semiOPBART_thresholdTermNames(nrow(summary_df))
  summary_df$architecture <- architecture
  summary_df[, c("architecture", "term", "mean", "sd", "ci_lo", "ci_hi")]
}


# ============================================================================
# 2. VARIABLE IMPORTANCE (BART split counts)
# ============================================================================

#' Architecture F's variable importance -- var_counts_mean was already
#' pooled n-weighted across sites by ds.semiOPBARTCombineF(), this just
#' reshapes it to match D/E's tidy format.
#' @export
semiOPBART_explainVarImportanceF <- function(fit_F) {
  if (is.null(fit_F) || is.null(fit_F$var_counts)) return(NULL)
  vc <- fit_F$var_counts
  data.frame(architecture = "F",
             term = if (!is.null(names(vc))) names(vc) else paste0("X", seq_along(vc)),
             importance = as.numeric(vc), stringsAsFactors = FALSE)
}

#' Architectures D/E's variable importance -- wraps the pooled numeric
#' vector ds.semiOPBARTVarImportance() (ds_semiopbart.R) already returns.
#' @param pooled_counts  output of ds.semiOPBARTVarImportance()
#' @export
semiOPBART_explainVarImportanceDE <- function(pooled_counts, architecture = c("D", "E")) {
  architecture <- match.arg(architecture)
  if (is.null(pooled_counts)) return(NULL)
  data.frame(architecture = architecture,
             term = if (!is.null(names(pooled_counts))) names(pooled_counts)
                    else paste0("X", seq_along(pooled_counts)),
             importance = as.numeric(pooled_counts), stringsAsFactors = FALSE)
}


# ============================================================================
# 3. PLOTS
# ============================================================================

#' Forest/coefficient plot for a semiOPBART_explainCoefficients()/
#' semiOPBART_explainThresholds() data.frame. Terms with sd = NA (F's
#' thresholds; E without theta_draws) plot as a point with no error bar
#' rather than being dropped.
#' @export
semiOPBART_plotCoefficients <- function(df, title = "Coefficient estimates") {
  if (is.null(df) || nrow(df) == 0) return(invisible(NULL))
  df$term <- factor(df$term, levels = df$term[order(df$mean)])
  ggplot2::ggplot(df, ggplot2::aes(x = term, y = mean)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::geom_pointrange(ggplot2::aes(ymin = ci_lo, ymax = ci_hi), na.rm = TRUE) +
    ggplot2::geom_point(data = df[is.na(df$sd), ], size = 2, na.rm = TRUE) +
    ggplot2::coord_flip() +
    ggplot2::labs(title = title, x = NULL, y = "Estimate (95% interval where available)") +
    ggplot2::theme_minimal(base_size = 12)
}

#' Bar chart for a semiOPBART_explainVarImportance{F,DE}() data.frame.
#' @export
semiOPBART_plotVarImportance <- function(df, title = "Variable importance (mean split count)") {
  if (is.null(df) || nrow(df) == 0) return(invisible(NULL))
  df$term <- factor(df$term, levels = df$term[order(df$importance)])
  ggplot2::ggplot(df, ggplot2::aes(x = term, y = importance)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::labs(title = title, x = NULL, y = "Mean split count") +
    ggplot2::theme_minimal(base_size = 12)
}


# ============================================================================
# 4. TOP-LEVEL ENTRYPOINT -- run for every architecture present in one
#    run_semiOPBART_orchestrator() result
# ============================================================================

#' Build (and optionally save) coefficient, threshold, and variable-
#' importance summaries/plots for every architecture (D/E/F) present in
#' one run_semiOPBART_orchestrator() result. No orchestrator changes
#' needed -- this consumes its return value.
#'
#' @param run_result  one run_semiOPBART_orchestrator() result (has fit_D/
#'   fit_E/fit_F for whichever architectures were run; others are NULL and
#'   silently skipped)
#' @param w_features  the SAME w_features passed to that orchestrator call
#' @param datasources  needed only if fit_D/fit_E is present (variable
#'   importance needs one fresh aggregate() round trip); NULL default
#'   finds every current connection
#' @param out_dir  if given, saves one CSV per summary and one PNG per
#'   plot here (created if missing); NULL (default) returns everything
#'   without writing to disk
#' @param weight_by_n  passed to ds.semiOPBARTVarImportance() for D/E
#' @export
semiOPBART_explainRun <- function(run_result, w_features = NULL, datasources = NULL,
                                   out_dir = NULL, weight_by_n = TRUE) {
  if (!is.null(out_dir) && !dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  coef_frames <- list(); thresh_frames <- list(); importance_frames <- list()

  if (!is.null(run_result$fit_D)) {
    coef_frames$D   <- semiOPBART_explainCoefficients(run_result$fit_D, w_features, "D")
    thresh_frames$D <- semiOPBART_explainThresholds(run_result$fit_D, "D")
    pooled_D <- tryCatch(
      ds.semiOPBARTVarImportance(run_result$fit_D, datasources, weight_by_n),
      error = function(e) {
        message("semiOPBART_explainRun(): Architecture D variable importance skipped: ",
                conditionMessage(e))
        NULL
      })
    importance_frames$D <- semiOPBART_explainVarImportanceDE(pooled_D, "D")
  }

  if (!is.null(run_result$fit_E)) {
    coef_frames$E   <- semiOPBART_explainCoefficients(run_result$fit_E, w_features, "E")
    thresh_frames$E <- semiOPBART_explainThresholds(run_result$fit_E, "E")
    pooled_E <- tryCatch(
      ds.semiOPBARTVarImportance(run_result$fit_E, datasources, weight_by_n),
      error = function(e) {
        message("semiOPBART_explainRun(): Architecture E variable importance skipped: ",
                conditionMessage(e))
        NULL
      })
    importance_frames$E <- semiOPBART_explainVarImportanceDE(pooled_E, "E")
  }

  if (!is.null(run_result$fit_F)) {
    coef_frames$F       <- semiOPBART_explainCoefficients(run_result$fit_F, w_features, "F")
    thresh_frames$F     <- semiOPBART_explainThresholds(run_result$fit_F, "F")
    importance_frames$F <- semiOPBART_explainVarImportanceF(run_result$fit_F)
  }

  coef_df       <- dplyr::bind_rows(coef_frames)
  thresh_df     <- dplyr::bind_rows(thresh_frames)
  importance_df <- dplyr::bind_rows(importance_frames)

  plots <- list()
  for (arch in names(coef_frames)) {
    if (!is.null(coef_frames[[arch]]))
      plots[[paste0("coef_", arch)]] <- semiOPBART_plotCoefficients(
        coef_frames[[arch]], paste0("Architecture ", arch, ": Linear Coefficients"))
    if (!is.null(thresh_frames[[arch]]))
      plots[[paste0("thresh_", arch)]] <- semiOPBART_plotCoefficients(
        thresh_frames[[arch]], paste0("Architecture ", arch, ": Thresholds"))
    if (!is.null(importance_frames[[arch]]))
      plots[[paste0("importance_", arch)]] <- semiOPBART_plotVarImportance(
        importance_frames[[arch]], paste0("Architecture ", arch, ": Variable Importance"))
  }

  if (!is.null(out_dir)) {
    if (!is.null(coef_df) && nrow(coef_df))
      utils::write.csv(coef_df, file.path(out_dir, "explain_coefficients.csv"), row.names = FALSE)
    if (!is.null(thresh_df) && nrow(thresh_df))
      utils::write.csv(thresh_df, file.path(out_dir, "explain_thresholds.csv"), row.names = FALSE)
    if (!is.null(importance_df) && nrow(importance_df))
      utils::write.csv(importance_df, file.path(out_dir, "explain_var_importance.csv"), row.names = FALSE)
    for (nm in names(plots)) {
      if (!is.null(plots[[nm]]))
        ggplot2::ggsave(file.path(out_dir, paste0(nm, ".png")), plots[[nm]],
                         width = 8, height = 6, dpi = 300)
    }
  }

  list(coefficients = coef_df, thresholds = thresh_df,
       importance = importance_df, plots = plots)
}
