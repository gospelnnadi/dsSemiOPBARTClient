# dsSemiOPBARTLocalCombineF.R
# ---------------------------------------------------------------------------
# Architecture F's combine -- pure client-side arithmetic, no datasources
# argument, no aggregate() call, no disclosure surface at all: it only
# ever touches what ds.semiOPBARTFitF() already returned. Pools ONLY
# theta (inverse-variance) and us (n- or equal-weighted) ONCE, after
# every site's independent smopbart() run finished -- see ds_localFit.R's
# header for how this differs from Architecture E.
# ---------------------------------------------------------------------------

#' @param site_fits      output of ds.semiOPBARTFitF()
#' @param combine_method "inverse_variance" (precision-weighted average --
#'   treats each site's local posterior as an independent noisy estimate of
#'   one shared truth) or "stack" (n- or equal-weighted mixture).
#' @param weight_by_n   used only for combine_method = "stack".
#' @export
ds.semiOPBARTCombineF <- function(site_fits,
                                  combine_method = c("inverse_variance", "stack"),
                                  weight_by_n = TRUE) {
  combine_method <- match.arg(combine_method)

  n <- vapply(site_fits, `[[`, numeric(1), "n")
  w <- if (weight_by_n) n / sum(n) else rep(1 / length(site_fits), length(site_fits))

  theta_combined <- if (combine_method == "inverse_variance") {
    precisions <- lapply(site_fits, function(s) solve(s$theta_cov))
    V <- solve(Reduce(`+`, precisions))
    as.numeric(V %*% Reduce(`+`, Map(function(s, P) P %*% s$theta_mean,
                                     site_fits, precisions)))
  } else {
    Reduce(`+`, Map(function(s, wi) wi * s$theta_mean, site_fits, w))
  }
  us_combined <- Reduce(`+`, Map(function(s, wi) wi * s$us_mean, site_fits, w))

  list(theta = theta_combined, us = us_combined, weights_used = w,
       per_site_fits = site_fits, combine_method = combine_method,
       linear_formula = attr(site_fits, "linear_formula"),
       state_name = attr(site_fits, "state_name"))
}
