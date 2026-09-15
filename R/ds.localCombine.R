
# dsSemiOPBARTLocalCombineF.R
# ---------------------------------------------------------------------------
# Architecture F's combine -- pure client-side arithmetic, no datasources
# argument, no aggregate() call, no disclosure surface at all: it only
# ever touches what ds.semiOPBARTFitF() already returned. Pools ONLY
# theta (inverse-variance) and us (n- or equal-weighted) ONCE, after
# every site's independent smopbart() run finished -- see ds_localFit.R's
# header for how this differs from Architecture E.
#
# ON A 1-SITE RUN UNDERPERFORMING THE NON-FEDERATED smopbart() RUN: if
# you're comparing federated performance against local smopbart() and
# seeing a drop even with only one site involved, the combine step below
# is NOT the cause -- see the single-site short-circuit inside
# ds.semiOPBARTCombineF(), which makes a 1-site run return that site's
# ds.semiOPBARTFitF() output completely unchanged (no solve(), no
# weighting, nothing that could introduce a difference). If a gap still
# shows up at n_sites == 1, look instead at: (a) for Architecture D/E,
# .semiOPBART_poolTheta()'s regularizing prior (prior_var = 100 by
# default), which applies on EVERY sync regardless of site count and is
# a real, deliberate difference from smopbart()'s flat/improper prior --
# not present in Architecture F at all, since F's single-site path now
# skips theta pooling entirely; (b) whether the production block-based
# MCMC (semiOPBARTLocalMCMC/Improved) is actually running the SAME sweep
# order as smopbart() -- see semiOPBARTLocalMCMC's own docstring for a
# confirmed, still-unpatched deviation there; (c) the data-prep/
# normalization pathway (ds.semiOPBARTPrepare) against whatever local
# pipeline smopbart() was run through directly; (d) train/test split
# RNG/stratification differences. Isolating and closing any 1-site gap
# FIRST is what makes a later multi-site comparison actually mean
# something -- right now it's confounded by whichever of these is live.
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

  # SINGLE-SITE SHORT-CIRCUIT: with only one contributor, every weighting
  # scheme below is mathematically a no-op anyway (weight = 1;
  # inverse-variance pooling of a single precision matrix returns that
  # same matrix, solve(solve(theta_cov)) == theta_cov up to floating-
  # point noise) -- but routing through solve()/Reduce() regardless adds
  # real risk for zero benefit: a single site's theta_cov can be poorly
  # conditioned in a way that's harmless as a plain value but produces a
  # numerically noisy round trip through two solve() calls. Skip
  # straight to returning that site's own values, unchanged.
  #
  # This also makes a 1-site "federated" run byte-for-byte identical to
  # that site's own ds.semiOPBARTFitF() output, which matters for
  # diagnosing federated-vs-local performance gaps -- see the extended
  # note in this file's header on why isolating this case is the right
  # FIRST step before concluding anything about whether combining
  # multiple sites helps or hurts.
  if (length(site_fits) == 1) {
    print("Number of sites is 1. No combine parameters")
    print("")
    s <- site_fits[[1]]
    return(list(theta = s$theta_mean, theta_cov = s$theta_cov,
                us = s$us_mean, var_counts = s$var_counts_mean,
                weights_used = 1,
                per_site_fits = site_fits, combine_method = combine_method,
                linear_formula = attr(site_fits, "linear_formula"),
                state_name = attr(site_fits, "state_name")))
  }

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
  # theta_cov: was computed above (as V) for inverse-variance but only ever
  # used internally and then discarded -- keep it, and give "stack" a
  # reasonable equivalent (weighted mixture of the site covariances; not
  # exact since it ignores between-site mean disagreement, but a real
  # uncertainty estimate is more useful to an explainability plot than
  # none at all). Used by semiOPBART_explainCoefficients() for F's
  # coefficient forest plot -- see dsSemiOPBARTExplain.R.
  theta_cov_combined <- if (combine_method == "inverse_variance") {
    V
  } else {
    Reduce(`+`, Map(function(s, wi) wi^2 * s$theta_cov, site_fits, w))
  }
  us_combined <- Reduce(`+`, Map(function(s, wi) wi * s$us_mean, site_fits, w))

  # var_counts_mean: only present if every site's semiOPBARTLocalFitFDS()
  # is new enough to return it -- fall back to NULL (not an error) so an
  # older/mixed-version deployment still combines theta/us fine, it just
  # won't get a variable-importance plot.
  has_var_counts <- all(vapply(site_fits, function(s) !is.null(s$var_counts_mean), logical(1)))
  var_counts_combined <- if (has_var_counts) {
    Reduce(`+`, Map(function(s, wi) wi * s$var_counts_mean, site_fits, w))
  } else NULL

  list(theta = theta_combined, theta_cov = theta_cov_combined,
       us = us_combined, var_counts = var_counts_combined,
       weights_used = w,
       per_site_fits = site_fits, combine_method = combine_method,
       linear_formula = attr(site_fits, "linear_formula"),
       state_name = attr(site_fits, "state_name"))
}






























# # dsSemiOPBARTLocalCombineF.R
# # ---------------------------------------------------------------------------
# # Architecture F's combine -- pure client-side arithmetic, no datasources
# # argument, no aggregate() call, no disclosure surface at all: it only
# # ever touches what ds.semiOPBARTFitF() already returned. Pools theta and
# # (as of us_cov being added to semiOPBARTLocalFitFDS()) us across sites
# # ONCE, after every site's independent smopbart() run finished -- see
# # ds_localFit.R's header for how this differs from Architecture E.
# #
# # THREE combine_method OPTIONS:
# #   "inverse_variance" (default, unchanged from before): FIXED-EFFECT
# #     meta-analysis -- assumes every site is estimating the exact SAME
# #     true parameter value, and all cross-site variation is just
# #     within-site sampling noise. Precision-weighted average.
# #   "random_effects": adds a between-site heterogeneity variance (tau^2,
# #     estimated via REML) on top of each site's own within-site
# #     variance -- assumes sites may genuinely differ, not just be noisy
# #     estimates of one shared truth. See .semiOPBART_remlTau2()'s
# #     docstring for the estimating equation and citation (Viechtbauer
# #     2010). Recommended over "inverse_variance" specifically when sites
# #     might have real between-site heterogeneity (different patient
# #     populations/practices) -- which real clinical federated sites
# #     usually do. With very few sites (2-3, typical here), this is also
# #     specifically preferred over the more commonly-taught DerSimonian-
# #     Laird estimator, which is known to be biased in that regime -- see
# #     this file's design discussion for the full citation trail.
# #   "stack": simple n- or equal-weighted mixture, no variance modeling
# #     at all (the only option that still works if some site is missing
# #     theta_cov/us_cov entirely, e.g. from an older package version).
# #
# # BOTH theta and us go through the SAME combine_method now -- previously
# # us was ALWAYS simple-averaged regardless of combine_method, purely
# # because no us_cov existed to weight by (fixed via localFitDS.R). us[1]
# # (fixed at exactly 0 by construction) is never combined -- it's always
# # just 0, with 0 variance, at every site.
# #
# # ON A 1-SITE RUN UNDERPERFORMING THE NON-FEDERATED smopbart() RUN: if
# # you're comparing federated performance against local smopbart() and
# # seeing a drop even with only one site involved, the combine step below
# # is NOT the cause -- see the single-site short-circuit inside
# # ds.semiOPBARTCombineF(), which makes a 1-site run return that site's
# # ds.semiOPBARTFitF() output completely unchanged (no solve(), no
# # weighting, nothing that could introduce a difference). If a gap still
# # shows up at n_sites == 1, look instead at: (a) for Architecture D/E,
# # .semiOPBART_poolTheta()'s regularizing prior (prior_var = 100 by
# # default), which applies on EVERY sync regardless of site count and is
# # a real, deliberate difference from smopbart()'s flat/improper prior --
# # not present in Architecture F at all, since F's single-site path now
# # skips theta pooling entirely; (b) whether the production block-based
# # MCMC (semiOPBARTLocalMCMC/Improved) is actually running the SAME sweep
# # order as smopbart() -- see semiOPBARTLocalMCMC's own docstring for a
# # confirmed, still-unpatched deviation there; (c) the data-prep/
# # normalization pathway (ds.semiOPBARTPrepare) against whatever local
# # pipeline smopbart() was run through directly; (d) train/test split
# # RNG/stratification differences. Isolating and closing any 1-site gap
# # FIRST is what makes a later multi-site comparison actually mean
# # something -- right now it's confounded by whichever of these is live.
# # ---------------------------------------------------------------------------

# #' REML estimate of the between-site heterogeneity variance (tau^2) for a
# #' single parameter, via the standard intercept-only random-effects
# #' meta-analysis REML estimating equation -- the same equation
# #' metafor::rma(method = "REML") solves (see Viechtbauer, W. (2010).
# #' "Conducting Meta-Analyses in R with the metafor Package." Journal of
# #' Statistical Software, 36(3), 1-48). Solved with uniroot() rather than
# #' a hand-derived Newton-Raphson step -- a well-tested, bulletproof base-R
# #' root finder, since a subtly wrong analytic derivative here would be a
# #' much worse failure mode (silently wrong tau^2, no error at all) than
# #' uniroot() occasionally failing loudly.
# #'
# #' @param y  per-site point estimates for this ONE parameter (length K)
# #' @param v  per-site WITHIN-site sampling variances for this same
# #'   parameter (length K) -- e.g. diag(theta_cov)[j] at each site
# #' @return REML tau^2 estimate (>= 0). Exactly 0 is a real, common
# #'   outcome -- it means REML found no evidence of between-site
# #'   heterogeneity beyond within-site sampling noise, in which case this
# #'   reduces to the ordinary fixed-effect estimate automatically.
# #' @export
# .semiOPBART_remlTau2 <- function(y, v) {
#   K <- length(y)
#   if (K < 2) return(0)   # nothing to estimate heterogeneity from with <2 sites

#   score <- function(tau2) {
#     w <- 1 / (v + tau2)
#     mu_hat <- sum(w * y) / sum(w)
#     sum(w) - sum(w^2 * (y - mu_hat)^2) - sum(w^2) / sum(w)
#   }

#   if (score(0) <= 0) return(0)   # boundary solution -- REML's true optimum
#                                   # is at or below 0, clamped to 0 (a
#                                   # variance can't be negative)

#   upper <- max(10 * mean(v), diff(range(y))^2, 1) * 100
#   tryCatch(
#     stats::uniroot(score, interval = c(0, upper))$root,
#     error = function(e) 0   # if uniroot genuinely can't bracket a root
#                              # (shouldn't happen given the score(0) check
#                              # above and a generous upper bound), fall
#                              # back to fixed-effect rather than
#                              # propagating an error into the whole combine
#   )
# }

# #' Inverse-variance combination of one parameter across K sites, with an
# #' optional between-site heterogeneity variance (tau2) folded into the
# #' weights. tau2 = 0 (default) is the ordinary FIXED-EFFECT combine
# #' (identical to what this file always did for theta); a REML-estimated
# #' tau2 (see .semiOPBART_remlTau2()) makes it RANDOM-EFFECTS instead --
# #' same formula either way, just a different tau2 input, which is why
# #' both combine_method options in ds.semiOPBARTCombineF() below share
# #' this one function.
# #' @return list(mu, var, tau2) -- var is the variance of the COMBINED
# #'   estimate itself (not a per-site variance), already properly wider
# #'   than a tau2 = 0 combine would give when tau2 > 0.
# #' @export
# .semiOPBART_ivCombine <- function(y, v, tau2 = 0) {
#   w <- 1 / (v + tau2)
#   list(mu = sum(w * y) / sum(w), var = 1 / sum(w), tau2 = tau2)
# }

# #' Combine one parameter (vectors of per-site means `y` and per-site
# #' variances `v`) across sites for a given combine_method -- shared by
# #' both theta and us below so they're treated identically rather than us
# #' silently falling back to something simpler.
# .semiOPBART_combineParam <- function(y, v, combine_method) {
#   tau2 <- if (combine_method == "random_effects") .semiOPBART_remlTau2(y, v) else 0
#   .semiOPBART_ivCombine(y, v, tau2)
# }

# #' @param site_fits      output of ds.semiOPBARTFitF()
# #' @param combine_method "inverse_variance" (fixed-effect, default) --
# #'   full-matrix precision weighting, correctly modeling cross-parameter
# #'   correlation via each site's full theta_cov (and, when every site has
# #'   it, us_cov). "random_effects" (REML -- see this file's header) --
# #'   adds a per-parameter between-site heterogeneity variance, but
# #'   deliberately PER-PARAMETER only (no cross-parameter correlation
# #'   modeled), since full multivariate random-effects estimation is
# #'   impractical with the typically 2-3 sites this combines. "stack" --
# #'   simple n-/equal-weighted mixture, no variance modeling at all (the
# #'   only option that still works if some site is missing theta_cov/
# #'   us_cov entirely, e.g. an older package version).
# #' @param weight_by_n   used only for combine_method = "stack".
# #' @export
# ds.semiOPBARTCombineF <- function(site_fits,
#                                   combine_method = c( "random_effects","inverse_variance", "stack"),
#                                   weight_by_n = TRUE) {
#   combine_method <- match.arg(combine_method)

#   # SINGLE-SITE SHORT-CIRCUIT: with only one contributor, every weighting
#   # scheme below is mathematically a no-op anyway (weight = 1;
#   # inverse-variance pooling of a single precision matrix returns that
#   # same matrix, solve(solve(theta_cov)) == theta_cov up to floating-
#   # point noise; tau2 can't be estimated from one site either, see
#   # .semiOPBART_remlTau2()) -- but routing through solve()/uniroot()
#   # regardless adds real risk for zero benefit. Skip straight to
#   # returning that site's own values, unchanged.
#   #
#   # This also makes a 1-site "federated" run byte-for-byte identical to
#   # that site's own ds.semiOPBARTFitF() output, which matters for
#   # diagnosing federated-vs-local performance gaps -- see the extended
#   # note in this file's header on why isolating this case is the right
#   # FIRST step before concluding anything about whether combining
#   # multiple sites helps or hurts.
#   if (length(site_fits) == 1) {
#     s <- site_fits[[1]]
#     return(list(theta = s$theta_mean, theta_cov = s$theta_cov,
#                 us = s$us_mean, us_cov = s$us_cov, var_counts = s$var_counts_mean,
#                 weights_used = 1, tau2 = NULL,
#                 per_site_fits = site_fits, combine_method = combine_method,
#                 linear_formula = attr(site_fits, "linear_formula"),
#                 state_name = attr(site_fits, "state_name")))
#   }

#   n <- vapply(site_fits, `[[`, numeric(1), "n")
#   w <- if (weight_by_n) n / sum(n) else rep(1 / length(site_fits), length(site_fits))

#   if (combine_method == "stack") {
#     print("using stack parameter combination")
#     theta_combined <- Reduce(`+`, Map(function(s, wi) wi * s$theta_mean, site_fits, w))
#     theta_cov_combined <- Reduce(`+`, Map(function(s, wi) wi^2 * s$theta_cov, site_fits, w))
#     us_combined <- Reduce(`+`, Map(function(s, wi) wi * s$us_mean, site_fits, w))
#     theta_tau2 <- us_tau2 <- NULL

#   } else if (combine_method == "inverse_variance") {
#     print("using inverse_variance parameter combination")
#     # FULL-MATRIX precision weighting -- unchanged from before this
#     # turn's changes. Correctly models cross-parameter correlation (via
#     # each site's full theta_cov, not just its diagonal), which the
#     # random_effects branch below deliberately does NOT do -- kept
#     # exactly as it was rather than weakened for code-sharing
#     # convenience with a method that has a real reason not to do this.
#     precisions <- lapply(site_fits, function(s) solve(s$theta_cov))
#     V <- solve(Reduce(`+`, precisions))
#     theta_combined <- as.numeric(V %*% Reduce(`+`, Map(function(s, P) P %*% s$theta_mean,
#                                                           site_fits, precisions)))
#     theta_cov_combined <- V

#     if (all(vapply(site_fits, function(s) !is.null(s$us_cov), logical(1))) &&
#         length(site_fits[[1]]$us_mean) >= 2) {
#       us_precisions <- lapply(site_fits, function(s) solve(s$us_cov))
#       Vu <- solve(Reduce(`+`, us_precisions))
#       us_free_combined <- as.numeric(Vu %*% Reduce(`+`, Map(function(s, P)
#         P %*% s$us_mean[-1], site_fits, us_precisions)))
#       us_combined <- c(0, us_free_combined)   # us[1] is always exactly 0
#     } else {
#       us_combined <- Reduce(`+`, Map(function(s, wi) wi * s$us_mean, site_fits, w))
#     }
#     theta_tau2 <- us_tau2 <- NULL

#   } else {
#     print("using random effects parameter combination")
#     # "random_effects" -- PER PARAMETER (i.e. independently for each
#     # component of theta, and each free threshold of us), via
#     # .semiOPBART_combineParam()/.semiOPBART_remlTau2(). Unlike
#     # "inverse_variance" above, this deliberately does NOT model
#     # cross-parameter correlation: a full multivariate random-effects
#     # fit needs to estimate a between-site heterogeneity COVARIANCE
#     # matrix, which is a much harder, much less stable estimation
#     # problem than a single tau^2 per parameter, and essentially
#     # impractical with only 2-3 "studies"/sites -- see this file's
#     # header. Still properly accounts for each parameter's OWN
#     # between-site heterogeneity, just not their joint structure.
#     p <- length(site_fits[[1]]$theta_mean)
#     theta_means_mat <- do.call(rbind, lapply(site_fits, `[[`, "theta_mean"))       # K x p
#     theta_vars_mat  <- do.call(rbind, lapply(site_fits, function(s) diag(s$theta_cov)))  # K x p

#     theta_results <- lapply(seq_len(p), function(j)
#       .semiOPBART_combineParam(theta_means_mat[, j], theta_vars_mat[, j], combine_method))
#     theta_combined <- vapply(theta_results, `[[`, numeric(1), "mu")
#     theta_var_combined <- vapply(theta_results, `[[`, numeric(1), "var")
#     theta_tau2 <- vapply(theta_results, `[[`, numeric(1), "tau2")
#     # Cross-parameter correlation not modeled here -- see the note above
#     # -- so the combined covariance is diagonal, not a full matrix, even
#     # though each site's OWN theta_cov (in per_site_fits) is full.
#     theta_cov_combined <- diag(theta_var_combined, p)

#     # us: only if EVERY site has us_cov (added to semiOPBARTLocalFitFDS()
#     # alongside theta_cov -- an older/mixed-version deployment might not
#     # have it yet). Falls back to simple n-weighted averaging for us
#     # specifically if not, rather than erroring the whole combine over a
#     # missing field on a secondary parameter.
#     has_us_cov <- all(vapply(site_fits, function(s) !is.null(s$us_cov), logical(1)))
#     if (has_us_cov && length(site_fits[[1]]$us_mean) >= 2) {
#       n_thresh <- length(site_fits[[1]]$us_mean)
#       us_means_mat <- do.call(rbind, lapply(site_fits, `[[`, "us_mean"))          # K x n_thresh
#       # us_cov only covers the FREE thresholds (columns 2..n_thresh) --
#       # see semiOPBARTLocalFitFDS()'s docstring for why us[1] is excluded
#       us_vars_mat <- do.call(rbind, lapply(site_fits, function(s) diag(s$us_cov)))  # K x (n_thresh-1)

#       us_results <- lapply(seq_len(n_thresh - 1), function(j)
#         .semiOPBART_combineParam(us_means_mat[, j + 1], us_vars_mat[, j], combine_method))
#       us_free_combined <- vapply(us_results, `[[`, numeric(1), "mu")
#       us_tau2 <- c(NA_real_, vapply(us_results, `[[`, numeric(1), "tau2"))  # NA for the fixed us[1]
#       us_combined <- c(0, us_free_combined)   # us[1] is always exactly 0
#     } else {
#       us_combined <- Reduce(`+`, Map(function(s, wi) wi * s$us_mean, site_fits, w))
#       us_tau2 <- NULL
#     }
#   }

#   # var_counts_mean: only present if every site's semiOPBARTLocalFitFDS()
#   # is new enough to return it -- fall back to NULL (not an error) so an
#   # older/mixed-version deployment still combines theta/us fine, it just
#   # won't get a variable-importance plot. Always simple n-weighted
#   # averaging -- a split-count has no natural per-site "variance" to
#   # weight by the way theta/us do, so combine_method doesn't apply to it.
#   has_var_counts <- all(vapply(site_fits, function(s) !is.null(s$var_counts_mean), logical(1)))
#   var_counts_combined <- if (has_var_counts) {
#     Reduce(`+`, Map(function(s, wi) wi * s$var_counts_mean, site_fits, w))
#   } else NULL

#   list(theta = theta_combined, theta_cov = theta_cov_combined,
#        us = us_combined, var_counts = var_counts_combined,
#        weights_used = w, tau2 = list(theta = theta_tau2, us = us_tau2),
#        per_site_fits = site_fits, combine_method = combine_method,
#        linear_formula = attr(site_fits, "linear_formula"),
#        state_name = attr(site_fits, "state_name"))
# }