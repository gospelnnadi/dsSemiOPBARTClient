# dsSemiOPBARTLocalFitF.R
# ---------------------------------------------------------------------------
# ARCHITECTURE F: the ORIGINAL one-shot design, kept alongside Architecture
# E rather than replaced by it. Each site runs its own COMPLETE,
# independent smopbart() Gibbs run to convergence; only AFTER every site
# finishes does ds.semiOPBARTCombineF() pool theta_mean/theta_cov via
# inverse-variance meta-analysis. Sites never inform each other's f(x)/
# theta estimates DURING training -- that's the real difference from
# Architecture E's ds.semiOPBARTTrainE() (ds_semiopbart.R), which pools
# theta EVERY sweep instead. F is cheaper (one server round-trip per site,
# vs. num_burn+num_save round-trips for E) and simpler to reason about,
# at the cost of not letting sites' data actually inform each other's fit.
# ---------------------------------------------------------------------------

#' Run the local fit at every site and return the raw per-site summaries.
#' Pass into ds.semiOPBARTCombineF() -- pure client-side arithmetic,
#' re-runnable with a different combine_method/weight_by_n without
#' re-fitting anywhere.
#'
#' @param data.name  the TRAIN object at each site (from ds.semiOPBARTSplit())
#' @param state_name  where this site's fit gets cached (read later by
#'   semiOPBARTLocalPredictFDS()) -- MUST be unique per run (pathology/
#'   model) if more than one is ever trained in the same session, or a
#'   later run's fit silently overwrites this one at every site
#' @export
# ds.semiOPBARTFitF <- function(formula, linear_formula,
#                               data.name = "semiOPBART_train",
#                               datasources = NULL,
#                               num_tree = 20, k = 1,
#                               num_burn = 1000, num_save = 1000,
#                               nfilter = 5, state_name = ".semiOPBART_local_fit_F", seed = 35) {
#   set.seed(seed)
#   if (is.null(datasources)) datasources <- DSI::datashield.connections_find()

#   site_fits <- DSI::datashield.aggregate(datasources,
#       call("semiOPBARTLocalFitFDS", semiOPBART_toSerialize(deparse1(formula)),
#            semiOPBART_toSerialize(deparse1(linear_formula)),
#            data.name, num_tree, k, num_burn, num_save, nfilter, state_name, seed = seed))

#   attr(site_fits, "linear_formula") <- linear_formula
#   attr(site_fits, "data.name") <- data.name
#   attr(site_fits, "state_name") <- state_name
#   site_fits
# }

ds.semiOPBARTFitF <- function(formula, linear_formula,
                              data.name = "semiOPBART_train",
                              datasources = NULL,
                              num_tree = 20, k = 1,
                              num_burn = 1000, num_save = 1000,
                              nfilter = 5, state_name = ".semiOPBART_local_fit_F", seed = 35) {
  set.seed(seed)
  
  # ---- FILTER: Only connections with train object ----
  if (is.null(datasources)) datasources <- DSI::datashield.connections_find()
  
  trainable_conns <- .filter_trainable_connections(
    conns = datasources,
    train.name = data.name,
    verbose = TRUE
  )
  
  if (length(trainable_conns) == 0) {
    stop("No trainable sites. All sites were auto-demoted.")
  }
  
  site_names <- names(trainable_conns)
  message(sprintf("[ds.semiOPBARTFitF] Fitting on %d site(s): %s", 
                  length(site_names), paste(site_names, collapse = ", ")))

  site_fits <- DSI::datashield.aggregate(trainable_conns,
      call("semiOPBARTLocalFitFDS", semiOPBART_toSerialize(deparse1(formula)),
           semiOPBART_toSerialize(deparse1(linear_formula)),
           data.name, num_tree, k, num_burn, num_save, nfilter, state_name, seed = seed))

  attr(site_fits, "linear_formula") <- linear_formula
  attr(site_fits, "data.name") <- data.name
  attr(site_fits, "state_name") <- state_name
  site_fits
}