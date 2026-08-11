`%||%` <- function(a, b) if (is.null(a)) b else a

# dsSemiOPBARTExternalPredict.R
# ---------------------------------------------------------------------------
# OPT-IN prediction at sites that never trained, Architecture D.
#
# UPDATED for the swap-based redesign (semiopbartDS.R/ds_semiopbart.R):
# there is no longer a single pooled "global" tree ensemble to reassemble
# -- every site keeps its OWN complete, permanent forest. So instead of
# pooling ALL sites' trees, the orchestrator/user now picks exactly ONE
# already-trained site whose CURRENT forest gets shipped as the reference
# model for never-trained sites (ds.semiOPBARTExportReferenceForest()
# below) -- same "pick a reference site" pattern already used for Option
# C earlier, just simpler here since there's no ecdf/dummyVars mismatch
# to work around: this site's forest was already grown on the shared
# federated_minmax/federated_ecdf scale (semiOPBARTLocalInitDS() enforces
# that hard requirement at training time), so it's directly usable
# anywhere else on that same scale, no site-specific normalization needs
# shipping alongside it.
#
# ds.semiOPBARTPredict() (dsSemiOPBARTTrainTest.R) is UNCHANGED: it still
# requires .semiOPBART_state (i.e. the site ran ds.semiOPBARTTrain()
# itself) and fails otherwise. That remains the default, lower-disclosure
# path for sites that trained.
#
# This file's function ships, to sites that never trained:
#   - ONE named site's own current tree ensemble (not everyone's, and not
#     the whole federation's aggregate num_tree * n_sites -- just that
#     one site's num_tree trees)
#   - the GLOBAL normalization reference ("global Z"): either the
#     global_min/global_max from ds.semiOPBARTComputeGlobalRange()
#     (normalize_method = "federated_minmax") or the pooled histogram-CDF
#     from ds.semiOPBARTComputeGlobalECDF() (normalize_method =
#     "federated_ecdf") -- whichever the original training run used. This
#     is what lets a site build a correctly-scaled X without ever having
#     trained.
#
# PRECONDITION: this only works if normalize_method was "federated_minmax"
# or "federated_ecdf" for the original ds.semiOPBARTPrepare() call.
# "local_ecdf" normalization is inherently site-specific -- there is no
# single global reference to ship, so this function refuses to run
# against it (enforced at training time already, by
# semiOPBARTLocalInitDS()'s own hard requirement).
# ---------------------------------------------------------------------------

#source("ds.serialize.R")

# ---- client side ------------------------------------------------------------

#' Pull ONE already-trained site's own current forest, to be shipped
#' elsewhere as the reference model for never-trained-site prediction.
#' Deliberately restricted to a single site -- this is what makes it a
#' "reference site" mechanism, not a second pooling step.
#'
#' @param fit             output of ds.semiOPBARTTrain() -- used only to
#'   confirm `reference_site` actually trained (fit$site_names)
#' @param reference_site  a single site name, must be one of `fit$site_names`
#' @export
ds.semiOPBARTExportReferenceForest <- function(fit, reference_site, datasources = NULL) {
  if (is.null(datasources)) datasources <- DSI::datashield.connections_find()
  if (length(reference_site) != 1)
    stop("ds.semiOPBARTExportReferenceForest(): reference_site must name exactly one site")
  if (!reference_site %in% fit$site_names)
    stop("ds.semiOPBARTExportReferenceForest(): '", reference_site,
         "' is not in fit$site_names (", paste(fit$site_names, collapse = ", "),
         ") -- it did not take part in the ds.semiOPBARTTrain() run `fit` came from")
  if (!reference_site %in% names(datasources))
    stop("ds.semiOPBARTExportReferenceForest(): '", reference_site,
         "' is not among `datasources`")

  message("ds.semiOPBARTExportReferenceForest(): exporting site '", reference_site,
          "'s own current forest (", fit$num_tree, " trees) as the reference ",
          "model for never-trained-site prediction. This ships that one ",
          "site's actual grown tree cutpoints -- more disclosive than the ",
          "default ds.semiOPBARTPredict(), which never leaves the training ",
          "sites at all. Confirm '", reference_site, "' has agreed to this.")

  out <- DSI::datashield.aggregate(datasources[reference_site],
      call("semiOPBARTLocalExportForestDS", fit$state_name %||% ".semiOPBART_state"))[[1]]
  list(trees_Serialize = out$trees_Serialize, num_tree = out$num_tree,
       reference_site = reference_site)
}

#' @param fit                output of ds.semiOPBARTTrain()
#' @param reference          output of ds.semiOPBARTExportReferenceForest()
#'   -- the ONE trained site's forest being shipped as the reference model
#' @param data.name_test     name of the object to predict at, at each of
#'   `datasources` -- see `already_prepared` for what this must contain
#' @param normalize_method   "federated_minmax" or "federated_ecdf" --
#'   MUST match whichever was used for the original training run's
#'   ds.semiOPBARTPrepare() call
#' @param global_range   REQUIRED if normalize_method = "federated_minmax":
#'   output of ds.semiOPBARTComputeGlobalRange() (or
#'   range_dict_to_global_range()) from the ORIGINAL training run
#' @param global_ecdf    REQUIRED if normalize_method = "federated_ecdf":
#'   output of ds.semiOPBARTComputeGlobalECDF() from the ORIGINAL training
#'   run. Either way, this is the "global Z" being shipped -- explicit,
#'   non-optional, because it's the ingredient that makes this function
#'   fundamentally different (and more disclosive) than the default
#'   ds.semiOPBARTPredict().
#' @param already_prepared   TRUE (default): `data.name_test` already
#'   points at a prepared list(X, W, Y, ...) object at every site in
#'   `datasources` (e.g. these sites were included in the original shared
#'   ds.semiOPBARTPrepare() call, just never trained). FALSE: these sites
#'   have never been prepared at all -- a genuinely new site -- and this
#'   function will run Transform + Normalize there first, using the
#'   shipped normalization reference (requires x_features, w_features,
#'   levels, transform_recipe to be supplied).
#' @param outcome_col  required in all cases, per the requirement that the
#'   orchestrator always names the label/data being used explicitly --
#'   used here only to build a placeholder Y column if a raw (unprepared)
#'   table is being prepared fresh; unlabeled rows are fine (see
#'   dsSemiOPBARTDataPrep.R's labeled-row handling).
#' @export
ds.semiOPBARTPredictExternal <- function(fit, reference, data.name_test, outcome_col,
                                          normalize_method = c("federated_minmax", "federated_ecdf"),
                                          global_range = NULL, global_ecdf = NULL,
                                          already_prepared = TRUE,
                                          x_features = NULL, w_features = NULL,
                                          levels = NULL, transform_recipe = "none",
                                          newobj_pred = "semiOPBART_pred",
                                          nfilter = 5, datasources = NULL, seed = 35) {
  set.seed(seed)
  normalize_method <- match.arg(normalize_method)
  if (normalize_method == "federated_minmax" && is.null(global_range))
    stop("global_range is required when normalize_method = 'federated_minmax' -- ",
         "pass the ds.semiOPBARTComputeGlobalRange()/range_dict_to_global_range() ",
         "result from the original training run")
  if (normalize_method == "federated_ecdf" && is.null(global_ecdf))
    stop("global_ecdf is required when normalize_method = 'federated_ecdf' -- ",
         "pass the ds.semiOPBARTComputeGlobalECDF() result from the original ",
         "training run")
  if (is.null(datasources)) datasources <- DSI::datashield.connections_find()

  if (!already_prepared) {
    if (is.null(x_features) || is.null(w_features) || is.null(levels))
      stop("x_features, w_features, and levels are required when ",
           "already_prepared = FALSE")
    ds.semiOPBARTTransform(data.name_test, outcome_col, levels, x_features,
                            w_features, transform_recipe,
                            newobj = "semiOPBART_external_transformed",
                            nfilter = nfilter, datasources = datasources)
    # NormalizeSplit(), not the old single-object Normalize() -- this is
    # an ad-hoc site with only a "test" object, no train/holdout, which is
    # exactly why normalize_method is hard-required to be one of the two
    # federated methods here: local_ecdf has nothing to fit ECDFs on at a
    # site that never trained. train.name/holdout.name = "null" (the
    # string sentinel, not bare NULL, per this codebase's wire-format
    # rule) tells NormalizeSplit those objects don't exist at this site.
    ds.semiOPBARTNormalizeSplit(train.name = "null",
                                 test.name = "semiOPBART_external_transformed",
                                 holdout.name = "null",
                                 normalize_method = normalize_method,
                                 global_range = global_range, global_ecdf = global_ecdf,
                                 nfilter = nfilter, datasources = datasources)
    data.name_test <- "semiOPBART_external_transformed"
  }

  message("ds.semiOPBARTPredictExternal(): shipping reference site '",
          reference$reference_site, "'s own current forest (",
          reference$num_tree, " trees) and the global normalization ",
          "reference (", normalize_method, ") to site(s) that never ",
          "trained: ", paste(names(datasources), collapse = ", "),
          ". This is strictly more information leaving that one training ",
          "site than the default ds.semiOPBARTPredict() ever exposes.")

  DSI::datashield.aggregate(datasources,
      call("semiOPBARTLocalPredictExternalDS", semiOPBART_toSerialize(fit$theta),
           semiOPBART_toSerialize(fit$us), reference$trees_Serialize, data.name_test,
           reference$num_tree, fit$k, newobj_pred, nfilter, seed = seed))
}
