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
      call("semiOPBARTLocalExportForestDS", fit$state_name))[[1]]

  if (is.null(out$local_ecdf_Serialize) &&
    is.null(out$federated_ecdf_Serialize)) {
  stop(
    "ds.semiOPBARTExportReferenceForest(): reference site exported ",
    "neither local_ecdf nor federated_ecdf."
  )
}
  list(
  trees_Serialize = out$trees_Serialize,
  num_tree = out$num_tree,

  local_ecdf_Serialize = out$local_ecdf_Serialize,
  federated_ecdf_Serialize = out$federated_ecdf_Serialize,

  has_local_ecdf = out$has_local_ecdf,
  has_federated_ecdf = out$has_federated_ecdf,

  reference_site = reference_site
)
}

#' Stage an externally exported forest and normalization reference
#'
#' Transfers the serialized reference forest and ECDF to destination sites
#' in chunks, avoiding embedding the complete serialized object in a single
#' DataSHIELD call.
#'
#' IMPORTANT:
#' - Never pads chunks with NA.
#' - Tree and ECDF streams are chunked independently.
#' - Empty strings are never sent as replacement data.
#' - The destination server assembles the exact original serialized strings.
#'
#' @param reference
#'   Output of ds.semiOPBARTExportReferenceForest().
#' @param reference_name
#'   Unique identifier for the staged reference.
#' @param datasources
#'   Destination DataSHIELD connections.
#' @param chunk_size
#'   Maximum number of characters per transfer.
#'
#' @return
#'   A list containing the destination server reference state name,
#'   normalization type, and basic diagnostics.
#'
#' @export
ds.semiOPBARTStageExternalReference <- function(
    reference,
    reference_name,
    datasources,
    chunk_size = 4000,
    quiet = TRUE
) {

  # ------------------------------------------------------------------
  # Validation
  # ------------------------------------------------------------------
   # PERFORMANCE: suppress progress bars if quiet
  if (quiet) {
    old_progress <- getOption("datashield.progress")
    options(datashield.progress = FALSE)
    on.exit(options(datashield.progress = old_progress), add = TRUE)
  }

  if (missing(reference)) {
    stop("reference is required")
  }

  if (missing(reference_name) ||
      length(reference_name) != 1L ||
      is.na(reference_name) ||
      !nzchar(reference_name)) {

    stop(
      "reference_name must be one non-empty, non-NA character string"
    )
  }

  if (missing(datasources) ||
      length(datasources) == 0L) {

    stop("datasources must contain at least one connection")
  }

  if (length(chunk_size) != 1L ||
      is.na(chunk_size) ||
      chunk_size <= 0) {

    stop(
      "chunk_size must be one positive integer"
    )
  }

  chunk_size <- as.integer(chunk_size)

  # ------------------------------------------------------------------
  # Validate forest serialization
  # ------------------------------------------------------------------

  trees_Serialize <- reference$trees_Serialize

  if (
    is.null(trees_Serialize) ||
    length(trees_Serialize) != 1L ||
    is.na(trees_Serialize)
  ) {

    stop(
      "reference$trees_Serialize must be one non-NA character string"
    )
  }

  if (!nzchar(trees_Serialize)) {
    stop(
      "reference$trees_Serialize is empty"
    )
  }

  # ------------------------------------------------------------------
  # Validate tree serialization header
  # ------------------------------------------------------------------

  if (!startsWith(
      trees_Serialize,
      "580a"
  )) {

    warning(
      "Tree serialization does not begin with expected '580a' header. ",
      "First bytes: ",
      substr(trees_Serialize, 1L, 32L)
    )
  }

  # ------------------------------------------------------------------
  # Select normalization reference
  #
  # For cross-site external prediction, federated_ecdf is preferred.
  # local_ecdf is only a fallback.
  # ------------------------------------------------------------------

  if (
    isTRUE(reference$has_federated_ecdf) &&
    !is.null(reference$federated_ecdf_Serialize)
  ) {

    ecdf_reference_Serialize <-
      reference$federated_ecdf_Serialize

    normalization_used <- "federated_ecdf"

  } else if (
    isTRUE(reference$has_local_ecdf) &&
    !is.null(reference$local_ecdf_Serialize)
  ) {

    ecdf_reference_Serialize <-
      reference$local_ecdf_Serialize

    normalization_used <- "local_ecdf"

  } else {

    stop(
      "Reference contains neither federated_ecdf_Serialize nor ",
      "local_ecdf_Serialize."
    )
  }

  if (
    length(ecdf_reference_Serialize) != 1L ||
    is.na(ecdf_reference_Serialize) ||
    !nzchar(ecdf_reference_Serialize)
  ) {

    stop(
      "Selected ECDF serialization is empty, NA, or not scalar"
    )
  }

  if (!startsWith(
      ecdf_reference_Serialize,
      "580a"
  )) {

    warning(
      "ECDF serialization does not begin with expected '580a' header. ",
      "First bytes: ",
      substr(ecdf_reference_Serialize, 1L, 32L)
    )
  }

  # ------------------------------------------------------------------
  # Tree chunk helper
  # ------------------------------------------------------------------

  make_chunks <- function(x, chunk_size) {

    n <- nchar(x)

    if (n == 0L) {
      stop("Cannot chunk an empty serialized object")
    }

    starts <- seq(
      from = 1L,
      to = n,
      by = chunk_size
    )

    chunks <- lapply(
      starts,
      function(start) {

        stop_at <- min(
          start + chunk_size - 1L,
          n
        )

        substring(
          x,
          start,
          stop_at
        )
      }
    )

    # Strong guarantee: no chunk is NULL or NA.
    bad <- vapply(
      chunks,
      function(z) {
        length(z) != 1L ||
        is.na(z) ||
        !nzchar(z)
      },
      logical(1)
    )

    if (any(bad)) {
      stop(
        "Internal chunking error: generated NULL/NA/empty chunk"
      )
    }

    chunks
  }

  # ------------------------------------------------------------------
  # Split independently
  # ------------------------------------------------------------------

  tree_chunks <- make_chunks(
    trees_Serialize,
    chunk_size
  )

  ecdf_chunks <- make_chunks(
    ecdf_reference_Serialize,
    chunk_size
  )

  n_tree_chunks <- length(tree_chunks)
  n_ecdf_chunks <- length(ecdf_chunks)

  # ------------------------------------------------------------------
  # Create a clean destination-side staging identifier
  # ------------------------------------------------------------------

  staging_id <- paste0(
    reference_name,
    "_",
    normalization_used
  )

  # ------------------------------------------------------------------
  # Store tree chunks
  # ------------------------------------------------------------------

  for (i in seq_along(tree_chunks)) {

    chunk <- tree_chunks[[i]]

    DSI::datashield.aggregate(
      datasources,
      call(
        "semiOPBARTLocalStoreExternalReferenceChunkDS",

        staging_id,
        "trees",

        chunk,

        as.integer(i),
        as.integer(n_tree_chunks)
      )
    )
  }

  # ------------------------------------------------------------------
  # Store ECDF chunks
  # ------------------------------------------------------------------

  for (i in seq_along(ecdf_chunks)) {

    chunk <- ecdf_chunks[[i]]

    DSI::datashield.aggregate(
      datasources,
      call(
        "semiOPBARTLocalStoreExternalReferenceChunkDS",

        staging_id,
        "ecdf",

        chunk,

        as.integer(i),
        as.integer(n_ecdf_chunks)
      )
    )
  }

  # ------------------------------------------------------------------
  # Assemble exact serialized strings on destination servers
  # ------------------------------------------------------------------

  assembled <-
    DSI::datashield.aggregate(
      datasources,
      call(
        "semiOPBARTLocalAssembleExternalReferenceDS",
        staging_id
      )
    )

  # ------------------------------------------------------------------
  # Destination state name
  # ------------------------------------------------------------------

  reference_state_name <- paste0(
    ".semiOPBART_external_reference_final_",
    staging_id
  )

  # ------------------------------------------------------------------
  # Diagnostics
  # ------------------------------------------------------------------

  if (length(assembled) != length(datasources)) {

    warning(
      "Number of assembly responses does not match number of destinations"
    )
  }

  message(
    "[ExternalReference] staged reference '",
    reference_name,
    "'"
  )

  message(
    "[ExternalReference] normalization = ",
    normalization_used
  )

  message(
    "[ExternalReference] tree serialization length = ",
    nchar(trees_Serialize)
  )

  message(
    "[ExternalReference] ECDF serialization length = ",
    nchar(ecdf_reference_Serialize)
  )

  message(
    "[ExternalReference] tree chunks = ",
    n_tree_chunks
  )

  message(
    "[ExternalReference] ECDF chunks = ",
    n_ecdf_chunks
  )

  message(
    "[ExternalReference] destination state = ",
    reference_state_name
  )

  list(
    reference_state_name = reference_state_name,

    normalization_used = normalization_used,

    num_tree = reference$num_tree,

    tree_serialized_length =
      nchar(trees_Serialize),

    ecdf_serialized_length =
      nchar(ecdf_reference_Serialize),

    n_tree_chunks = n_tree_chunks,
    n_ecdf_chunks = n_ecdf_chunks,

    assemble_result = assembled
  )
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
#' @param federated_ecdf    REQUIRED if normalize_method = "federated_ecdf":
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
ds.semiOPBARTPredictExternal <- function(
    fit,
    reference,
    data.name_test,
    outcome_col,
    levels = 0:4,
    newobj_pred = "semiOPBART_pred",
    nfilter = 5,
    datasources = NULL,
    seed = 35,
    quiet = TRUE
) {

  set.seed(seed)
   # PERFORMANCE: suppress progress bars if quiet
  if (quiet) {
    old_progress <- getOption("datashield.progress")
    options(datashield.progress = FALSE)
    on.exit(options(datashield.progress = old_progress), add = TRUE)
  }

  if (is.null(datasources)) {
    datasources <- DSI::datashield.connections_find()
  }
  # print("fit")
  # print(fit)

  # ---------------------------------------------------------------
  # Normalization reference MUST have been exported
  # ---------------------------------------------------------------
   normalization_used <- "unknown"
  if (
    !is.null(reference$federated_ecdf_Serialize) &&
    isTRUE(reference$has_federated_ecdf)
  ) {

    ecdf_reference_Serialize <-
      reference$federated_ecdf_Serialize

    normalization_used <- "federated_ecdf"
    print("ds.semiOPBARTPredictExternal(): using federated_ecdf normalization")

  } else if (
    !is.null(reference$local_ecdf_Serialize) &&
    isTRUE(reference$has_local_ecdf)
  ) {

    ecdf_reference_Serialize <-
      reference$local_ecdf_Serialize

    normalization_used <- "local_ecdf"
    print("ds.semiOPBARTPredictExternal(): using local_ecdf normalization")

  } else {

    stop(
      "ds.semiOPBARTPredictExternal(): reference contains neither ",
      "federated_ecdf nor local_ecdf."
    )
  }

  if (
    is.null(reference$num_tree) ||
    reference$num_tree <= 0
  ) {
    stop(
      "ds.semiOPBARTPredictExternal(): invalid reference$num_tree: ",
      reference$num_tree
    )
  }

  staged <- ds.semiOPBARTStageExternalReference(
  reference = reference,
  reference_name = paste0(
    "config_",
    reference$reference_site
  ),
  datasources = datasources
)

DSI::datashield.aggregate(
  datasources,
  call(
    "semiOPBARTLocalPredictExternalStoredDS",

    semiOPBART_toSerialize(fit$theta),
    semiOPBART_toSerialize(fit$us),

    staged$reference_state_name,

    data.name_test,
    staged$normalization_used,
     

    as.integer(reference$num_tree),
    as.numeric(fit$k),

    newobj_pred,

    semiOPBART_toSerialize(levels),

    as.integer(nfilter),

    seed = as.integer(seed)
  )
)

}

#' @export
ds.semiOPBARTsetExternalPredictedHoldout <- function(
    data.name,
    holdout.name,
    predicted_col = "predicted",
    newobj = "semiOPBART_dfPred",
    nfilter = 5,
    datasources = NULL
) {

  if (is.null(datasources)) {
    datasources <- datashield.connections_find()
  }

  res <- DSI::datashield.aggregate(
    conns = datasources,
    expr = call(
      "semiOPBARTsetExternalPredictedHoldoutDS",
      as.symbol(data.name),
      as.symbol(holdout.name),
      predicted_col,
      newobj,
      nfilter
    )
  )

  invisible(res)
}


# #' Pull ONE already-trained site's own current forest, to be shipped
# #' elsewhere as the reference model for never-trained-site prediction.
# #' Deliberately restricted to a single site -- this is what makes it a
# #' "reference site" mechanism, not a second pooling step.
# #'
# #' @param fit             output of ds.semiOPBARTTrain() -- used only to
# #'   confirm `reference_site` actually trained (fit$site_names)
# #' @param reference_site  a single site name, must be one of `fit$site_names`
# #' @export
# ds.semiOPBARTExportReferenceForest <- function(fit, reference_site, datasources = NULL) {
#   if (is.null(datasources)) datasources <- DSI::datashield.connections_find()
#   if (length(reference_site) != 1)
#     stop("ds.semiOPBARTExportReferenceForest(): reference_site must name exactly one site")
#   if (!reference_site %in% fit$site_names)
#     stop("ds.semiOPBARTExportReferenceForest(): '", reference_site,
#          "' is not in fit$site_names (", paste(fit$site_names, collapse = ", "),
#          ") -- it did not take part in the ds.semiOPBARTTrain() run `fit` came from")
#   if (!reference_site %in% names(datasources))
#     stop("ds.semiOPBARTExportReferenceForest(): '", reference_site,
#          "' is not among `datasources`")

#   message("ds.semiOPBARTExportReferenceForest(): exporting site '", reference_site,
#           "'s own current forest (", fit$num_tree, " trees) as the reference ",
#           "model for never-trained-site prediction. This ships that one ",
#           "site's actual grown tree cutpoints -- more disclosive than the ",
#           "default ds.semiOPBARTPredict(), which never leaves the training ",
#           "sites at all. Confirm '", reference_site, "' has agreed to this.")

#   out <- DSI::datashield.aggregate(datasources[reference_site],
#       call("semiOPBARTLocalExportForestDS", fit$state_name %||% ".semiOPBART_state"))[[1]]
#   list(trees_Serialize = out$trees_Serialize, num_tree = out$num_tree,
#        reference_site = reference_site)
# }

# #' @param fit                output of ds.semiOPBARTTrain()
# #' @param reference          output of ds.semiOPBARTExportReferenceForest()
# #'   -- the ONE trained site's forest being shipped as the reference model
# #' @param data.name_test     name of the object to predict at, at each of
# #'   `datasources` -- see `already_prepared` for what this must contain
# #' @param normalize_method   "federated_minmax" or "federated_ecdf" --
# #'   MUST match whichever was used for the original training run's
# #'   ds.semiOPBARTPrepare() call
# #' @param global_range   REQUIRED if normalize_method = "federated_minmax":
# #'   output of ds.semiOPBARTComputeGlobalRange() (or
# #'   range_dict_to_global_range()) from the ORIGINAL training run
# #' @param federated_ecdf    REQUIRED if normalize_method = "federated_ecdf":
# #'   output of ds.semiOPBARTComputeGlobalECDF() from the ORIGINAL training
# #'   run. Either way, this is the "global Z" being shipped -- explicit,
# #'   non-optional, because it's the ingredient that makes this function
# #'   fundamentally different (and more disclosive) than the default
# #'   ds.semiOPBARTPredict().
# #' @param already_prepared   TRUE (default): `data.name_test` already
# #'   points at a prepared list(X, W, Y, ...) object at every site in
# #'   `datasources` (e.g. these sites were included in the original shared
# #'   ds.semiOPBARTPrepare() call, just never trained). FALSE: these sites
# #'   have never been prepared at all -- a genuinely new site -- and this
# #'   function will run Transform + Normalize there first, using the
# #'   shipped normalization reference (requires x_features, w_features,
# #'   levels, transform_recipe to be supplied).
# #' @param outcome_col  required in all cases, per the requirement that the
# #'   orchestrator always names the label/data being used explicitly --
# #'   used here only to build a placeholder Y column if a raw (unprepared)
# #'   table is being prepared fresh; unlabeled rows are fine (see
# #'   dsSemiOPBARTDataPrep.R's labeled-row handling).
# #' @export
# ds.semiOPBARTPredictExternal <- function(fit, reference, data.name_test, outcome_col,
#                                           normalize_method = c("federated_minmax", "federated_ecdf"),
#                                           global_range = NULL, federated_ecdf = NULL,
#                                           already_prepared = TRUE,
#                                           x_features = NULL, w_features = NULL,
#                                           levels = NULL, transform_recipe = "none",
#                                           newobj_pred = "semiOPBART_pred",
#                                           nfilter = 5, datasources = NULL, seed = 35) {
#   set.seed(seed)
#   normalize_method <- match.arg(normalize_method)
#   if (normalize_method == "federated_minmax" && is.null(global_range))
#     stop("global_range is required when normalize_method = 'federated_minmax' -- ",
#          "pass the ds.semiOPBARTComputeGlobalRange()/range_dict_to_global_range() ",
#          "result from the original training run")
#   if (normalize_method == "federated_ecdf" && is.null(federated_ecdf))
#     stop("federated_ecdf is required when normalize_method = 'federated_ecdf' -- ",
#          "pass the ds.semiOPBARTComputeGlobalECDF() result from the original ",
#          "training run")
#   if (is.null(datasources)) datasources <- DSI::datashield.connections_find()

#   if (!already_prepared) {
#     if (is.null(x_features) || is.null(w_features) || is.null(levels))
#       stop("x_features, w_features, and levels are required when ",
#            "already_prepared = FALSE")
#     ds.semiOPBARTTransform(data.name_test, outcome_col, levels, x_features,
#                             w_features, transform_recipe,
#                             newobj = "semiOPBART_external_transformed",
#                             nfilter = nfilter, datasources = datasources)
#     # NormalizeSplit(), not the old single-object Normalize() -- this is
#     # an ad-hoc site with only a "test" object, no train/holdout, which is
#     # exactly why normalize_method is hard-required to be one of the two
#     # federated methods here: local_ecdf has nothing to fit ECDFs on at a
#     # site that never trained. train.name/holdout.name = "null" (the
#     # string sentinel, not bare NULL, per this codebase's wire-format
#     # rule) tells NormalizeSplit those objects don't exist at this site.
#     ds.semiOPBARTNormalizeSplit(train.name = "null",
#                                  test.name = "semiOPBART_external_transformed",
#                                  holdout.name = "null",
#                                  normalize_method = normalize_method,
#                                  global_range = global_range, federated_ecdf = federated_ecdf,
#                                  nfilter = nfilter, datasources = datasources)
#     data.name_test <- "semiOPBART_external_transformed"
#   }

#   message("ds.semiOPBARTPredictExternal(): shipping reference site '",
#           reference$reference_site, "'s own current forest (",
#           reference$num_tree, " trees) and the global normalization ",
#           "reference (", normalize_method, ") to site(s) that never ",
#           "trained: ", paste(names(datasources), collapse = ", "),
#           ". This is strictly more information leaving that one training ",
#           "site than the default ds.semiOPBARTPredict() ever exposes.")

#   DSI::datashield.aggregate(datasources,
#       call("semiOPBARTLocalPredictExternalDS", semiOPBART_toSerialize(fit$theta),
#            semiOPBART_toSerialize(fit$us), reference$trees_Serialize, data.name_test,
#            reference$num_tree, fit$k, newobj_pred, nfilter, seed = seed))
# }
