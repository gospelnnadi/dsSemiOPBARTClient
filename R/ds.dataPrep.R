# dsSemiOPBARTDataPrep.R
# ---------------------------------------------------------------------------
# Data preparation is now three explicit steps rather than one, so that
# normalization can be swapped independently of feature engineering:
#
#   1. TRANSFORM  -- vetted feature engineering + dummy encoding + outcome
#                    factor coercion. Produces UNNORMALISED X, plus W and Y.
#   2. RANGE       -- (only if normalize_method = "federated_minmax") each
#                    site reports its local per-column min/max on the
#                    transformed (not yet normalised) X; client aggregates
#                    to a single GLOBAL min/max per column and broadcasts it
#                    back. Skipped entirely for "local_ecdf".
#   3. NORMALIZE   -- applies either the ORIGINAL per-site local ECDF
#                    (matches the non-federated smopbart() exactly, but is
#                    only a same-site approximation of the pooled
#                    distribution -- the bias flagged early in this
#                    conversation) or the newly-broadcast federated
#                    min/max (every site normalises against the SAME global
#                    min/max, removing that bias, at a much smaller
#                    disclosure cost than the quantile-grid alternative
#                    floated earlier: two numbers per column instead of a
#                    whole grid).
#
# ds.semiOPBARTPrepare() below still offers the single-call convenience
# interface from before; it just now orchestrates three round trips instead
# of one when normalize_method = "federated_minmax".
#
# DISCLOSURE NOTE ON federated_minmax, READ BEFORE USING IT:
# the global min and global max ARE, by construction, an exact covariate
# value belonging to one real patient at one real site (whichever site
# happens to hold the true extreme). That is a smaller disclosure surface
# than publishing a full quantile grid (two numbers per column instead of
# 20-50), but each of those two numbers is maximally identifying FOR THAT
# ONE RECORD -- it is by definition the most unusual value for that
# covariate in the whole federation. If any covariate is itself sensitive
# or rare (e.g. an unusually high/low lab value that could re-identify a
# patient in combination with other public knowledge), don't default to
# this option for that covariate. It is offered as opt-in, per this
# request, not as the new default.
# ---------------------------------------------------------------------------

#source("ds.serialize.R")

# ---- vetted transform registry (extend by adding named entries here) ------

semiOPBARTTransformRegistry <- list(
  rheumatoid_ip_ep = function(d) {
    within(d, {
      bclass_IL6i  <- as.integer(bDMARD == 2)
      ESR_GAP      <- abs(ESR + SJC28)
      CRP_collapsed <- as.integer(10 * CRP < ESR)
      CSR_ESR_GAP  <- abs(ESR - SJC28)
      CRP_noscale  <- CRP + 1e-3 +
        (0.2 * CSR_ESR_GAP * CRP_collapsed) +
        (0.6 * CSR_ESR_GAP * bclass_IL6i * as.integer(CRP == 0))
      CRP <- CRP_noscale
    })
  },
  none = function(d) d
)

#' Canonical object/state names for one run_tag. Used by BOTH
#' ds.semiOPBARTClearState() AND run_semiOPBART_orchestrator() itself, so
#' the two can never drift out of sync -- the earlier hardcoded ClearState
#' list was wrong precisely because it was maintained BY HAND, separately
#' from what the orchestrator actually creates.
#'
#' @param run_tag  NULL for the legacy untagged convention (one run only,
#'   ever, per session). A short string (e.g. "ip", "ep") appends `_tag`
#'   to every object AND internal state name, so two runs (different
#'   pathology/model/architecture) never collide at the same site --
#'   this is what was missing when a second orchestrator call silently
#'   overwrote the first one's D_train/D_test/D_holdout AND its
#'   .semiOPBART_state, .semiOPBART_local_fit_F.
#' @export
semiOPBART_run_names <- function(run_tag = NULL) {
  suf <- function(base) if (is.null(run_tag)) base else paste0(base, "_", run_tag)
  list(
    transformed    = suf("semiOPBART_transformed"),
    train          = suf("D_train"),
    test           = suf("D_test"),
    holdout        = suf("D_holdout"),
    state_D        = suf(".semiOPBART_state_D"),
    state_E        = suf(".semiOPBART_state_E"),
    state_F        = suf(".semiOPBART_local_fit_F"),
    pred_D_test    = suf("pred_D_test"),
    pred_D_holdout = suf("pred_D_holdout"),
    pred_E_test    = suf("pred_E_test"),
    pred_E_holdout = suf("pred_E_holdout"),
    pred_F_test    = suf("pred_F_test"),
    pred_F_holdout = suf("pred_F_holdout")
  )
}

# ---- client side ------------------------------------------------------------
#' Clear semiOPBART state on every site.
#'
#' @param run_tag  NULL (default): clears the legacy untagged name set
#'   (semiOPBART_toSerializes.semiOPBART_run_names(NULL)). Pass the SAME
#'   run_tag a previous run_semiOPBART_orchestrator() call used to clear
#'   exactly that run's objects and internal state, without touching any
#'   OTHER run_tag's state that might still be in use at the same sites.
#' @export
ds.semiOPBARTClearState <- function(run_tag = NULL, datasources = NULL) {
  if (is.null(datasources)) {
    datasources <- DSI::datashield.connections_find()
  }
  names_Serialize <- semiOPBART_toSerialize(unlist(semiOPBART_run_names(run_tag),
                                                     use.names = FALSE))

  DSI::datashield.aggregate(
    datasources,
    call("semiOPBARTLocalClearStateDS", names_Serialize)
  )
}
# =====================  1. TRANSFORM  ========================================

# ---- client side ------------------------------------------------------------
#' @export
ds.semiOPBARTTransform <- function(data.name, outcome_col, levels,
                                    x_features, w_features,
                                    transform_recipe = "none",
                                    newobj = "semiOPBART_transformed",
                                    nfilter = 5, datasources = NULL) {
  if (is.null(datasources)) datasources <- DSI::datashield.connections_find()

  out <- DSI::datashield.aggregate(datasources,
      call("semiOPBARTLocalTransformDS", data.name, outcome_col,
           semiOPBART_toSerialize(levels), semiOPBART_toSerialize(x_features),
           semiOPBART_toSerialize(w_features), transform_recipe, newobj, nfilter))

  x_col_sets <- lapply(out, `[[`, "x_cols")
  if (length(unique(x_col_sets)) > 1)
    warning("sites produced DIFFERENT dummy-encoded X columns after transform ",
            "-- likely a categorical feature with different levels present ",
            "at different sites. Training/prediction downstream assumes a ",
            "common column set; fix by specifying factor levels explicitly ",
            "in the transform recipe rather than letting dummyVars infer them.")
  out
}

# =====================  2. FEDERATED RANGE  ==================================
# Only needed for normalize_method = "federated_minmax". Skip entirely for
# "local_ecdf".
# ---- client side ------------------------------------------------------------

#' Aggregate every site's local min/max into one global min/max per column
#' (global_min = min of local mins, global_max = max of local maxs), then
#' return it for the caller to pass into ds.semiOPBARTNormalizeSplit(). Does NOT
#' broadcast automatically -- kept explicit so the disclosure note is seen
#' at the call site, not buried inside a longer pipeline.
#'
#' @param data.name  should be the TRAIN object name (e.g. "D_train"),
#'   called against `datasources` restricted to train_sites only -- see
#'   ds.semiOPBARTPrepare() for how this is normally wired up
#' @export
ds.semiOPBARTComputeGlobalRange <- function(data.name = "semiOPBART_train",
                                             datasources = NULL) {
  if (is.null(datasources)) datasources <- DSI::datashield.connections_find()
  site_ranges <- DSI::datashield.aggregate(datasources,
      call("semiOPBARTLocalFeatureRangeDS", data.name))

  cols <- site_ranges[[1]]$cols
  if (!all(vapply(site_ranges, function(r) identical(r$cols, cols), logical(1))))
    stop("sites report different columns for feature range -- run ",
         "ds.semiOPBARTTransform() consistently at every site first")

  global_min <- do.call(pmin, lapply(site_ranges, `[[`, "min"))
  global_max <- do.call(pmax, lapply(site_ranges, `[[`, "max"))
  names(global_min) <- names(global_max) <- cols

  message("ds.semiOPBARTComputeGlobalRange(): global min/max computed. ",
          "Each of these ", 2 * length(cols), " numbers is an exact ",
          "covariate value from some real patient at some site -- see the ",
          "disclosure note in dsSemiOPBARTDataPrep.R before broadcasting ",
          "these onward via ds.semiOPBARTNormalizeSplit().")

  list(cols = cols, global_min = global_min, global_max = global_max)
}

#' Build a global_range object (same shape ds.semiOPBARTComputeGlobalRange()
#' returns) directly from a KNOWN, analyst-specified range dictionary --
#' NOT derived from the data at all. This is preferable whenever the
#' variables have well-established clinical/domain bounds (as is common
#' for structured EHR/registry variables): it eliminates the range-
#' computation round trip AND the exposure that comes with it entirely
#' (ds.semiOPBARTComputeGlobalRange()'s global_min/global_max are exact
#' real patient values; a known clinical bound is not data-derived at
#' all). Only the histogram bin counts (ds.semiOPBARTComputeGlobalECDF())
#' or, for federated_minmax, nothing further, ever needs to touch a site.
#'
#' TRADE-OFF, not a defect: if the actual data doesn't span the full
#' known range (e.g. this cohort's real Age_diagnosis is 20-75 within a
#' stated valid range of 18-120), histogram bins near the edges will be
#' sparse or empty -- use more bins, or clinically-informed non-uniform
#' edges (denser where data plausibly concentrates), to compensate.
#'
#' @param known_range  a named list, one c(min, max) per column, e.g.
#'   list(Age_diagnosis = c(1, 120), DAS28 = c(0, 10), Sex = c(0, 1), ...)
#' @export
range_dict_to_global_range <- function(known_range) {
  cols <- names(known_range)
  global_min <- vapply(known_range, `[`, numeric(1), 1)
  global_max <- vapply(known_range, `[`, numeric(1), 2)
  names(global_min) <- names(global_max) <- cols
  list(cols = cols, global_min = global_min, global_max = global_max)
}

# =====================  2b. FEDERATED HISTOGRAM (for federated_ecdf) =========
#
# federated_minmax anchors only the two ENDPOINTS of the normalized scale
# -- it says nothing about the SHAPE of the pooled distribution in
# between. If sites have very different marginal distributions of a
# covariate, min-max normalization keeps tree-sharing coherent (any
# consistent global rescaling does) but doesn't give the roughly-uniform
# spread that the real algorithm's local ecdf() normalization is designed
# to produce. This section builds a genuinely federated approximation to
# the POOLED ecdf, without ever moving a raw value: sites agree on shared
# bin edges (from the already-computed global range), report LOCAL COUNTS
# per bin, and those counts are summed across sites -- an aggregate
# statistic, structurally identical to any ds.glm-style sufficient
# statistic. A bin count is a much smaller disclosure surface than an
# exact min/max: it summarizes a whole interval of values, not one
# identifiable extreme.

# ---- client side ------------------------------------------------------------

#' Build shared bin edges from the global range, aggregate every site's
#' local bin counts (summed -- an aggregate, not any one site's raw
#' contribution), and turn the pooled cumulative counts into a piecewise-
#' linear approximation of the FEDERATION-WIDE ecdf per column.
#'
#' @param global_range  output of ds.semiOPBARTComputeGlobalRange() --
#'   reused as the outer bounds of the histogram, so this needs no
#'   separate min/max round trip of its own
#' @param num_bins  more bins = closer to the true pooled ecdf, at the
#'   cost of a smaller (and therefore more disclosure-sensitive) count
#'   per bin, and a larger payload broadcast to every site. 20-50 is a
#'   reasonable range for most covariates; the pooled (not per-site) bin
#'   count is what gets checked against `nfilter`.
#' @param nfilter  minimum POOLED (summed-across-sites) count per bin --
#'   mirrors the same "check the aggregate, not each site's raw
#'   contribution" pattern already used for confusion matrices in
#'   dsSemiOPBARTEvaluate.R. Per-site contributions to a bin (including
#'   zero) DO cross the wire raw, same disclosure caveat as that pattern.
#' @param custom_bin_edges  NULL (default): uniform `seq(a, b, length.out =
#'   num_bins + 1)` edges for every column, as before. Otherwise a named
#'   list, one numeric edge vector per column that needs NON-uniform
#'   edges (e.g. denser near where a right-skewed lab value like CRP/ESR
#'   actually concentrates) -- columns NOT named here still fall back to
#'   the uniform `num_bins` rule. Every supplied vector's first/last edge
#'   MUST equal that column's `global_range$global_min`/`global_max`
#'   exactly (checked below) so the pooled ecdf still spans the full
#'   analyst-declared range, matching what `federated_minmax` and any
#'   never-trained site's normalization reference both assume.
#' @export
ds.semiOPBARTComputeGlobalECDF <- function(global_range, num_bins = 30,
                                            data.name = "semiOPBART_train",
                                            nfilter = 5, custom_bin_edges = NULL,
                                            datasources = NULL) {
  if (is.null(datasources)) datasources <- DSI::datashield.connections_find()
  cols <- global_range$cols

  bin_edges <- setNames(lapply(cols, function(cn) {
    a <- global_range$global_min[cn]; b <- global_range$global_max[cn]
    if (a == b) return(c(a, a + 1))   # degenerate column -- single bin,
                                       # avoids a zero-width edge sequence
    if (!is.null(custom_bin_edges) && !is.null(custom_bin_edges[[cn]])) {
      edges <- custom_bin_edges[[cn]]
      if (edges[1] != a || edges[length(edges)] != b)
        stop("ds.semiOPBARTComputeGlobalECDF(): custom_bin_edges[['", cn,
             "']] must start at global_min (", a, ") and end at global_max (",
             b, ") -- got [", edges[1], ", ", edges[length(edges)], "]. ",
             "A partial-range histogram would silently under-cover the ",
             "column relative to what federated_minmax/never-trained-site ",
             "normalization both assume.")
      if (is.unsorted(edges, strictly = TRUE))
        stop("ds.semiOPBARTComputeGlobalECDF(): custom_bin_edges[['", cn,
             "']] must be strictly increasing.")
      return(edges)
    }
    seq(a, b, length.out = num_bins + 1)
  }), cols)

  site_hist <- DSI::datashield.aggregate(datasources,
      call("semiOPBARTLocalHistogramDS", data.name, semiOPBART_toSerialize(bin_edges)))

  pooled_counts <- setNames(lapply(cols, function(cn) {
    Reduce(`+`, lapply(site_hist, function(h) h$counts[[cn]]))
  }), cols)

  small_bins <- vapply(pooled_counts, function(cnts) any(cnts > 0 & cnts < nfilter), logical(1))
  if (any(small_bins))
    warning("ds.semiOPBARTComputeGlobalECDF(): pooled bin count below nfilter=",
            nfilter, " for column(s): ", paste(cols[small_bins], collapse = ", "),
            " -- consider fewer bins (coarser resolution) for these columns")

  cum_frac <- setNames(lapply(cols, function(cn) {
    cs <- cumsum(pooled_counts[[cn]])
    c(0, cs / cs[length(cs)])   # length num_bins + 1, aligned with bin_edges
  }), cols)

  message("ds.semiOPBARTComputeGlobalECDF(): federated ecdf computed from ",
          "pooled bin counts across all sites -- no individual raw value ",
          "left any site; only aggregated counts per bin did.")

  list(cols = cols, bin_edges = bin_edges, cum_frac = cum_frac, num_bins = num_bins)
}

# =====================  3. NORMALIZE (runs AFTER split, not before) =========
#
# CONFIRMED against the real smopbart() source: ECDFs are fit on X_train
# ONLY, then applied out-of-sample to X_test. This section replaces an
# earlier version that normalised BEFORE splitting -- which let a site's
# own test/holdout rows influence that SAME site's local_ecdf scale, a
# real (if mild) fidelity deviation, not a privacy issue. Normalize now
# takes the ALREADY-SPLIT train/test/holdout object names together, fits
# on train only (local_ecdf) or applies the shipped global range
# (federated_minmax) to whichever of the three objects exist at this site.

# ---- client side ------------------------------------------------------------
#' @export
ds.semiOPBARTNormalizeSplit <- function(train.name = "semiOPBART_train",
                                         test.name = "semiOPBART_test",
                                         holdout.name = "semiOPBART_holdout",
                                         normalize_method = c("local_ecdf", "federated_ecdf"),
                                         global_range = NULL, global_ecdf = NULL,
                                         nfilter = 5, datasources = NULL) {
  normalize_method <- match.arg(normalize_method)
  if (is.null(datasources)) datasources <- DSI::datashield.connections_find()

  gmin_Serialize <- gmax_Serialize <- gecdf_Serialize <- "null"
  if (normalize_method == "federated_ecdf") {
    if (is.null(global_ecdf))
      stop("normalize_method = 'federated_ecdf' requires `global_ecdf` ",
           "from ds.semiOPBARTComputeGlobalECDF()")
    gecdf_Serialize <- semiOPBART_toSerialize(global_ecdf)
  }

  DSI::datashield.aggregate(datasources,
      call("semiOPBARTLocalNormalizeSplitDS", train.name, test.name, holdout.name,
           normalize_method, gecdf_Serialize, nfilter))
}

# ---- convenience wrapper: the original single-call interface --------------

#' Orchestrates Transform -> Split -> (optional) ComputeGlobalRange/ECDF ->
#' NormalizeSplit, in that order -- normalization now correctly happens
#' AFTER splitting, fitting on train only where applicable (see the
#' section header above for why this order matters).
#'
#' @param normalize_method  "local_ecdf" (default), 
#'   "federated_ecdf" (a genuinely pooled empirical CDF from aggregated
#'   histogram bin counts -- see the section 2b header for why this is
#'   preferable to federated_minmax when distribution SHAPE, not just
#'   endpoints, matters, e.g. for Architecture D's tree-sharing).
#' @param known_range  NULL (default): the global range needed for
#'   "federated_ecdf" is computed FROM THE DATA via
#'   ds.semiOPBARTComputeGlobalRange() (exposes exact per-site extreme
#'   values -- see that function's disclosure note). A named list of
#'   c(min, max) per column (e.g. from a known clinical/domain range
#'   table): skips that round trip AND its exposure entirely -- the range
#'   becomes a priori knowledge, not data-derived. See
#'   range_dict_to_global_range()'s docstring for the resolution
#'   trade-off if the real data doesn't span the full known range.
#' @param num_ecdf_bins  only used for normalize_method = "federated_ecdf"
#'   -- see ds.semiOPBARTComputeGlobalECDF()'s `num_bins`
#' @param custom_bin_edges  only used for normalize_method = "federated_ecdf"
#'   -- passed straight through to ds.semiOPBARTComputeGlobalECDF()'s
#'   argument of the same name; a named list of non-uniform edge vectors
#'   for specific (typically right-skewed) columns, uniform `num_ecdf_bins`
#'   still applies to any column not named here
#' @param site_roles, train_ratio, seed  see ds.semiOPBARTSplit()
#' @export
ds.semiOPBARTPrepare <- function(data.name, outcome_col, levels,
                                  x_features, w_features,
                                  transform_recipe = "none",
                                  normalize_method = c("local_ecdf", "federated_ecdf"),
                                  known_range = NULL, num_ecdf_bins = 30,
                                  custom_bin_edges = NULL,
                                  site_roles = NULL, train_ratio = NULL, seed = NULL,
                                  newobj_train = "semiOPBART_train",
                                  newobj_test = "semiOPBART_test",
                                  newobj_holdout = "semiOPBART_holdout",
                                  nfilter = 5, datasources = NULL,    
                                  remove_rare_classes = TRUE,#FALSE,
                                  remove_classes_below = 15,# NULL,
                                  removed_newobj = "semiOPBART_filtered",
                                  balance_classes = FALSE,
                                  balance_target = 5
                                  ) {

  normalize_method <- match.arg(normalize_method)
  if (is.null(datasources)) datasources <- DSI::datashield.connections_find()
  if (nfilter < 1) stop("nfilter must be >= 1; it is the minimum POOLED count per bin.")
  if (remove_rare_classes) {

  if (is.null(remove_classes_below))
    stop(
      "`remove_classes_below` must be supplied when ",
      "`remove_rare_classes = TRUE`."
    )

  rare_result <- ds.semiOPBARTRemoveRareClasses( data.name = data.name, outcome_col = outcome_col, 
  levels = levels, remove_classes_below = remove_classes_below,
     newobj = removed_newobj,  datasources = datasources
  )
  if (rare_result$per_site != NULL){
    transform_data_name <- data.name
  }
  else{
      transform_data_name <- removed_newobj
      # use the levels that survived the pooled filtering
      levels <- rare_result$updated_levels
  }
  
} else {
  rare_result <- NULL
  transform_data_name <- data.name
  }
  # ds.semiOPBARTTransform(data.name, outcome_col, levels, x_features, w_features,
  #                         transform_recipe, newobj = "semiOPBART_transformed",
  #                         nfilter = nfilter, datasources = datasources)
  ds.semiOPBARTTransform(transform_data_name, outcome_col, levels, x_features, w_features,
                          transform_recipe, newobj = "semiOPBART_transformed",
                          nfilter = nfilter, datasources = datasources)
  print("ds.semiOPBARTPrepare(): transform complete, now splitting train/test/holdout...")
  split_plan <- ds.semiOPBARTSplit(
      data.name = "semiOPBART_transformed", outcome_col = outcome_col,
      site_roles = site_roles, train_ratio = train_ratio, seed = seed,
      newobj_train = newobj_train, newobj_test = newobj_test,
      newobj_holdout = newobj_holdout, nfilter = nfilter, datasources = datasources,levels,
      balance_classes = balance_classes, balance_target = balance_target,  balance_seed = seed)

  global_range <- NULL; global_ecdf <- NULL
  if (normalize_method %in% c("federated_ecdf") &&
      length(split_plan$train_sites)) {
    global_range <- if (!is.null(known_range)) {
      message("ds.semiOPBARTPrepare(): using a KNOWN range dictionary -- no ",
              "data-derived range round trip, no per-site extreme value ",
              "exposure for the range itself.")
      range_dict_to_global_range(known_range)
    } else {
      ds.semiOPBARTComputeGlobalRange(data.name = newobj_train,
          datasources = datasources[split_plan$train_sites])
    }
    if (normalize_method == "federated_ecdf")
      global_ecdf <- ds.semiOPBARTComputeGlobalECDF(global_range = global_range,
          num_bins = num_ecdf_bins, data.name = newobj_train,
          nfilter = nfilter, custom_bin_edges = custom_bin_edges,
          datasources = datasources[split_plan$train_sites])
  }
  print("ds.semiOPBARTPrepare(): normalization complete...")
  norm_result <- ds.semiOPBARTNormalizeSplit(
      train.name = newobj_train, test.name = newobj_test, holdout.name = newobj_holdout,
      normalize_method = normalize_method, global_range = global_range,
      global_ecdf = global_ecdf, nfilter = nfilter, datasources = datasources)
print("ds.semiOPBARTPrepare(): all steps complete.")
  # global_range/global_ecdf are returned (not just used internally) so the
  # caller -- typically the orchestrator -- can hang onto them for
  # ds.semiOPBARTPredictExternal() later, which needs the EXACT same
  # normalization reference a never-trained site is meant to be evaluated
  # against.
  list(split_plan = split_plan, per_site = norm_result,
       normalize_method = normalize_method,
       global_range = global_range, global_ecdf = global_ecdf,  levels=levels)
}

# =====================  SPLIT (unchanged logic, Serialize-safe args) =============
# ---- client side ------------------------------------------------------------

#' @param site_roles  named character vector ("train"/"test"/"split"), names
#'   matching names(datasources). NULL -> every site is "split" (or "train"
#'   if train_ratio is also NULL).
#' @param train_ratio single numeric (applied to every "split" site) or a
#'   named numeric vector (per-"split"-site ratio).
#' @param newobj_test     Group A: held-out rows WITH outcome_col --
#'   standard, label-based evaluation reads this population only.
#' @param newobj_holdout  Group B: Group A's rows again, unioned with every
#'   row that never had outcome_col -- evaluation against a separate
#'   validation_col reads this population only. Deliberately a different
#'   object from Group A; see semiOPBARTLocalSplitDS() for why they must
#'   stay separate rather than being merged into one test object.
#' @export
ds.semiOPBARTSplit <- function(data.name = "semiOPBART_transformed", outcome_col,
                                site_roles = NULL, train_ratio = NULL, seed = NULL,
                                newobj_train = "semiOPBART_train",
                                newobj_test  = "semiOPBART_test",
                                newobj_holdout = "semiOPBART_holdout",
                                nfilter = 5, datasources = NULL,
                                levels=0:4,
                                balance_classes = FALSE,
                                balance_target = 8,
                                balance_seed = NULL
) {
  if (is.null(datasources)) datasources <- DSI::datashield.connections_find()
  site_names <- names(datasources)
  if (is.null(site_names)) stop("datasources must be named so sites can be ",
                                 "matched to site_roles")

  if (is.null(site_roles))
    site_roles <- setNames(rep(if (is.null(train_ratio)) "train" else "split",
                                length(site_names)), site_names)
  if (!all(site_names %in% names(site_roles)))
    stop("site_roles missing an entry for: ",
         paste(setdiff(site_names, names(site_roles)), collapse = ", "))
  if (!all(site_roles %in% c("train", "test", "split")))
    stop("site_roles values must be 'train', 'test', or 'split'")

  ratio_for <- function(site) {
    if (is.null(train_ratio)) stop("train_ratio required for 'split'-role site: ", site)
    if (length(train_ratio) == 1 && is.null(names(train_ratio))) return(train_ratio)
    if (!site %in% names(train_ratio))
      stop("train_ratio has no entry for 'split'-role site: ", site)
    train_ratio[[site]]
  }

  raw_results <- lapply(site_names, function(site) {
    ds1 <- datasources[site]; role <- site_roles[[site]]
    if (role == "train") {
      DSI::datashield.aggregate(ds1, call("semiOPBARTLocalSplitDS", data.name,
          outcome_col, NULL, seed, newobj_train, newobj_test, newobj_holdout,
          nfilter,levels,balance_classes,balance_target, balance_seed))[[1]]
    } else if (role == "test") {
      c(DSI::datashield.aggregate(ds1,
          call("semiOPBARTLocalAsTestDS", data.name, newobj_test, newobj_holdout,
               nfilter))[[1]],
        auto_demoted = FALSE)
    } else {
      DSI::datashield.aggregate(ds1, call("semiOPBARTLocalSplitDS", data.name,
          outcome_col, ratio_for(site), seed, newobj_train, newobj_test,
          newobj_holdout, nfilter,levels, balance_classes,balance_target, balance_seed))[[1]]
    }
  })
  names(raw_results) <- site_names

  # AUTO-DEMOTION: a "train" or "split" site with too few labeled rows
  # comes back with auto_demoted = TRUE from semiOPBARTLocalSplitDS() --
  # move it from train_sites to holdout-only in the role bookkeeping here.
  demoted <- site_names[vapply(raw_results, function(r) isTRUE(r$auto_demoted), logical(1))]
  effective_roles <- site_roles
  if (length(demoted)) {
    effective_roles[demoted] <- "test"
    message("ds.semiOPBARTSplit(): the following sites had too few labeled ",
            "rows to train and were auto-demoted to predict-only: ",
            paste(demoted, collapse = ", "))
  }

  # test_sites: sites with a Group A object (n_test > 0) -- for standard
  # evaluation. holdout_sites: sites with a Group B object (n_holdout > 0)
  # -- for validation_col evaluation. These are DELIBERATELY tracked
  # separately, not derived from site_roles alone, since Group A can exist
  # even at role="test" sites (see semiOPBARTLocalAsTestDS()), and Group B
  # can be empty at a "train"-role site with no unlabeled rows at all.
  test_sites    <- site_names[vapply(raw_results, function(r) isTRUE(r$n_test > 0), logical(1))]
  holdout_sites <- site_names[vapply(raw_results, function(r) isTRUE(r$n_holdout > 0), logical(1))]

  list(per_site = raw_results,
       train_sites = site_names[effective_roles %in% c("train", "split")],
       test_sites = test_sites, holdout_sites = holdout_sites,
       site_roles = site_roles, effective_roles = effective_roles)
}


# ============================================================
# CLIENT-SIDE POOLED RARE-CLASS FILTERING
# ============================================================

#' Compute class counts pooled across all supplied sites.
#'
#' Only aggregate class counts are returned from each site.
#' Classes whose pooled count is below `remove_classes_below`
#' can subsequently be removed from every site.
#'
#' @export
ds.semiOPBARTComputeGlobalClassCounts <- function(
    data.name,
    outcome_col,
    levels,
    datasources = NULL) {

  if (is.null(datasources)) {
    datasources <- DSI::datashield.connections_find()
  }

  requested_levels <- as.character(levels)

  results <- DSI::datashield.aggregate(
    datasources,
    call(
      "semiOPBARTLocalClassCountsDS",
      data.name,
      outcome_col,
      requested_levels
    )
  )

  # Do NOT compare site-returned factor levels.
  # Every site is now guaranteed to return counts
  # in exactly requested_levels order.
  counts <- Reduce(
    `+`,
    lapply(results, function(x) {
      if (length(x$counts) != length(requested_levels)) {
        stop("Site returned an invalid number of class counts.")
      }

      as.numeric(x$counts)
    })
  )

  names(counts) <- requested_levels

  list(
    levels = requested_levels,
    counts = counts
  )
}

#' Remove outcome classes whose pooled count across all sites is
#' below a specified threshold.
#'
#' The pooled class counts are computed first. The corresponding
#' class labels are then sent to every site, and each site removes
#' those rows locally.
#'
#' Unlabeled rows are retained.
#'
#' @export
ds.semiOPBARTRemoveRareClasses <- function(
    data.name,
    outcome_col,
    levels,
    remove_classes_below,
    newobj = "semiOPBART_filtered",
    datasources = NULL) {

  if (is.null(datasources)) {
    datasources <- DSI::datashield.connections_find()
  }

  if (length(remove_classes_below) != 1L ||
      !is.numeric(remove_classes_below) ||
      is.na(remove_classes_below) ||
      remove_classes_below < 0) {
    stop("remove_classes_below must be a single non-negative number.")
  }

  levels <- as.character(levels)

  # Get pooled counts using the original requested levels
  global_counts <- ds.semiOPBARTComputeGlobalClassCounts(
    data.name = data.name,
    outcome_col = outcome_col,
    levels = levels,
    datasources = datasources
  )

  # Classes to remove
  remove_levels <- global_counts$levels[
    global_counts$counts < remove_classes_below
  ]

  # Levels remaining after rare-class removal
  updated_levels <- global_counts$levels[
    global_counts$counts >= remove_classes_below
  ]

  # Remove rare classes at every site
  # per_site <- DSI::datashield.aggregate(
  #   datasources,
  #   call(
  #     "semiOPBARTLocalRemoveClassesDS",
  #     data.name,
  #     outcome_col,
  #     remove_levels,
  #     newobj
  #   )
  # )
  # Remove rare classes at every site
if (length(remove_levels) == 0L) {

  # Nothing to remove. Keep the original object.
  per_site <- NULL


} else {

  per_site <- DSI::datashield.aggregate(
    datasources,
    call(
      "semiOPBARTLocalRemoveClassesDS",
      data.name,
      outcome_col,
      semiOPBART_toSerialize(remove_levels),
      newobj
    )
  )
}

  list(
    global_counts = global_counts,
    removed_classes = remove_levels,
    updated_levels = updated_levels,
    threshold = remove_classes_below,
    per_site = per_site
  )
}