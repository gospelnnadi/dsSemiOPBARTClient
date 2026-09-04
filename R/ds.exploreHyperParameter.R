# dsSemiOPBARTExploreHyperParameter.R
# ---------------------------------------------------------------------------
# FEDERATED HYPERPARAMETER EXPLORATION -- Architectures D, E, F
#
# WHY THIS FILE WAS REWRITTEN
# ---------------------------------------------------------------------------
# The previous version of this file called seven server-side aggregate
# functions that do not exist anywhere in the dsSemiOPBART package:
#   semiOPBARTLocalFitFHyper, semiOPBARTLocalPredictFHyper,
#   semiOPBARTLocalInitEHyper, semiOPBARTLocalGibbsEHyper,
#   semiOPBARTLocalPredictEHyper, semiOPBARTLocalInitDHyper,
#   semiOPBARTLocalGibbsDHyper
# (grep across every server file in this package confirms none of these
# names are ever `<- function(...)`-defined -- only *called*, from here.)
# Every `datashield.aggregate()` call using one of those names would fail
# on a real Opal/DSLite server. On top of that, even where a call somehow
# returned a result, the code assumed the *raw predictions*
# (`p$pred`, `p$probs`) come back to the client -- they never do, by
# design (see trainTestDS.R's header): predict functions store
# `list(prob, map, truth)` LOCALLY at each site and return only a row
# count. Comparing `p$pred` to `p$probs[,1]` was therefore comparing two
# things that were never both present, which is why test_acc was always
# exactly 0 and test_mae was some meaningless leftover NA-arithmetic
# artifact -- not an evaluation of the model at all.
#
# THE FIX: don't reimplement the Gibbs loop / prediction / evaluation a
# second time for "exploration" purposes. Every one of D, E, F already has
# a complete, tested client entrypoint:
#   Architecture F: ds.semiOPBARTFitF() + ds.semiOPBARTCombineF()
#                   + ds.semiOPBARTPredictF()      (ds_localFit.R,
#                     ds_localCombine.R, ds_trainTest.R)
#   Architecture E: ds.semiOPBARTTrainE() + ds.semiOPBARTPredictE()
#                     (ds_semiopbart.R, ds_trainTest.R)
#   Architecture D: ds.semiOPBARTTrain()  + ds.semiOPBARTPredict()
#                     (ds_semiopbart.R, ds_trainTest.R)
# and evaluation for all three funnels through the SAME disclosure-safe,
# confusion-matrix-based ds.semiOPBARTEvaluate() (ds_evaluate.R). This
# file now does nothing but grid-search over that existing machinery,
# catch per-config failures, and tabulate the results.
# ---------------------------------------------------------------------------

# ============================================================================
# SHARED HELPERS
# ============================================================================

#' Cap a parameter grid at `max_configs` rows via a reproducible random
#' subsample, same behaviour the old code had (kept as a shared helper
#' instead of copy-pasted three times).
.semiOPBART_capConfigs <- function(param_grid, max_configs, seed = 42) {
  if (!is.null(max_configs) && nrow(param_grid) > max_configs) {
    set.seed(seed)
    param_grid <- param_grid[sample(nrow(param_grid), max_configs), , drop = FALSE]
    rownames(param_grid) <- NULL
  }
  param_grid
}



#' Turn one ds.semiOPBARTEvaluate() result (or NULL, if evaluation wasn't
#' run/failed for this config) into the flat metrics used in every results
#' row below.
#'
#' Deliberately does NOT surface raw ("Accuracy") as a headline metric.
#' These outcomes (IP.probability_0_4_BM / EP.probability_0_4_BM, levels
#' 0:4) are ordinal and typically class-imbalanced -- a model that always
#' predicts the majority class can score a high raw accuracy while being
#' useless for the minority categories. Balanced accuracy (mean per-class
#' recall) and macro F1 (mean per-class F1) both weight every class
#' equally regardless of its size, so they don't reward that failure mode
#' the way raw accuracy does. Both are already computed by
#' ds.semiOPBARTEvaluate() (ds_evaluate.R) from the same confusion matrix
#' -- this just picks them out instead of "Accuracy".
.semiOPBART_evalToMetrics <- function(eval_res, levels) {
  if (is.null(eval_res)) {
    return(list(test_mae = NA_real_, test_acc = NA_real_, test_balacc = NA_real_,
                test_f1 = NA_real_, test_kappa = NA_real_,
                n_test = NA_integer_))
  }
  mae <- eval_res$metrics$value[eval_res$metrics$metric == "MAE"]
  acc <- eval_res$metrics$value[eval_res$metrics$metric == "Accuracy"]
  balacc <- eval_res$metrics$value[eval_res$metrics$metric == "Balanced Accuracy"]
  f1     <- eval_res$metrics$value[eval_res$metrics$metric == "Macro F1"]
  kappa  <- eval_res$metrics$value[eval_res$metrics$metric == "Kappa"]

  list(
    test_mae    = if (length(mae)) mae else NA_real_, 
    test_acc = if (length(acc)) acc else NA_real_,
    test_balacc = if (length(balacc)) balacc else NA_real_,
    test_f1     = if (length(f1))     f1     else NA_real_,
    test_kappa  = if (length(kappa))  kappa  else NA_real_,
    n_test      = eval_res$n_total
  )
}

#' Run `expr_fun()`, printing any error immediately via `message()` (which
#' always flushes right away, unlike `warning()` -- warnings can be
#' buffered by R and never surface if the whole exploration run later
#' errors out, which is exactly what was happening here: every config was
#' failing silently and only the final "No configurations completed"
#' summary error was visible). Returns NULL on failure instead of
#' aborting the whole exploration run.
.semiOPBART_safe <- function(expr_fun, config_id, step_label) {
  tryCatch(expr_fun(), error = function(e) {
    message(sprintf("[Config %d] %s FAILED: %s", config_id, step_label, conditionMessage(e)))
    NULL
  })
}

#' One empty metrics row, used whenever test_data.name is NULL or every
#' step up to evaluation failed -- keeps every architecture's results data
#' frame the same shape.
.semiOPBART_naMetrics <- function() {
  list(test_mae = NA_real_, test_balacc = NA_real_,
       test_f1 = NA_real_, test_kappa = NA_real_, n_test = NA_integer_)
}


# ============================================================================
# ARCHITECTURE F HYPERPARAMETER EXPLORATION
# ============================================================================

#' Client: Explore hyperparameters for Architecture F, using the real
#' ds.semiOPBARTFitF() / ds.semiOPBARTCombineF() / ds.semiOPBARTPredictF()
#' / ds.semiOPBARTEvaluate() pipeline for every grid point.
#'
#' @param conns  list of DataSHIELD connections to the sites
#' @param data.name  the TRAIN object at each site (from ds.semiOPBARTPrepare()/Split())
#' @param outcome_col  name of the outcome column in the data
#' @param levels  factor levels of the outcome variable (analyst-specified,
#'   same set passed to ds.semiOPBARTPrepare())
#' @param x_features  vector of names of the predictor columns
#' @param w_features  vector of names of the linear predictor columns
#' @param test_data.name  name of the TEST object at each site (optional;
#'   without it, only n_sites/config metadata is reported, no test metrics)
#' @param n_samp_range,n_burn_range,n_tree_range,seed_range,k_values  grids
#'   to sweep, forwarded to ds.semiOPBARTFitF()'s num_save/num_burn/num_tree/seed/k
#' @param max_configs  cap the grid at this many randomly-sampled rows (optional)
#' @param nfilter  disclosure floor forwarded to fit/predict/evaluate
#' @param combine_method  "inverse_variance" or "stack" -- forwarded to ds.semiOPBARTCombineF()
#' @param verbose  print progress
#' @param state_name  base name for this run's per-config server-side state
#' @export
ds.semiOPBARTExploreHyperF <- function(conns,
                                        data.name,
                                        outcome_col,
                                        levels,
                                        x_features,
                                        w_features,
                                        test_data.name = NULL,
                                        n_samp_range = c(500, 1000, 2000),
                                        n_burn_range = c(500, 1000, 2000),
                                        n_tree_range = c(50, 100, 200),
                                        seed_range = c(42, 123, 456),
                                        k_values = c(1, 2, 3),
                                        max_configs = NULL,
                                        nfilter = 5,
                                        combine_method = c("random_effects", "inverse_variance", "stack"),
                                        verbose = TRUE,
                                        state_name = ".semiOPBART_hyper") {

  combine_method <- match.arg(combine_method)
  if (is.null(conns)) conns <- DSI::datashield.connections_find()

   # ---- FILTER: Only connections with train object ----
  trainable_conns <- .filter_trainable_connections(
    conns = conns,
    train.name = data.name,
    verbose = verbose
  )
  
  site_names <- names(trainable_conns)
  #site_names <- names(conns)

  param_grid <- expand.grid(
    n_samp = n_samp_range, n_burn = n_burn_range, n_tree = n_tree_range,
    seed = seed_range, k = k_values, stringsAsFactors = FALSE
  )
  param_grid <- .semiOPBART_capConfigs(param_grid, max_configs)

  formula_obj        <- stats::as.formula(paste(outcome_col, "~", paste(x_features, collapse = " + ")))
  linear_formula_obj <- stats::as.formula(paste("~", paste(w_features, collapse = " + ")))

  if (verbose) message(sprintf("Testing %d configurations for Architecture F across %d sites",
                                nrow(param_grid), length(site_names)))

  all_results <- vector("list", nrow(param_grid))
  if (verbose) pb <- txtProgressBar(min = 0, max = nrow(param_grid), style = 3)

  for (i in seq_len(nrow(param_grid))) {
    params <- param_grid[i, ]
    if (verbose) setTxtProgressBar(pb, i)
    print( paste0(state_name, "_F_config_", i))
    config_state <- paste0(state_name, "_config_")

    site_fits <- .semiOPBART_safe(function() {
      ds.semiOPBARTFitF(
        formula = formula_obj, linear_formula = linear_formula_obj,
        data.name = data.name, datasources = trainable_conns,
        num_tree = as.integer(params$n_tree), k = params$k,
        num_burn = as.integer(params$n_burn), num_save = as.integer(params$n_samp),
        nfilter = nfilter, state_name = config_state, seed = as.integer(params$seed))
    }, i, "fit (ds.semiOPBARTFitF)")
    if (is.null(site_fits)) next

    combined <- .semiOPBART_safe(function() {
      ds.semiOPBARTCombineF(site_fits, combine_method = combine_method)
    }, i, "combine (ds.semiOPBARTCombineF)")
    if (is.null(combined)) next

    metrics <- .semiOPBART_naMetrics()
    if (!is.null(test_data.name)) {
      pred_obj <- paste0(config_state, "_pred")
      pred_ok <- .semiOPBART_safe(function() {
        ds.semiOPBARTPredictF(
          combined, data.name_test = test_data.name, outcome_col = outcome_col,
          newobj_pred = pred_obj, nfilter = nfilter, datasources = trainable_conns,
          seed = as.integer(params$seed))
      }, i, "predict (ds.semiOPBARTPredictF)")

      if (!is.null(pred_ok)) {
        eval_res <- .semiOPBART_safe(function() {
          ds.semiOPBARTEvaluate(
            pred_obj = pred_obj, levels = levels, nfilter = nfilter,
            datasources = trainable_conns, label = "F_config_")
        }, i, "evaluate (ds.semiOPBARTEvaluate)")
        metrics <- .semiOPBART_evalToMetrics(eval_res, levels)
      }
    }




    all_results[[i]] <- data.frame(
      config_id = i, n_samp = as.integer(params$n_samp), n_burn = as.integer(params$n_burn),
      n_tree = as.integer(params$n_tree), seed = as.integer(params$seed), k = params$k,
      n_swap = NA_integer_, sync_every = NA_integer_, architecture = "F", n_sites = length(site_fits),
      test_mae = metrics$test_mae, test_acc = metrics$test_acc, test_balacc = metrics$test_balacc,
      test_f1 = metrics$test_f1, test_kappa = metrics$test_kappa,
      n_test = metrics$n_test, stringsAsFactors = FALSE
    )
  }
  if (verbose) close(pb)

  results_df <- dplyr::bind_rows(all_results)
  if (nrow(results_df) == 0) stop("No configurations completed successfully for Architecture F")

  list(results = results_df, architecture = "F",
       best_config = results_df[which.min(results_df$test_mae), ],
       param_grid = param_grid, site_names = site_names)
}


# ============================================================================
# ARCHITECTURE E HYPERPARAMETER EXPLORATION
# ============================================================================

#' Client: Explore hyperparameters for Architecture E, using the real
#' ds.semiOPBARTTrainE() / ds.semiOPBARTPredictE() / ds.semiOPBARTEvaluate()
#' pipeline (periodic theta/us sync, trees always local) for every grid point.
#'
#' @param conns  list of DataSHIELD connections to the sites
#' @param data.name  the TRAIN object at each site
#' @param outcome_col,x_features,w_features  unused by Architecture E's own
#'   training call (ds.semiOPBARTTrainE() reads directly from the already-
#'   prepared `data.name` object) but kept for a consistent signature across
#'   ds.semiOPBARTExploreHyperF/E/D() and for evaluate()'s `levels`
#' @param levels  factor levels of the outcome (forwarded to ds.semiOPBARTEvaluate())
#' @param test_data.name  name of the TEST object at each site (optional)
#' @param n_samp_range,n_burn_range,n_tree_range,seed_range,k_values  grids to sweep
#' @param warmup_burn  warm-start iterations before the main loop, forwarded to
#'   ds.semiOPBARTTrainE(); default `min(50, n_burn)` per config, kept short
#'   because hyperparameter exploration is already many configs x many sweeps
#' @param max_configs  cap the grid at this many randomly-sampled rows (optional)
#' @param nfilter  disclosure floor forwarded to train/predict/evaluate
#' @param verbose  print progress
#' @param state_name  base name for this run's per-config server-side state
#' @export
ds.semiOPBARTExploreHyperE <- function(conns,
                                        data.name,
                                        outcome_col,
                                        levels,
                                        x_features,
                                        w_features,
                                        test_data.name = NULL,
                                        n_samp_range = c(500, 1000, 2000),
                                        n_burn_range = c(500, 1000, 2000),
                                        n_tree_range = c(50, 100, 200),
                                        seed_range = c(42, 123, 456),
                                       sync_every= c(10,15, 20,50,100),
                                        k_values = c(1, 2, 3),
                                        warmup_burn =1, #  NULL, 
                                        max_configs = NULL,
                                        nfilter = 5,
                                        verbose = TRUE,
                                        state_name = ".semiOPBART_hyper") {

  if (is.null(conns)) conns <- DSI::datashield.connections_find()
  # ---- FILTER: Only connections with train object ----
  trainable_conns <- .filter_trainable_connections(
    conns = conns,
    train.name = data.name,
    verbose = verbose
  )
  
  site_names <- names(trainable_conns)
  #site_names <- names(conns)

  param_grid <- expand.grid(
    n_samp = n_samp_range, n_burn = n_burn_range, n_tree = n_tree_range,
    seed = seed_range,sync_every=sync_every, k = k_values, stringsAsFactors = FALSE
  )
  param_grid <- .semiOPBART_capConfigs(param_grid, max_configs)

  if (verbose) message(sprintf("Testing %d configurations for Architecture E across %d sites",
                                nrow(param_grid), length(site_names)))

  all_results <- vector("list", nrow(param_grid))
  if (verbose) pb <- txtProgressBar(min = 0, max = nrow(param_grid), style = 3)

  for (i in seq_len(nrow(param_grid))) {
    print(paste0("i  : ", i))
    params <- param_grid[i, ]
    if (verbose) setTxtProgressBar(pb, i)
    print( paste0(state_name, "_E_config_", i))
    config_state <- paste0(state_name, "_config_")
    wb <- if (is.null(warmup_burn)) min(50, as.integer(params$n_burn)) else warmup_burn

    fit <- .semiOPBART_safe(function() {
      ds.semiOPBARTTrainE(
        data.name = data.name, datasources = trainable_conns,
        num_tree = as.integer(params$n_tree), k = params$k,
        num_burn = as.integer(params$n_burn), num_save = as.integer(params$n_samp),
        sync_every=as.integer(params$sync_every),
        warmup_burn = wb, nfilter_threshold = nfilter,
        state_name = config_state, seed = as.integer(params$seed))
      }, i, "train (ds.semiOPBARTTrainE)")

    if (is.null(fit)) next

    metrics <- .semiOPBART_naMetrics()
    if (!is.null(test_data.name)) {
      pred_obj <- paste0(config_state, "_pred")
      pred_ok <- .semiOPBART_safe(function() {
        ds.semiOPBARTPredictE(
          fit, data.name_test = test_data.name, newobj_pred = pred_obj,
          #prediction_method = "point",
          nfilter = nfilter, datasources = trainable_conns,
          seed = as.integer(params$seed))
      }, i, "predict (ds.semiOPBARTPredictE)")

      if (!is.null(pred_ok)) {
        eval_res <- .semiOPBART_safe(function() {
          ds.semiOPBARTEvaluate(
            pred_obj = pred_obj, levels = levels, nfilter = nfilter,
            datasources = trainable_conns, label = "E_config_" )
        }, i, "evaluate (ds.semiOPBARTEvaluate)")
        metrics <- .semiOPBART_evalToMetrics(eval_res, levels)
      }
    }

    all_results[[i]] <- data.frame(
      config_id = i, n_samp = as.integer(params$n_samp), n_burn = as.integer(params$n_burn),
      n_tree = as.integer(params$n_tree),
      sync_every = as.integer(params$sync_every),
      seed = as.integer(params$seed), k = params$k,
      n_swap = NA_integer_, architecture = "E", n_sites = length(fit$site_names),
      test_mae = metrics$test_mae,test_acc = metrics$test_acc, test_balacc = metrics$test_balacc,
      test_f1 = metrics$test_f1, test_kappa = metrics$test_kappa,
      n_test = metrics$n_test, stringsAsFactors = FALSE
    )
  }
  if (verbose) close(pb)

  results_df <- dplyr::bind_rows(all_results)
  if (nrow(results_df) == 0) stop("No configurations completed successfully for Architecture E")

  list(results = results_df, architecture = "E",
       best_config = results_df[which.min(results_df$test_mae), ],
       param_grid = param_grid, site_names = site_names)
}


# ============================================================================
# CONVENIENCE WRAPPER -- run D, E, F for one outcome model in one call
# ============================================================================

#' Run Architecture F, E, and D hyperparameter exploration for a single
#' outcome model (e.g. IP or EP) and return one tidy, combined results
#' data frame plus the raw per-architecture outputs -- replaces the
#' repeated F/E/D block previously copy-pasted once per model in the
#' orchestrator script.
#'
#' @param model_label  short label stamped into the `model` column of the
#'   combined results (e.g. "IP", "EP")
#' @param ...  all other arguments are shared across F/E/D and forwarded
#'   as-is to ds.semiOPBARTExploreHyperF/E/D() -- see those for details
#'   (conns, data.name, outcome_col, levels, x_features, w_features,
#'   test_data.name, n_samp_range, n_burn_range, n_tree_range, seed_range,
#'   k_values, max_configs, nfilter, verbose)
#' @param n_swap_range  Architecture D only
#' @param combine_method  Architecture F only
#' @export
ds.semiOPBARTExploreHyperparameters <- function(model_label,
                                                 conns, data.name, outcome_col, levels,
                                                 x_features, w_features,
                                                 test_data.name = NULL,
                                                 n_samp_range = c(500, 1000, 2000),
                                                 n_burn_range = c(500, 1000, 2000),
                                                 n_tree_range = c(50, 100, 200),
                                                 seed_range = c(42, 123, 456),
                                                 k_values = c(1, 2, 3),
                                                 n_swap_range = c(2, 4, 8),
                                                sync_every= c(10,15),
                                                 combine_method = c(  "inverse_variance", "random_effects", "stack"), #c( "stack", "inverse_variance", "random_effects" ), #c( "random_effects", "stack", "inverse_variance"), #
                                                 max_configs = NULL,
                                                 nfilter =0,# 5,
                                                 verbose = TRUE) {

  combine_method <- match.arg(combine_method)

  if (verbose) message(sprintf(">>> %s: Architecture F (independent local fits)", model_label))
  # res_F <- ds.semiOPBARTExploreHyperF(
  #   conns = conns, data.name = data.name, outcome_col = outcome_col, levels = levels,
  #   x_features = x_features, w_features = w_features, test_data.name = test_data.name,
  #   n_samp_range = n_samp_range, n_burn_range = n_burn_range, n_tree_range = n_tree_range,
  #   seed_range = seed_range, k_values = k_values, max_configs = max_configs,
  #   nfilter = nfilter, combine_method = combine_method, verbose = verbose,
  #   state_name = ".semiOPBART_hyper")
  res_F <- NULL

#   combined <- dplyr::bind_rows(res_F$results, NULL, NULL)
#   combined$model <- model_label

#   list(
#     combined_results = combined,
#     by_architecture = list(F = res_F),
#     best_overall = combined[which.min(combined$test_mae), ],
#     best_by_architecture = do.call(rbind, lapply(split(combined, combined$architecture), function(d) {
#       d[which.min(d$test_mae), ]
#     }))
#   )
# }



if (verbose) message(sprintf(">>> %s: Architecture E (theta/us pooling, local trees)", model_label))
  res_E <- ds.semiOPBARTExploreHyperE(
    conns = conns, data.name = data.name, outcome_col = outcome_col, levels = levels,
    x_features = x_features, w_features = w_features, test_data.name = test_data.name,
    n_samp_range = n_samp_range, n_burn_range = n_burn_range, n_tree_range = n_tree_range,
    seed_range = seed_range,sync_every= sync_every, k_values = k_values, max_configs = max_configs,
    nfilter = nfilter, verbose = verbose,
    state_name = ".semiOPBART_hyper")

  # if (verbose) message(sprintf(">>> %s: Architecture D (tree swapping + theta/us pooling)", model_label))
  # res_D <- ds.semiOPBARTExploreHyperD(
  #   conns = conns, data.name = data.name, outcome_col = outcome_col, levels = levels,
  #   x_features = x_features, w_features = w_features, test_data.name = test_data.name,
  #   n_samp_range = n_samp_range, n_burn_range = n_burn_range, n_tree_range = n_tree_range,
  #   seed_range = seed_range, k_values = k_values, n_swap_range = n_swap_range,
  #   max_configs = max_configs, nfilter = nfilter, verbose = verbose,
  #   state_name = ".semiOPBART_hyper")
  #
  # combined <- dplyr::bind_rows(res_F$results, res_E$results, res_D$results)
  # combined$model <- model_label
  #
  # list(
  #   combined_results = combined,
  #   by_architecture = list(F = res_F, E = res_E, D = res_D),
  #   best_overall = combined[which.min(combined$test_mae), ],
  #   best_by_architecture = do.call(rbind, lapply(split(combined, combined$architecture), function(d) {
  #     d[which.min(d$test_mae), ]
  #   }))
  # )

  combined <- dplyr::bind_rows(res_F$results, res_E$results, NULL)
  combined$model <- model_label

  list(
    combined_results = combined,
    by_architecture = list(F = res_F, E = res_E, D = NULL),
    best_overall = combined[which.min(combined$test_mae), ],
    best_by_architecture = do.call(rbind, lapply(split(combined, combined$architecture), function(d) {
      d[which.min(d$test_mae), ]
    }))
  )
}




#' Fit the best Architecture E/F configuration and perform external prediction
#'
#' Uses the best configuration selected by the orchestrator, trains that
#' architecture on the trainable sites, chooses the trained site with the
#' largest number of observations as the reference site, exports that site's
#' forest together with its normalization reference, and performs prediction
#' at non-trained sites.
#'
#' Prediction is always returned. If validation_col is supplied, the function
#' also runs ds.semiOPBARTEvaluate() and returns the evaluation metrics.
#'
#' @param best_config One row from the orchestrator results, e.g. best_ip.
#' @param architecture "E" or "F". If NULL, uses best_config$architecture.
#' @param conns DataSHIELD connections.
#' @param data.name Training object name.
#' @param test_data.name Test/validation object name.
#' @param outcome_col Outcome column name.
#' @param validation_col Optional validation outcome column. If supplied,
#'   evaluation is performed.
#' @param levels Outcome levels.
#' @param x_features Predictor columns for Architecture F.
#' @param w_features Linear predictor columns for Architecture F.
#' @param combine_method Combination method for Architecture F.
#' @param nfilter Disclosure threshold.
#' @param warmup_burn Warm-up iterations for Architecture E.
#' @param state_prefix Server-side state prefix.
#' @param prediction_name Server-side prediction object name.
#' @param verbose Logical.
#'
#' @export
ds.semiOPBARTFitBestAndExternalPredict <- function(
    best_config,
    architecture = NULL,
    conns = NULL,
    data.name,
    test_data.name,
    outcome_col,
    validation_col = NULL,
    levels = 0:4,
    x_features = NULL,
    w_features = NULL,
    combine_method = "inverse_variance",
    nfilter = 5,
    warmup_burn = 1,
    state_prefix = ".semiOPBART_best_external",
    prediction_name = NULL,
    verbose = TRUE
) {

  if (is.null(conns)) {
    conns <- DSI::datashield.connections_find()
  }

  if (nrow(best_config) != 1L) {
    stop("best_config must contain exactly one configuration row")
  }

  if (is.null(architecture)) {
    if (!"architecture" %in% names(best_config)) {
      stop("best_config does not contain an 'architecture' column")
    }

    architecture <- as.character(best_config$architecture[[1]])
  }

  architecture <- toupper(architecture)

  if (!architecture %in% c("E", "F")) {
    stop(
      "Only Architecture E or F is supported, got '",
      architecture, "'"
    )
  }

  required_params <- c(
    "n_tree",
    "k",
    "n_burn",
    "n_samp",
    "seed"
  )

  missing_params <- setdiff(
    required_params,
    names(best_config)
  )

  if (length(missing_params)) {
    stop(
      "best_config is missing required parameters: ",
      paste(missing_params, collapse = ", ")
    )
  }

  # ---------------------------------------------------------------
  # 1. Restrict training to sites that actually have data.name
  # ---------------------------------------------------------------

  trainable_conns <- .filter_trainable_connections(
    conns = conns,
    train.name = data.name,
    verbose = verbose
  )

  if (!length(trainable_conns)) {
    stop("No trainable sites were found")
  }

  train_site_names <- names(trainable_conns)

  if (verbose) {
    message(
      "[BestExternal] Trainable sites: ",
      paste(train_site_names, collapse = ", ")
    )
  }

  # ---------------------------------------------------------------
  # 2. State name for this single final fit
  # ---------------------------------------------------------------

  config_id <- if ("config_id" %in% names(best_config)) {
    as.integer(best_config$config_id[[1]])
  } else {
    1L
  }

  config_state <- paste0(
    state_prefix,
    "_",
    architecture,
    "_config_",
    config_id
  )

  pred_obj <- if (is.null(prediction_name)) {
    paste0(config_state, "_external_pred")
  } else {
    prediction_name
  }

  # ---------------------------------------------------------------
  # 3. Extract selected hyperparameters
  # ---------------------------------------------------------------

  n_tree <- as.integer(best_config$n_tree[[1]])
  k <- as.numeric(best_config$k[[1]])
  n_burn <- as.integer(best_config$n_burn[[1]])
  n_samp <- as.integer(best_config$n_samp[[1]])
  seed <- as.integer(best_config$seed[[1]])

  # Architecture E has sync_every; F does not.
  sync_every <- if ("sync_every" %in% names(best_config)) {
    as.integer(best_config$sync_every[[1]])
  } else {
    1L
  }

  # ---------------------------------------------------------------
  # 4. Fit selected architecture
  # ---------------------------------------------------------------

  if (verbose) {
    message(
      "[BestExternal] Fitting Architecture ", architecture,
      " with n_tree=", n_tree,
      ", k=", k,
      ", burn=", n_burn,
      ", save=", n_samp,
      if (architecture == "E")
        paste0(", sync_every=", sync_every)
      else
        ""
    )
  }

  if (architecture == "E") {

    fit <- ds.semiOPBARTTrainE(
      data.name = data.name,
      datasources = trainable_conns,

      num_tree = n_tree,
      k = k,
      num_burn = n_burn,
      num_save = n_samp,

      sync_every = sync_every,
      warmup_burn = warmup_burn,

      nfilter_threshold = nfilter,
      state_name = config_state,
      seed = seed
    )

  } else {

    if (is.null(x_features) || is.null(w_features)) {
      stop(
        "x_features and w_features are required for Architecture F"
      )
    }

    formula_obj <- stats::as.formula(
      paste(
        outcome_col,
        "~",
        paste(x_features, collapse = " + ")
      )
    )

    linear_formula_obj <- stats::as.formula(
      paste(
        "~",
        paste(w_features, collapse = " + ")
      )
    )

    site_fits <- ds.semiOPBARTFitF(
      formula = formula_obj,
      linear_formula = linear_formula_obj,

      data.name = data.name,
      datasources = trainable_conns,

      num_tree = n_tree,
      k = k,
      num_burn = n_burn,
      num_save = n_samp,

      nfilter = nfilter,
      state_name = config_state,
      seed = seed
    )

    fit <- ds.semiOPBARTCombineF(
      site_fits,
      combine_method = combine_method
    )

    # Preserve the metadata needed below.
    attr(fit, "site_fits") <- site_fits
    attr(fit, "site_names") <- train_site_names
    attr(fit, "state_name") <- config_state
    attr(fit, "num_tree") <- n_tree
    attr(fit, "k") <- k
  }


print("fit")
print(fit)
  # ---------------------------------------------------------------
  # 5. Determine how many rows each trained site has
  # ---------------------------------------------------------------
site_n <- fit$site_n
  # if (architecture == "E") {

  #   if (is.null(fit$site_n)) {
  #     stop(
  #       "Architecture E fit does not contain site_n. ",
  #       "Modify ds.semiOPBARTTrainE() to return site_n for selecting ",
  #       "the largest trained reference site."
  #     )
  #   }

  #   site_n <- fit$site_n

  # } else {

  #   site_fits <- attr(fit, "site_fits")

  #   if (is.null(site_fits)) {
  #     stop(
  #       "Architecture F fit does not retain site_fits, so the largest ",
  #       "reference site cannot be identified."
  #     )
  #   }

  #   site_n <- vapply(
  #     site_fits,
  #     function(z) {

  #       if (!is.null(z$n)) {
  #         return(as.integer(z$n))
  #       }

  #       if (!is.null(z$n_total)) {
  #         return(as.integer(z$n_total))
  #       }

  #       NA_integer_
  #     },
  #     integer(1)
  #   )
    print("site_n")
    print(site_n)

    #names(site_n) <- names(trainable_conns)

    if (all(is.na(site_n))) {
      stop(
        "Could not determine the number of observations at Architecture F ",
        "training sites."
      )
    }
  #}

  # ---------------------------------------------------------------
  # 6. Largest trained site becomes reference site
  # ---------------------------------------------------------------

  # reference_site <- names(
  #   site_n[which.max(site_n)]
  # )
   reference_site <- names(
    trainable_conns[which.max(site_n)]
  )
  if (verbose) {
    message(
      "[BestExternal] Reference site: ",
      reference_site,
      " (n=", site_n[[reference_site]], ")"
    )
  }

  # ---------------------------------------------------------------
  # 7. Export forest + normalization from reference site
  # ---------------------------------------------------------------

  reference <- ds.semiOPBARTExportReferenceForest(
    fit = if (architecture == "E") {
      fit
    } else {
      list(
        site_names = train_site_names,
        state_name = config_state,
        num_tree = n_tree
      )
    },
    reference_site = reference_site,
    datasources = trainable_conns
  )
# print("reference")
# print(reference)
  # ---------------------------------------------------------------
  # 8. Select normalization reference
  #
  # For cross-site external prediction, global ECDF is preferred.
  # local_ecdf is only a fallback.
  # ---------------------------------------------------------------

  if (
    !is.null(reference$global_ecdf_Serialize) &&
    isTRUE(reference$has_global_ecdf)
  ) {

    normalization_used <- "global_ecdf"

  } else if (
    !is.null(reference$local_ecdf_Serialize) &&
    isTRUE(reference$has_local_ecdf)
  ) {

    normalization_used <- "local_ecdf"

  } else {

    stop(
      "Reference site exported neither global_ecdf nor local_ecdf. ",
      "External prediction cannot proceed."
    )
  }

  # ---------------------------------------------------------------
  # 9. Determine sites that did NOT train
  # ---------------------------------------------------------------

  external_site_names <- setdiff(
    names(conns),
    train_site_names
  )

  if (!length(external_site_names)) {

    if (verbose) {
      message(
        "[BestExternal] All sites trained; no external prediction sites."
      )
    }

    return(list(
      fit = fit,
      architecture = architecture,
      reference = reference,
      reference_site = reference_site,
      site_n = site_n,
      normalization_used = normalization_used,
      external_sites = character(0),
      predictions = NULL,
      evaluation = NULL
    ))
  }

  external_conns <- conns[external_site_names]

  # ---------------------------------------------------------------
  # 10. Obtain posterior draws
  # ---------------------------------------------------------------
  if (architecture == "E") {

    theta_draws <- fit$theta_draws
    us_draws <- fit$us_draws

  } else {

    # The combined F object should contain posterior draws.
    theta_draws <- fit$theta_draws
    us_draws <- fit$us_draws

    if (is.null(theta_draws) || is.null(us_draws)) {

      # Allow alternative naming used by some CombineF versions.
      theta_draws <- fit$theta
      us_draws <- fit$us
    }

    if (is.null(theta_draws) || is.null(us_draws)) {
      stop(
        "Architecture F combined fit does not contain theta/us posterior ",
        "draws required for external prediction."
      )
    }
  }

  # ---------------------------------------------------------------
  # 11. External prediction
  # ---------------------------------------------------------------

  if (verbose) {
    message(
      "[BestExternal] External prediction on: ",
      paste(external_site_names, collapse = ", ")
    )
  }

  prediction_result <- ds.semiOPBARTPredictExternal(
    fit = list(
      theta = theta_draws,
      us = us_draws,
      k = k
    ),
    reference = reference,
    data.name_test = test_data.name,
    outcome_col = outcome_col,

    newobj_pred = pred_obj,
    nfilter = nfilter,

    datasources = external_conns,
    seed = seed
  )
  print("prediction_result")
  print(prediction_result)
  # ---------------------------------------------------------------
  # 12. Evaluate only when validation_col is supplied
  # ---------------------------------------------------------------

  evaluation <- NULL
  #TODO: This is a bit of a hack. The prediction_result object should contain the truth column if it was supplied, but we also want to check if validation_col was supplied to this function. If both are present, we run evaluation.
  if (!is.null(validation_col) && !is.null(prediction_result$truth))  {

    # The server-side prediction object already contains truth.
    # validation_col is retained here as the explicit indication that
    # validation/evaluation was requested.
    evaluation <- ds.semiOPBARTEvaluate(
      pred_obj = pred_obj,
      levels = levels,
      nfilter = nfilter,
      datasources = external_conns,
      label = paste0(
        architecture,
        "_best_external_"
      )
    )
  }
  else {
    if (verbose) {
      message(
        "[BestExternal] No validation_col supplied; skipping evaluation."
      )
    }
  }

  # ---------------------------------------------------------------
  # 13. Always return predictions
  # ---------------------------------------------------------------

  list(
    fit = fit,

    architecture = architecture,

    config = best_config,

    reference = reference,
    reference_site = reference_site,
    reference_n = site_n[[reference_site]],

    site_n = site_n,

    train_sites = train_site_names,
    external_sites = external_site_names,

    normalization_used = normalization_used,

    prediction_object = pred_obj,
    prediction_result = prediction_result,

    evaluation = evaluation
  )
}